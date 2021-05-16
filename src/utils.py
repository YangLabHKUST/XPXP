import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import logging
import sys
from functools import reduce
from scipy.sparse.linalg import cg
from scipy import linalg
import os
import scipy.stats as st
import tempfile
import seaborn as sns
import matplotlib.pyplot as plt
import glob as gb


def configure_logging(logger_name): 
    LOG_LEVEL = logging.INFO
    log_filename = logger_name+'.log'
    importer_logger = logging.getLogger('importer_logger')
    importer_logger.setLevel(LOG_LEVEL)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')

    fh = logging.FileHandler(filename=log_filename)
    fh.setLevel(LOG_LEVEL)
    fh.setFormatter(formatter)
    importer_logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(LOG_LEVEL)
    sh.setFormatter(formatter)
    importer_logger.addHandler(sh)
    return importer_logger


def ReadPlink(plink_file,bim,dtype=np.float32):
    Genotype = read_plink1_bin(plink_file+".bed", plink_file+".bim", plink_file+".fam", verbose=False)
    Genotype = Genotype.where(Genotype.snp.isin(Genotype.snp.values[bim['index'].values]), drop=True)
    Genotype = Genotype.astype(np.int8)
    G_geno = Genotype.values
    G_geno[np.isnan(G_geno)] = 2
    G_geno = 2-G_geno
    return G_geno.astype(dtype)


def ScaleGenotype(G_geno,dtype=np.float32):
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    G_geno = ((G_geno-G_geno_mean)/G_geno_sd)/np.sqrt(G_geno.shape[1])
    return G_geno.astype(dtype), G_geno_mean/2, G_geno_sd


def GenotypeMoment(G_geno):
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    return G_geno_mean, G_geno_sd


def CrossCorrSS(z1,z2,K1,K2,K12,n1,n2,logger,Z1=None,Z2=None,group=None):
    m1 = K1.shape[0]
    m2 = K2.shape[0]
    p = z1.shape[0]
    ngroup = np.unique(group).shape[0]
    
    # calculate h1^2 for y1
    if Z1 is None:
        Z1 = np.ones((m1,1))
        M1 = np.diag(np.ones(m1))-np.full((m1,m1),1/m1)
        MK1 = K1
        MK12 = K12
    else:
        Z1 = np.hstack((np.ones((m1,1)),Z1))
        M1 = np.diag(np.ones(m1))-Z1.dot(np.linalg.inv(Z1.T.dot(Z1))).dot(Z1.T)
        MK1 = M1.dot(K1)
        MK12 = M1.dot(K12)
    q1 = Z1.shape[1]
    trK1 = MK1.trace()
    trSQK1 = trK1**2
    trK1SQ = (MK1**2).sum()
    S1 = (trK1SQ-trSQK1/(m1-q1)) / (m1-q1)**2
    zz1 = z1**2/n1-1/n1
    c1 = zz1.sum()/p
    h1 = c1/S1
    
    # calculate h2^2 for y2
    if Z2 is None:
        Z2 = np.ones((m2,1))
        M2 = np.diag(np.ones(m2))-np.full((m2,m2),1/m2)
        MK2 = K2
        K12M = K12
    else:
        Z2 = np.hstack((np.ones((m2,1)),Z2))
        M2 = np.diag(np.ones(m2))-Z2.dot(np.linalg.inv(Z2.T.dot(Z2))).dot(Z2.T)
        MK2 = M2.dot(K2)
        K12M = K12.dot(M2)
    q2 = Z2.shape[1]
    trK2 = MK2.trace()
    trSQK2 = trK2**2
    trK2SQ = (MK2**2).sum()
    S2 = (trK2SQ-trSQK2/(m2-q2)) / (m2-q2)**2
    zz2 = z2**2/n2-1/n2
    c2 = zz2.sum()/p
    h2 = c2/S2

    # safe guard 
    h1 = h1 if h1>1e-6 else 1e-6
    h2 = h2 if h2>1e-6 else 1e-6
    h1 = h1 if h1<1 else 0.99
    h2 = h2 if h2<1 else 0.99

    # calculate co-heritability
    if np.array_equal(K1,K2):
        S3 = S1
    else:
        trK12SQ = (MK12*K12M).sum()
        S3 = trK12SQ/(m1-q1)/(m2-q2)
    zz12 = z1*z2/np.sqrt(n1)/np.sqrt(n2)
    c3 = (zz12).sum()/p
    h12 = c3/S3
    rho = h12/np.sqrt(h1*h2)
    
    # safe guard 
    rho = rho if rho <= .95*(h12/(h1*h2)**.5) else .95*(h12/(h1*h2)**.5)
    rho = rho if rho >= -.95*(h12/(h1*h2)**.5) else -.95*(h12/(h1*h2)**.5)
    h12 = rho*np.sqrt(h1*h2)

    # calculate variance
    zj = np.array([[zz1[group==_].sum(),zz2[group==_].sum(),zz12[group==_].sum(),(group==_).sum()] for _ in range(ngroup)])
    c1_jf = (zz1.sum()-zj[:,0])/(p-zj[:,3])
    sd_h1 = (c1_jf.var()/S1/S1*(ngroup-1))**.5
    
    c2_jf = (zz2.sum()-zj[:,1])/(p-zj[:,3])
    sd_h2 = (c2_jf.var()/S2/S2*(ngroup-1))**.5
    
    c3_jf = (zz12.sum()-zj[:,2])/(p-zj[:,3])
    sd_h12 = (c3_jf.var()/S3/S3*(ngroup-1))**.5
    
    sd_rho = ((c3_jf/np.sqrt(c1_jf*c2_jf)).var() * S1*S2/S3/S3 * (ngroup-1))**.5
    res = pd.DataFrame([[h1,h2,h12,rho],[sd_h1,sd_h2,sd_h12,sd_rho]],index=['value','se'],columns=['h1','h2','h12','rho'])
    return res


def ld_clump(df,file_ref1,p1,r2=0.1):
    if 'P' not in df.columns:
        df['P'] = df['Z'].map(lambda x:st.norm.sf(abs(x))*2)
    ftmp = tempfile.NamedTemporaryFile()
    ftmp_res = tempfile.NamedTemporaryFile()
    df[['SNP','P']].to_csv(ftmp.name,sep='\t',index=None)
    clump_com = 'plink --bfile {} --clump {} --clump-p1 {} --clump-r2 {} --clump-kb 1000 --out {}'.format(file_ref1,
        ftmp.name,p1,r2,ftmp_res.name)
    res = os.system(clump_com)
    if res != 0:
        print('error in clumping')
        sys.exit()
    snps_clump = pd.read_csv(ftmp_res.name+'.clumped',delim_whitespace=True)['SNP'].values
    ftmp.close()
    ftmp_res.close()
    [os.remove(_) for _ in gb.glob(ftmp_res.name+'.*')]
    return snps_clump


def ImputeZscore(p,X,zsb):
    X = X/np.sqrt(X.shape[0])
    L = X.T.dot(X)
    Lp = L*p
    Lp[np.diag_indices_from(Lp)] += 0.1
    Sigmatt = Lp[~np.isnan(zsb)][:,~np.isnan(zsb)]
    Sigmait = Lp[np.isnan(zsb)][:,~np.isnan(zsb)]
    try:
        SigmattInvZ = cg(Sigmatt,zsb[~np.isnan(zsb)])[0]
    except:
        SigmattInv = linalg.cho_solve(linalg.cho_factor(Sigmatt,lower=True), np.eye(Sigmatt.shape[0]))
        SigmattInvZ = SigmattInv.dot(zsb[~np.isnan(zsb)])
    zsb[np.isnan(zsb)] = Sigmait.dot(SigmattInvZ)
    return zsb 


def ComputePMXPXP(bi, sums_names, sumst_ref1, ns, p, L, L2, Sigma_beta, Sigma_e, zsb, logger, use_cho=True, use_cg=True, return_LDpredinf=True):
    pb = L.shape[0] 
    ref1_n = len(sumst_ref1)

    # construct matrix A 
    Sigma_beta_inv = np.linalg.inv(Sigma_beta)
    ns_mat = np.zeros((len(sums_names),len(sums_names)))
    for i,sums_n in enumerate(sums_names):
        for j,sums_n2 in enumerate(sums_names):
            ns_mat[i,j] = (ns[sums_n]*ns[sums_n2])**.5
    Sigma_e_inv = np.linalg.inv(Sigma_e)
    Sigma_e_inv_ns = Sigma_e_inv*ns_mat
    A = np.zeros((pb*Sigma_beta.shape[0],pb*Sigma_beta.shape[0]))
    #A[:pb*ref1_n,:pb*ref1_n] = FastKron(Sigma_e_inv_ns[:ref1_n,:ref1_n],L)#np.kron(Sigma_e_inv_ns[:ref1_n,:ref1_n], L)
    #A[pb*ref1_n:,pb*ref1_n:] = FastKron(Sigma_e_inv_ns[ref1_n:,ref1_n:],L2)#np.kron(Sigma_e_inv_ns[ref1_n:,ref1_n:], L2)
    for i in range(ref1_n):
        for j in range(ref1_n):
            A[i*pb:(i+1)*pb,j*pb:(j+1)*pb] = np.multiply(L,Sigma_e_inv_ns[i,j]) #Sigma_e_inv_ns[i,j]*L
    
    for i in range(ref1_n,len(sums_names)):
        for j in range(ref1_n,len(sums_names)):
            A[i*pb:(i+1)*pb,j*pb:(j+1)*pb] = np.multiply(L2,Sigma_e_inv_ns[i,j]) #Sigma_e_inv_ns[i,j]*L2
            
    for i in range(len(sums_names)):
        for j in range(len(sums_names)):
            A[i*pb:(i+1)*pb,j*pb:(j+1)*pb][np.diag_indices(pb)] += Sigma_beta_inv[i,j]
    
    # Construct Z-score vector
    #zb = np.concatenate(tuple([zsb[_]*np.sqrt(ns[_]/p) for _ in sums_names]))
    zb = np.zeros(pb*len(sums_names))
    for i in range(len(sums_names)):
        tmp = np.zeros(pb)
        for j,sums_n in enumerate(sums_names):
            tmp += Sigma_e_inv[i,j]*zsb[sums_n]
        zb[pb*i:pb*(i+1)] = tmp
    
    mu_xpss, status = cg(A,zb,maxiter=1000)
    if status != 0:
        logger.warn("Does not converge at {} block, info: {}, SNPs num: {}. Swich to matrix inverse version (note: check the positive definiteness of SigmaBeta and SigamE; check the allele frequency of reference panel)".format(bi, status, pb))
        #print("Does not converge at {} block, info: {}, SNPs num: {}. Swich to matrix inverse version".format(bi, status, pb))
        mu_xpss = np.linalg.inv(A).dot(zb)
            
    # mu, mu_xpass for X1 
    if return_LDpredinf:
        mu = np.zeros((pb,len(sums_names)))
        for i,sums_n in enumerate(sums_names[:ref1_n]):
            S1 = L*ns[sums_n]
            h1 = Sigma_beta[i,i]
            S1[np.diag_indices_from(S1)] += (1-h1)/h1
            if use_cg:
                S1invz1 = cg(S1,zsb[sums_n])[0]
            else:
                S1invz1 = linalg.cho_solve(linalg.cho_factor(S1,lower=True), np.eye(pb)).dot(zsb[sums_n])
            mu[:,i] = np.sqrt(ns[sums_n]/p)*S1invz1
        for i,sums_n in enumerate(sums_names[ref1_n:]):
            S1 = L2*ns[sums_n]
            h1 = Sigma_beta[i+ref1_n,i+ref1_n]
            S1[np.diag_indices_from(S1)] += (1-h1)/h1
            if use_cg:
                S1invz1 = cg(S1,zsb[sums_n])[0]
            else:
                S1invz1 = linalg.cho_solve(linalg.cho_factor(S1,lower=True), np.eye(pb)).dot(zsb[sums_n])
            mu[:,i+ref1_n] = np.sqrt(ns[sums_n]/p)*S1invz1
    else:
        mu = np.zeros((pb,len(sums_names)))
    
    return mu, mu_xpss.reshape(len(sums_names),pb).T


def PlotCorrelatin(tme_corr,tme_corr_pvalue,label,ismask=False):
    tme_corr_anno = tme_corr.round(3).astype('str').values
    star_index = tme_corr_pvalue<=(0.05)
    tme_corr_anno[star_index] += ' *'
    star_star_index = tme_corr_pvalue<=(0.05/(tme_corr_pvalue.shape[0]))
    tme_corr_anno[star_star_index] += '*'
    sns.set_context('paper',font_scale=1.6) 
    plt.figure(figsize=(18,18),dpi=200)
    if ismask:
        mask = np.tril(np.ones_like(tme_corr, dtype=np.bool))
        ax = sns.heatmap(tme_corr.T, mask=mask.T, square = True, linewidths = 1,
                              cmap='coolwarm', center=0,alpha=0.9,
            cbar_kws={"shrink": 0.7,'label':label,"orientation": "vertical"},
                         annot=tme_corr_anno.T, fmt='',annot_kws={"size": 16})
    else:
        ax = sns.heatmap(tme_corr.T,square = True, linewidths = 1,
                              cmap='coolwarm', center=0,alpha=0.9,
            cbar_kws={"shrink": 0.7,'label':label,"orientation": "vertical"},
                         annot=tme_corr_anno.T, fmt='',annot_kws={"size": 16})
    ax.tick_params(axis='y',rotation=0)
    ax.tick_params(axis='x',rotation=90)
    ax.tick_params(axis='both', which='both', length=0)
    plt.show()
    

def cov_to_corr(df):
    columns = df.columns
    indexs = df.index
    values = df.values
    corr = np.zeros((df.shape[0],df.shape[1]))
    for i in range(df.shape[0]):
        for j in range(i+1,df.shape[0]):
            corr[i,j] = values[i,j]/(values[i,i]*values[j,j])**.5
    corr = corr.T+corr
    return pd.DataFrame(corr,columns=columns,index=indexs)
