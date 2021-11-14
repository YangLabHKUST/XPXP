import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import logging
import sys
sys.path.append('../XPXP/src')
from utils import *
from functools import reduce
import argparse
import os
from utils import *
from multiprocessing import Pool
import glob as gb


def worker(x):
    return os.system(x)


def get_correlated_traits_EAS(pheno_bbj,eas_cohri_path,corr_bbj_threshold):
    pheno_bbj_include = []
    for p in pheno_bbj:
        ps = gb.glob(eas_cohri_path+'*{}*.log'.format(p))
        for i in ps:
            if 'EAS' in i.split('/')[-1]:
                continue
            corr = float(open(i,'r').readlines()[-1].split()[-5])
            pval = float(open(i,'r').readlines()[-1].split()[-7])
            if abs(corr)>corr_bbj_threshold and pval<0.0001:
                #print(corr,pval,i)
                try:
                    pheno_bbj_include.append(list(set([i.split('/')[-1].split('_')[0],i.split('/')[-1].split('_')[-2]]).difference(pheno_bbj))[0])
                except:
                    pass
    pheno_bbj_include = sorted(list(set(pheno_bbj_include)))
    return pheno_bbj_include


def get_correlated_traits_EUR(pheno_ukb,eur_cohri_path,corr_ukb_threshold):
    pheno_ukb_include = []
    for p in pheno_ukb:
        ps = gb.glob(eur_cohri_path+'*{}*.log'.format(p))
        for i in ps:
            if 'EUR' in i.split('/')[-1] or 'Locke' in i.split('/')[-1] or 'finngen' in i.split('/')[-1]:
                continue
            corr = float(open(i,'r').readlines()[-1].split()[-5])
            pval = float(open(i,'r').readlines()[-1].split()[-7])
            if abs(corr)>corr_ukb_threshold and pval<0.0001:
                #print(corr,pval,i)
                try:
                    pheno_ukb_include.append(list(set([i.split('/')[-1].split('_')[0],i.split('/')[-1].split('_')[-2]]).difference(pheno_ukb))[0])
                except:
                    pass
    pheno_ukb_include = sorted(list(set(pheno_ukb_include)))
    return pheno_ukb_include


if __name__ == '__main__':
    pool=Pool(processes=4)
    eas_cohri_path = './gera/EAS_coheritability/'
    eur_cohri_path = './gera/EUR_coheritability/'
    ldsc_format_path = './gera/ldsc_format/'
    xpxp_path = '../XPXP/src'
    use_snps = './gera/used_SNPs_3M.txt'
    out_path = './gera/t2d-bbj-metabolicOnly'
    genetic_corr_path = out_path+'/genetic_corr'
    pheno_bbj = ['T2D']
    pheno_ukb = ['T2D']
    add_pheno = ['BMI']
    pheno_bbj_include = []
    pheno_ukb_include = []
    metabolic_pheno = ['TC','TG','HDL','LDL']
    
    BBJ_include_pheno = pheno_bbj + add_pheno + pheno_bbj_include + metabolic_pheno
    UKB_include_pheno = pheno_ukb + add_pheno + pheno_ukb_include
    EAS_pheno = [_+'-EAS' for _ in BBJ_include_pheno]
    EUR_pheno = [_+'-EUR' for _ in UKB_include_pheno]
    
    print(BBJ_include_pheno)
    print(UKB_include_pheno)
    
    f_eas = []
    for bs in pheno_bbj + add_pheno + pheno_bbj_include:
        f = gb.glob(ldsc_format_path+'{}_BBJ_summary_format_ldsc.txt'.format(bs))[0]
        f_eas.append(f)
    for bs in metabolic_pheno:
        f = gb.glob(ldsc_format_path+'{}_BBJ_summary_format_ldsc.txt'.format(bs))[0]
        f_eas.append(f)        
    
    f_eur = []
    for bs in UKB_include_pheno:
        f = gb.glob(ldsc_format_path+'{}_UKB_summary_format_ldsc.txt'.format(bs))[0]
        f_eur.append(f)        
    
    print(f_eas)
    print(f_eur)
    
    
    # Common SNPs
    
    
    ns_smt = {}
    ns_smt2 = {}
    snps_common = pd.read_csv(use_snps,header=None)
    snps_ref = pd.read_csv('/home/share/UKB/ld_ref_2k/height_affy_ldpred_ref_2000_noMHC.bim',sep='\t',header=None)
    snps_common = snps_common.merge(snps_ref[[1]],left_on=0,right_on=1,how='inner')
    for i,f in enumerate(f_eas):
        df_tmp = pd.read_csv(f,sep='\t')
        ns_smt[EAS_pheno[i]] = df_tmp['N'].median()
        ns_smt2[EAS_pheno[i]] = df_tmp['N'].median()
        snps_common = snps_common[[0]].merge(df_tmp[['SNP']],left_on=0,right_on='SNP',how='inner')
    for i,f in enumerate(f_eur):
        df_tmp = pd.read_csv(f,sep='\t')
        ns_smt2[EUR_pheno[i]] = df_tmp['N'].median()
        snps_common = snps_common[[0]].merge(df_tmp[['SNP']],left_on=0,right_on='SNP',how='inner')    
    
    if not os.path.exists(out_path):
        os.mkdir(out_path)
    use_snps = out_path+'/used_SNPs_{}.txt'.format(pheno_bbj[0])
    snps_common[[0]].to_csv(use_snps,index=None,header=None)
    
    # Genetic corr
    if not os.path.exists(genetic_corr_path):
        os.mkdir(genetic_corr_path)
        
    # Parameters: Trans
    if not os.path.exists(genetic_corr_path+'/trans'):
        os.mkdir(genetic_corr_path+'/trans')
    com = 'python {}/TranGC.py --use_snp {} --save {}/trans/{}_EAS-{}_EUR --sumst_files {},{}'
    coms = []
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(UKB_include_pheno):
            coms.append(com.format(xpxp_path,use_snps,genetic_corr_path,peas,peur,f_eas[i],f_eur[j]))
    output = pool.map(worker,coms)
    
   
    # Parameters: Target
    ldsc_com_heri = 'conda run -n ldsc python ldsc.py --h2 {0} --ref-ld-chr {1} --w-ld-chr {1} --out {2}'
    ldsc_com_rg = 'conda run -n ldsc python ldsc.py --rg {0} --ref-ld-chr {1} --w-ld-chr {1} --out {2}'
    if not os.path.exists(genetic_corr_path+'/target'):
        os.mkdir(genetic_corr_path+'/target')
    for i,peas in enumerate(BBJ_include_pheno):
        os.system(ldsc_com_heri.format(f_eas[i],./XPXP_demo/eas_ldscores/,genetic_corr_path+'/target/'+peas))    
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(BBJ_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(f_eas[i]+','+f_eas[j+i+1],./XPXP_demo/eas_ldscores/,genetic_corr_path+'/target/'+peas+'-'+peur))

    # Parameters: Auxiliary
    if not os.path.exists(genetic_corr_path+'/auxiliary'):
        os.mkdir(genetic_corr_path+'/auxiliary')
    for i,peas in enumerate(UKB_include_pheno):
        os.system(ldsc_com_heri.format(f_eur[i],./XPXP_demo/eur_ldscores/,genetic_corr_path+'/auxiliary/'+peas))
    
    for i,peas in enumerate(UKB_include_pheno):
        for j,peur in enumerate(UKB_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(f_eur[i]+','+f_eur[j+i+1],./XPXP_demo/eur_ldscores/,genetic_corr_path+'/auxiliary/'+peas+'-'+peur))


    # Organize result: heri
    eas_heri = np.zeros(len(BBJ_include_pheno))
    eur_heri = np.zeros(len(UKB_include_pheno))

    for i,peas in enumerate(BBJ_include_pheno):
        fc = '{}/target/{}.log'.format(genetic_corr_path,peas)
        eas_heri[i] = float(open(fc).readlines()[-7].split()[-2])
    for i,peur in enumerate(UKB_include_pheno):
        fc = '{}/auxiliary/{}.log'.format(genetic_corr_path,peur)
        eur_heri[i] = float(open(fc).readlines()[-7].split()[-2])
        
    # Organize result: trans-cov
    eas_eur_cov = np.zeros((len(BBJ_include_pheno),len(UKB_include_pheno)))
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(UKB_include_pheno):
            fc = '{}/trans/{}-{}.log'.format(genetic_corr_path,peas,peur)
            corr = float(open(fc).readlines()[-2].split()[-1])
            eas_eur_cov[i,j] = corr*(eas_heri[i]*eur_heri[j])**.5
            
    eas_eur_cov_df = pd.DataFrame(eas_eur_cov,index=BBJ_include_pheno,columns=UKB_include_pheno)
    eas_eur_cov_e_df = pd.DataFrame(np.zeros((len(BBJ_include_pheno),len(UKB_include_pheno))),index=BBJ_include_pheno,columns=UKB_include_pheno)

    # Organize result: eas-cov
    eas_eas_cov = np.zeros((len(BBJ_include_pheno),len(BBJ_include_pheno)))
    eas_eas_cov_e = np.zeros((len(BBJ_include_pheno),len(BBJ_include_pheno)))
    for i,peas in enumerate(BBJ_include_pheno):
        eas_eas_cov[i,i] = eas_heri[i]
        eas_eas_cov_e[i,i] = 1-eas_heri[i]
        for j,peas2 in enumerate(BBJ_include_pheno[i+1:]):
            fc = '{}/target/{}-{}.log'.format(genetic_corr_path,peas,peas2)
            for ln,line in enumerate(open(fc).readlines()):
                if 'gencov:' in line:
                    eas_eas_cov[i,i+1+j] = float(line.split()[-2])
                    break
            eas_eas_cov_e[i,i+1+j] = float(open(fc).readlines()[ln+2].split()[-2])-eas_eas_cov[i,i+1+j]
    eas_eas_cov = eas_eas_cov+eas_eas_cov.T-np.diag(np.diag(eas_eas_cov))
    eas_eas_cov_e = eas_eas_cov_e+eas_eas_cov_e.T-np.diag(np.diag(eas_eas_cov_e))
    eas_eas_cov_df = pd.DataFrame(eas_eas_cov,index=BBJ_include_pheno,columns=BBJ_include_pheno)
    eas_eas_cov_e_df = pd.DataFrame(eas_eas_cov_e,index=BBJ_include_pheno,columns=BBJ_include_pheno)   

    # Organize result: eur-cov
    eur_eur_cov = np.zeros((len(UKB_include_pheno),len(UKB_include_pheno)))
    eur_eur_cov_e = np.zeros((len(UKB_include_pheno),len(UKB_include_pheno)))
    for i,peas in enumerate(UKB_include_pheno):
        eur_eur_cov[i,i] = eur_heri[i]
        eur_eur_cov_e[i,i] = 1-eur_heri[i]
        for j,peas2 in enumerate(UKB_include_pheno[i+1:]):
            fc = '{}/auxiliary/{}-{}.log'.format(genetic_corr_path,peas,peas2)
            for ln,line in enumerate(open(fc).readlines()):
                if 'gencov:' in line:
                    eur_eur_cov[i,i+1+j] = float(line.split()[-2])
                    break
            eur_eur_cov_e[i,i+1+j] = float(open(fc).readlines()[ln+2].split()[-2])-eur_eur_cov[i,i+1+j]
    eur_eur_cov = eur_eur_cov+eur_eur_cov.T-np.diag(np.diag(eur_eur_cov))
    eur_eur_cov_e = eur_eur_cov_e+eur_eur_cov_e.T-np.diag(np.diag(eur_eur_cov_e))
    eur_eur_cov_df = pd.DataFrame(eur_eur_cov,index=UKB_include_pheno,columns=UKB_include_pheno)
    eur_eur_cov_e_df = pd.DataFrame(eur_eur_cov_e,index=UKB_include_pheno,columns=UKB_include_pheno)        

    cov_df = pd.concat((pd.concat((eas_eas_cov_df,eas_eur_cov_df),axis=1),pd.concat((eas_eur_cov_df.T,eur_eur_cov_df),axis=1)),axis=0)
    cov_e_df = pd.concat((pd.concat((eas_eas_cov_e_df,eas_eur_cov_e_df),axis=1),pd.concat((eas_eur_cov_e_df.T,eur_eur_cov_e_df),axis=1)),axis=0)

    # guarantee the positive definiteness of covariance matrix
    for c1 in BBJ_include_pheno+UKB_include_pheno:
        for c2 in BBJ_include_pheno+UKB_include_pheno:
            if c1==c2:
                continue
            else:
                if cov_df.loc[c1,c2] > 0.95*(cov_df.loc[c1,c1]*cov_df.loc[c2,c2])**.5:
                    cov_df.loc[c1,c2] = 0.95*(cov_df.loc[c1,c1]*cov_df.loc[c2,c2])**.5
                elif cov_df.loc[c1,c2] < -0.95*(cov_df.loc[c1,c1]*cov_df.loc[c2,c2])**.5:
                    cov_df.loc[c1,c2] = -0.95*(cov_df.loc[c1,c1]*cov_df.loc[c2,c2])**.5
                else:
                    pass
                if cov_e_df.loc[c1,c2] > 0.95*(cov_e_df.loc[c1,c1]*cov_e_df.loc[c2,c2])**.5:
                    cov_e_df.loc[c1,c2] = 0.95*(cov_e_df.loc[c1,c1]*cov_e_df.loc[c2,c2])**.5
                elif cov_e_df.loc[c1,c2] < -0.95*(cov_e_df.loc[c1,c1]*cov_e_df.loc[c2,c2])**.5:
                    cov_e_df.loc[c1,c2] = -0.95*(cov_e_df.loc[c1,c1]*cov_e_df.loc[c2,c2])**.5
                else:
                    pass
    
    corr_df = cov_to_corr(cov_df)
    corr_df.to_csv(genetic_corr_path+'/gcorr_{}.csv'.format(pheno_bbj[0]))
    cov_df.to_csv(genetic_corr_path+'/gcov_{}.csv'.format(pheno_bbj[0]))
    cov_e_df.to_csv(genetic_corr_path+'/ecov_{}.csv'.format(pheno_bbj[0]))
    

    # XPXP
    
    if not os.path.exists(out_path+'/xpxp_res'):
        os.mkdir(out_path+'/xpxp_res')
    if not os.path.exists(out_path+'/xpxp_res/predict'):
        os.mkdir(out_path+'/xpxp_res/predict')
    com_xpxp = 'python ../XPXP/src/XPXP.py \
    --num_threads 40 \
    --return_LDpredinf \
    --save {} \
    --gc_file {} \
    --ec_file {} \
    --sumst_files {} \
    --sumst_ss {} \
    --use_snp {} \
    --ld_block_file /home/share/xiaojs/database/prs/EAS_fourier_ls-all.bed'
    com_xpxp1 = com_xpxp.format(out_path+'/xpxp_res/MT-PM-{}.csv'.format(pheno_bbj[0]),
                                genetic_corr_path+'/gcov_{}.csv'.format(pheno_bbj[0]),
                                genetic_corr_path+'/ecov_{}.csv'.format(pheno_bbj[0]),
                                ','.join(f_eas)+','+','.join(f_eur),
                                ','.join(EAS_pheno)+'+'+','.join(EUR_pheno),
                                use_snps)
    print(com_xpxp1)
    os.system(com_xpxp1)
                 
    com_xpxp = 'python ../XPXP/src/XPXP.py \
    --num_threads 40 \
    --save {} \
    --gc_file {} \
    --ec_file {} \
    --sumst_files {} \
    --sumst_ss {} \
    --use_snp {} \
    --ld_block_file /home/share/xiaojs/database/prs/EAS_fourier_ls-all.bed'
    com_xpxp2 = com_xpxp.format(out_path+'/xpxp_res/MT-PM-{}-EASonly.csv'.format(pheno_bbj[0]),
                                genetic_corr_path+'/gcov_{}.csv'.format(pheno_bbj[0]),
                                genetic_corr_path+'/ecov_{}.csv'.format(pheno_bbj[0]),
                                ','.join(f_eas),
                                ','.join(EAS_pheno)+'+',
                                use_snps)
    print(com_xpxp2)
    os.system(com_xpxp2)
    
    com_xpxp = 'python ../XPXP/src/XPXP.py \
    --num_threads 40 \
    --save {} \
    --gc_file {} \
    --sumst_files {} \
    --sumst_ss {} \
    --use_snp {} \
    --ld_block_file /home/share/xiaojs/database/prs/EAS_fourier_ls-all.bed'
    com_xpxp3 = com_xpxp.format(out_path+'/xpxp_res/MT-PM-{}-only.csv'.format(pheno_bbj[0]),
                                genetic_corr_path+'/gcov_{}.csv'.format(pheno_bbj[0]),
                                f_eas[0]+','+f_eur[0],
                                EAS_pheno[0]+'+'+EUR_pheno[0],
                                use_snps)
    print(com_xpxp3)
    os.system(com_xpxp3)
    
    
    # XPXP predict
    
    xpxp1 = out_path+'/xpxp_res/MT-PM-{}.csv'.format(pheno_bbj[0])
    xpxp2 = out_path+'/xpxp_res/MT-PM-{}-EASonly.csv'.format(pheno_bbj[0])
    xpxp3 = out_path+'/xpxp_res/MT-PM-{}-only.csv'.format(pheno_bbj[0])
    
    com = 'plink \
    --bfile ./gera/geno/eas_qc3_3M \
    --score {0} 2 4 {1} header sum \
    --out {2}/xpxp_res/predict/{3}'
    coms = []
    header = open(xpxp1,'r').readline().strip().split('\t')
    header2 = open(xpxp2,'r').readline().strip().split('\t')
    header3 = open(xpxp3,'r').readline().strip().split('\t')
    for p in EAS_pheno+EUR_pheno:
        col = 1+header.index('{}-mu'.format(p))
        coms.append(com.format(xpxp1,col,out_path,p))
        col = 1+header.index('{}-muxpxp'.format(p))
        coms.append(com.format(xpxp1,col,out_path,p+'-xpxpall'))
    for p in EAS_pheno[:1]:
        col = 1+header2.index('{}-muxpxp'.format(p))
        coms.append(com.format(xpxp2,col,out_path,p+'-xpxpeas'))
        col = 1+header3.index('{}-muxpxp'.format(p))
        coms.append(com.format(xpxp3,col,out_path,p+'-xpxp'))
        
        
    pool=Pool(processes=10)
    output = pool.map(worker,coms)
    
    
    # smtpred EAS
    
    if not os.path.exists(out_path+'/smtpred'):
        os.mkdir(out_path+'/smtpred')
    if not os.path.exists(out_path+'/smtpred/predict'):
        os.mkdir(out_path+'/smtpred/predict')
    heri_smt = {}
    rgs_smt = {}
    
    for i,p1 in enumerate(EAS_pheno):
        for j,p2 in enumerate((EAS_pheno)[i+1:]):
            rgs_smt[p1+'+'+p2] = corr_df.loc[p1,p2]
            heri_smt[p1] = cov_df.loc[p1,p1]
            heri_smt[p2] = cov_df.loc[p2,p2]
    
    with open(out_path+'/smtpred/ns-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for p in EAS_pheno:
            f.write(p+'\t'+str(ns_smt[p])+'\n')
        f.close()
    with open(out_path+'/smtpred/h2s-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for p in EAS_pheno:
            f.write(p+'\t'+str(heri_smt[p])+'\n')
        f.close()
    with open(out_path+'/smtpred/rgs-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for i,p in enumerate(EAS_pheno):
            for p2 in EAS_pheno[i+1:]:
                f.write(p+'\t'+p2+'\t'+str(rgs_smt[p+'+'+p2])+'\n')
        f.close()
        
    scorefiles = []
    for p in EAS_pheno:
        scorefiles.append(out_path+'/xpxp_res/predict/'+p+'.profile')

    com = 'conda run -n ldsc python /home/share/xiaojs/software/smtpred/smtpred.py \
        --h2file {} \
        --rgfile {} \
        --nfile {} \
        --scorefiles {} \
        --out {}/smtpred/predict/{}-wMT-SBLUP \
        --alltraits \
        --blup'.format(out_path+'/smtpred/h2s-{}.txt'.format(pheno_bbj[0]),
                       out_path+'/smtpred/rgs-{}.txt'.format(pheno_bbj[0]),
                       out_path+'/smtpred/ns-{}.txt'.format(pheno_bbj[0]),
                       ' '.join(scorefiles),
                       out_path,pheno_bbj[0])
    
    print(com)
    os.system(com)
    
    # smtpred EAS+EUR
    
    if not os.path.exists(out_path+'/smtpred_all'):
        os.mkdir(out_path+'/smtpred_all')
    if not os.path.exists(out_path+'/smtpred_all/predict'):
        os.mkdir(out_path+'/smtpred_all/predict')
    heri_smt2 = {}
    rgs_smt2 = {}
    for i,p1 in enumerate(EAS_pheno+EUR_pheno):
        for j,p2 in enumerate((EAS_pheno+EUR_pheno)[i+1:]):
            rgs_smt2[p1+'+'+p2] = corr_df.loc[p1,p2]
            heri_smt2[p1] = cov_df.loc[p1,p1]
            heri_smt2[p2] = cov_df.loc[p2,p2]

    
    with open(out_path+'/smtpred_all/ns-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for p in EAS_pheno+EUR_pheno:
            f.write(p+'\t'+str(ns_smt2[p])+'\n')
        f.close()
    with open(out_path+'/smtpred_all/h2s-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for p in EAS_pheno+EUR_pheno:
            f.write(p+'\t'+str(heri_smt2[p])+'\n')
        f.close()
    with open(out_path+'/smtpred_all/rgs-{}.txt'.format(pheno_bbj[0]), 'w') as f:
        for i,p in enumerate(EAS_pheno+EUR_pheno):
            for p2 in (EAS_pheno+EUR_pheno)[i+1:]:
                f.write(p+'\t'+p2+'\t'+str(rgs_smt2[p+'+'+p2])+'\n')
        f.close()
        
    
    scorefiles2 = []
    for p in EAS_pheno+EUR_pheno:
        scorefiles2.append(out_path+'/xpxp_res/predict/'+p+'.profile')
        
    com = 'conda run -n ldsc python /home/share/xiaojs/software/smtpred/smtpred.py \
        --h2file {} \
        --rgfile {} \
        --nfile {} \
        --scorefiles {} \
        --out {}/smtpred_all/predict/{}-wMT-SBLUP \
        --alltraits \
        --blup'.format(out_path+'/smtpred_all/h2s-{}.txt'.format(pheno_bbj[0]),
                       out_path+'/smtpred_all/rgs-{}.txt'.format(pheno_bbj[0]),
                       out_path+'/smtpred_all/ns-{}.txt'.format(pheno_bbj[0]),
                       ' '.join(scorefiles2),
                       out_path,pheno_bbj[0])
    
    print(com)
    os.system(com)
    
    
