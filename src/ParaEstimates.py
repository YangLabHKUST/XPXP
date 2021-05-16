import pandas as pd
import numpy as np
import argparse
import sys
import os 
from utils import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='estimate parameters using GWAS summary')
    parser.add_argument('--save', type=str, help='output path', require=True)
    parser.add_argument('--sumst_files', type=str, \
        help='summary statisitc files, separated by comma',required=True)
    parser.add_argument('--sumst_names', type=str, \
        help='summary statisitc names, separated by comma, corresponding to the --sumst_files, different populations are separated by "+", .e.g. Target+Auxiliary',required=True)
    parser.add_argument('--ref_files', type=str, \
        help='LD reference files path, plink1 file version, seperated by comma.',\
        default='/home/share/UKB/1kg_ref/1000G.EAS.QC.hm3.ind,/home/share/UKB/1kg_ref/1000G.EUR.QC.hm3.ind')
    parser.add_argument('--covar_files', type=str, \
        help='LD reference covariate files path, seperated by comma.',\
        default='/home/share/UKB/1kg_ref/1000G.EAS.QC.hm3.ind.pc5.txt,/home/share/UKB/1kg_ref/1000G.EUR.QC.hm3.ind.pc20.txt')
    parser.add_argument('--ldsc_files', type=str, \
        help='LDscore files path, seperated by comma.',\
        default='/home/share/xiaojs/database/1kg_EAS_ldsc/eas_ldscores/,/home/share/xiaojs/database/1kg_EUR/eur_w_ld_chr/')    
    parser.add_argument('--ld_block_file', type=str, \
        help='ld block file path',default='/home/share/xiaojs/database/prs/EAS_fourier_ls-all.bed')
    parser.add_argument('--num_threads', type=str, help='number of threads', default="22")

    args = parser.parse_args()

    BBJ_include_pheno = args.sumst_names.split('+')[0].split(',')
    UKB_include_pheno = args.sumst_names.split('+')[1].split(',')
    f_eas = args.sumst_files.split(',')[:len(BBJ_include_pheno)]
    f_eur = args.sumst_files.split(',')[len(BBJ_include_pheno):]
    genetic_corr_path = args.save

    if not os.path.exists(genetic_corr_path):
        os.mkdir(genetic_corr_path)
    if not os.path.exists(genetic_corr_path+'/trans'):
        os.mkdir(genetic_corr_path+'/trans')
    com = 'python TransGC.py --save {}/trans/{}-{} --sumst_files {},{}'
    coms = []
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(UKB_include_pheno):
            os.system(com.format(genetic_corr_path,peas,peur,f_eas[i],f_eur[j]))

    ldsc_com_heri = 'conda run -n ldsc python /home/share/xiaojs/software/ldsc/ldsc.py --h2 {0} --ref-ld-chr {1} --w-ld-chr {1} --out {2}'
    ldsc_com_rg = 'conda run -n ldsc python /home/share/xiaojs/software/ldsc/ldsc.py --rg {0} --ref-ld-chr {1} --w-ld-chr {1} --out {2}'
    if not os.path.exists(genetic_corr_path+'/target'):
        os.mkdir(genetic_corr_path+'/target')
    for i,peas in enumerate(BBJ_include_pheno):
        os.system(ldsc_com_heri.format(f_eas[i],args.ldsc_files.split(',')[0],genetic_corr_path+'/target/'+peas))
    
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(BBJ_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(f_eas[i]+','+f_eas[j+i+1],args.ldsc_files.split(',')[0],genetic_corr_path+'/target/'+peas+'-'+peur))

    if not os.path.exists(genetic_corr_path+'/auxiliary'):
        os.mkdir(genetic_corr_path+'/auxiliary')
    for i,peas in enumerate(UKB_include_pheno):
        os.system(ldsc_com_heri.format(f_eur[i],args.ldsc_files.split(',')[1],genetic_corr_path+'/auxiliary/'+peas))
    
    for i,peas in enumerate(UKB_include_pheno):
        for j,peur in enumerate(UKB_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(f_eur[i]+','+f_eur[j+i+1],args.ldsc_files.split(',')[1],genetic_corr_path+'/auxiliary/'+peas+'-'+peur))

    
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
            tmp_corr = eas_eas_cov[i,i+1+j]/(eas_heri[i]*eas_heri[i+1+j])**.5
            eas_eas_cov_e[i,i+1+j] = (float(open(fc).readlines()[ln+2].split()[-2])-\
                                    tmp_corr)*((1-eas_heri[i])*(1-eas_heri[i+1+j]))**.5
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
            for line in open(fc).readlines():
                if 'Genetic Correlation:' in line:
                    tmp_corr = float(line.split()[-2])
                    break
            eur_eur_cov_e[i,i+1+j] = (float(open(fc).readlines()[ln+2].split()[-2])-\
                                    tmp_corr)/((1-eur_heri[i])*(1-eur_heri[i+1+j]))**.5
    eur_eur_cov = eur_eur_cov+eur_eur_cov.T-np.diag(np.diag(eur_eur_cov))
    eur_eur_cov_e = eur_eur_cov_e+eur_eur_cov_e.T-np.diag(np.diag(eur_eur_cov_e))
    eur_eur_cov_df = pd.DataFrame(eur_eur_cov,index=UKB_include_pheno,columns=UKB_include_pheno)
    eur_eur_cov_e_df = pd.DataFrame(eur_eur_cov_e,index=UKB_include_pheno,columns=UKB_include_pheno)        

    cov_df = pd.concat((pd.concat((eas_eas_cov_df,eas_eur_cov_df),axis=1),pd.concat((eas_eur_cov_df.T,eur_eur_cov_df),axis=1)),axis=0)
    cov_e_df = pd.concat((pd.concat((eas_eas_cov_e_df,eas_eur_cov_e_df),axis=1),pd.concat((eas_eur_cov_e_df.T,eur_eur_cov_e_df),axis=1)),axis=0)

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
    corr_df.to_csv(genetic_corr_path+'/gcorr.csv')
    cov_df.to_csv(genetic_corr_path+'/gcov.csv')
    cov_e_df.to_csv(genetic_corr_path+'/ecov.csv')









        