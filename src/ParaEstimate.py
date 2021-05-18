import pandas as pd
import numpy as np
import argparse
import sys
import os 
from utils import *


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='estimate parameters using GWAS summary')
    parser.add_argument('--save_dir', type=str, help='output path', required=True)
    parser.add_argument('--ldsc_path', type=str, help='LDSC path', required=True)
    parser.add_argument('--ldsc_files', type=str, \
        help='LDscore files path, seperated by comma.',required=True)    
    parser.add_argument('--merge_alleles', type=str, help='file used for matching alleles', required=True)
    parser.add_argument('--sumst_files', type=str, \
        help='summary statisitc files, separated by comma',required=True)
    parser.add_argument('--sumst_names', type=str, \
        help='summary statisitc names, separated by comma, the order is corresponds to the --sumst_files, different populations are separated by "+", .e.g. Target+Auxiliary',required=True)
    parser.add_argument('--ref_files', type=str, \
        help='LD reference files path, plink1 file version, seperated by comma.',required=True)
    parser.add_argument('--covar_files', type=str, \
        help='LD reference covariate files path, seperated by comma.',required=True)
    parser.add_argument('--ld_block_file', type=str,required=True)
    parser.add_argument('--num_threads', type=str, help='number of threads', default="22")

    args = parser.parse_args()

    BBJ_include_pheno = args.sumst_names.split('+')[0].split(',')
    UKB_include_pheno = args.sumst_names.split('+')[1].split(',')

    f_eas = args.sumst_files.split(',')[:len(BBJ_include_pheno)]
    f_eas_ldsc = []
    f_eur = args.sumst_files.split(',')[len(BBJ_include_pheno):]
    f_eur_ldsc = []

    # mounge SNPs
    ldsc_munge_com = 'conda run -n ldsc python {}/munge_sumstats.py \
        --chunksize 500000 \
        --merge-alleles {} \
        --sumstats {} \
        --out {}'
    for f in f_eas:
        f_eas_tmp = os.path.splitext(f)[0]
        f_eas_ldsc.append(f_eas_tmp+'.sumstats.gz')
        #os.system(ldsc_munge_com.format(args.ldsc_path,args.merge_alleles,f,f_eas_tmp))
    for f in f_eur:
        f_eur_tmp = os.path.splitext(f)[0]
        f_eur_ldsc.append(f_eur_tmp+'.sumstats.gz')
        #os.system(ldsc_munge_com.format(args.ldsc_path,args.merge_alleles,f,f_eur_tmp))


    genetic_corr_path = args.save_dir
    
    # Parameters: Trans
    if not os.path.exists(genetic_corr_path):
        os.mkdir(genetic_corr_path)
    if not os.path.exists(genetic_corr_path+'/trans'):
        os.mkdir(genetic_corr_path+'/trans')
    com = 'python {}/TransGC.py --save {}/trans/{}-{} --sumst_files {},{} --ref_files {} --covar_files {} --ld_block_file {}'
    coms = []
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(UKB_include_pheno):
            os.system(com.format(os.path.dirname(os.path.realpath(__file__)),genetic_corr_path,peas,peur,f_eas[i],f_eur[j],args.ref_files,args.covar_files,args.ld_block_file))
    
    # Parameters: Target
    ldsc_com_heri = 'conda run -n ldsc python {0}/ldsc.py --h2 {1} --ref-ld-chr {2} --w-ld-chr {2} --out {3}'
    ldsc_com_rg = 'conda run -n ldsc python {0}/ldsc.py --rg {1} --ref-ld-chr {2} --w-ld-chr {2} --out {3}'
    if not os.path.exists(genetic_corr_path+'/target'):
        os.mkdir(genetic_corr_path+'/target')
    for i,peas in enumerate(BBJ_include_pheno):
        os.system(ldsc_com_heri.format(args.ldsc_path,f_eas_ldsc[i],args.ldsc_files.split(',')[0],genetic_corr_path+'/target/'+peas))    
    for i,peas in enumerate(BBJ_include_pheno):
        for j,peur in enumerate(BBJ_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(args.ldsc_path,f_eas_ldsc[i]+','+f_eas_ldsc[j+i+1],args.ldsc_files.split(',')[0],genetic_corr_path+'/target/'+peas+'-'+peur))

    # Parameters: Auxiliary
    if not os.path.exists(genetic_corr_path+'/auxiliary'):
        os.mkdir(genetic_corr_path+'/auxiliary')
    for i,peas in enumerate(UKB_include_pheno):
        os.system(ldsc_com_heri.format(args.ldsc_path,f_eur_ldsc[i],args.ldsc_files.split(',')[1],genetic_corr_path+'/auxiliary/'+peas))
    
    for i,peas in enumerate(UKB_include_pheno):
        for j,peur in enumerate(UKB_include_pheno[(i+1):]):
            os.system(ldsc_com_rg.format(args.ldsc_path,f_eur_ldsc[i]+','+f_eur_ldsc[j+i+1],args.ldsc_files.split(',')[1],genetic_corr_path+'/auxiliary/'+peas+'-'+peur))


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
    corr_df.to_csv(genetic_corr_path+'/gcorr.csv')
    cov_df.to_csv(genetic_corr_path+'/gcov.csv')
    cov_e_df.to_csv(genetic_corr_path+'/ecov.csv')









        