import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import sys
import copy
import argparse
import os
from utils import *

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='impute Z scores for untyped SNPs')
    parser.add_argument('--sumst_file', type=str, help='summary statisitc file',required=True)
    parser.add_argument('--use_snp', type=str, help='SNPs need to be imputed')  
    parser.add_argument('--ref_file', type=str, help='LD reference files path, plink1 file version',required=True)
    parser.add_argument('--ld_block_file', type=str, help='ld block file path',\
        default='/home/jxiaoae/database/EAS_fourier_ls-all.bed')
    parser.add_argument('--snp', type=str, help='Name of SNP column', default='SNP')
    parser.add_argument('--N', type=str, help='Name of N column', default='N')
    parser.add_argument('--assignN', type=int, help='assign N for GWAS summary', default=0)
    parser.add_argument('--a1', type=str, help='Name of effect allele column', default='A1')
    parser.add_argument('--a2', type=str, help='Name of non-effect allele column', default='A2')
    parser.add_argument('--Z', type=str, help='Name of Z-score column', default='Z')
    parser.add_argument('--beta', type=str, help='Name of beta column', default='BETA')
    parser.add_argument('--se', type=str, help='Name of se column', default='SE')

    '''
    usage:
    python /home/jxiaoae/cross-popu/src/Impute_Zscore.py \
        --a1 ALT --a2 REF \
        --sumst_file /import/home/share/xiaojs/cross-trait/sums-stats/bbj_sex_stratified/BMI/Female_2017_BMI_BBJ_autosome.txt.gz 
    '''

    args = parser.parse_args() 

    save = '.'.join(args.sumst_file.split('.')[:-1])+'.3MImp'
    logger = configure_logging(save)
    logger.info("GWAS summary statistic file path: {} ".format(args.sumst_file))
    if args.use_snp is None:
        logger.warning('Common SNPs are not provided, using merged SNPs')
    else:
        logger.info('load common SNPs from {}'.format(args.use_snp))
    logger.info("Plink reference file path: {} ".format(args.ref_file))
    logger.info("LD block file path: {} ".format(args.ld_block_file))
    logger.info("output path: {} ".format(save))

    if args.sumst_file.endswith('gz') or args.sumst_file.endswith('zip'):
        df = pd.read_csv(args.sumst_file,compression='infer',sep='\t')
    else:
        df = pd.read_csv(args.sumst_file,sep='\t')
    if not args.Z in df.columns:
        df[args.Z] = df.apply(lambda x:x[args.beta]/x[args.se], axis=1)
    if not args.N in df.columns:
        if args.assignN == 0:
            logger.error("Need assign N")
        logger.info("Assigned N: {}".format(args.assignN))
        df['N'] = args.assignN
    df = df[[args.snp, args.N, args.Z, args.a1, args.a2]]
    df.columns = ['SNP','N','Z','A1','A2']

    df_columns = df.columns.values.tolist()
    ref1_info = pd.read_csv(args.ref_file+'.bim',sep='\t',header=None)
    ref1_info.columns = ['chr','SNP','cm','bp','A1','A2']
    logger.info("Ref bim shape: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    logger.info("GWAS summary statistic shape: {}x{} ".format(df.shape[0],df.shape[1])) 

    if args.use_snp is None:
        snp_common = np.intersect1d(df['SNP'].values,ref1_info['SNP'].values)
    else:
        snp_common = np.loadtxt(args.use_snp,dtype=str)
    
    df = df.loc[df['SNP'].isin(snp_common)]
    ref1_info = ref1_info.loc[ref1_info['SNP'].isin(snp_common)]
    ref1_info = ref1_info.reset_index()
    logger.info("Ref bim shape after matching: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    logger.info("GWAS summary statistic after matching: {}x{}".format(df.shape[0],df.shape[1]))

    df = ref1_info[['SNP','A1','A2']].merge(df,left_on='SNP',right_on='SNP',how='left')
    df.loc[df['A1_y'].isnull(),'A1'] = df.loc[df['A1_y'].isnull(),'A1_x']
    df.loc[~df['A1_y'].isnull(),'A1'] = df.loc[~df['A1_y'].isnull(),'A1_y']
    df.loc[df['A2_y'].isnull(),'A2'] = df.loc[df['A2_y'].isnull(),'A2_x']
    df.loc[~df['A2_y'].isnull(),'A2'] = df.loc[~df['A2_y'].isnull(),'A2_y']
    df['N']= df['N'].median()
    df = df[df_columns]
    logger.info("GWAS summary statistic after filling with NAN: {}x{}".format(df.shape[0],df.shape[1]))

    logger.info('Based on the minor allel of ref panel, correct for the z score...')
    # based on the minor allel of reference panel, correct for the z score in summary statistics
    sums_inverse = (ref1_info['A1'].values != df['A1'].values)
    df.loc[sums_inverse,'Z'] = -df.loc[sums_inverse,'Z']
    df_A1_store = df.loc[sums_inverse,'A1']
    df.loc[sums_inverse,'A1'] = df.loc[sums_inverse,'A2']
    df.loc[sums_inverse,'A2'] = df_A1_store
    logger.info("SNPs need aligned for: {}".format(sums_inverse.sum()))

    if not df['Z'].isnull().any():
        logger.info('No need to impute Z-score, All Z-score exist!')
        df.to_csv(save+'.txt',sep='\t',index=None)
        sys.exit()

    logger.info('Loading reference genotype...')
    X1 = ReadPlink(args.ref_file,ref1_info,np.int8)
    logger.info("Reference genotype shape: {}x{} ".format(X1.shape[0],X1.shape[1])) 
    #X1, X1_maf, X1_sd = ScaleGenotype(X1)
    X1_mean, X1_sd = GenotypeMoment(X1)

    p = ref1_info.shape[0]

    block = pd.read_csv(args.ld_block_file,delim_whitespace=True)
    block['chr'] = block['chr'].map(lambda x:int(x.replace('chr','')))
    ngroup = block.shape[0]
    logger.info("{} blocks loaded from {}".format(ngroup, args.ld_block_file))
    ref1_info_block = np.zeros(ref1_info.shape[0],dtype=int)
    for b_ in block.iterrows():
        ref1_info_block[ref1_info.loc[(ref1_info['chr']==b_[1]['chr'])&(ref1_info['bp']>=b_[1]['start'])&(ref1_info['bp']<=b_[1]['stop'])].index.values]= b_[0]
    ref1_info['block'] = ref1_info_block

    for i in range(ngroup):
        if i % 100 == 0:
            logger.info("Fininshing  {} block".format(i))
        idx_i = ref1_info.loc[ref1_info['block']==i,].index.values
        zs_b = df.loc[idx_i,'Z'].values
        X1_b = ((X1[:,idx_i]-X1_mean[idx_i])/X1_sd[idx_i])/np.sqrt(X1.shape[1])
        if np.isnan(zs_b).any():
            zs_b_imp = ImputeZscore(p,X1_b,zs_b)
            df.loc[idx_i,'Z'] = zs_b_imp
        else:
            continue
    df.to_csv(save+'.txt',sep='\t',index=None)
    