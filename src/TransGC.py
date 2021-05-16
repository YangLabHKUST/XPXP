import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import sys
import argparse
import os
from utils import *


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='estimate trans-ethinic genetic correlation using GWAS summary')
    parser.add_argument('--save', type=str, help='output log path')
    parser.add_argument('--sumst_files', type=str, help='summary statisitc files, seperated by comma')
    parser.add_argument('--use_snps', type=str, help='use assigned SNPs, skip SNPs matching step')
    parser.add_argument('--ref_files', type=str,\
        help='LD reference files path, plink1 file version, seperated by comma.',\
        default='/home/share/UKB/1kg_ref/1000G.EAS.QC.hm3.ind,/home/share/UKB/1kg_ref/1000G.EUR.QC.hm3.ind')
    parser.add_argument('--covar_files', type=str,\
        help='LD reference covariate files path, seperated by comma.',\
        default='/home/share/UKB/1kg_ref/1000G.EAS.QC.hm3.ind.pc5.txt,/home/share/UKB/1kg_ref/1000G.EUR.QC.hm3.ind.pc20.txt')
    parser.add_argument('--ld_block_file', type=str, help='ld block file path',default='/home/share/xiaojs/database/prs/EUR_fourier_ls-all.bed')
    parser.add_argument('--num_threads', type=str, help='number of threads', default="22")


    args = parser.parse_args()

    os.environ["OMP_NUM_THREADS"] = args.num_threads
    os.environ["OPENBLAS_NUM_THREADS"] = args.num_threads
    os.environ["MKL_NUM_THREADS"] = args.num_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = args.num_threads
    os.environ["NUMEXPR_NUM_THREADS"] = args.num_threads

    logger = configure_logging(args.save)
    sumst_files = args.sumst_files.split(',')
    file_ref1, file_ref2 = args.ref_files.split(',')
    file_cov1, file_cov2 = args.covar_files.split(',')
    ld_block_file = args.ld_block_file
    
    logger.info("Plink reference 1 file path: {} ".format(file_ref1))
    logger.info("Plink reference 2 file path: {} ".format(file_ref2))
    logger.info("Reference covariate 1 file path: {} ".format(file_cov1))
    logger.info("Reference covariate 2 file path: {} ".format(file_cov2))
    logger.info("LD block file path: {} ".format(ld_block_file))
    logger.info("log path: {} ".format(args.save+'.log'))
    for sumsf in sumst_files:
        logger.info("GWAS summary statistic file path: {} ".format(sumsf)) 

    zf1 = pd.read_csv(sumst_files[0],sep='\t')
    zf1 = zf1.drop_duplicates('SNP')
    zf2 = pd.read_csv(sumst_files[1],sep='\t')
    zf2 = zf2.drop_duplicates('SNP')
    logger.info("GWAS summary statistic 1 shape: {}x{} ".format(zf1.shape[0],zf1.shape[1])) 
    logger.info("GWAS summary statistic 2 shape: {}x{} ".format(zf2.shape[0],zf2.shape[1])) 

    ref1_info = pd.read_csv(file_ref1+'.bim',sep='\t',header=None)
    ref1_info.columns = ['chr','SNP','cm','bp','A1','A2']
    ref2_info = pd.read_csv(file_ref2+'.bim',sep='\t',header=None)
    ref2_info.columns = ['chr','SNP','cm','bp','A1','A2']
    logger.info("Ref1 bim shape: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    logger.info("Ref2 bim shape: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1])) 

    # matching SNPs
    if args.use_snps is None:
        logger.info('Matching SNPs from LD reference and GWAS summary...')
        snps_lst = [set(ref1_info['SNP'].values), set(ref2_info['SNP'].values),\
            set(zf1['SNP'].values), set(zf2['SNP'].values)]
        snps_common = list(snps_lst[0].intersection(*snps_lst))
    else:
        logger.info('load common SNPs from {}'.format(args.use_snps))
        snps_common = np.loadtxt(args.use_snps,dtype=str)

    ref1_info = ref1_info.loc[ref1_info['SNP'].isin(snps_common)]
    ref2_info = ref2_info.loc[ref2_info['SNP'].isin(snps_common)]
    zf1 = ref1_info[['SNP']].merge(zf1,left_on='SNP',right_on='SNP',how='left')
    zf2 = ref1_info[['SNP']].merge(zf2,left_on='SNP',right_on='SNP',how='left')
    logger.info("Sumstats 1 shape after matching: {}x{}".format(zf1.shape[0],zf1.shape[1]))
    logger.info("Sumstats 2 shape after matching: {}x{}".format(zf2.shape[0],zf2.shape[1]))
    logger.info("Ref1 bim shape after matching: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    logger.info("Ref2 bim shape after matching: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1])) 

    # replace T with A, replace G with C; A=1, C=2
    ref1_info['A1_int'] = ref1_info['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    zf1['A1_int'] = zf1['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref1_info['A2_int'] = ref1_info['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    zf1['A2_int'] = zf1['A2'].map(lambda x: 1 if x in ['A','T'] else 2)

    ref2_info['A1_int'] = ref2_info['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    zf2['A1_int'] = zf2['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref2_info['A2_int'] = ref2_info['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    zf2['A2_int'] = zf2['A2'].map(lambda x: 1 if x in ['A','T'] else 2)

    snps_rm = (
    (ref1_info['A1_int'].values+ref1_info['A2_int'].values != zf1['A1_int'].values+zf1['A2_int'].values) | \
    (ref2_info['A1_int'].values+ref2_info['A2_int'].values != zf2['A1_int'].values+zf2['A2_int'].values) | \
    (ref1_info['A1_int'].values+ref1_info['A2_int'].values != ref2_info['A1_int'].values+ref2_info['A2_int'].values)) 

    # remove ambiguous SNPs
    ref1_info = ref1_info.reset_index()
    ref1_info = ref1_info.loc[~snps_rm,:]
    ref2_info = ref2_info.reset_index()
    ref2_info = ref2_info.loc[~snps_rm,:]
    zf1 = zf1.reset_index()
    zf1 = zf1.loc[~snps_rm,:]
    zf2 = zf2.reset_index()
    zf2 = zf2.loc[~snps_rm,:]
    logger.info("Sumstats 1 shape after removing ambiguous SNPs: {}x{}".format(zf1.shape[0],zf1.shape[1]))
    logger.info("Sumstats 2 shape after removing ambiguous SNPs: {}x{}".format(zf2.shape[0],zf2.shape[1]))
    logger.info("Ref1 bim shape after removing ambiguous SNPs: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    logger.info("Ref2 bim shape after removing ambiguous SNPs: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1])) 

    ref1_info = ref1_info.reset_index()
    ref2_info = ref2_info.reset_index()
    zf1 = zf1.reset_index()
    zf2 = zf2.reset_index()

    # based on the minor allel of reference panel, correct for the z score in summary statistics
    logger.info('Based on the minor allel of ref1 panel, correct for the z score...')
    logger.info("Z score inner product before correction: {}".format(zf1['Z'].values@zf2['Z'].values/zf1.shape[0]))
    ind1 = (ref1_info['A1'].values != zf1['A1'].values)
    ind2 = (ref1_info['A1'].values != zf2['A1'].values)
    z_score1 = copy.deepcopy(zf1['Z'].values)
    z_score2 = copy.deepcopy(zf2['Z'].values)
    logger.info("SNPs need to be aligned for z1: {}".format(ind1.sum()))
    logger.info("SNPs need to be aligned for z2: {}".format(ind2.sum()))
    z_score1[ind1] = -z_score1[ind1]
    z_score2[ind2] = -z_score2[ind2]
    logger.info("Z score inner product after correction: {}".format(z_score1.dot(z_score2)/z_score1.shape[0]))

    logger.info('Loading first reference genotype...')
    X1 = ReadPlink(file_ref1,ref1_info)
    logger.info("Reference 1 genotype shape: {}x{} ".format(X1.shape[0],X1.shape[1]))
    logger.info('Loading second reference genotype...')
    X2 = ReadPlink(file_ref2,ref2_info)
    logger.info("Reference 2 genotype shape: {}x{} ".format(X2.shape[0],X2.shape[1]))

    ind_ref = (ref1_info['A1']!=ref2_info['A1']).values
    logger.info("{} SNPs need to be flipped in ref2 to match ref1".format(ind_ref.sum()))
    if ind_ref.sum()>0:
        X2[:,ind_ref] = 2-X2[:,ind_ref]
    
    logger.info('Normalizing (zero mean, 1/p variance) reference genotype...')
    X1, X1_maf, X1_sd = ScaleGenotype(X1)
    X2, X2_maf, X2_sd = ScaleGenotype(X2)

    K1 = X1.dot(X1.T)
    K2 = X2.dot(X2.T)
    K12 = X1.dot(X2.T)

    block = pd.read_csv(ld_block_file,delim_whitespace=True)
    block['chr'] = block['chr'].map(lambda x:int(x.replace('chr','')))
    ngroup = block.shape[0]
    logger.info("{} blocks loaded from {}".format(ngroup, ld_block_file))
    ref1_info_block = np.zeros(ref1_info.shape[0],dtype=int)
    for b_ in block.iterrows():
        ref1_info_block[ref1_info.loc[(ref1_info['chr']==b_[1]['chr'])&(ref1_info['bp']>=b_[1]['start'])&(ref1_info['bp']<=b_[1]['stop'])].index.values]= b_[0]
    ref1_info['block'] = ref1_info_block

    if file_cov1 == '':
        cov1 = None
    else:
        cov1 = pd.read_csv(file_cov1,delim_whitespace=True,header=None).values
    if file_cov2 == '':
        cov2 = None
    else:
        cov2 = pd.read_csv(file_cov2,delim_whitespace=True,header=None).values

    df_corr = CrossCorrSS(z_score1,z_score2,K1,K2,K12,zf1['N'].median(),zf2['N'].median(),logger,Z1=cov1,Z2=cov2,group=ref1_info['block'].values)
    logger.info(df_corr.to_string().replace('\n', '\n\t'))
    df_corr.to_csv(args.save+'.csv',sep='\t')



