import pandas as pd
import numpy as np
import copy
import argparse
from utils import *

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='Predict PRS with summary statistic')
    parser.add_argument('--ref_file', type=str, help='LD reference files path, plink1 file version')
    parser.add_argument('--sumst_file', type=str, help='input GWAS summmary statistic file', required=True)
    parser.add_argument('--beta', type=str, help='input esitmated SNP effect size file', required=True)
    parser.add_argument('--col_name', type=str, help='input col_name of beta file', required=True)
    parser.add_argument('--ld_block_file', type=str, help='ld block file path',\
        default='./XPXP_demo/EAS_fourier_ls-all.bed') # EUR_fourier_ls-all.bed 

    args = parser.parse_args()

    file_ref1 = args.ref_file 
    test_gwas_file = args.sumst_file
    ref1_info = pd.read_csv(file_ref1+'.bim',sep='\t',header=None)
    ref1_info.columns = ['chr','SNP','cm','bp','A1','A2']
    df_gwas = pd.read_csv(test_gwas_file,sep='\t')
    df_gwas = df_gwas.drop_duplicates('SNP') 
    df_pm = pd.read_csv(args.beta,sep='\t')

    snp_common = np.intersect1d(df_pm['SNP'].values,df_gwas['SNP'].values)
    print('Overlapping SNPs number: {}'.format(snp_common.shape[0]))
    df_gwas = df_gwas.loc[df_gwas['SNP'].isin(snp_common)]
    df_gwas = df_gwas.reset_index()
    df_pm = df_pm.loc[df_pm['SNP'].isin(snp_common)]
    df_pm = df_pm.reset_index()

    ref1_info = ref1_info.loc[ref1_info['SNP'].isin(snp_common)]
    ref1_info = ref1_info.reset_index()
    X1 = ReadPlink(file_ref1,ref1_info)
    X1, X1_maf, X1_sd = ScaleGenotype(X1)
    X1_stan = X1*X1.shape[1]**.5

    ld_block_file = args.ld_block_file
    block = pd.read_csv(ld_block_file,delim_whitespace=True)
    block['chr'] = block['chr'].map(lambda x:int(x.replace('chr','')))
    ngroup = block.shape[0]
    ref1_info_block = np.zeros(ref1_info.shape[0],dtype=int)
    for b_ in block.iterrows():
        ref1_info_block[ref1_info.loc[(ref1_info['chr']==b_[1]['chr'])&(ref1_info['bp']>=b_[1]['start'])&(ref1_info['bp']<=b_[1]['stop'])].index.values]= b_[0]
    ref1_info['block'] = ref1_info_block

    sums_inverse = (ref1_info['A1'].values != df_gwas['A1'].values)
    df_gwas.loc[sums_inverse,'Z'] = -df_gwas.loc[sums_inverse,'Z'] 
    n_test = df_gwas['N'].median()
    zs = copy.deepcopy(df_gwas['Z'].values)  

    beta1 = df_pm[args.col_name].values*X1_sd
    denominator = 0
    for i in range(ngroup):
        idx_i = ref1_info.loc[ref1_info['block']==i,].index.values
        X_stan_b = X1_stan[:,idx_i]
        L_b = X_stan_b.T@X_stan_b/X1_stan.shape[0]
        denominator += beta1[idx_i]@L_b@beta1[idx_i]
    corr = (beta1.dot(zs)/(n_test**.5))/((denominator)**.5)
    print('R2 for {}: {}'.format(args.col_name,corr**2))


