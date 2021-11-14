import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.stats import norm
import argparse


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='reformat ss')
    parser.add_argument('--inp', type=str, help='output file path', required=True)
    args = parser.parse_args()
    
    df = pd.read_csv(args.inp,sep='\t')
    
    bim = pd.read_csv('/import/home/share/xiaojs/cross-trait/simulation/genodata/wegene_merge_qc2_330k_noMHC_2kSample.bim',sep='\t',header=None)
    
    df = df.merge(bim[[0,1,3]],left_on='SNP',right_on=1,how='inner')
    
    se = 1/df['N']**.5
    beta = df['Z']*se
    df['beta'] = beta
    df['beta_se'] = se
    df = df[['SNP',0,3,'A2','A1','beta','beta_se','N']]
    df.columns = 'rsid,chr,pos,a0,a1,beta,beta_se,N'.split(',')
    df.to_csv(args.inp+'.ldpred2_format.txt',index=None)  