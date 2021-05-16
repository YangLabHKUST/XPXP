import pandas as pd
import sys
import tempfile
import argparse
import os
from scipy.stats.stats import pearsonr


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='Predict PRS from individual level data with plink')
    parser.add_argument('--save', type=str, help='output file path', required=True)
    parser.add_argument('--geno', type=str, help='input genotype file, plink1 version', required=True)
    parser.add_argument('--pheno', type=str, help='input true phenotype file, tsv format', required=False)
    parser.add_argument('--col_name', type=str, help='column of input true phenotype file, tsv format', required=False)
    parser.add_argument('--beta', type=str, help='input esitmated SNP effect size file', required=True)
    args = parser.parse_args()

    if args.pheno is not None:
        pheno = pd.read_csv(args.pheno,sep='\t')

    plink_com = 'plink \
        --bfile {0} \
        --score {1} 2 4 {2} header sum \
        --out {3}'

    header = open(args.beta,'r').readline().strip().split('\t')
    ftmp = tempfile.NamedTemporaryFile()
    res = os.system(plink_com.format(args.geno,args.beta,6,ftmp.name))
    if res != 0:
        print('error in predicting')
        sys.exit()
    predict_res = pd.read_csv(ftmp.name+'.profile',delim_whitespace=True)
    predict_res = predict_res[['FID','IID','SCORESUM']]
    predict_res.columns = ['FID','IID',header[5]]
    ftmp.close()
    if args.pheno is not None:
            predict_res_tmp = predict_res.merge(pheno[['IID',args.col_name]],left_on='IID',right_on='IID',how='inner')
            corr = pearsonr(predict_res_tmp[header[5]],predict_res_tmp[args.col_name])[0]
            print('R2 for {}: {}'.format(header[5],corr**2))
    for col,trait in enumerate(header[6:]):
        ftmp = tempfile.NamedTemporaryFile()
        res = os.system(plink_com.format(args.geno,args.beta,col+7,ftmp.name))
        predict_res_tmp = pd.read_csv(ftmp.name+'.profile',delim_whitespace=True)
        predict_res[trait] = predict_res_tmp['SCORESUM'].values
        if args.pheno is not None:
            predict_res_tmp = predict_res_tmp.merge(pheno[['IID',args.col_name]],left_on='IID',right_on='IID',how='inner')
            corr = pearsonr(predict_res_tmp['SCORESUM'],predict_res_tmp[args.col_name])[0]
            print('R2 for {}: {}'.format(trait,corr**2))
        ftmp.close()
    if args.save.endswith('.csv'):
        predict_res.to_csv(args.save,sep='\t',index=None)
    else:
        predict_res.to_csv(args.save+'.csv',sep='\t',index=None)
