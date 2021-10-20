import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin
import logging
import copy
from scipy.sparse.linalg import cg
from scipy import linalg
import argparse
import sys
import os
from utils import * 

__version__ = '1.0.1'
SOFTWARE_CORRESPONDENCE_EMAIL1 = 'jxiaoae@connect.ust.hk'
SOFTWARE_CORRESPONDENCE_EMAIL2 = 'mcaiad@connect.ust.hk' 

HEADER = """
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
<> 
<> XPXP: Improving polygenic prediction by cross-population and cross-phenotype analysis
<> Version: %s
<> MIT License
<>
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
<> Software-related correspondence: %s or %s
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>       
""" % (__version__, SOFTWARE_CORRESPONDENCE_EMAIL1, SOFTWARE_CORRESPONDENCE_EMAIL2)
  
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='posterior mean of efffect size for multiple GWAS summary')
    parser.add_argument('--save', type=str, help='output path',required=True)
    parser.add_argument('--gc_file', type=str, help='genetic covariance file path',required=True)
    parser.add_argument('--ec_file', type=str, help='environment covariance file path',required=False)
    parser.add_argument('--sumst_files', type=str, help='summary statisitc files, separated by comma, ordered by population: Target+Auxiliary',required=True)
    parser.add_argument('--sumst_names', type=str, help='summary statisitc names, separated by comma, the order corresponds to the summary statisitc files, different populations are separated by "+"',required=True)
    parser.add_argument('--use_snps', type=str, help='use assigned SNPs, skip SNPs matching step',required=False)
    parser.add_argument('--ref_files', type=str, help='LD reference files path, plink1 file version, seperated by comma.',required=False)
    parser.add_argument('--ld_block_file', type=str, help='ld block file path',required=False)
    parser.add_argument('--return_LDpredinf', help='return LDpredinf for each traits', action="store_true")
    parser.add_argument('--num_threads', type=str, help='number of threads', default="22")
    parser.add_argument('--fix_effect_traits', type=str, help='traits to incorporate fix large genetic effect, seperated by comma',required=False)
    parser.add_argument('--pvalue', type=float, help='p-value threshold for fix effect', default=1e-6)
    parser.add_argument('--alpha', type=float, help='selection-related parameter: [2f(1 −f)]^α', default=-0.25)


    args = parser.parse_args()

    os.environ["OMP_NUM_THREADS"] = args.num_threads
    os.environ["OPENBLAS_NUM_THREADS"] = args.num_threads
    os.environ["MKL_NUM_THREADS"] = args.num_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = args.num_threads
    os.environ["NUMEXPR_NUM_THREADS"] = args.num_threads

    if not args.save.endswith('.csv'):
        args.save = args.save+'.csv'

    logger = configure_logging('.'.join(args.save.split('.')[:-1]))

    logger.info(HEADER)
    logger.info("See full log at: %s\n", os.path.abspath('.'.join(args.save.split('.')[:-1])+'.log'))
    logger.info("\nProgram executed via:\n%s\n", ' '.join(sys.argv).replace("--", " \\ \n\t--"))

    res_save_file = args.save
    corr_file = args.gc_file
    if len(args.ref_files.split(',')) == 2:
        n_pop = 2
        file_ref1, file_ref2 = args.ref_files.split(',')
    elif len(args.ref_files.split(',')) == 3:
        n_pop = 3
        file_ref1, file_ref2, file_ref3 = args.ref_files.split(',')
    elif len(args.ref_files.split(',')) == 1:
        n_pop = 1
        file_ref1 = args.ref_files
    else:
        logger.error('detect more than three LD reference files')
        sys.exit()

    logger.info("Find {} populations".format(n_pop))

    ld_block_file = args.ld_block_file
    sumst_files = args.sumst_files.split(',')

    if n_pop != len(args.sumst_names.split('+')):
        logger.error('the number of populations in LD reference files dose not match with number of populations of summary statistics files')
        sys.exit()

    sumst_names = args.sumst_names.replace('+',',').split(',')
    sumst_ref1 = args.sumst_names.split('+')[0].split(',')
    if n_pop>1:
        sumst_ref2 = args.sumst_names.split('+')[1].split(',')
    else:
        sumst_ref2 = []

    logger.info("Plink reference 1 file path: {} ".format(file_ref1))
    if n_pop>1:
        logger.info("Plink reference 2 file path: {} ".format(file_ref2))
    if n_pop>2:
        logger.info("Plink reference 3 file path: {} ".format(file_ref3))

    logger.info("LD block file path: {} ".format(ld_block_file))
    if args.return_LDpredinf:
        logger.info("return LDpredinf for each GWAS summary")
    logger.info("Genetic covariance file path: {} ".format(corr_file))
    if args.ec_file is not None:
        logger.info("Environmental covariance file path: {} ".format(args.ec_file))

    logger.info("output path: {} ".format(args.save))

    for sumsf in sumst_files:
        logger.info("GWAS summary statistic file path: {} ".format(sumsf))

    if len(sumst_names) != len(sumst_files):
        logger.error("number of summary statistic files doesn't equals to length of summary statistic name list")
        sys.exit()
    
    try:
        height_corr = pd.read_csv(corr_file,index_col=0).loc[sumst_names,sumst_names]
    except:
        logger.error("summary statistic name dosen't match index/column name of genetic covariance file")
        sys.exit()

    for i,sums_n in enumerate(sumst_names):
        if sums_n in sumst_ref1:
            logger.info("GWAS summary statistic from population 1, name {}: {} ".format(i+1,sums_n))
        elif sums_n in sumst_ref2:
            logger.info("GWAS summary statistic from population 2, name {}: {} ".format(i+1,sums_n))
        else:
            logger.info("GWAS summary statistic from population 3, name {}: {} ".format(i+1,sums_n))

    sumstats_dict = {}
    for i,sumsf in enumerate(sumst_files):
        try:
            tmp_df = pd.read_csv(sumsf,sep='\t',compression='infer')  
        except:
            tmp_df = pd.read_csv(sumsf,sep='\t')    
        sumstats_dict[sumst_names[i]] = tmp_df.dropna().drop_duplicates('SNP') 
        logger.info("Sumstats {} shape: {}".format(sumst_names[i], sumstats_dict[sumst_names[i]].shape)) 

    if args.fix_effect_traits is None:
        le_traits = []
        logger.info('No fix effect applied!')
    else:
        le_traits = args.fix_effect_traits.split(',')
        for le_trait in le_traits:
            if le_trait not in sumst_names:
                logger.error("fix effect trait name dosen't match with summary statistic names: {}".format(','.join(sumst_names)))
                sys.exit()
            else:
                if le_trait in sumst_ref1:
                    logger.info('Estimate fix effect from GWAS summary: {} in population 1'.format(le_trait))
                elif le_trait in sumst_ref2:
                    logger.info('Estimate fix effect from GWAS summary: {} in population 2'.format(le_trait))
                else:
                    logger.info('Estimate fix effect from GWAS summary: {} in population 3'.format(le_trait))

    ref1_info = pd.read_csv(file_ref1+'.bim',sep='\t',header=None)
    ref1_info.columns = ['chr','SNP','cm','bp','A1','A2']
    logger.info("Ref1 bim shape: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    if n_pop>1:
        ref2_info = pd.read_csv(file_ref2+'.bim',sep='\t',header=None)
        ref2_info.columns = ['chr','SNP','cm','bp','A1','A2']
        logger.info("Ref2 bim shape: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1]))     
    else:
        ref2_info = ref1_info.copy()
    if n_pop>2:
        ref3_info = pd.read_csv(file_ref3+'.bim',sep='\t',header=None)
        ref3_info.columns = ['chr','SNP','cm','bp','A1','A2'] 
        logger.info("Ref3 bim shape: {}x{} ".format(ref3_info.shape[0],ref3_info.shape[1]))
    else:
        ref3_info = ref1_info.copy()

    # matching SNPs
    snps_lst = [set(ref1_info['SNP'].values),set(ref2_info['SNP'].values),set(ref3_info['SNP'].values)]+\
            [set(sumstats_dict[_]['SNP'].values) for _ in sumst_names]
    if args.use_snps is None:
        logger.info('Matching SNPs from LD reference and {} GWAS summary...'.format(','.join(sumst_names)))
        snps_common = list(snps_lst[0].intersection(*snps_lst))
    else:
        logger.info('load user supplied SNPs from {}'.format(args.use_snps))
        snps_lst += [set(np.loadtxt(args.use_snps,dtype=str))]
        snps_common = list(snps_lst[0].intersection(*snps_lst))
        
    ref1_info = ref1_info.loc[ref1_info['SNP'].isin(snps_common)]
    ref2_info = ref2_info.loc[ref2_info['SNP'].isin(snps_common)]
    ref3_info = ref3_info.loc[ref3_info['SNP'].isin(snps_common)]
    logger.info("Ref1 bim shape after matching: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    if n_pop>1:
        logger.info("Ref2 bim shape after matching: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1])) 
    if n_pop>2:
        logger.info("Ref3 bim shape after matching: {}x{} ".format(ref3_info.shape[0],ref3_info.shape[1])) 

    for i,sums_n in enumerate(sumst_names):
        sumstats_dict[sums_n] = sumstats_dict[sums_n].loc[sumstats_dict[sums_n]['SNP'].isin(snps_common)]

    sumstats_df_columns = ['SNP','N','Z','A1','A2']
    for sums_n in sumst_names:
        sumstats_dict[sums_n] = ref1_info[['SNP','A1','A2']].merge(sumstats_dict[sums_n],left_on='SNP',right_on='SNP',how='left')
        sumstats_dict[sums_n].loc[sumstats_dict[sums_n]['A1_y'].isnull(),'A1'] = sumstats_dict[sums_n].loc[sumstats_dict[sums_n]['A1_y'].isnull(),'A1_x']
        sumstats_dict[sums_n].loc[~sumstats_dict[sums_n]['A1_y'].isnull(),'A1'] = sumstats_dict[sums_n].loc[~sumstats_dict[sums_n]['A1_y'].isnull(),'A1_y']
        sumstats_dict[sums_n].loc[sumstats_dict[sums_n]['A2_y'].isnull(),'A2'] = sumstats_dict[sums_n].loc[sumstats_dict[sums_n]['A2_y'].isnull(),'A2_x']
        sumstats_dict[sums_n].loc[~sumstats_dict[sums_n]['A2_y'].isnull(),'A2'] = sumstats_dict[sums_n].loc[~sumstats_dict[sums_n]['A2_y'].isnull(),'A2_y']
        sumstats_dict[sums_n]['N']= sumstats_dict[sums_n]['N'].median()
        sumstats_dict[sums_n] = sumstats_dict[sums_n][sumstats_df_columns]
    
    for i,sums_n in enumerate(sumst_names):
        logger.info("Sumstats {} shape after matching: {}x{}".format(sums_n, sumstats_dict[sums_n].shape[0],\
            sumstats_dict[sums_n].shape[1]))
    

    logger.info('Removing ambiguous SNPs...')
    # replace T with A, replace G with C; A=1, C=2
    ref1_info['A1_int'] = ref1_info['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref1_info['A2_int'] = ref1_info['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref2_info['A1_int'] = ref2_info['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref2_info['A2_int'] = ref2_info['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref3_info['A1_int'] = ref3_info['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
    ref3_info['A2_int'] = ref3_info['A2'].map(lambda x: 1 if x in ['A','T'] else 2)

    for i,sums_n in enumerate(sumst_names):
        sumstats_dict[sums_n]['A1_int'] = sumstats_dict[sums_n]['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
        sumstats_dict[sums_n]['A2_int'] = sumstats_dict[sums_n]['A2'].map(lambda x: 1 if x in ['A','T'] else 2)

    snps_rm_idx = (np.array([(ref1_info['A1_int'].values+ref1_info['A2_int'].values != ref2_info['A1_int'].values+ref2_info['A2_int'].values)]+\
        [(ref1_info['A1_int'].values+ref1_info['A2_int'].values != ref3_info['A1_int'].values+ref3_info['A2_int'].values)] +\
            [(ref2_info['A1_int'].values+ref2_info['A2_int'].values != ref3_info['A1_int'].values+ref3_info['A2_int'].values)] +\
        [(ref1_info['A1_int'].values+ref1_info['A2_int'].values != sumstats_dict[sums_n]['A1_int'].values+\
          sumstats_dict[sums_n]['A2_int'].values) for sums_n in sumst_names])).sum(axis=0)>0

    ref1_info = ref1_info.reset_index()
    ref1_info = ref1_info.loc[~snps_rm_idx,:]
    ref1_info = ref1_info.reset_index()
    ref2_info = ref2_info.reset_index()
    ref2_info = ref2_info.loc[~snps_rm_idx,:]
    ref2_info = ref2_info.reset_index()
    ref3_info = ref3_info.reset_index()
    ref3_info = ref3_info.loc[~snps_rm_idx,:]
    ref3_info = ref3_info.reset_index()
    logger.info("Ref1 bim shape after removing ambiguous SNPs: {}x{} ".format(ref1_info.shape[0],ref1_info.shape[1])) 
    if n_pop>1:
        logger.info("Ref2 bim shape after removing ambiguous SNPs: {}x{} ".format(ref2_info.shape[0],ref2_info.shape[1])) 
    if n_pop>2:
        logger.info("Ref3 bim shape after removing ambiguous SNPs: {}x{} ".format(ref3_info.shape[0],ref3_info.shape[1])) 

    for i,sums_n in enumerate(sumst_names):
        sumstats_dict[sums_n] = sumstats_dict[sums_n].reset_index()
        sumstats_dict[sums_n] = sumstats_dict[sums_n].loc[~snps_rm_idx,:]
        sumstats_dict[sums_n] = sumstats_dict[sums_n].reset_index()
        logger.info("Sumstats {} shape after removing ambiguous SNPs: {}x{}".format(sums_n, sumstats_dict[sums_n].shape[0],\
            sumstats_dict[sums_n].shape[1]))
    
    logger.info('Loading first reference genotype...')
    X1 = ReadPlink(file_ref1,ref1_info,np.int8)
    logger.info("Reference 1 genotype shape: {}x{} ".format(X1.shape[0],X1.shape[1]))
    if n_pop>1:
        logger.info('Loading second reference genotype...')
        X2 = ReadPlink(file_ref2,ref2_info,np.int8)
        logger.info("Reference 2 genotype shape: {}x{} ".format(X2.shape[0],X2.shape[1]))
        ind_ref = (ref1_info['A1']!=ref2_info['A1']).values
        logger.info("{} SNPs need to be flipped in ref2 to match ref1".format(ind_ref.sum()))
        if ind_ref.sum()>0:
            X2[:,ind_ref] = 2-X2[:,ind_ref]
    else:
        X2 = X1
    if n_pop>2:
        logger.info('Loading third reference genotype...')
        X3 = ReadPlink(file_ref3,ref3_info,np.int8)
        logger.info("Reference 3 genotype shape: {}x{} ".format(X3.shape[0],X3.shape[1]))
        ind_ref = (ref1_info['A1']!=ref3_info['A1']).values
        logger.info("{} SNPs need to be flipped in ref3 to match ref1".format(ind_ref.sum()))
        if ind_ref.sum()>0:
            X3[:,ind_ref] = 2-X3[:,ind_ref]
    else:
        X3 = X1
    
    
    #X1, X1_maf, X1_sd = ScaleGenotype(X1)
    #X2, X2_maf, X2_sd = ScaleGenotype(X2)
    X1_mean, X1_sd = GenotypeMoment(X1)
    X2_mean, X2_sd = GenotypeMoment(X2)
    X3_mean, X3_sd = GenotypeMoment(X3)

    logger.info('Based on the minor allel of ref1 panel, correct for the z score...')
    # based on the minor allel of reference panel, correct for the z score in summary statistics
    for i,sums_n in enumerate(sumst_names):
        sums_inverse = (ref1_info['A1'].values != sumstats_dict[sums_n]['A1'].values)
        sumstats_dict[sums_n].loc[sums_inverse,'Z'] = -sumstats_dict[sums_n].loc[sums_inverse,'Z'] 
        logger.info("SNPs need aligned for {}: {}".format(sums_n,sums_inverse.sum()))

    p = ref1_info.shape[0]
    logger.info("GWAS summary SNP size: {}".format(p))
    ns = {}
    for i,sums_n in enumerate(sumst_names):
        ns[sums_n] = sumstats_dict[sums_n]['N'].median()
        logger.info("{} GWAS summary sample size: {}".format(sums_n, ns[sums_n]))
    
    Sigma_beta = height_corr.values
    if args.ec_file is None:
        Sigma_e = np.diag(1-np.diag(height_corr))
        logger.info("Assume all GWASs are mutual independent...")
    else:
        logger.info("Load environment covariance from {}".format(args.ec_file))
        try:
            Sigma_e = pd.read_csv(args.ec_file,index_col=0).loc[sumst_names,sumst_names].values
        except:
            logger.error("summary statistic name dosen't match index/column name of environment covariance file")
            sys.exit()

    block = pd.read_csv(ld_block_file,delim_whitespace=True,dtype={'start':int,'stop':int})
    block['chr'] = block['chr'].map(lambda x:int(x.replace('chr','')))
    ngroup = block.shape[0]
    logger.info("{} blocks loaded from {}".format(ngroup, ld_block_file))
    ref1_info_block = np.zeros(ref1_info.shape[0],dtype=int)
    for b_ in block.iterrows():
        ref1_info_block[ref1_info.loc[(ref1_info['chr']==b_[1]['chr'])&(ref1_info['bp']>=b_[1]['start'])&(ref1_info['bp']<=b_[1]['stop'])].index.values]= b_[0]
    ref1_info['block'] = ref1_info_block

    # fix effects

    if args.fix_effect_traits is None:
        pass
    else:
        snps_clump = {}
        idx_snps_clump = {}
        beta_snps_clump = {}
        for le_trait in le_traits:
            if le_trait in sumst_ref1:
                snps_clump[le_trait] = ld_clump(sumstats_dict[le_trait],file_ref1,args.pvalue)
                idx_snps_clump[le_trait] = ref1_info.loc[ref1_info['SNP'].isin(snps_clump[le_trait])].index.values
                beta_snps_clump[le_trait] = np.zeros(idx_snps_clump[le_trait].shape[0])
            elif le_trait in sumst_ref2:
                snps_clump[le_trait] = ld_clump(sumstats_dict[le_trait],file_ref2,args.pvalue)
                idx_snps_clump[le_trait] = ref2_info.loc[ref2_info['SNP'].isin(snps_clump[le_trait])].index.values
                beta_snps_clump[le_trait] = np.zeros(idx_snps_clump[le_trait].shape[0])
            else:
                snps_clump[le_trait] = ld_clump(sumstats_dict[le_trait],file_ref3,args.pvalue)
                idx_snps_clump[le_trait] = ref3_info.loc[ref3_info['SNP'].isin(snps_clump[le_trait])].index.values
                beta_snps_clump[le_trait] = np.zeros(idx_snps_clump[le_trait].shape[0])

    beta_estimate_mu = np.zeros((p,len(sumst_names)))
    beta_estimate_muxpxp = np.zeros((p,len(sumst_names)))
    logger.info('Computing posterior mean...')
    for i in range(ngroup):
        if i % 10 == 0:
            logger.info("Fininshing {} block".format(i))
        idx_i = ref1_info.loc[ref1_info['block']==i,].index.values
        zs_b = {}
        for i_,sums_n in enumerate(sumst_names):
            zs_b[sums_n] = np.sqrt(ns[sums_n]/p)*sumstats_dict[sums_n]['Z'].values[idx_i]
        X1_b = ((X1[:,idx_i]-X1_mean[idx_i])/X1_sd[idx_i])/np.sqrt(X1.shape[1])
        L1 = X1_b.T.dot(X1_b)/X1_b.shape[0]
        if n_pop>1:
            X2_b = ((X2[:,idx_i]-X2_mean[idx_i])/X2_sd[idx_i])/np.sqrt(X2.shape[1])
            L2 = X2_b.T.dot(X2_b)/X2_b.shape[0]
        else:
            X2_b = X1_b
            L2 = L1
        if n_pop>2:
            X3_b = ((X3[:,idx_i]-X3_mean[idx_i])/X3_sd[idx_i])/np.sqrt(X3.shape[1])
            L3 = X3_b.T.dot(X3_b)/X3_b.shape[0]
        else:
            X3_b = X1_b
            L3 = L1
        for le_trait in le_traits:
            isle = np.intersect1d(idx_i,idx_snps_clump[le_trait])
            if isle.shape[0]==0:
                pass
            else:
                logger.info("Finding {} snps with large effect in {} block for trait {}".format(isle.shape[0],i,le_trait))
                idx_le = [idx_snps_clump[le_trait].tolist().index(_) for _ in isle]
                if le_trait in sumst_ref1:
                    Xl1_b = ((X1[:,isle]-X1_mean[isle])/X1_sd[isle])/np.sqrt(X1.shape[1])
                    L_sl1 = X1_b.T@Xl1_b*ns[le_trait]/X1_b.shape[0]
                    L_ss = L1*ns[le_trait]
                elif le_trait in sumst_ref2:
                    Xl1_b = ((X2[:,isle]-X2_mean[isle])/X2_sd[isle])/np.sqrt(X2.shape[1])
                    L_sl1 = X2_b.T@Xl1_b*ns[le_trait]/X2_b.shape[0]
                    L_ss = L2*ns[le_trait]
                else:
                    Xl1_b = ((X3[:,isle]-X3_mean[isle])/X3_sd[isle])/np.sqrt(X3.shape[1])
                    L_sl1 = X3_b.T@Xl1_b*ns[le_trait]/X3_b.shape[0]
                    L_ss = L3*ns[le_trait]
                L_l1 = Xl1_b.T.dot(Xl1_b)*ns[le_trait]/Xl1_b.shape[0]
                zl = (np.sqrt(ns[le_trait]/p)*(sumstats_dict[le_trait]['Z'].values[isle])).reshape(-1,)
                zs = zs_b[le_trait].reshape(-1,)
                beta_le_b = compute_fe(L_l1,L_sl1,L_ss,zs,zl,height_corr.loc[le_trait,le_trait])
                #beta_le_b = linalg.inv(L_l1) @ (np.sqrt(ns[le_trait]/p)*(sumstats_dict[le_trait]['Z'].values[isle])).reshape(-1,)
                beta_snps_clump[le_trait][idx_le] = beta_le_b
                zs_b[le_trait] -= L_sl1@beta_le_b
        
        beta_estimate_mu[idx_i,:], beta_estimate_muxpxp[idx_i,:] = ComputePMXPXP(i, sumst_names, \
                sumst_ref1, sumst_ref2, ns, p, L1, L2, L3, Sigma_beta, Sigma_e, zs_b, logger, return_LDpredinf=args.return_LDpredinf)

    beta_estimate_mu[:,:len(sumst_ref1)] = (beta_estimate_mu[:,:len(sumst_ref1)]*X1_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)
    beta_estimate_mu[:,len(sumst_ref1):len(sumst_ref1+sumst_ref2)] = (beta_estimate_mu[:,len(sumst_ref1):len(sumst_ref1+sumst_ref2)]*X2_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)
    beta_estimate_mu[:,len(sumst_ref1+sumst_ref2):] = (beta_estimate_mu[:,len(sumst_ref1+sumst_ref2):]*X3_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)

    beta_estimate_muxpxp[:,:len(sumst_ref1)] = (beta_estimate_muxpxp[:,:len(sumst_ref1)]*X1_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)
    beta_estimate_muxpxp[:,len(sumst_ref1):len(sumst_ref1+sumst_ref2)] = (beta_estimate_muxpxp[:,len(sumst_ref1):len(sumst_ref1+sumst_ref2)]*X2_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)
    beta_estimate_muxpxp[:,len(sumst_ref1+sumst_ref2):] = (beta_estimate_muxpxp[:,len(sumst_ref1+sumst_ref2):]*X3_sd.reshape(-1,1)**args.alpha)/np.sqrt(p)
    
    if args.fix_effect_traits is not None:
        beta_estimate_fixeffect = ref1_info[['chr','SNP','bp','A1','A2']]
        for le_trait in le_traits:
            if le_trait in sumst_ref1:
                beta_snps_clump[le_trait] = beta_snps_clump[le_trait]*X1_sd[idx_snps_clump[le_trait]]**args.alpha/np.sqrt(p)
            elif le_trait in sumst_ref2:
                beta_snps_clump[le_trait] = beta_snps_clump[le_trait]*X2_sd[idx_snps_clump[le_trait]]**args.alpha/np.sqrt(p)
            else:
                beta_snps_clump[le_trait] = beta_snps_clump[le_trait]*X3_sd[idx_snps_clump[le_trait]]**args.alpha/np.sqrt(p)
            beta_estimate_fixeffect.loc[idx_snps_clump[le_trait],le_trait] = beta_snps_clump[le_trait]
            le_trait_idx = sumst_names.index(le_trait)
            beta_estimate_muxpxp[idx_snps_clump[le_trait],le_trait_idx] = beta_estimate_muxpxp[idx_snps_clump[le_trait],le_trait_idx]+beta_snps_clump[le_trait]
        beta_estimate_fixeffect = beta_estimate_fixeffect.dropna(subset=le_traits, how='all')
        beta_estimate_fixeffect = beta_estimate_fixeffect.fillna(0)
        beta_estimate_fixeffect.to_csv(res_save_file.replace('.csv','')+'-large-effect.csv',index=None,sep='\t')

    bim = ref1_info[['chr','SNP','bp','A1','A2']]
    for i,sums_n in enumerate(sumst_names):
        if args.return_LDpredinf:
            bim.loc[:,sums_n+'-mu'] = beta_estimate_mu[:,i]
        bim.loc[:,sums_n+'-muxpxp'] = beta_estimate_muxpxp[:,i]
    bim.to_csv(res_save_file,index=None,sep='\t')
