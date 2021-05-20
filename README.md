# XPXP
The XPXP package for improving polygenic risk score (PRS) prediction by cross-population and cross-phenotype analysis

# Installation
```bash
git clone https://github.com/YangLabHKUST/XPXP.git
cd XPXP
conda env create -f environment.yml
conda activate xpxp
```
check the installation status
```bash
python ./src/XPXP.py -h
python ./src/ParaEstimate.py -h
```

# Quick start

We illustrate the usage of XPXP using the GWAS summary statistics of height from BBJ, UKBB European and GIANT cohorts. For demonstration, we use the easily accessible 1000 Genomes project genotypes as reference panels. The datasets involved in the following example is availabel from [here](https://www.dropbox.com/s/b7g1438pl78c9q9/XPXP_demo.tar.gz?dl=0)

## Data preparation

Input files of XPXP include:

- summay statistics files of the target and auxiliary population
- summay statistics names of the target and auxiliary population, population is separated by '+'
- reference panel of the target and auxiliary population in plink 1 format
- genetic covariance file
- environmental covariance file if sample overlap exist
- LD blocks file
- summay statistics names for incorporating population-specific of phenotype specific large genetic effects

The XPXP format GWAS summary statistics file has at least 6 fields:

- SNP: SNP rsid
- N: sample size
- Z: Z-scores
- A1: effect allele
- A2: other allele. 
- P: p-value 


e.g. 
```bash
$ head ./XPXP_demo/height_BBJ_summary_hm3.txt
SNP     N       Z       A1      A2      P
rs1048488       159095  0.659462022101568       T       C       0.509599125893493
rs3115850       159095  0.658089413463925       C       T       0.510480678080961
rs4040617       159095  -0.37360256722158497    G       A       0.708700023669772
```


The genetic covariance file with trait names as index and header:
```bash
$ cat ./XPXP_demo/gcov_height_BBJ.csv
,height_BBJ-EAS,height_Wood-EUR,height_UKB-EUR
height_BBJ-EAS,0.3894,0.1866552482397845,0.20741027144622678
height_Wood-EUR,0.1866552482397845,0.2596,0.2525
height_UKB-EUR,0.20741027144622678,0.2525,0.369
```


## Run XPXP to compute SNPs effect size

Once the input files are formatted, XPXP will automatically process the datasets, including SNPs overlapping and allele matching.
Run XPXP with the following comand (**delete "#" and comments when run it**):
```bash
$ python [INSTALL PATH]/XPXP/src/XPXP.py \
--num_threads 40 \
--save ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv \ # output file path 
--gc_file ./XPXP_demo/gcov_height_BBJ.csv \ # genetic covariance file
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_GIANT_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \ # summary statistics files, Target+Auxiliary 
--sumst_names height_BBJ-EAS+height_Wood-EUR,height_UKB-EUR \ # summary statistics names, the order corresponds to the summary statistics files, population is seprated by "+", e.g., Target+Auxiliary
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \ # reference panels, Target+Auxiliary
--fix_effect_traits height_BBJ-EAS # traits to incorporate fixed large genetic effect
```

XPXP returns the estimated SNPs effect size in ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv
```bash
$ head ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv
chr     SNP     bp      A1      A2      height_BBJ-EAS-muxpxp   height_Wood-EUR-muxpxp  height_UKB-EUR-muxpxp
1       rs3934834       1005806 T       C       0.001490031305500368    0.0008587749623481179   0.001348694388142856
1       rs3766191       1017587 T       C       0.001294311742987717    0.0007056128453579005   0.0011593703907873875
1       rs9442372       1018704 A       G       0.001354169642801807    0.0007882917460030307   0.0009947411386104607
```

where A1 is the effect allele, A2 is the other allele. \<TraitName-muxpxp\> is the estimated SNPs effect size of \<TraitName\> computed by XPXP. If argument "--return_LDpredinf" is given, then XPXP will also output the estimated SNPs effect size computed by LDpred-inf (\<TraitName-mu\>) using the GWAS summary statistic of \<TraitName\> only, see the exmaple below. 

```bash
$ python [INSTALL PATH]/XPXP/src/XPXP.py \
--num_threads 40 \
--save ./XPXP_demo/PM_height_BBJ-GIANT-UKB_AddLDpredinf.csv \
--gc_file ./XPXP_demo/gcov_height_BBJ.csv \
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_GIANT_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \
--sumst_names height_BBJ-EAS+height_Wood-EUR,height_UKB-EUR \
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \
--fix_effect_traits height_BBJ-EAS \
--return_LDpredinf

$ head ./XPXP_demo/PM_height_BBJ-GIANT-UKB_AddLDpredinf.csv
chr     SNP     bp      A1      A2      height_BBJ-EAS-mu       height_BBJ-EAS-muxpxp   height_Wood-EUR-mu      height_Wood-EUR-muxpxp  height_UKB-EUR-mu       height_UKB-EUR-muxpxp
1       rs3934834       1005806 T       C       0.0005643129160874588   0.001490031305500368    -6.687588995674362e-05  0.0008587749623481179   0.001128698464435024    0.001348694388142856
1       rs3766191       1017587 T       C       0.0005275501273343171   0.001294311742987717    -9.677694792123072e-05  0.0007056128453579005   0.0010440824756778855   0.0011593703907873875
1       rs9442372       1018704 A       G       0.0005044304401368068   0.001354169642801807    0.00018283924802650126  0.0007882917460030307   0.0006365824088935744   0.0009947411386104607
```


## Evaluate the prediction performance using individual-level GWAS data

Input files:

- ```--geno``` genotype file of testing data (UKBB Chinese, n=1,439), plink1 version
- ```--beta``` the estimated SNPs effect size returned by XPXP

```bash
$ python [INSTALL PATH]/XPXP/src/Predict.py \
--save ./XPXP_demo_NotAvailable/predict_ukb_chn.csv \
--geno ./XPXP_demo_NotAvailable/ukb_chn_qc1 \
--beta ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv 
```
The predicted risk scores are returned in ./XPXP_demo_NotAvailable/predict_ukb_chn.csv
```bash
$ head ./XPXP_demo_NotAvailable/predict_ukb_chn.csv
FID     IID     height_BBJ-EAS-muxpxp   height_Wood-EUR-muxpxp  height_UKB-EUR-muxpxp
1002529 1002529 -0.472046       -1.1332799999999998     -0.8280420000000001
1006646 1006646 -0.57514        -0.7420140000000001     -0.48337600000000003
1006809 1006809 -0.772656       -0.966025       -0.8342870000000001
```
where column \<TraitName-muxpxp\> is the predicted PRS using the estimated SNPs effect size (\<TraitName-muxpxp\>) returned by XPXP

However, the individual-level GWAS data of UKBB is not availabel due to the data sharing restriction. We therefore show how to use the height GWAS of UKBB Chinese as validation dataset

## Evaluate the prediction performance using GWAS summary statistics

We follow [XPASS](https://github.com/YangLabHKUST/XPASS) to use the following equation:

<img src="https://latex.codecogs.com/svg.image?R^2=corr(y,\hat{y})^2=\left(\frac{cov(y,\hat{y})}{\sqrt{var(y)var(\hat{y})}}\right)^2=\left(\frac{z^T\tilde{\mu}/\sqrt{n}}{\sqrt{\tilde{\mu}^T\Sigma\tilde{\mu}}}\right)^2," title="R^2=corr(y,\hat{y})^2=\left(\frac{cov(y,\hat{y})}{\sqrt{var(y)var(\hat{y})}}\right)^2=\left(\frac{z^T\tilde{\mu}/\sqrt{n}}{\sqrt{\tilde{\mu}^T\Sigma\tilde{\mu}}}\right)^2," />

where z is the z-score of external summsry statistics, n is its sample size, <img src="https://latex.codecogs.com/svg.image?\tilde{\mu}" title="\tilde{\mu}" /> is the posterior mean of effect size at the standardized genotype scale (mean 0 and variance 1), <img src="https://latex.codecogs.com/svg.image?\Sigma" title="\Sigma" /> is the SNPs correlation matrix computed from a reference panel.

Input files:

- ```--sumst_file``` GWAS summary statistics of UKBB Chinese height
- ```--beta``` the estimated SNPs effect size returned by XPXP
- ```--col_name``` specify the column name of SNPs effect size file

```bash
$ python [INSTALL PATH]/XPXP/src/PredictSS.py \
--ref_file ./XPXP_demo/1000G.EAS.QC.hm3.ind \
--sumst_file ./XPXP_demo/UKB_Chinese_height_GWAS_summary.txt \
--beta ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv  \
--col_name height_BBJ-EAS-muxpxp

Output: R2 for height_BBJ-EAS-muxpxp: 0.13039907941388484
```

Compared to XPXP trained on the BBJ and UKBB datasets only:
```bash
$ python [INSTALL PATH]/XPXP/src/XPXP.py \
--num_threads 40 \
--save ./XPXP_demo/PM_height_BBJ-UKB.csv \
--gc_file ./XPXP_demo/gcov_height_BBJ.csv \
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \
--sumst_names height_BBJ-EAS+height_UKB-EUR \
--use_snps ./XPXP_demo/snps_overlap.txt \ # use the same set of SNPs for fairness
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \
--fix_effect_traits height_BBJ-EAS

$ python [INSTALL PATH]/XPXP/src/PredictSS.py \
--ref_file ./XPXP_demo/1000G.EAS.QC.hm3.ind \
--sumst_file ./XPXP_demo/UKB_Chinese_height_GWAS_summary.txt \
--beta ./XPXP_demo/PM_height_BBJ-UKB.csv  \
--col_name height_BBJ-EAS-muxpxp

Output: R2 for height_BBJ-EAS-muxpxp: 0.1256874019449264
```
The predicted R2 declined a little bit due to the removing of GIANT training data

## Generate genetic and environmental covariance matrix
XPXP requires genetic and environmental covariance matrix estimates for computing the posterior mean of SNPs effect size. For parameters within a population (e.g., SNP-heritability, genetic covariance, and environmental covariance for pair of traits with substiantial sample overlap), we apply LD score regression ([LDSC](https://github.com/bulik/ldsc)) to obtain their estimates. For parameters cross two populations (e.g., trans-ancestry genetic covariance), we follow [XPASS](https://github.com/YangLabHKUST/XPASS) to estimate the trans-ancestry genetic covariance using fully LD matrix rather than the LD information from local SNPs utilized in LDSC.

Here we provide a helper script (```ParaEstimate.py```, a wrapper of LDSC and XPASS) to conveniently obtain the input parameters for  ```XPXP.py```. 
First of all, we need install the LDSC v1.0.1 using conda:
```bash
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda env create -f environment.yml
# no need to run 'conda activate ldsc'
```
please note that we <ins>**do not need to activate the ldsc environment**</ins>

then we run ```ParaEstimate.py``` in ```xpxp``` environment as following:
```bash
$ python [INSTALL PATH]/XPXP/src/ParaEstimate.py \
--save_dir ./XPXP_demo/params \
--ldsc_path [LDSC PATH] \
--ldsc_files ./XPXP_demo/eas_ldscores/,./XPXP_demo/eur_ldscores/ \
--merge_alleles ./XPXP_demo/w_hm3.snplist \
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_GIANT_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \
--sumst_names height_BBJ-EAS+height_Wood-EUR,height_UKB-EUR \
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \
--covar_files ./XPXP_demo/1000G.EAS.QC.hm3.ind.pc5.txt,./XPXP_demo/1000G.EUR.QC.hm3.ind.pc20.txt
```
Inputs:
- ```--save_dir``` output dir path
- ```--ldsc_path``` LDSC install path
- ```--ldsc_files``` LDscore files
- ```--merge_alleles``` file used for matching alleles
- ```--sumst_files``` summary statisitc files, separated by comma
- ```--sumst_names``` summary statisitc names, separated by comma, the order is corresponds to the --sumst_files, different populations are separated by "+", .e.g. Target+Auxiliary'
- ```--ref_files``` LD reference files path, plink1 file version, seperated by comma
- ```--covar_files``` LD reference covariate files path, seperated by comma
- ```--ld_block_file``` LD block file

Outputs:
- Genetic covariance file
```bash
$ cat ./XPXP_demo/params/gcov.csv
,height_BBJ-EAS,height_Wood-EUR,height_UKB-EUR
height_BBJ-EAS,0.3894,0.1940299242269646,0.21439033369134217
height_Wood-EUR,0.1940299242269646,0.2711,0.2606
height_UKB-EUR,0.21439033369134217,0.2606,0.369
```
- Genetic correlation file
```bash
$ cat ./XPXP_demo/params/gcorr.csv
,height_BBJ-EAS,height_Wood-EUR,height_UKB-EUR
height_BBJ-EAS,0.0,0.597181,0.56558
height_Wood-EUR,0.597181,0.0,0.8239416739808213
height_UKB-EUR,0.56558,0.8239416739808213,0.0
```
- Environmental covariance file
```bash
$ cat ./XPXP_demo/params/ecov.csv
,height_BBJ-EAS,height_Wood-EUR,height_UKB-EUR
height_BBJ-EAS,0.6106,0.0,0.0
height_Wood-EUR,0.0,0.7289,0.08950000000000002
height_UKB-EUR,0.0,0.08950000000000002,0.631
```


# Development

The XPXP package is developed by Jiashun Xiao (jxiaoae@connect.ust.hk).

# Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.












