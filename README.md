# XPXP
The XPXP package for improving polygenic risk score (PRS) prediction by cross-population and cross-phenotype analysis

# Quick start

We illustrate the usage of XPXP using the GWAS summary statistics of height from BBJ, UKBB European and GIANT cohorts. For demonstration, we use the easily accessible 1000 Genomes project genotypes as reference panels.

## Data preparation

Input files of XPXP include:

- summay statistics files of the target population
- summay statistics files of the auxiliary population
- reference panel of the target population in plink 1 format
- reference panel of the auxiliary population in plink 1 format
- genetic covariance file
- environmental covariance file if sample overlap exist
- LD blocks file

The XPXP format GWAS summary statistics file has at least 5 fields:

- SNP: SNP rsid
- N: sample size
- Z: Z-scores
- A1: effect allele
- A2: other allele. 


e.g. 
```bash
$ head ./XPXP_demo/height_BBJ_summary_hm3.txt
SNP     N       Z       A1      A2      P
rs1048488       159095  0.659462022101568       T       C       0.509599125893493
rs3115850       159095  0.658089413463925       C       T       0.510480678080961
rs4040617       159095  -0.37360256722158497    G       A       0.708700023669772
```


The genetic covariance file is a covariance matrix with trait names as index and header:
```bash
$ cat ./XPXP_demo/gcov_height_BBJ.csv
,height_BBJ-EAS,height_Wood-EUR,height_UKB-EUR
height_BBJ-EAS,0.3894,0.1866552482397845,0.20741027144622678
height_Wood-EUR,0.1866552482397845,0.2596,0.2525
height_UKB-EUR,0.20741027144622678,0.2525,0.369
```


## Run XPXP to compute SNPs effect size

Once the imput files are formatted, XPXP will automatically process the datasets, including SNPs overlapping and allele matching.
Run XPXP with the following comand:
```bash
python [INSTALL PATH]/XPXP/src/XPXP.py \
--num_threads 40 \
--save ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv \
--gc_file ./XPXP_demo/gcov_height_BBJ.csv \
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_GIANT_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \
--sumst_names height_BBJ-EAS+height_Wood-EUR,height_UKB-EUR \
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \
--fix_effect_traits height_BBJ-EAS
```

XPXP returns the estimated SNPs effect size in ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv
```bash
$ head ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv
chr     SNP     bp      A1      A2      height_BBJ-EAS-muxpxp   height_Wood-EUR-muxpxp  height_UKB-EUR-muxpxp
1       rs3934834       1005806 T       C       0.001490031305500368    0.0008587749623481179   0.001348694388142856
1       rs3766191       1017587 T       C       0.001294311742987717    0.0007056128453579005   0.0011593703907873875
1       rs9442372       1018704 A       G       0.001354169642801807    0.0007882917460030307   0.0009947411386104607
```

## Evaluate the prediction performance using individual-level GWAS data

Input files:

- --geno genotype file of testing data (UKBB Chinese, n=1,439), plink1 version
- --beta the estimated SNPs effect size returned by XPXP

```bash
$ python [INSTALL PATH]/XPXP/src/Predict.py \
--save ./XPXP_demo_NotAvailable/predict_ukb_chn \
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
However, the individual-level GWAS data of UKBB is availabel due to the data sharing restriction. We therefore show how to use the height GWAS of UKBB Chinese as validation dataset

## Evaluate the prediction performance using GWAS summary statistics

We follow [XPASS](!https://github.com/YangLabHKUST/XPASS) to use the following equation:

<img src="https://latex.codecogs.com/svg.image?R^2=corr(y,\hat{y})^2=\left(\frac{cov(y,\hat{y})}{\sqrt{var(y)var(\hat{y})}}\right)^2=\left(\frac{z^T\tilde{\mu}/\sqrt{n}}{\sqrt{\tilde{\mu}^T\Sigma\tilde{\mu}}}\right)^2," title="R^2=corr(y,\hat{y})^2=\left(\frac{cov(y,\hat{y})}{\sqrt{var(y)var(\hat{y})}}\right)^2=\left(\frac{z^T\tilde{\mu}/\sqrt{n}}{\sqrt{\tilde{\mu}^T\Sigma\tilde{\mu}}}\right)^2," />

where z is the z-score of external summsry statistics, n is its sample size, <img src="https://latex.codecogs.com/svg.image?\tilde{\mu}" title="\tilde{\mu}" /> is the posterior mean of effect size at the standardized genotype scale, <img src="https://latex.codecogs.com/svg.image?\Sigma" title="\Sigma" /> is the SNPs correlation matrix.

Input files:

- --sumst_file GWAS summary statistics of UKBB Chinese height
- --beta the estimated SNPs effect size returned by XPXP
- --col_name specify the column name of SNPs effect size file

```bash
$ python [INSTALL PATH]/XPXP/src/PredictSS.py \
--ref_file ./XPXP_demo/1000G.EAS.QC.hm3.ind \
--sumst_file ./XPXP_demo/UKB_Chinese_height_GWAS_summary.txt \
--beta ./XPXP_demo/PM_height_BBJ-GIANT-UKB.csv  \
--col_name height_BBJ-EAS-muxpxp

Output: R2 for height_BBJ-EAS-muxpxp: 0.13039907941388484
```

Compared to XPXP trained on the BBJ and UKB only:
```bash
python [INSTALL PATH]/XPXP/src/XPXP.py \
--num_threads 40 \
--save ./XPXP_demo/PM_height_BBJ-UKB.csv \
--gc_file ./XPXP_demo/gcov_height_BBJ.csv \
--sumst_files ./XPXP_demo/height_BBJ_summary_hm3.txt,./XPXP_demo/height_UKB_summary_hm3.txt \
--sumst_names height_BBJ-EAS+height_UKB-EUR \
--use_snps ./XPXP_demo/snps_overlap.txt \ # use the same set of SNPs
--ld_block_file ./XPXP_demo/EAS_fourier_ls-all.bed \
--ref_files ./XPXP_demo/1000G.EAS.QC.hm3.ind,./XPXP_demo/1000G.EUR.QC.hm3.ind \
--fix_effect_traits height_BBJ-EAS

python [INSTALL PATH]/XPXP/src/PredictSS.py \
--ref_file ./XPXP_demo/1000G.EAS.QC.hm3.ind \
--sumst_file ./XPXP_demo/UKB_Chinese_height_GWAS_summary.txt \
--beta ./XPXP_demo/PM_height_BBJ-UKB.csv  \
--trait height_BBJ-EAS-muxpxp

Output: R2 for height_BBJ-EAS-muxpxp: 0.1256874019449264
```


# Development

The XPASS package is developed by Jiashun Xiao (jxiaoae@connect.ust.hk).

# Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.












