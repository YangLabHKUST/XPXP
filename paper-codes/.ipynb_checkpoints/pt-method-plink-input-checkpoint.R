#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(mvtnorm)
library(XPASS)
library(data.table)
library(RhpcBLASctl)
library(ieugwasr)
blas_set_num_threads(20)

thr_all <- c(5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.2,0.5,1)

ref_wegene <- paste0("/home/share/UKB/ld_ref_2k/ukb_AFR_ldpred_ref_2000_noMHC")

test_geno <- read_data("/home/share/xiaojs/xpa/African/combine_ipm_ukb/Afr_merge_test")
inp_sumss = args[1]
out_res = args[2]
#print(inp_sumss)
#print(out_res)

bim <- fread(paste0(ref_wegene,".bim"))
sumstats_wegene <- fread(inp_sumss)  # read raw sumstats with beta

bim <- bim[match(sumstats_wegene$ID,bim$V2),]

idx_snp <- match(sumstats_wegene$ID,test_geno$snps)
X <- test_geno$X[,idx_snp]
rm(test_geno)
gc()

out_mtag <- list()
PRS_mtag <- matrix(0,nrow(X),length(thr_all))

idx_flip <- which(bim$V5!=sumstats_wegene$ALT)     # find flipped allele
sumstats_wegene$BETA[idx_flip] <- -sumstats_wegene$BETA[idx_flip]    # flip beta
pval_wegene <- data.frame(rsid=sumstats_wegene$ID,pval=sumstats_wegene$P)

for(l in 1:length(thr_all)){
  if(all(pval_wegene$pval>thr_all[l])){   # if no significant SNPs passing thr
    clp_wegene <- NULL
  } else{
    clp_wegene <- ld_clump(pval_wegene, clump_kb=1000, clump_r2=0.1, clump_p=thr_all[l],
                           bfile=ref_wegene,
                           plink_bin="/home/share/xiaojs/software/plink")
    idx_clp <- match(clp_wegene$rsid,sumstats_wegene$ID)
    clp_wegene$beta <- sumstats_wegene$BETA[idx_clp]
    PRS_mtag[,l] <- matrix(X[,idx_clp],ncol = length(idx_clp)) %*% clp_wegene$beta
  }
  out_mtag[[as.character(thr_all[l])]] <- clp_wegene
  cat(l,"-th thr finished. \n")
}
colnames(PRS_mtag) <- as.character(thr_all)
#save(out_mtag,PRS_mtag,file = "/home/share/mingxuan/wegene_qc2/height_results/height_affy_MTAG_C+T.RData")
write.table(PRS_mtag,file=out_res,col.names = T,row.names = F,quote=F,sep="\t")

