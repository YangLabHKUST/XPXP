args <- commandArgs(trailingOnly = TRUE)

library(bigsnpr)
library(dplyr)

obj.bigSNP <- snp_attach(args[1])
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores()
# Read external summary statistics
sumstats <- bigreadr::fread2(args[2])
        
sumstats$n_eff <- sumstats$N
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
df_beta <- snp_match(sumstats, map)
df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos

POS2 <- obj.bigSNP$map$genetic.dist
tmp <- tempfile(tmpdir = ".")

for (chr in 1:22) {
  
  # print(chr)
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)
  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3  # file size in GB

(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))

h2_est <- ldsc[["h2"]]

multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(0.01,0.5, length.out = 5),
                               ncores = NCORES)
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)

df_beta$beta_auto = beta_auto  
                    
df = select(df_beta, 'chr', 'rsid', 'pos', 'a1', 'a0', 'beta_auto')
write.csv(df,args[3],row.names = FALSE,quote=F)