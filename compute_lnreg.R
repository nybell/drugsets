##### ----- SETWD, LOAD DATA, PACKAGES ----- #####

args = commandArgs(trailingOnly = TRUE)  # setcorrs.rdata  (1), metadata (2), group tested (3), setsize (4), out (5), outdir (6) 
# {corrdata} {metaRdata} {group} {nsize} {out_name}_{group.upper()} {out_dir} {correct}

# library
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

# assign input args
corr.data.path    = as.character(args[1])
meta.data.path    = as.character(args[2])
test.group        = as.character(args[3])
setsize           = as.numeric(args[4])
out.name          = as.character(args[5])
out.dir           = as.character(args[6])
correct.type      = as.character(args[7])

# # testing
# corr.data.path = "/Users/nyb/phd_code/drugsets-may2026/output/enrich-out-magma-z-no-mhc/dgidb2-median/mp0/cond/adhd.no.mhc.pgc2022.magZ.cond.cov.TRUE_setcorrs.rdata"
# meta.data.path = "/Users/nyb/drugsets/DATA/meta.data.2026.RData"
# test.group     = "atc3"
# setsize        = 5
# out.name       = "adhd.no.mhc.pgc2022"
# out.dir        = "/Users/nyb/phd_code/drugsets-may2026/output/test"
# correct.type   = "bonf"

# load data from compute_corrs.r
load(corr.data.path)

# assign vars
set.info = corr.data[["set.info"]]
cov.matrix = corr.data[["cov.matrix"]]

# load gene set meta 
load(meta.data.path)
# load("/Users/nyb/drugsets/DATA/meta.data.2026.RData")
# load("/Users/nyb/drugsets/DATA/meta-data-old/metadata.RData")
# meta = master.save

##### ----- SPECIFY PARAMETERS ----- ######

# merge with set.info to remove any drugs not in analysis
merged = merge(set.info, meta.data, by = 'set.name') # drug_name, DRUG = old

# specifiy codes and create analysis column
if (tolower(test.group) == 'atc1') {
  codes = codes.atc1
  merged$data <- merged$atc1
} else if (tolower(test.group) == 'atc2') {
  codes = codes.atc2
  merged$data = merged$atc2
} else if (tolower(test.group) == 'atc3') {
  codes = codes.atc3
  merged$data = merged$atc3
} else if (tolower(test.group) == 'atc4') {
  codes = codes.atc4
  merged$data = merged$atc4
} else if (tolower(test.group) == 'moa_clue') {
  codes = codes.moa.clue
  merged$data = merged$moa_clue
} else if (tolower(test.group) == 'moa_targ_clue') {
  codes = codes.moa.targ.clue
  merged$data = merged$moa_targ_clue
} else if (tolower(test.group) == 'moa_chembl') {
  codes = codes.moa.chembl
  merged$data = merged$moa_chembl
} else if (tolower(test.group) == 'moa_targ_chembl') {
  codes = codes.moa.targ.chembl
  merged$data = merged$moa_targ_chembl
} else {
  message('Invalid input for enrichment group.')
}

##### ----- FUNCTIONS ----- #####

# function to run multiple linear regression
lnreg_dep <- function(X, cov.matrix, y) {
  
  N = length(y); K = ncol(X)                                               # define N & K
  df <- (N - K)                                                            # calculate degrees of freedom 
  Xt <- t(X)                                                               # tranpose X 
  W <- solve(Xt %*% cov.matrix %*% X)                                      # compute W matrix 
  B <- W %*% Xt %*% cov.matrix %*% y                                       # compute beta 
  sigma2 <- (t(y - (X %*% B)) %*% cov.matrix %*% (y - (X %*% B))) / df     # compute sigma squared
  t <- B[2] / sqrt(sigma2 * W[2,2])                                        # compute test statistic 
  p <- pt(t, df, lower.tail = F)                                           # compute two tailed p-value for test stat: add abs(), *2 for two-tailed
  se <- sqrt(sigma2 * W[2,2])                                                 # /sqrt(N) calculate standard error
  error <- qt(0.975,df=(N-K))*se                                              # calculate error margin
  left <- (B[2] - error)                                                      # calculate lower bound of ci
  right <- (B[2] + error)                                                     # calculate upper bound of ci
  ci <- c(left,right)                                                         # confidence interval              

  # save results 
  results <- c(B[2],se,ci,t,df,p)
  
  return(results)
  
}

##### ----- ANALYSIS ----- #####

# y vector 
y <- merged$stat

# 1's vector
ones <- as.vector(rep(c(1),each=length(merged$data)))

# gene set size vector
size <- as.vector(as.numeric(merged$size))

# define empty data frame 
results <- data.frame(matrix(NA, nrow = 0, ncol = 9))

# loop through all unique drug group codes in data and run linear regression
for (code in codes) {

  # get indices of group 
  ind <- which(sapply(merged$data, function(y) code %in% y))         # get indices of drugs in group
  
  if (length(ind) >= setsize) {   # args[4]

    # binarize group
    merged$bin <- 0                                                  # set all values equal to 0
    merged$bin[ind] <- 1                                             # set drugs in group equal to 1
    s <- merged$bin                                                  # s vector 
    
    # X matrix
    X <- cbind(ones,s,size,log(size))
    
    # get results 
    outcome <- lnreg_dep(X, cov.matrix, y)
    
    # paste results with code 
    outcome <- c(code, outcome, length(ind))
    
    # add to dataframe 
    results <- rbind(outcome, results)

  } else {
    next
  }

}

# add column names to results
colnames(results) <- c('GROUP','BETA','SE','CI_lower','CI_upper','T','DF','P', 'GROUP_N')
results <- results %>% mutate_at(c('BETA','SE','CI_lower','CI_upper','T','DF','P'), as.numeric)
results <- results[order(results$P),]

if (correct.type == 'bonf') {
  
  # Bonferroni correction
  bonf <- 0.05 / nrow(results)
  results_bonf <- subset(results, P < bonf)
  
} else if (correct.type == 'fdr') {
  
  # FDR correction
  results$fdrP <- p.adjust(results$P, method = 'BH')
  results_fdr <- subset(results, fdrP < 0.05)
  
}

# save results 
write.table(results, file.path(out.dir, sprintf('%s.groups.out', out.name)), row.names = FALSE, sep='\t', quote = FALSE)

# save corrected results
if (correct.type == 'bonf') {
  if (nrow(results_bonf) > 0) {
    write.table(results_bonf, file.path(out.dir, sprintf('%s.%s.groups.out', out.name, correct.type)), row.names = FALSE, sep='\t', quote = FALSE)
  }
} else if (correct.type == 'fdr') {
  if (nrow(results_fdr) > 0) {
    write.table(results_fdr, file.path(out.dir, sprintf('%s.%s.groups.out', out.name, correct.type)), row.names = FALSE, sep='\t', quote = FALSE)
  }
}


# args: 1:{corrdata} 2:{metaRdata} 3:{group} 4:{nsize} 5:{out_name}_{group.upper()} 6:{out_dir} 7:{correct}
##### ----- END ----- #####



