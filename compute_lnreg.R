##### ----- SETWD, LOAD DATA, PACKAGES ----- #####

args = commandArgs(trailingOnly = TRUE)  # setcorrs.rdata  (1), metadata (2), group tested (3), setsize (4), out (5), outdir (6) 

# # for testing purposes
# setwd('~/drugsets/')
# load('test_t2d.rdata')
# load('/Users/nyb/drugsets/DATA/metadata.rdata')

# library
library(tidyr)
library(dplyr)
library(FSA)

# load data from compute_corrs.r
load(args[1])

# load gene set meta 
load(args[2])

##### ----- SPECIFY PARAMETERS ----- ######

# merge with set.info to remove any drugs not in analysis
merged <- merge(set.info, meta, by.x = 'set.name', by.y = 'DRUG')

# specifiy codes and create analysis column
if (tolower(args[3]) == 'atc') {
  codes <- codesAtc
  merged$data <- merged$atc
} else if (tolower(args[3]) == 'moa') {
  codes <- codesMoa
  merged$data <- merged$moa
} else if (tolower(args[3]) == 'ind') {
  codes <- codesInd
  merged$data <- merged$indication
} else {
  print('Invalid input for enrichment group.')
}

##### ----- FUNCTIONS ----- #####

# function to run multiple linear regression
lnreg_dep <- function(X, set.corrs.inv, y) {
  
  N = length(y); K = ncol(X) # define N & K
  Xt <- t(X)     # tranpose X 
  W <- solve(Xt %*% set.corrs.inv %*% X)     # compute W matrix 
  B <- W %*% Xt %*% set.corrs.inv %*% y      # compute beta 
  sigma <- t(y - (X %*% B)) %*% set.corrs.inv %*% (y - (X %*% B)) / (N - K)     # compute sigma squared
  t <- B[2] / sqrt(sigma * W[2,2])     # compute test statistic 
  p <- pt(t, (N - K), lower.tail = F)            # compute two tailed p-value for test stat: add abs(), *2 for two-tailed
  
  results <- c(B[2],sigma,t,p)    # save results
  
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
results <- data.frame(matrix(NA, nrow = 0, ncol = 5))

# loop through all unique ATC codes in data set and run linear regression
for (code in codes) {

  # get indices of group 
  ind <- which(sapply(merged$data, function(y) code %in% y))            # get group indices
  
  if (length(ind) >= as.numeric(5)) {   # args[4]

    # binarize group
    merged$bin <- 0                                                  # set all values equal to 0
    merged$bin[ind] <- 1                                             # set drugs in group equal to 1
    s <- merged$bin                                                  # s vector 
    
    # X matrix
    X <- cbind(ones,s,size,log(size))
    
    # get results 
    outcome <- lnreg_dep(X, set.corrs.inv, y)
    
    # paste results with code 
    outcome <- append(code,outcome)
    
    # add to dataframe 
    results <- rbind(outcome, results)
    
  } else {
    next
  }
  
  
}

# add column names to results
colnames(results) <- c('GROUP','BETA','SIGMA','T','P')
results <- results %>% mutate_at(c('BETA','SIGMA','T','P'), as.numeric)
results <- results[order(results$P),]

if (args[7] == 'bonf') {
  
  # Bonferroni correction
  bonf <- 0.05 / nrow(results)
  results_bonf <- subset(results, P < bonf)
  
} else if (args[7] == 'fdr') {
  
  # FDR correction
  results$fdrP <- p.adjust(results$P, method = 'BH')
  results_fdr <- subset(results, fdrP < 0.05)
  
}

# save results 
write.table(results, sprintf('%s/lnreg.%s.%s.csv', args[6], args[3], args[5]), row.names = FALSE, sep=',', quote = FALSE)   # raw 

# save corrected results
if (args[7] == 'bonf') {
  if (nrow(results_bonf) > 0) {
    write.table(results_bonf, sprintf('%s/lnreg.%s.%s.%s.csv',args[6], args[7], args[3], args[5]), row.names = FALSE, sep=',', quote = FALSE)  # Bonf
  }
} else if (args[7] == 'fdr') {
  if (nrow(results_fdr) > 0) {
    write.table(results_fdr, sprintf('%s/lnreg.%s.%s.%s.csv',args[6], args[7], args[3], args[5]), row.names = FALSE, sep=',', quote = FALSE)  # Bonf
  }
}

##### ----- END ----- #####



