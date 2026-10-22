##### ----- COMPUTE CORRELATION MATRIX ----- #####

args = commandArgs(trailingOnly = TRUE)     # input arguments in order: .raw file, gene set file, MAGMA output file, pheno name, working directory

# # for testing purposes
# # setwd('~/drugsets/')
# raw.file = "/Users/nyb/phd_code/drugsets-may2026/output/magma-v110-s1s2-pops-no-mhc/adhd.no.mhc.pgc2022.step2.genes.raw"
# set.file = '/Users/nyb/drugsets/DATA/GENESETS/drug_mp0_dgidb_median_ensembl_cond_sets.txt'
# res.file = '/Users/nyb/phd_code/drugsets-may2026/output/test/adhd.no.mhc.pgc2022.magZ.cond.drug_mp0.gsa.out'
# setwd('/Users/nyb/phd_code/drugsets-may2026/output/test/')
# correct = as.logical(TRUE)

# set working directory 
setwd(args[5])

# set file names to load
raw.file = args[1]                          # "SCZ_SAMPLE.genes.raw"
set.file = args[2]                          # "entrez_cond_sets.txt"
res.file = args[3]                          # "cond.gsa.out"
correct = as.logical(args[6])

if (correct == TRUE) {
  prune.thresh = 0.1  #used for gene-set correlation matrix as well
  
  raw.data = strsplit(scan(raw.file, what="", comment.char="#", sep="\n"), " ")
  info.length = length(raw.data[[1]])
  
  
  raw.info = data.frame(matrix(sapply(raw.data, head, info.length), ncol=info.length, byrow=T), stringsAsFactors=F)[,-(3:4)]
  names(raw.info) = c("gene", "chr", "nsnps", "nparam", "nsamp", "mac", "zstat")
  for (i in 3:ncol(raw.info)) raw.info[,i] = as.numeric(raw.info[,i])
  
  
  raw.corrs = lapply(raw.data, tail, -info.length)
  chr.names = unique(raw.info$chr)
  projection = list(); proj.index = list(); no.proj = 0
  for (chr in chr.names) {
  	curr.raw = raw.corrs[raw.info$chr == chr]
  	curr.size = length(curr.raw)
  	curr.corr = matrix(0, ncol=curr.size, nrow=curr.size)
  	for (i in 2:curr.size) {
  		len = length(curr.raw[[i]])
  		if (len > 0) curr.corr[i,1:len+(i-len-1)] = as.numeric(curr.raw[[i]])
  	}
  	curr.corr = curr.corr + t(curr.corr)
  	diag(curr.corr) = 1
  
  	eig = eigen(curr.corr)
  	use = which(eig$values >= prune.thresh)
    projection[[chr]] = eig$vectors[,use] %*% diag(1/sqrt(eig$values[use]))
  	proj.index[[chr]] = 1:length(use) + no.proj
  	no.proj = no.proj + length(use)
  }
  
  project.data = function(data) {
  	out = matrix(NA, nrow=no.proj, ncol=ncol(data))
  	for (chr in chr.names) {
  		out.index = proj.index[[chr]]
  		out[out.index,] = t(projection[[chr]]) %*% data[raw.info$chr == chr,]
  	}
  	return(out)
  }
  
  
  residualize = cbind(data.matrix(raw.info[,c("nsnps", "nparam", "nsamp")]), 1/raw.info$mac)
  residualize = residualize[,apply(residualize, 2, var) > 0]
  residualize = cbind(1, scale(cbind(residualize, log(residualize))))
  residualize = project.data(residualize)
  
  outcome = project.data(matrix(raw.info$zstat, ncol=1))
  outcome = lm(outcome ~ residualize - 1)$residuals
  
  
  gsa.res = read.table(res.file, header=T, stringsAsFactors=F)
  sets = strsplit(scan(set.file, what="", sep="\n"), "\t")
  sets.use = c(gsa.res$FULL_NAME, "druggable")
  
  set.names = sapply(sets, head, 1)
  set.genes = lapply(sets[set.names %in% sets.use], tail, -1)
  set.names = set.names[set.names %in% sets.use]
  
  sets = matrix(0, nrow=nrow(raw.info), ncol=length(set.names))
  for (i in 1:length(set.names)) sets[raw.info$gene %in% set.genes[[i]],i] = 1
  sets = project.data(sets)
  
  #assuming conditioning on 'druggable'
  covar = cbind(residualize[,1], sets[,set.names == "druggable"])
  sets = sets[,set.names != "druggable"]
  
  ctc.inv = solve(t(covar) %*% covar)
  cts = t(covar) %*% sets
  det = apply(sets^2, 2, sum) - apply(cts * (ctc.inv %*% cts), 2, sum)
  V = sweep(sets - covar %*% ctc.inv %*% cts, 2, det, FUN="/")
  
  beta = t(V) %*% outcome
  set.corrs = cov2cor(t(V) %*% V)
  
  eig = eigen(set.corrs)
  use = which(eig$values >= prune.thresh)
  set.corrs.inv = eig$vectors[,use] %*% diag(1/eig$values[use]) %*% t(eig$vectors[,use])
  
  set.info = data.frame(
    set.name = gsa.res$FULL_NAME,
    size = gsa.res$NGENES,
    stat = gsa.res$BETA / gsa.res$SE,
    stringsAsFactors=F
  )
  
  # create list to save
  corr.data = list(set.info = set.info, cov.matrix = set.corrs.inv)
  # save
  save(corr.data, file=file.path(args[5], paste0(args[4], "_setcorrs.rdata")))            
  
  
} else if (correct == FALSE) {
  
  # load it 
  gsa.res = read.table(res.file, header=T, stringsAsFactors=F)
  
  # remove "druggable" from the data
  gsa.res = subset(gsa.res, FULL_NAME != "druggable")
  
  # Create an identity matrix
  identity.matrix = diag(nrow(gsa.res))
  
  set.info = data.frame(
    set.name = gsa.res$FULL_NAME,
    size = gsa.res$NGENES,
    stat = gsa.res$BETA / gsa.res$SE,
    stringsAsFactors=F
  )
  
  # create list to save
  corr.data = list(set.info = set.info, cov.matrix = identity.matrix)
  # save
  save(corr.data, file=file.path(args[5], paste0(args[4], "_setcorrs.rdata")))          
  
}
  

