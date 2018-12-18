require("minfi", quietly = TRUE)
options(warn = -1)

GRSet <- get(load("GRSet.rdata"))

pheno <- read.table("pheno.txt",skip = 1)

group <- pheno$V2

pair <- factor(pheno$V3)

design.matrix <- model.matrix(~ group + pair)

maxGap <- 250

if(is.null(GRSet$cluster)){
  cluster = NULL
  maxGap = maxGap
} else {
  cluster = GRSet$cluster
  maxGap = NULL
}

cutoff <- 0.1
B <- 0 #as.numeric('$number_of_resamples')
nullMethod <- "permutation"
coef <- 2

dmrs <- bumphunter(GRSet,
                   design = design.matrix, 
                   cluster = cluster,
                   maxGap = maxGap,
                   cutoff = cutoff, 
                   nullMethod = nullMethod,
                   B = B)

dmrGR <- with(dmrs$table,GRanges(chr,IRanges(start,end),area=area,value=value))

dmrGR <- as.data.frame(dmrGR)

colnames(dmrGR) <- c("seqnames","start","end","width","strand","area","value")

dmrGR$strand <- NULL

dmrGR$area <- NULL

write.table(dmrGR, file= 'dmr.bedgraph', quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")

