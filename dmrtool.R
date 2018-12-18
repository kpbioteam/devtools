require("minfi", quietly = TRUE)

set <- get(load("GRSet.rdata"))

genomeranges <- as.data.frame(ranges(set))

beta <- getBeta(set)

pheno <- read.table("pheno.txt",skip=1)

type <- "categorical"

qCutoff <- 1

shrinkVar <- FALSE

tab <- read.table("ucsc.gtf")

tab <- tab[,-(11:14),drop=FALSE] 

colnames(tab) <- c("seqname","source","feature","start","end","score","strand", "frame","attributes", "names")

tab[,"source"] <- NULL

tab[,"frame"] <- NULL

tab[,"attributes"] <- NULL

dmp <- dmpFinder(beta, pheno[,"V2"], type = type, qCutoff = qCutoff, shrinkVar = shrinkVar)

dmp[,"names"] <- rownames(dmp)

data <- merge(dmp, tab, by="names",sort = TRUE)

data <- data[c("seqname","start","end","names","score","strand", "feature","intercept", "f", "pval","qval")]

write.table(data, file= "dmp.bedgraph", quote = FALSE,col.names = FALSE, row.names = FALSE, sep = "\t")
