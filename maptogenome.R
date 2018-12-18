
require("minfi", quietly = TRUE)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly = TRUE)

RSet <- get(load("RSet.rdata"))

GRSet <- mapToGenome(RSet)

save(GRSet,file = "GRSet.rdata")

