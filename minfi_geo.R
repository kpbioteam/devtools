
require("minfi", quietly = TRUE)
require("GEOquery", quietly = TRUE)
require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)

options(warn = -1)
gset <- getGEO('GSE42752')
class(gr)

if(length(gset)==0) stop("Empty list retrieved from GEO.")
if(length(gset)>1){
  warning("More than one ExpressionSet found:\n",names(gset),"\nUsing entry ",1)
  gset <- gset[[1]]
} else gset <- gset[[1]]

platform <- annotation(gset)

gr <-  getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19,orderByLocation = TRUE)

locusNames <- names(gr)

sampleNames(gset) <- gset$title

common <- intersect(locusNames, featureNames(gset))
if(length(common)==0)
  stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")

ind1 <- match(common,fData(gset)$Name)
ind2 <- match(common,locusNames)

preprocessing <- c(rg.norm=paste0('See GEO ','$input_geoid',' for details'))

what <- "Beta"
if(what=="Beta"){
  out <- GenomicRatioSet(gr=gr[ind2,],
                         Beta=gset[ind1,,drop=FALSE],
                         annotation=c(array = "IlluminaHumanMethylation450k"),
                         preprocessMethod=preprocessing)
} else {
  
  out <- GenomicRatioSet(gr=gr[ind2,],
                         M=mat[ind1,,drop=FALSE],
                         annotation=c(array = "IlluminaHumanMethylation450k"),
                         preprocessMethod=preprocessing)
}
save(out, file = "GRSetfromGEO.rdata")

