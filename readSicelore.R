library(argparse)
library(stringr)
library(Seurat)

#############################
parser <- ArgumentParser(description='Loading Sicelore Pipeline Result')
parser$add_argument("--isog",
                    type="character",
                    default=NULL,
                    help="path of *_genematrix.txt")


parser$add_argument("--iso",
                    type="character",
                    default=NULL,
                    help="path of *_isomatrix.txt")

parser$add_argument("--junc",
                   type="character",
                   default=NULL,
                   help="path of *_juncmatrix.txt")

parser$add_argument("--illuminaObj",
                    type="character",
                    help="path of illuminaO Seurat Object or matrix path ,barcode must be same in Sicelore Result")

parser$add_argument("--outdir",
                    type="character",
                    default="Sicelore_result")

args <- parser$parse_args()

outdir<-args$outdir
options(stringsAsFactors=FALSE)
if(!dir.exists(outdir)){
        dir.create(outdir,recursive=TRUE)
}


############################
message("INFO : Loading illumina Object..")
seurat<-readRDS(args$illuminaObj)

############################
message("INFO : Loading Sicelore Result ..")
message("INFO : Loading isoform genematrix ..")
all = read.delim(args$isog, stringsAsFactors = F)
ISOG = all[,2:ncol(all)]
colnames(ISOG) <- paste(colnames(ISOG),"-1", sep="")
rownames(ISOG) <- all$geneId
ISOG <- ISOG[,colnames(seurat)]
seurat[["ISOG"]] <- CreateAssayObject(counts = ISOG)
seurat <- NormalizeData(object = seurat, assay = "ISOG")
seurat <- ScaleData(object = seurat, assay = "ISOG")

###########################
message("INFO : Loading isoform matrix ..")
all = read.delim(args$iso, stringsAsFactors = F)
ISO = all[,4:ncol(all)]
colnames(ISO) <- paste(colnames(ISO),"-1", sep="")
rownames(ISO) <- paste(all$geneId, all$transcriptId, sep="..")
idx <- grep("undef", rownames(ISO), invert = TRUE)
ISO <- ISO[idx,colnames(seurat)]
seurat[["ISO"]] <- CreateAssayObject(counts = ISO)
seurat <- NormalizeData(object = seurat, assay = "ISO")
seurat <- ScaleData(seurat, assay="ISO")

############################
message("INFO : Loading junc matrix ..")
all = read.delim(args$junc, stringsAsFactors = F)
JUNC = all[,2:ncol(all)]
colnames(JUNC) <- paste(colnames(JUNC),"-1", sep="")
rownames(JUNC) <- str_replace(all$junctionId, ":", "..")
JUNC <- JUNC[,colnames(seurat)]
seurat[["JUNC"]] <- CreateAssayObject(counts = JUNC)
seurat <- NormalizeData(object = seurat, assay = "JUNC")
seurat <- ScaleData(seurat, assay="JUNC")


message("INFO : Save ..")
saveRDS(seurat,file.path(outdir,"SiceloreSeuratObj.rds"))
message("INFO : Done!")



