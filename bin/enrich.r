#!/usr/bin/env Rscript

#######################################################################
#######################################################################
## Created on Nov. 10, 2020 enrichment analysis
## Copyright (c) 2020 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################
c("optparse", "ChIPpeakAnno", "clusterProfiler", "pathview", "biomaRt")
# set libPath to pwd
pwd <- getwd()
pwd <- file.path(pwd, "lib")
dir.create(pwd)
.libPaths(c(pwd, .libPaths()))

library(ChIPpeakAnno)
library(clusterProfiler)
library(pathview)
library(biomaRt)
library(optparse)
library(ggplot2)
library(scales)
writeLines(as.character(packageVersion("ChIPpeakAnno")), "ChIPpeakAnno.version.txt")
writeLines(as.character(packageVersion("clusterProfiler")), "clusterProfiler.version.txt")
writeLines(as.character(packageVersion("pathview")), "pathview.version.txt")
writeLines(as.character(packageVersion("biomaRt")), "biomaRt.version.txt")


option_list <- list(make_option(c("-s", "--genome"), type="character", default=NULL, help="genome assembly name", metavar="string"),
                    make_option(c("-n", "--ucscname"), type="character", default=NULL, help="ucscname", metavar="string"),
                    make_option(c("-q", "--fdr"), type="double", default=0.05, help="false discovery rate cutoff value", metavar="double"),
                    make_option(c("-o", "--output"), type="character", default=".", help="output folder", metavar="string"),
                    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of cores", metavar="integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$genome)){
    print_help(opt_parser)
    stop("Please provide genome name.", call.=FALSE)
}
if (!is.null(opt$output)){
    pf <- opt$output
}else{
    pf <- "."
}
FDRcutoff <- 0.05
if (!is.null(opt$fdr)){
    FDRcutoff <- as.numeric(opt$fdr)[1]
    if(FDRcutoff<0 | FDRcutoff>1){
        stop("cutoff fdr should be a number between 0 and 1.")
    }
}

dir.create(pf, recursive = TRUE, showWarnings = FALSE)

attr <- c("hsapiens_homolog_associated_gene_name")
scientificName <- c("GRCh37"="Homo sapiens",
                    "GRCh38"="Homo sapiens",
                    "GRCm38"="Mus musculus",
                    "TAIR10"="Arabidopsis thaliana",
                    "EB2"="Bacillus subtilis 168",
                    "UMD3.1"="Bos taurus",
                    "WBcel235"="Caenorhabditis elegans",
                    "CanFam3.1"="Canis familiaris",
                    "GRCz10"="Danio rerio",
                    "BDGP6"="Drosophila melanogaster",
                    "EquCab2"="Equus caballus",
                    "EB1"="Escherichia coli str. K12a",
                    "Galgal4"="Gallus gallus",
                    "Gm01"="Glycine max",
                    "Mmul_1"="Macaca mulatta",
                    "IRGSP-1.0"="Oryza sativa japonica",
                    "CHIMP2.1.4"="Pan troglodytes",
                    "Rnor_6.0"="Rattus norvegicus",
                    "R64-1-1"="Saccharomyces cerevisiae",
                    "EF2"="Schizosaccharomyces pombe",
                    "Sbi1"="Sorghum bicolor",
                    "Sscrofa10.2"="Sus scrofa",
                    "AGPv3"="Zea mays",
                    "hg38"="Homo sapiens",
                    "hg19"="Homo sapiens",
                    "hg38"="Homo sapiens",
                    "mm10"="Mus musculus",
                    "bosTau8"="Bos taurus",
                    "ce10"="Caenorhabditis elegans",
                    "canFam3"="Canis familiaris",
                    "danRer10"="Danio rerio",
                    "dm6"="Drosophila melanogaster",
                    "equCab2"="Equus caballus",
                    "galGal4"="Gallus gallus",
                    "panTro4"="Pan troglodytes",
                    "rn6"="Rattus norvegicus",
                    "sacCer3"="Saccharomyces cerevisiae",
                    "susScr3"="Sus scrofa")
if(is.null(opt$ucscname)){
    opt$ucscname <- opt$genome
}
if(opt$ucscname %in% names(scientificName)){
    scientificName <- scientificName[opt$ucscname]
}else{
    if(opt$genome %in% names(scientificName)){
        scientificName <- scientificName[opt$genome]
    }else{
        stop("Not a valid genome for enrichment analysis.")
    }
}

scientificName2martTab <- function(.ele){
    .ele <- tolower(strsplit(.ele)[[1]])
    paste0(substring(.ele[1], 1, 1), .ele[2])
}
org <- egOrgMap(scientificName)
lib <- .libPaths()
if(file.access(lib[1], mode=2)!=0){
    pwd <- getwd()
    pwd <- file.path(pwd, "lib")
    dir.create(pwd)
    .libPaths(c(pwd, lib))
}
while(!require(org, character.only = TRUE)) BiocManager::install(org, update = FALSE, ask = FALSE)
library(org, character.only = TRUE)
organism <- ChIPpeakAnno:::.findKEGGRESTOrganismName(org)
org <- get(org)

shortStrs <- function(strs, len=60){
    if(length(strs)==0) return(strs)
    strs <- as.character(strs)
    shortStr <- function(str, len=60){
        stopifnot(length(str)==1)
        stopifnot(is.character(str))
        if(nchar(str)<=len) return(str)
        strs <- strsplit(str, " ")[[1]]
        nc <- nchar(strs)
        nclast <- nc[length(nc)] + 3
        paste0(substring(str, first = 1, last = len-nclast), "...", strs[length(strs)])
    }
    strs <- sapply(strs, shortStr, len=len)
    make.unique(strs)
}

files <- dir(".", "DEtable.*.anno.csv", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("padj", files)]

gmt <- "ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v7.2.symbols.gmt"
for(file in files){
    data <- read.csv(file)
    if(length(data$gene)!=nrow(data)){
        data$gene <- data$symbol
    }
    if(length(data$gene)!=nrow(data)){ ## interactive data
        data$gene <- ifelse(data$shortestDistance.anchor1<data$shortestDistance.anchor2,
                            data$symbol.anchor1, data$symbol.anchor2)
    }
    data.s <- data[data$FDR<FDRcutoff, , drop=FALSE]
    pff <- file.path(pf, basename(dirname(file)))
    dir.create(pff, recursive = TRUE)
    if(nrow(data.s)>1){
        gene.df <- bitr(data.s$gene,
                        fromType=ifelse(grepl("^ENS", as.character(data.s$gene)[1]),
                                        "ENSEMBL", "SYMBOL"),
                        toType = c("ENTREZID", "SYMBOL"), OrgDb = org)

        ego <- sapply(c("BP", "MF", "CC"), function(.onto){
            enrichGO(gene = gene.df$ENTREZID,
                    OrgDb = org,
                    ont = .onto,
                    readable = TRUE
            )
        })

        null <- mapply(ego, names(ego), FUN=function(.ele, .name){
            write.csv(.ele, file.path(pff, paste0("GO.", .name, ".enrichment.for.padj0.05.csv")))

            .ele <- as.data.frame(.ele)
            if(nrow(.ele)>1){
                colnames(.ele) <- tolower(colnames(.ele))
                .ele$qvalue <- -log10(.ele$pvalue)
                plotdata <- .ele[!is.na(.ele$qvalue), c("description", "qvalue", "count")]
                if(nrow(plotdata)>20) plotdata <- plotdata[1:20, ]
                plotdata$description <- shortStrs(plotdata$description)
                ggplot(plotdata, aes(x=reorder(description, -qvalue), y=qvalue, fill=count, label=count)) +
                    scale_fill_gradient2(low = muted("blue"), high = muted("red"), oob = scales::squish) +
                    geom_bar(stat="identity") + scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
                    geom_text(vjust=-.1) +
                    xlab("") + ylab("-log10(p-value)") +
                    theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
                ggsave(file.path(pff, paste("GO.", .name, ".enrichment.for.padj0.05.top.pdf", sep = ".")), width = 6, height = 6)
            }
        })

        kk <- enrichKEGG(gene = gene.df$ENTREZID, organism = organism)
        kk <- as.data.frame(kk)
        eid <- strsplit(kk$geneID, "\\/")
        symbol <- lapply(eid, function(.ele) gene.df[match(.ele, gene.df$ENTREZID), "SYMBOL"])
        symbol <- sapply(symbol, paste, collapse="/")
        kk$geneSYM <- symbol
        write.csv(kk, file.path(pff, "KEGGenrichment.for.padj0.05.csv"))

        ds <- data$log2FoldChange
        if(length(ds)!=nrow(data)){
            ds <- data$logFC
        }
        gene.df <- bitr(data$gene,
                        fromType=ifelse(grepl("^ENS", as.character(data.s$gene)[1]),
                                        "ENSEMBL", "SYMBOL"),
                        toType = c("ENTREZID", "SYMBOL"), OrgDb = org)
        names(ds) <-
            gene.df[match(data$gene,
                        gene.df[, ifelse(grepl("^ENS",
                                            as.character(data.s$gene)[1]),
                                        "ENSEMBL", "SYMBOL")]), "ENTREZID"]
        ds <- ds[!is.na(names(ds))]
        ds <- ds[!is.na(ds)]
        if(length(ds)>2){
            p <- file.path(pff, "pathview")
            dir.create(p)
            for (pid in kk$ID[-seq.int(45)]) {
                tryCatch(pv.out <- pathview(gene.data = ds, pathway.id = pid,
                                            species=organism, kegg.dir = p,
                                            kegg.native=TRUE),
                            error=function(.e) message(.e))
            }
        }

        pngs <- dir(".", "pathview.png")
        file.rename(pngs, file.path(p, pngs))

        rnk <- file.path(pff, sub(".csv", ".preranked.rnk", basename(file)))
        stat <- ifelse("stat" %in% colnames(data), "stat", "F")
        if(scientificName!="Homo sapiens"){
            mart <- useMart("ensembl", paste0(scientificName2martTab(scientificName), "_gene_ensembl"))
            bm <- getBM(values = unique(data$ensembl_id),
                        attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
                        filters = "ensembl_gene_id",
                        mart = mart)
            data$hsapiens_homolog_associated_gene_name <-
                bm[match(data$ensembl_id, bm$ensembl_gene_id),
                "hsapiens_homolog_associated_gene_name"]
            data.rnk <- data[, c("hsapiens_homolog_associated_gene_name", stat)]
        }else{
            data.rnk <- data[, c("gene", stat)]
        }
        colnames(data.rnk)[1] <- c("hsapiens_gene_name")
        data.rnk <- data.rnk[order(data[[stat]], decreasing=TRUE), ]
        data.rnk <- data.rnk[!is.na(data.rnk[, 1]), , drop=FALSE]
        data.rnk <- data.rnk[data.rnk[, 1]!="", , drop=FALSE]
        data.rnk <- data.rnk[!duplicated(data.rnk[, 1]), , drop=FALSE]
        write.table(data.rnk, file = rnk, quote=FALSE, row.names = FALSE, sep="\t")
        rpt_label <- "c2.all.v7.2"
        outfolder <- file.path(pff, "GSEA")
        cmd <- paste("gsea-cli.sh GSEAPreranked -gmx", gmt, "-norm meandiv -nperm 1000 -rnk",
                    rnk, "-scoring_scheme weighted -rpt_label", rpt_label,
                    "-create_svgs true -make_sets true -plot_top_x 100 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out",
                    outfolder)
        tryCatch(system(cmd), error=function(e){message(e)})
    }else{
        writeLines(paste("No available enrichment results by FDR <", FDRcutoff), file.path(pff, sub(".csv", ".txt", basename(file))))
    }
}
