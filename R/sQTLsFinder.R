sQTLsFinder <- function(ASdb=NULL,Total.snpdata=NULL,Total.snplocus=NULL,GroupSam=NULL,method="lm",CalIndex=NULL,Ncor=1,out.dir=NULL){
    cal.met <- method
    total.result <- NULL
    Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(Exon.ratio.mat) <- c("ES","ASS","IR")
    total.list <- Exon.ratio.mat
    if (ASdb@"Ratio"[["ES"]][1,1] != "NA"){
        Exon.ratio.mat$"ES" <- ASdb@"Ratio"[["ES"]]
        if (length(CalIndex) != 0)    Exon.ratio.mat$"ES" <- rbind(Exon.ratio.mat$"ES"[is.element(Exon.ratio.mat$"ES"[,"Index"],CalIndex),])
    }
    if (ASdb@"Ratio"[["ASS"]][1,1] != "NA"){
        Exon.ratio.mat$"ASS" <- ASdb@"Ratio"[["ASS"]]
        if (length(CalIndex) != 0)    Exon.ratio.mat$"ASS" <- rbind(Exon.ratio.mat$"ASS"[is.element(Exon.ratio.mat$"ASS"[,"Index"],CalIndex),])
    }
    if (ASdb@"Ratio"[["IR"]][1,1] != "NA"){
        Exon.ratio.mat$"IR" <- ASdb@"Ratio"[["IR"]]
        if (length(CalIndex) != 0)    Exon.ratio.mat$"IR" <- rbind(Exon.ratio.mat$"IR"[is.element(Exon.ratio.mat$"IR"[,"Index"],CalIndex),])
    }
    
    subtypes <- names(Exon.ratio.mat[Exon.ratio.mat != "NA" & lengths(Exon.ratio.mat) != 0])
    registerDoParallel(cores=Ncor)
    Total.snplocus <- gsub(" ","",as.matrix(Total.snplocus))
    len.chr <- grep("chr",(Exon.ratio.mat[[subtypes[1]]][,"Nchr"]))
    if (length(len.chr) != 0){
        Total.snplocus[,2] <- paste("chr",gsub("chr","",Total.snplocus[,2]),sep="")
    }
    else if (length(len.chr) == 0){
        Total.snplocus[,2] <- gsub("chr","",Total.snplocus[,2])
    }
    
    
    if (length(Total.snpdata) != 0){
        Total.snpdata <- gsub(" ","",as.matrix(Total.snpdata))
        total.result <- NULL
        final.result <- lapply(paste(subtypes,method,sep="-"),function(each.type){
            each.type <- unlist(strsplit(each.type,"-"))
            cal.met <- each.type[2]
            each.type <- each.type[1]
            sub.exon.ratio.mat <- gsub(" ","",as.matrix(Exon.ratio.mat[[each.type]]))
            sub.snplocus <- rbind(Total.snplocus[is.element(Total.snplocus[,"CHR"],unique(sub.exon.ratio.mat[,"Nchr"])),])
            inter.snp <- intersect(rownames(Total.snpdata),sub.snplocus[,"SNP"])
            sub.snpdata <- rbind(Total.snpdata[inter.snp,])
            rownames(sub.snpdata) <- inter.snp
            if (length(sub.exon.ratio.mat) != 0 & length(sub.snplocus) != 0 & length(sub.snpdata) != 0){
                over.samples <- intersect(colnames(sub.exon.ratio.mat),colnames(Total.snpdata))
                sub.exon.ratio <- sub.exon.ratio.mat
                called.packages <- c("lme4","GenomicRanges","GenomicFeatures")
                Total.chr <- unique(as.matrix(sub.snplocus[,"CHR"]))
                Total.chr <- Total.chr[order(as.integer(Total.chr))]
                for (j in 1:length(Total.chr)){
                    if (cal.met != "boxplot"){
                        print (paste("-------------------Processing : chr",Total.chr[j]," (",each.type,") -------------------",sep=""))
                    }
                    ch.sub.exon.ratio <- rbind(sub.exon.ratio[sub.exon.ratio[,"Nchr"] == Total.chr[j],])
                    ch.snp.locus <- rbind(sub.snplocus[sub.snplocus[,"CHR"] == Total.chr[j],])
                    inter.snp <- intersect(rownames(sub.snpdata),ch.snp.locus[,"SNP"])
                    ch.snp.data <- rbind(sub.snpdata[inter.snp,])
                    rownames(ch.snp.data) <- inter.snp
                    i = NULL
                    pa.result <- foreach(i=1:nrow(ch.sub.exon.ratio),.packages=called.packages,.combine=rbind) %dopar% {
                        each.sub.exon.ratio <- rbind(ch.sub.exon.ratio[i,])
                        test.expdata <- rbind(each.sub.exon.ratio[,over.samples])
                        ex.regions <- do.call(rbind,strsplit(each.sub.exon.ratio[,is.element(colnames(each.sub.exon.ratio),c("DownEX","UpEX","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX"))],"-"))
                        ex.regions <- unlist(strsplit(ex.regions,","))
                        ex.regions <- ex.regions[ex.regions != "NA" & ex.regions != "NaN" & !is.na(ex.regions)]
                        ex.regions <- cbind(min(as.integer(ex.regions)),max(as.integer(ex.regions)))
                        colnames(ex.regions) <- c("start","end")
                        EX.region.range <- GRanges(seqnames=Rle(each.sub.exon.ratio[,"Nchr"]),ranges=IRanges(start=as.integer(ex.regions[,"start"]),end=as.integer(ex.regions[,"end"])))
                        EX.region.range <- list(EX.region.range)
                        names(EX.region.range) <- "alterIntron"
                        SNPragne <- GRanges(seqnames=Rle(ch.snp.locus[,"CHR"]),ranges=IRanges(start=as.integer(ch.snp.locus[,"locus"]),end=as.integer(ch.snp.locus[,"locus"])),metadata=ch.snp.locus[,"SNP"])
                        overlapsnp <- findOversnp(EX.region.range,SNPragne)
                        if (length(overlapsnp) != 0 & length(test.expdata[test.expdata!="NA"]) != 0){
                            inter.snp <- intersect(rownames(ch.snp.data),overlapsnp[,"snp"])
                            if (length(inter.snp) != 0){
                                test.snpdata <- rbind(ch.snp.data[inter.snp,over.samples])
                                rownames(test.snpdata) <- inter.snp
                                test.snplocus <- rbind(ch.snp.locus[is.element(ch.snp.locus[,"SNP"],overlapsnp[,"snp"]),])
                                sig.result <- CalSigSNP(ratio.mat=test.expdata,snp.mat=test.snpdata,overlapsnp=overlapsnp,
                                    each.snplocus=test.snplocus,chr=Total.chr[j],each.gene=each.sub.exon.ratio[,"EnsID"],GroupSam,method=cal.met)
                                if (cal.met == "boxplot")    sig.result
                                else if (cal.met != "boxplot" & length(sig.result) != 0){
                                    inter.cn <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX","DownEX","UpEX","Types","Diff.P","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX","RetainEX")
                                    inter.cn <- inter.cn[is.element(inter.cn,colnames(each.sub.exon.ratio))]
                                    pre.inf <- rep(rbind(each.sub.exon.ratio[,inter.cn]),nrow(sig.result))
                                    pre.inf <- matrix(pre.inf,nrow=nrow(sig.result),byrow=TRUE)
                                    colnames(pre.inf) <- inter.cn
                                    cbind(sig.result[,"snpid"],pre.inf,rbind(sig.result[,is.element(colnames(sig.result),c("pByGeno","diff","pByGroups","OR","lowCI","highCI","met"))]))
                                }
                                else {NULL}
                            }
                            else {NULL}
                        }
                        else {NULL}
                    }
                    if (cal.met == "boxplot")    total.result <- pa.result
                    else {
                        total.result <- rbind(total.result,pa.result)
                        if (length(total.result) != 0)    colnames(total.result)[1] <- "SNP"
                    }
                }
                if (cal.met == "boxplot") total.result
                else    unique(total.result)
            }
            else    NULL
        })
        if (cal.met == "boxplot")    return (final.result[[1]])
    }
    names(final.result) <- subtypes
    if (length(final.result$"ES") != 0){
        each.result <- unique(final.result$"ES")
        rownames(each.result) <- 1:nrow(each.result)
        p.num.geno <- which(colnames(each.result)=="pByGeno")
        each.result <- cbind(rbind(each.result[,1:p.num.geno]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num.geno+1):ncol(each.result)]))
        if (length(GroupSam) != 0){
            p.num.group <- which(colnames(each.result)=="pByGroups")
            each.result <- cbind(rbind(each.result[,1:p.num.group]),FdrByGroups=p.adjust(as.double(each.result[,"pByGroups"]),"fdr"),rbind(each.result[,as.integer(p.num.group+1):ncol(each.result)]))
        }
        total.list$"ES" <- each.result
    }
    if (length(final.result$"ASS") != 0){
        each.result <- unique(final.result$"ASS")
        rownames(each.result) <- 1:nrow(each.result)
        p.num.geno <- which(colnames(each.result)=="pByGeno")
        each.result <- cbind(rbind(each.result[,1:p.num.geno]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num.geno+1):ncol(each.result)]))
        if (length(GroupSam) != 0){
            p.num.group <- which(colnames(each.result)=="pByGroups")
            each.result <- cbind(rbind(each.result[,1:p.num.group]),FdrByGroups=p.adjust(as.double(each.result[,"pByGroups"]),"fdr"),rbind(each.result[,as.integer(p.num.group+1):ncol(each.result)]))
        }
        total.list$"ASS" <- each.result
    }
    if (length(final.result$"IR") != 0){
        each.result <- unique(final.result$"IR")
        rownames(each.result) <- 1:nrow(each.result)
        p.num.geno <- which(colnames(each.result)=="pByGeno")
        each.result <- cbind(rbind(each.result[,1:p.num.geno]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num.geno+1):ncol(each.result)]))
        if (length(GroupSam) != 0){
            p.num.group <- which(colnames(each.result)=="pByGroups")
            each.result <- cbind(rbind(each.result[,1:p.num.group]),FdrByGroups=p.adjust(as.double(each.result[,"pByGroups"]),"fdr"),rbind(each.result[,as.integer(p.num.group+1):ncol(each.result)]))
        }
        total.list$"IR" <- each.result
    }
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",GroupDiff=ASdb@"GroupDiff",sQTLs=total.list,Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_sQTLs",sep=""))
        write.table(final.result[["ES"]],paste(out.dir,"/AS_sQTLs/ES_sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.result[["ASS"]],paste(out.dir,"/AS_sQTLs/ASS_sQTLs.txt",sep=""),sep='\t',quote=FALSE)
        write.table(final.result[["IR"]],paste(out.dir,"/AS_sQTLs/IR_sQTLs.txt",sep=""),sep='\t',quote=FALSE)
    }
    return (ASdb)
}

