sQTLsFinder <- function(ASdb=NULL,Total.snpdata=NULL,Total.snplocus=NULL,GroupSam=NULL,method="lm",CalIndex=NULL,Ncor=1,out.dir=NULL){
    testSNP <- function(ex.r.mat,MP){
        if (!length(ex.r.mat))    return (NULL)
        te.chr <- intersect(Total.snplocus[,"CHR"],unique(ex.r.mat[,"Nchr"]))
        ex.r.mat <- rbind(ex.r.mat[is.element(ex.r.mat[,"Nchr"],te.chr),])
        snplocus <- rbind(Total.snplocus[is.element(Total.snplocus[,"CHR"],
            te.chr),])
        o.snp <- intersect(rownames(Total.snpdata),snplocus[,"SNP"])
        snpdata <- rbind(Total.snpdata[o.snp,])
        rownames(snpdata) <- o.snp
        if (!length(ex.r.mat) | !length(snplocus) | !length(snpdata)){
            return (NULL)
        }
        ov.sam <- intersect(colnames(ex.r.mat),colnames(snpdata))
        mulsq <- function(i){
            ea.ex.ra <- rbind(ex.r.mat[i,])
            ea.cn <- colnames(ea.ex.ra)
            te.exp <- rbind(ea.ex.ra[,ov.sam])
            te.ex.renge <- ea.ex.ra[,is.element(ea.cn,ex.cns)]
            ex.res <- do.call(rbind,strsplit(rbind(te.ex.renge),"-"))
            ex.res <- unlist(strsplit(ex.res,","))
            nonNA <- ex.res != "NA" & ex.res != "NaN" & !is.na(ex.res)
            ex.res <- ex.res[nonNA]
            ex.res <- cbind(min(as.integer(ex.res)),max(as.integer(ex.res)))
            colnames(ex.res) <- c("start","end")
            Irang <- IRanges(ex.res[,"start"],ex.res[,"end"])
            EX.range <- GRanges(Rle(ea.ex.ra[,"Nchr"]),ranges=Irang)
            EX.range <- list(EX.range)
            names(EX.range) <- "alterIntron"
            SnpIran <- IRanges(as.integer(snplocus[,"locus"]),
                as.integer(snplocus[,"locus"]))
            SNPragne <- GRanges(Rle(snplocus[,"CHR"]),ranges=SnpIran,
                metadata=snplocus[,"SNP"])
            overlapsnp <- findOversnp(EX.range,SNPragne)
            inter.snp <- intersect(rownames(snpdata),overlapsnp[,"snp"])
            if (length(overlapsnp) & length(which(te.exp != "NA")) & 
                length(inter.snp)){
                te.snp <- rbind(snpdata[inter.snp,ov.sam])
                rownames(te.snp) <- inter.snp
                snp.lo <- is.element(snplocus[,"SNP"],overlapsnp[,"snp"])
                te.snplo <- rbind(snplocus[snp.lo,])
                sig.re <- CalSigSNP(ratio.mat=te.exp,snp.mat=te.snp,
                    overlapsnp=overlapsnp,each.snplocus=
                    te.snplo,chr=ea.ex.ra[,"Nchr"],
                    each.gene=ea.ex.ra[,"EnsID"],
                    GroupSam=GroupSam,method=method)
                if (method == "boxplot")    sig.re
                else if (method != "boxplot" & length(sig.re) != 0){
                    ov.cn <- inter.cns[is.element(inter.cns,ea.cn)]
                    pre.inf <- rep(rbind(ea.ex.ra[,ov.cn]),nrow(sig.re))
                    pre.inf <- matrix(pre.inf,nrow=nrow(sig.re),byrow=TRUE)
                    colnames(pre.inf) <- ov.cn
                    ov.p <- is.element(colnames(sig.re),p.cns)
                    cbind(sig.re[,"snpid"],pre.inf,rbind(sig.re[,ov.p]))
                }
                else    NULL
            }
            else    NULL
        }
        pa.result <- bplapply(seq_len(nrow(ex.r.mat)),mulsq,BPPARAM=MP)
        if (method  != "boxplot"){
            pa.result <- do.call(rbind,pa.result)
            if (length(pa.result))    colnames(pa.result)[1] <- "SNP"
        }
        return (pa.result)
    }
    FDR.cal <- function(each.result){
        if (!length(each.result))    return (NULL)
        rownames(each.result) <- 1:nrow(each.result)
        p.ge <- which(colnames(each.result)=="pByGeno")
        fdrva <- p.adjust(as.double(each.result[,"pByGeno"]),"fdr")
        fi.ea <- rbind(each.result[,1:p.ge])
        se.ea <- rbind(each.result[,as.integer(p.ge+1):ncol(each.result)])
        each.result <- cbind(fi.ea,FdrByGeno=fdrva,se.ea)
        p.gr <- which(colnames(each.result)=="pByGroups")
        if (length(p.gr)){
            fdrva <- p.adjust(as.double(each.result[,"pByGroups"]),"fdr")
            fi.ea <- rbind(each.result[,1:p.gr])
            se.ea <- each.result[,as.integer(p.gr+1):ncol(each.result)]
            se.ea <- rbind(se.ea)
            each.result <- cbind(fi.ea,FdrByGroups=fdrva,se.ea)
        }
        return (each.result)
    }
    ex.cns <- c("DownEX","UpEX","ShortEX","LongEX","NeighborEX",
        "ShortNeighborEX","LongNeighborEX")
    inter.cns <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX","DownEX",
        "UpEX","Types","Diff.P","ShortEX","LongEX","NeighborEX",
        "ShortNeighborEX","LongNeighborEX","RetainEX")
    p.cns <- c("pByGeno","diff","pByGroups","OR","lowCI","highCI","met")
    i <- NULL
    pa.result <- NULL
    sig.re <- NULL
    FdrByGroups <- NULL
    FdrByGeno <- NULL
    MP <- SnowParam(workers=Ncor,type="SOCK")
    Total.snplocus <- gsub(" ","",as.matrix(Total.snplocus))
    Total.snpdata <- gsub(" ","",as.matrix(Total.snpdata))
    Ex.f.re <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
    names(Exon.ratio.mat) <- c("ES","ASS","IR")
    ratio.mat <- ASdb@"Ratio"
    if (ncol(ratio.mat$ES) != 1){
        ES.mat <- ratio.mat$ES
        if (length(CalIndex)){
            ES.mat <- rbind(ES.mat[is.element(ES.mat[,"Index"],CalIndex),])
        }
        if (length(ES.mat)){
            ES.re <- testSNP(ES.mat,MP)
            if (method == "boxplot")    return (ES.re[[1]])
            ES.re <- FDR.cal(ES.re)
            if (length(ES.re))    Exon.ratio.mat$ES <- ES.re
        }
    }
    if (ncol(ratio.mat$ASS) != 1){
        ASS.mat <- ratio.mat$ASS
        if (length(CalIndex)){
            ASS.mat <- rbind(ASS.mat[is.element(ASS.mat[,"Index"],CalIndex),])
        }
        if (length(ASS.mat)){
            ASS.re <- testSNP(ASS.mat,MP)
            if (method == "boxplot")    return (ASS.re[[1]])
            ASS.re <- FDR.cal(ASS.re)
            if (length(ASS.re))    Exon.ratio.mat$ASS <- ASS.re
        }
    }
    if (ncol(ratio.mat$IR) != 1){
        IR.mat <- ratio.mat$IR
        if (length(CalIndex)){
            IR.mat <- rbind(IR.mat[is.element(IR.mat[,"Index"],CalIndex),])
        }
        if (length(IR.mat)){
            IR.re <- testSNP(IR.mat,MP)
            if (method == "boxplot")    return (IR.re[[1]])
            IR.re <- FDR.cal(IR.re)
            if (length(IR.re))    Exon.ratio.mat$IR <- IR.re
        }
    }
    ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",
        GroupDiff=ASdb@"GroupDiff",sQTLs=Exon.ratio.mat,
        Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
    if (length(out.dir) != 0){
        p.out <- paste(out.dir,"/AS_sQTLs/",sep="")
        system(paste("mkdir -p ",p.out,sep=""))
        write.table(ES.re,paste(p.out,"/ES_sQTLs.txt",
            sep=""),sep='\t',quote=FALSE)
        write.table(ASS.re,paste(p.out,"/ASS_sQTLs.txt",
            sep=""),sep='\t',quote=FALSE)
        write.table(IR.re,paste(p.out,"/IR_sQTLs.txt",sep=""),
            sep='\t',quote=FALSE)
    }
    return (ASdb)
}

