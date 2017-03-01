sQTLsFinder <- function(ASdb=NULL,Total.snpdata=NULL,Total.snplocus=NULL,GroupSam=NULL,method="lm",CalIndex=NULL,Ncor=1,out.dir=NULL){
  CalsigSNP_fre <- function(test.Fre){
    A.groups <- GroupSam$"GroupA"
    B.groups <- GroupSam$"GroupB"
    test.result <- lapply(1:nrow(test.Fre),function(each.nums){
      each.Fre <- rbind(test.Fre[each.nums,])
      A.fre <- each.Fre[,"FreA"]
      B.fre <- each.Fre[,"FreB"]
      A.odd <- A.fre/(1-A.fre)
      B.odd <- B.fre/(1-B.fre)
      odds.ratio <- A.odd/B.odd
      each.samples <- matrix(c(length(A.groups)*A.fre,length(A.groups)*(1-A.fre),length(B.groups)*B.fre,length(B.groups)*(1-B.fre)),ncol=2)
      sqr.value <- sum(1/each.samples)
      low.CI <- log(odds.ratio) - 1.96*sqrt(sqr.value)
      high.CI <- log(odds.ratio) + 1.96*sqrt(sqr.value)
      chi.p <- chisq.test(each.samples)$"p.value"
      c(odds.ratio,chi.p,low.CI,high.CI)
    })
    test.result <- do.call(rbind,test.result)
    colnames(test.result) <- c("Odds","Pvalue","low","high")
    return (test.result)
  }
  CalSigSNP <- function(ratio.mat,snp.mat,overlapsnp,each.snplocus,chr,each.gene){
    colnames(overlapsnp) <- c("snp","locus")
    realNA <- colnames(ratio.mat)[ratio.mat != "NA" & ratio.mat != "NaN"]
    test.exp <- rbind(as.double(ratio.mat[,realNA]))
    colnames(test.exp) <- realNA
    snprn <- rownames(snp.mat)
    test.snp <- rbind(snp.mat[,realNA])
    rownames(test.snp) <- snprn
    each.row.ratio <- test.exp
    diff.value <- 0.1
    numsamp <- 10
    if (length(realNA) > numsamp){
      stac.result <- lapply(1:nrow(test.snp),function(test.each.num){
        pre.result.glm <- NULL
        pre.result.lm <- NULL
        snpid <- rownames(test.snp)[test.each.num]
        test.each.snp <- rbind(snp.mat[test.each.num,realNA])
        skEX.ratio <- test.exp*100
        inEX.ratio <- (1-test.exp)*100
        genoform <- as.matrix(test.snp[test.each.num,])
        lm.geno <- unique(unlist(strsplit(rownames(table(genoform)),"")))
        if (length(lm.geno)>1){
          if(length(lm.geno)==2){
            pregeno <- gsub(lm.geno[1],0,genoform)
            pregeno <- gsub(lm.geno[2],1,pregeno)
          }
          else if(length(lm.geno)==3){
            pregeno <- gsub(lm.geno[1],0,genoform)
            pregeno <- gsub(lm.geno[2],1,pregeno)
            pregeno <- gsub(lm.geno[3],1,pregeno)
          }
          matrix.geno <- matrix(as.integer(gsub(10,1,pregeno)),ncol=length(pregeno))
          matrix.geno <- matrix(as.integer(gsub(01,1,matrix.geno)),ncol=length(matrix.geno))
          matrix.geno <- matrix(as.integer(gsub(11,2,matrix.geno)),ncol=length(matrix.geno))
          pregeno <- matrix(pregeno,ncol=length(matrix.geno))
          colnames(matrix.geno) <- rownames(genoform)
          lm.genoform <- matrix.geno
          Ratios <- NULL
          if (length(GroupSam) != 0){
            Asamples <- intersect(colnames(matrix.geno),GroupSam$"GroupA")
            Bsamples <- intersect(colnames(matrix.geno),GroupSam$"GroupB")
            Afre <- sum(matrix.geno[,Asamples])/(length(matrix.geno[,Asamples])*2)
            Bfre <- sum(matrix.geno[,Bsamples])/(length(matrix.geno[,Bsamples])*2)
            if (Afre != "NaN" & Bfre != "NaN"){
              TotalFre <- rbind(c(FreA=Afre,FreB=Bfre))
              rownames(TotalFre) <- "SNP"
              Ratios <- CalsigSNP_fre(TotalFre)
            }
          }
          auo.matrix <- data.frame(ratio=c(each.row.ratio),group=factor(lm.genoform))
          each.row.ratio.vec <- each.row.ratio[1,]
          split.g <- split(each.row.ratio.vec,auo.matrix$group)
          mamean <- sapply(split.g,median)
          if (method == "boxplot"){
            genos <- c(paste(lm.geno[1],lm.geno[1],sep=""),paste(lm.geno[1],lm.geno[2],sep=""),paste(lm.geno[2],lm.geno[2],sep=""))
            names(genos) <- c(0,1,2)
            names(split.g) <- genos[names(split.g)]
            list(exp=split.g,Ratios=Ratios)
          }
          else {
            dis.test <- "NonDiff"
            if(length(mamean) == 2){
              if(abs(mamean[1]-mamean[2])>diff.value){
                dis.test <- "diff"
              }
            }
            else if (length(mamean) == 3){
              if((mamean[1]-mamean[2]>diff.value & mamean[2]-mamean[3]>diff.value) | (mamean[1]-mamean[2]< -1*diff.value & mamean[2]-mamean[3]< -1*diff.value)){
                dis.test <- "diff"
              }
            }
            lm.auo.pvalue <- "NaN"
            glm.auo.pvalue <- "NaN"
            if (method == "lm" | method == "both"){
              auo <- summary(lm(formula = ratio ~ group, data = auo.matrix))
              if (length(auo$fstatistic[1])>0 &length(auo$fstatistic[2])>0 & length(auo$fstatistic[3])>0){
                lm.auo.pvalue <- pf(auo$fstatistic[1],auo$fstatistic[2],auo$fstatistic[3],lower.tail=FALSE)
                if (lm.auo.pvalue != "NaN" & lm.auo.pvalue != "NA"){
                  pre.result.lm <- cbind(snpid,lm.auo.pvalue,dis.test,Ratios,"lm")
                }
              }
            }
            if (method == "glm" | method == "both"){
              calmatrix <- cbind(round(inEX.ratio),round(skEX.ratio))
              calmatrix <- matrix(calmatrix,ncol=2)
              rownames(calmatrix) <- colnames(inEX.ratio)
              obs <- c(1:length(lm.genoform))
              test1 <- try(glmer(calmatrix ~ as.character(as.matrix(lm.genoform)) + (1|obs),na.action = na.exclude,family=binomial),silent=TRUE)
              test2 <- try(glmer(calmatrix ~ 1 +(1|obs),na.action = na.exclude,family=binomial),silent=TRUE)
              if (!( inherits(test1,"try-error") | inherits(test2,"try-error"))){
                glm.auo.pvalue <- anova(test1,test2)$"Pr(>Chisq)"[2]
                if (glm.auo.pvalue != "NaN" & glm.auo.pvalue != "NA"){
                  pre.result.glm <- cbind(snpid,glm.auo.pvalue,dis.test,Ratios,"glm")
                }
              }
            }
            final.result <- rbind(pre.result.lm,pre.result.glm)
            final.result
          }
        }
      })
    }
    if (method == "boxplot")  return (stac.result[[1]])
    stac.result <- do.call(rbind,stac.result)
    if (is.element("Odds",colnames(stac.result)))  colnames(stac.result) <- c("snpid","pByGeno","diff","pByGroups","OR","lowCI","highCI","met")
    else colnames(stac.result) <- c("snpid","pByGeno","diff","met")
    return (stac.result)
  }
  
  
  total.result <- NULL
  
  
  Exon.ratio.mat <- list(as.matrix("NA"),as.matrix("NA"),as.matrix("NA"))
  names(Exon.ratio.mat) <- c("ES","ASS","IR")
  total.list <- Exon.ratio.mat
  if (ASdb@"Ratio"[["ES"]][1,1] != "NA"){
    Exon.ratio.mat$"ES" <- ASdb@"Ratio"[["ES"]]
    if (length(CalIndex) != 0)  Exon.ratio.mat$"ES" <- rbind(Exon.ratio.mat$"ES"[is.element(Exon.ratio.mat$"ES"[,"Index"],CalIndex),])
  }
  if (ASdb@"Ratio"[["ASS"]][1,1] != "NA"){
    Exon.ratio.mat$"ASS" <- ASdb@"Ratio"[["ASS"]]
    if (length(CalIndex) != 0)  Exon.ratio.mat$"ASS" <- rbind(Exon.ratio.mat$"ASS"[is.element(Exon.ratio.mat$"ASS"[,"Index"],CalIndex),])
  }
  if (ASdb@"Ratio"[["IR"]][1,1] != "NA"){
    Exon.ratio.mat$"IR" <- ASdb@"Ratio"[["IR"]]
    if (length(CalIndex) != 0)  Exon.ratio.mat$"IR" <- rbind(Exon.ratio.mat$"IR"[is.element(Exon.ratio.mat$"IR"[,"Index"],CalIndex),])
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
    final.result <- lapply(subtypes,function(each.type){
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
          if (method != "boxplot"){
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
                sig.result <- CalSigSNP(test.expdata,test.snpdata,overlapsnp,test.snplocus,Total.chr[j],each.sub.exon.ratio[,"EnsID"])
                if (method == "boxplot")  sig.result
                else if (method != "boxplot"){
                  inter.cn <- c("Index","EnsID","Strand","Nchr","1stEX","2ndEX","DownEX","UpEX","Types","Diff.P","ShortEX","LongEX","NeighborEX","ShortNeighborEX","LongNeighborEX","RetainEX")
                  inter.cn <- inter.cn[is.element(inter.cn,colnames(each.sub.exon.ratio))]
                  pre.inf <- rep(rbind(each.sub.exon.ratio[,inter.cn]),nrow(sig.result))
                  pre.inf <- matrix(pre.inf,nrow=nrow(sig.result),byrow=T)
                  colnames(pre.inf) <- inter.cn
                  cbind(sig.result[,"snpid"],pre.inf,rbind(sig.result[,is.element(colnames(sig.result),c("pByGeno","diff","pByGroups","OR","lowCI","highCI","met"))]))
                }
              }
              else {NULL}
            }
            else {NULL}
          }
          if (method == "boxplot")  total.result <- pa.result
          else {
            total.result <- rbind(total.result,pa.result)
            if (length(total.result) != 0)  colnames(total.result)[1] <- "SNP"
          }
        }
        if (method == "boxplot") total.result
        else  unique(total.result)
      }
      else  NULL
    })
    if (method == "boxplot")  return (final.result[[1]])
  }
  names(final.result) <- subtypes
  if (length(final.result$"ES") != 0){
    each.result <- unique(final.result$"ES")
    rownames(each.result) <- 1:nrow(each.result)
    p.num <- which(colnames(each.result)=="pByGeno")
    each.result <- cbind(rbind(each.result[,1:p.num]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num+1):ncol(each.result)]))
    total.list$"ES" <- each.result
  }
  if (length(final.result$"ASS") != 0){
    each.result <- unique(final.result$"ASS")
    rownames(each.result) <- 1:nrow(each.result)
    p.num <- which(colnames(each.result)=="pByGeno")
    each.result <- cbind(rbind(each.result[,1:p.num]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num+1):ncol(each.result)]))
    total.list$"ASS" <- each.result
  }
  if (length(final.result$"IR") != 0){
    each.result <- unique(final.result$"IR")
    rownames(each.result) <- 1:nrow(each.result)
    p.num <- which(colnames(each.result)=="pByGeno")
    each.result <- cbind(rbind(each.result[,1:p.num]),FdrByGeno=p.adjust(as.double(each.result[,"pByGeno"]),"fdr"),rbind(each.result[,as.integer(p.num+1):ncol(each.result)]))
    total.list$"IR" <- each.result
  }
  ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=ASdb@"Ratio",GroupDiff=ASdb@"GroupDiff",sQTLs=total.list,Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
  if (length(out.dir) != 0){
    system(paste("mkdir -p ",out.dir,"/AS_sQTLs",sep=""))
    write.table(final.result[["ES"]],paste(out.dir,"/AS_sQTLs/ES_sQTLs.txt",sep=""),sep='\t',quote=F)
    write.table(final.result[["ASS"]],paste(out.dir,"/AS_sQTLs/ASS_sQTLs.txt",sep=""),sep='\t',quote=F)
    write.table(final.result[["IR"]],paste(out.dir,"/AS_sQTLs/IR_sQTLs.txt",sep=""),sep='\t',quote=F)
  }
  return (ASdb)
}