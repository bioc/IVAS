MsqtlFinder <- function(expdata=NULL,snpdata=NULL,snplocus=NULL,GTFdata=NULL,met=NULL,Ncor=1,bplotout=NULL){
    snp.cn <- colnames(snpdata)
    exp.cn <- colnames(expdata)
    over.exp <- expdata[,is.element(colnames(expdata),snp.cn)]
    over.snp <- snpdata[,is.element(colnames(snpdata),exp.cn)]
    over.snp <- over.snp[,names(over.exp)]
    expdata <- over.exp
    snpdata <- over.snp
    cal.snp <- is.element(snplocus[,"SNP"],rownames(snpdata))
    sub.chr <- snplocus[cal.snp,"CHR"]
    chr <- unique(as.matrix(sub.chr))
    predictSQTL <- NULL
    registerDoParallel(cores=Ncor)
    for (j in 1:length(chr)){
        pa.result <- NULL
        print (paste("chr",chr[j],":processing",sep=""))
        ch.snp.num <- snplocus[,"CHR"] == chr[j]
        ch.snp.locus <- as.matrix(snplocus[ch.snp.num,])
        over.snp <- is.element(ch.snp.locus[,"SNP"],rownames(snpdata))
        ch.snps <- matrix(ch.snp.locus[over.snp,],ncol=3)
        colnames(ch.snps) <- colnames(ch.snp.locus)
        irange <- IRanges(start=as.integer(ch.snps[,"locus"]),end=as.integer(ch.snps[,"locus"]))
        ch.snps.range <- GRanges(seqnames=Rle(chr[j]),ranges=irange,metadata=ch.snps[,"SNP"])
        sub.snpdata <- snpdata[is.element(rownames(snpdata),ch.snps[,"SNP"]),]
        transdb <- chrseparate(GTFdata,chr[j])
        if (length(transdb) == 0){break}
        trans.exon.range <- exonsBy(transdb,by="tx")
        trans.intron.range <- intronsByTranscript(transdb)
        transnum <- names(trans.exon.range)
        txTable <- select(transdb, keys=transnum, columns=c("TXID","TXNAME","GENEID","TXSTART","TXEND"), keytype="TXID")
        irange <- IRanges(start=as.integer(txTable[,"TXSTART"]),end=as.integer(txTable[,"TXEND"]))
        generange <- GRanges(seqnames=Rle(chr[j]),ranges=irange,metadata=txTable[,"GENEID"])
        over.gene <- as.matrix(findOverlaps(ch.snps.range,generange,select="all"))[,"subjectHits"]
        uniqueID1 <- unique(generange[over.gene]$metadata)
        over.tx <- is.element(txTable[,"TXNAME"],rownames(expdata))
        uniqueID2 <- unique(txTable[over.tx,"GENEID"])
        tx.num <- table(txTable[,"GENEID"])
        not.one.gene <- which(tx.num>1)
        up2gene <- names(not.one.gene)
        up2gene <- up2gene[is.element(up2gene,uniqueID1) & is.element(up2gene,uniqueID2)]
        i <- NULL
        pa.result <- foreach(i=1:length(up2gene),.combine=rbind) %dopar% {
            per <- length(up2gene)/10
            if (i == as.integer(per)*5){print ("50%")}
            Altvalue <- findAlternative(up2gene[i],txTable,trans.exon.range,trans.intron.range,chr[j])
            overlapsnp <- findOversnp(Altvalue,ch.snps.range)
            sqtlfinder(Altvalue,overlapsnp,expdata,sub.snpdata,met)
            }
        if (length(pa.result)>0){
            predictSQTL <- unique(rbind(predictSQTL,pa.result))
        }}
    lm.sig.sqtl <- NULL
    glm.sig.sqtl <- NULL
    if (met == "lm" | met == "both"){
        lm.predictSQTL <- predictSQTL[which(predictSQTL[,"method"] == "lm"),]
        fdr.p <- p.adjust(lm.predictSQTL[,"P.value"],method="fdr",n=length(lm.predictSQTL[,"P.value"]))
        fdr.sig <- which(as.double(fdr.p) < 0.01)
        lm.sig.sqtl <- matrix(lm.predictSQTL[fdr.sig,],ncol=ncol(lm.predictSQTL))
        colnames(lm.sig.sqtl) <- colnames(predictSQTL)
        lm.sig.sqtl <- unique(lm.sig.sqtl[which(lm.sig.sqtl[,"per.P.value"]=="sig"),])
        lm.sig.sqtl <- matrix(lm.sig.sqtl,ncol=10)
        }
    if (met == "glm" | met == "both"){
        glm.predictSQTL <- predictSQTL[which(predictSQTL[,"method"] == "glm"),]
        fdr.p <- p.adjust(glm.predictSQTL[,"P.value"],method="fdr",n=length(glm.predictSQTL[,"P.value"]))
        fdr.sig <- which(as.double(fdr.p) < 0.01)
        glm.sig.sqtl <- matrix(glm.predictSQTL[fdr.sig,],ncol=ncol(glm.predictSQTL))
        colnames(glm.sig.sqtl) <- colnames(predictSQTL)
        glm.sig.sqtl <- unique(glm.sig.sqtl[which(glm.sig.sqtl[,"per.P.value"]=="sig"),])
        glm.sig.sqtl <- matrix(glm.sig.sqtl,ncol=10)
        }
    if (length(bplotout) > 0){
        print ("save the boxplot image")
        if (met == "lm" | met == "both"){saveBplot(lm.sig.sqtl,expdata,snpdata,snplocus,GTFdata,bplotout)}
        if (met == "glm" | met == "both"){saveBplot(glm.sig.sqtl,expdata,snpdata,snplocus,GTFdata,bplotout)}
        }
    total.sqtl <- rbind(lm.sig.sqtl,glm.sig.sqtl)
    colnames(total.sqtl) <- colnames(predictSQTL)
    return (total.sqtl)
    }

