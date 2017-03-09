saveBplot <- function(ASdb=ASdb,Total.snpdata=NULL,Total.snplocus=NULL,CalIndex=NULL,out.dir=NULL){
    BoxsQTLs <- function(snp.result=NULL,test.snp.mat=NULL){
        snpplots.re <- lapply(rownames(test.snp.mat),function(snpid){
            each.result <- rbind(snp.result[snp.result[,"SNP"] == snpid])
            if (length(each.result) != 0){
                colnames(each.result) <- colnames(snp.result)
                each.snp.mat <- rbind(test.snp.mat[snpid,])
                rownames(each.snp.mat) <- snpid
                box.list <- sQTLsFinder(ASdb,each.snp.mat,snplocus,method="boxplot",CalIndex=CalIndex)
                box.exp <- box.list$"exp"
                box.OR <- box.list$"Ratios"
                text.box <- rbind(each.result[,!is.element(colnames(each.result),c("diff","Index","Strand"))])
                text.box[,"pByGeno"] <- round(as.double(text.box[,"pByGeno"]),digits=4)
                text.box <- gsub("-","-\n",text.box)            
                TotalRatio <- lapply(names(box.exp),function(each.geno){
                    each.genoBox <- cbind(exp=box.exp[[each.geno]],each.geno)
                    each.genoBox
                })
                TotalRatio <- do.call(rbind,TotalRatio)
                TotalRatio <- data.frame(as.double(TotalRatio[,"exp"]),TotalRatio[,"each.geno"])
                colnames(TotalRatio) <- c("Ratio","Genotype")
                gplot.result <- ggplot(data=TotalRatio,aes(x=Genotype,y=Ratio,fill=Genotype))+geom_boxplot(width=0.3,outlier.shape=NA)+geom_jitter(width = 0.03)
                gplot.result <- list(plot=gplot.result,text=text.box)
                gplot.result
            }
            else    NULL
        })
        names(snpplots.re) <- rownames(test.snp.mat)
        return (snpplots.re)
    }
    Genotype <- NULL
    Ratio <- NULL
    snplocus <- Total.snplocus
    if (length(grep("ES",CalIndex)) != 0)    testType <- "ES"
    if (length(grep("ASS",CalIndex)) != 0)    testType <- "ASS"
    if (length(grep("IR",CalIndex)) != 0)    testType <- "IR"
    
    testASmodel <- ASdb@"SplicingModel"[[testType]]
    if (ncol(testASmodel) == 1) return (NULL)
    exonRatio <- ASdb@"Ratio"[[testType]]
    sQTLs <- ASdb@"sQTLs"[[testType]]
    testRatio <- NULL
    testsQTLs <- NULL
    
    testASmodel <- rbind(testASmodel[testASmodel[,"Index"] == CalIndex,])
    if (length(testASmodel) == 0) return (NULL)
    if (ncol(exonRatio) != 1)    testRatio <- rbind(exonRatio[exonRatio[,"Index"] == CalIndex,])
    if (ncol(sQTLs) != 1)    testsQTLs <- rbind(sQTLs[sQTLs[,"Index"] == CalIndex,])
    
    
    test.snp <- Total.snpdata[is.element(rownames(Total.snpdata),testsQTLs[,"SNP"]),]
    sQTLs.plot <- BoxsQTLs(testsQTLs,test.snp)
    system(paste("mkdir -p ",out.dir,sep=""))
    processing <- lapply(sQTLs.plot,function(each.plot){
        png(paste(out.dir,"/",each.plot$"text"[,"SNP"],"_",testsQTLs[,"Index"],"_",each.plot$"text"[,"pByGeno"],"_",each.plot$"text"[,"met"],".png",sep=""))
        print (each.plot$"plot")
        dev.off()
    })
    return (NULL)
}
