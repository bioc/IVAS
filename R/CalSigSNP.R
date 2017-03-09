CalSigSNP <- function(ratio.mat=NULL,snp.mat=NULL,overlapsnp=NULL,each.snplocus=NULL,chr=NULL,each.gene=NULL,GroupSam=NULL,method="lm"){
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
    stac.result <- NULL
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
    else if (length(stac.result) == 0) return (NULL)
    if (method == "boxplot")    return (stac.result[[1]])
    stac.result <- do.call(rbind,stac.result)
    if (is.element("Odds",colnames(stac.result)))    colnames(stac.result) <- c("snpid","pByGeno","diff","pByGroups","OR","lowCI","highCI","met")
    else colnames(stac.result) <- c("snpid","pByGeno","diff","met")
    return (stac.result)
}