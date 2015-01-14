calSignificant <- function(tx.gene=NULL,total.locus=NULL,exon.locus=NULL,intron.locus=NULL,info.strand=NULL,overapvalue=NULL,chrnum=NULL,expdata=NULL,snpdata=NULL,method=NULL){
    if (total.locus[2,ncol(total.locus)] == "AS"){
        Nmax.num <- ncol(total.locus)
        max.intron.num <- ncol(total.locus)-1
        Nmax=as.integer(total.locus[1,Nmax.num])
        max.intron <- total.locus[1,max.intron.num]
        total.locus <- total.locus[,-(max.intron.num:Nmax.num)]
        }
    if (ncol(total.locus) == 0){return (NULL)}
    over.info <- is.element(overapvalue[,"range"],total.locus["bigrange",])
    if (!is.element("TRUE",over.info)){return (NULL)}
    oversnp.mat <- NULL
    oversnp.mat <- matrix(overapvalue[over.info,],ncol=2)
    colnames(oversnp.mat) <- c("range","snp")
    result <- NULL
    sp.ex.locus <- strsplit(exon.locus,"-")
    exinfo <- do.call(rbind,sp.ex.locus)
    colnames(exinfo) <- c("start","end")
    ex.start <- exinfo[,"start"]
    ex.end <- exinfo[,"end"]
    samp.snp <- NULL
    for (z in 1:ncol(total.locus)){
        cal.start <- as.integer(strsplit(total.locus["bigrange",z],"-")[[1]][1])
        cal.end <- as.integer(strsplit(total.locus["bigrange",z],"-")[[1]][2])
        fi.tx.eval <- (cal.start >= tx.gene[,"TXSTART"] & cal.end <= tx.gene[,"TXEND"])
        se.tx.eval <- (cal.end <= tx.gene[,"TXSTART"] & cal.start >= tx.gene[,"TXEND"])
        tx.gene.start <- tx.gene[fi.tx.eval | se.tx.eval,]
        target.snp <- oversnp.mat[which(oversnp.mat[,"range"] == total.locus["bigrange",z]),"snp"]
        samp.snp <- snpdata[is.element(rownames(snpdata),target.snp),]
        if (nrow(samp.snp) == 0){next}
        if (total.locus["type",z]=="AS"){
            tx.gene.name <- NULL
            int.cal.start <- strsplit(intron.locus[grep(cal.start,intron.locus)],"-")
            int.cal.end <- strsplit(intron.locus[grep(cal.end,intron.locus)],"-")
            over.start <- unique(do.call(rbind,int.cal.start)[,2])
            over.end <- unique(do.call(rbind,int.cal.end)[,1])
            mat.end <- lapply(over.start,function(o){
                int.over <- strsplit(intron.locus[grep(o,intron.locus)],"-")
                do.call(rbind,int.over)[,1]
                })
            over.end <- unique(c(over.end,unlist(mat.end)))
            mat.start <- lapply(over.end,function(o){
                int.over <- strsplit(intron.locus[grep(o,intron.locus)],"-")
                do.call(rbind,int.over)[,2]
                })
            over.start <- unique(c(over.start,unlist(mat.end)))
            over.st.name <- lapply(over.start,function(o){
                names(intron.locus[grep(o,intron.locus)])
                })
            over.en.name <- lapply(over.end,function(o){
                names(intron.locus[grep(o,intron.locus)])
                })
            tx.gene.name <- unlist(c(over.st.name,over.en.name))
            tx.gene.start <- na.omit(tx.gene[unique(tx.gene.name),])
            nextstart <- "not"
            if (z < Nmax){
                tx.gene.start <- tx.gene.start[(cal.end >= tx.gene.start[,"TXSTART"] & cal.end <= tx.gene.start[,"TXEND"]),]
                grep.intron <- grep(cal.end,intron.locus)
                if (length(grep.intron) != 0){
                    nextstart <- grep.intron[names(intron.locus[grep.intron]) == names(intron.locus[grep.intron+1])]
                    if (length(nextstart) != 0){
                        next.intron <- strsplit(intron.locus[as.integer(nextstart)+1],"-")
                        nextstart <- unique(do.call(rbind,next.intron)[,1])
                        }
                    else {nextstart <- "not"}
                    }
                cal.intron <- names(intron.locus[grep.intron])
                next.intron <- names(intron.locus[grep(nextstart[1],intron.locus)])
                cal.next <- unique(c(cal.intron,next.intron))
                have.int <- list(tx.gene.start[cal.next,"TXNAME"])
                not.int <- !is.element(tx.gene.start[,"TXNAME"],have.int[[1]])
                not.int <- list(tx.gene.start[not.int,"TXNAME"])
                }
            else if (z >= Nmax){
                tx.gene.start <- tx.gene.start[(cal.start >= tx.gene.start[,"TXSTART"] & cal.start <= tx.gene.start[,"TXEND"]),]
                grep.intron <- grep(cal.start,intron.locus)
                if (length(grep.intron) != 0){
                    nextstart <- grep.intron[names(intron.locus[grep.intron]) == names(intron.locus[grep.intron-1])]
                    if (length(nextstart) != 0){
                        next.intron <- strsplit(intron.locus[as.integer(nextstart)-1],"-")
                        nextstart <- unique(do.call(rbind,next.intron)[,2])
                        }
                    else {nextstart <- "not"}
                    }
                cal.intron <- names(intron.locus[grep.intron])
                next.intron <- names(intron.locus[grep(nextstart[1],intron.locus)])
                cal.next <- unique(c(cal.intron,next.intron))
                have.int <- list(tx.gene.start[cal.next,"TXNAME"])
                not.int <- !is.element(tx.gene.start[,"TXNAME"],have.int[[1]])
                not.int <- list(tx.gene.start[not.int,"TXNAME"])
                }
            in.exon <- which(as.integer(exinfo[,"start"])>cal.start & as.integer(exinfo[,"end"])<cal.end)
            in.exon <- in.exon[in.exon != 1 & in.exon != length(exon.locus)]
            if (length(in.exon) != 0){
                eval1 <- names(exon.locus[in.exon])!=names(exon.locus[in.exon-1])
                eval2 <- names(exon.locus[in.exon])!=names(exon.locus[in.exon+1])
                eval3 <- length(in.exon)==1 & in.exon==1
                eval4 <- length(in.exon)==1 & in.exon==length(exon.locus)
                rm.in.exon <- which((eval1 | eval2 | eval3 | eval4) == "TRUE")
                in.exon <- in.exon[-rm.in.exon]
                }
            in.trans <- unique(names(exon.locus[in.exon]))
            over.trans.int <- is.element(tx.gene.start[in.trans,"TXNAME"],not.int[[1]])
            over.in.trans <- tx.gene.start[in.trans,][over.trans.int,"TXID"]
            over.in.exon.ranges <- matrix(unique(exinfo[in.exon,]),ncol=2)
            colnames(over.in.exon.ranges) <- c("start","end")
            splice.type <- "SE"
            if (length(over.in.exon.ranges[,"start"]) !=0 & length(na.omit(have.int[[1]])) != 0 & length(na.omit(not.int[[1]])) !=0){
                sum.over.in.start <- paste(as.character(as.integer(over.in.exon.ranges[,"start"])-1),collapse="|")
                sum.over.in.end <- paste(as.character(as.integer(over.in.exon.ranges[,"end"])+1),collapse="|")
                test.intron <- unique(intron.locus[grep(sum.over.in.start,intron.locus)])
                test.intron <- unique(c(test.intron,intron.locus[grep(sum.over.in.end,intron.locus)]))
                if (length(test.intron) != 0){
                    test.intron <- strsplit(test.intron,"-")
                    test.intron <- do.call(rbind,test.intron)
                    over.tx.id <- tx.gene.start[is.element(tx.gene.start[,"TXNAME"],have.int[[1]]),"TXID"]
                    over.tx.id <- as.character(unique(over.tx.id))
                    test.exon <- exon.locus[is.element(names(exon.locus),over.tx.id)]
                    test.exon <- strsplit(test.exon,"-")
                    test.exon <- do.call(rbind,test.exon)
                    test.result <- NULL
                    for (v in 1:length(test.intron[,1])){
                        high.ex.start <- as.integer(test.intron[v,1])<as.integer(test.exon[,1])
                        low.ex.end <- as.integer(test.intron[v,2])>as.integer(test.exon[,2])
                        test.ex.name <- names(exon.locus[high.ex.start & low.ex.end])
                        test.result <- c(test.result,test.ex.name)
                        }
                    test.tx <- tx.gene.start[unique(test.result),"TXNAME"]
                    ex.int.tx <- is.element(test.tx,have.int[[1]])
                    have.test.int <- nrow(test.intron) > 1
                    if (length(grep("TRUE",(ex.int.tx & have.test.int))) != 0){
                        splice.type <- "MXE"
                        }
                    else {
                        splice.type <- "SE"
                        }
                    }
                }
            }
        else if(total.locus["type",z] == "tri"){
            cal.start.int <- names(intron.locus[grep(cal.start,intron.locus)])
            cal.end.int <- names(intron.locus[grep(cal.end,intron.locus)])
            have.int <- tx.gene.start[cal.start.int,"TXNAME"]
            not.int <- tx.gene.start[!is.element(tx.gene.start[,"TXNAME"],have.int),"TXNAME"]
            have.int <- list(have.int,tx.gene.start[cal.end.int,"TXNAME"])
            not.int <- list(not.int,tx.gene.start[!is.element(tx.gene.start[,"TXNAME"],have.int[[2]]),"TXNAME"])
            }
        else if(total.locus["type",z] != "tri" & total.locus["type",z]!="AS"){
            total.int <- names(which(total.locus["bigrange",z] ==intron.locus))
            have.int <- list(tx.gene.start[total.int,"TXNAME"])
            not.int <- list(tx.gene.start[!is.element(tx.gene.start[,"TXNAME"],have.int[[1]]),"TXNAME"])
            }
        for (y in 1:length(have.int)){
            have.exp <- na.omit(expdata[have.int[[y]],])
            not.exp <- na.omit(expdata[not.int[[y]],])
            sum.have.exp <- apply(have.exp,2,sum)
            sum.not.exp <- apply(not.exp,2,sum)
            ratio <- sum.have.exp/(sum.have.exp+sum.not.exp)
            notratio <- sum.not.exp/(sum.have.exp+sum.not.exp)*100
            haveratio <- sum.have.exp/(sum.have.exp+sum.not.exp)*100
            realNA <- names(haveratio[haveratio!="NaN" & notratio!="NaN"])
            ratio <- ratio[realNA]
            notratio <- notratio[realNA]
            haveratio <- haveratio[realNA]
            resnp <- samp.snp[,realNA]
            ratio <- ratio[colnames(resnp)]
            notratio <- notratio[colnames(resnp)]
            haveratio <- haveratio[colnames(resnp)]
            for (x in 1:nrow(samp.snp)){
                samp.form <- as.character(as.matrix(samp.snp[x,realNA]))
                if (length(unique(samp.form))>1 & nrow(have.exp) >0 & nrow(not.exp)>0){
                    genoform <- as.matrix(resnp[x,])
                    if (method == "boxplot"){
                        return (t(rbind(ratio,genoform)))
                        }
                    lm.geno <- unique(unlist(strsplit(rownames(table(genoform)),"")))
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
                    if(length(table(matrix.geno))<2){next}
                    colnames(matrix.geno) <- colnames(pregeno)
                    rownames(matrix.geno) <- rownames(pregeno)
                    lm.genoform <- matrix.geno
                    auo.matrix <- data.frame(ratio=ratio,group=factor(lm.genoform))
                    split.g <- split(auo.matrix$ratio,auo.matrix$group)
                    mamean <- sapply(split.g,median)
                    lm.auo.pvalue <- "NaN"
                    glm.auo.pvalue <- "NaN"
                    if (method == "lm" | method == "both"){
                        auo <- summary(lm(formula = ratio ~ as.integer(group), data = auo.matrix))
                        if (length(auo$fstatistic[1])>0 &length(auo$fstatistic[2])>0 & length(auo$fstatistic[3])>0){
                            lm.auo.pvalue <- pf(auo$fstatistic[1],auo$fstatistic[2],auo$fstatistic[3],lower.tail=FALSE)
                            }
                        }
                    if (method == "glm" | method == "both"){
                        calmatrix <- cbind(round(haveratio),round(notratio))
                        obs <- c(1:length(genoform))
                        test1 <- try(glmer(calmatrix ~ as.character(as.matrix(lm.genoform)) +(1|obs),na.action = na.exclude,family=binomial),silent=TRUE)
                        test2 <- try(glmer(calmatrix ~ 1 +(1|obs),na.action = na.exclude,family=binomial),silent=TRUE)
                        if (!( inherits(test1,"try-error") | inherits(test2,"try-error"))){
                               glm.auo.pvalue <- anova(test1,test2)$"Pr(>Chisq)"[2]
                               }
                        }
                    if ((lm.auo.pvalue !="NaN" & lm.auo.pvalue !="NA") | (glm.auo.pvalue !="NaN" & glm.auo.pvalue !="NA")){
                        perP <- "no sig"
                        if(length(mamean) == 2){
                            if(abs(mamean[1]-mamean[2])>0.1){
                                perP <- "sig"
                                }
                            }
                        else if (length(mamean) == 3){
                            if((mamean[1]-mamean[2]>0.1 & mamean[2]-mamean[3]>0.1) | (mamean[1]-mamean[2]< -1*0.1 & mamean[2]-mamean[3]< -1*0.1)){
                                perP <- "sig"
                                }
                            }
                        intron.range <- as.integer(strsplit(total.locus["bigrange",z],"-")[[1]])
                        junc1 <- sort(intron.range[2])+1
                        junc2 <- sort(intron.range[1])-1
                        if (total.locus["type",z]!="tri"){
                            intron.range <- total.locus["bigrange",z]
                            if (total.locus["type",z]=="A3SS"){
                                tar.ex <- paste(ex.start[ex.start == junc1],ex.end[ex.start == junc1],sep="-")[1]
                                if (info.strand == "+"){splice.type <- "A3SS"}
                                else {splice.type <- "A5SS"}
                                    }
                            else if (total.locus["type",z]=="A5SS"){
                                tar.ex <- paste(ex.start[ex.end == junc2],ex.end[ex.end == junc2],sep="-")[1]
                                if (info.strand == "+"){splice.type <- "A5SS"}
                                else {splice.type <- "A3SS"}
                                }
                            else if (total.locus["type",z]=="NI"){
                                tar.ex <- paste(ex.start[ex.end == junc2],ex.end[ex.start == junc1][1],sep="-")[1]
                                splice.type <- "IR"
                                }
                            else if (total.locus["type",z]=="AS"){
                                intron.range <- max.intron
                                over.tx <- which(tx.gene.start[,"TXNAME"]==have.int[[y]][1])
                                ASexon <- exon.locus[which(names(exon.locus)==tx.gene.start[over.tx,"TXID"])]
                                ASexon <- sort(ASexon)
                                if (z < Nmax){
                                    tar.ex <- ASexon[as.integer(grep(junc1,ASexon))][1]
                                    }
                                else if (z >= Nmax){
                                    tar.ex <- ASexon[as.integer(grep(junc2,ASexon))][1]
                                    }
                                }
                            if (lm.auo.pvalue !="NaN" & lm.auo.pvalue !="NA"){
                                each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,intron.range,splice.type,lm.auo.pvalue,perP,tx.gene[1,"GENEID"],"lm")
                                result <- rbind(result,each.result)
                                }
                            if (glm.auo.pvalue !="NaN" & glm.auo.pvalue !="NA"){
                                each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,intron.range,splice.type,glm.auo.pvalue,perP,tx.gene[1,"GENEID"],"glm")
                                result <- rbind(result,each.result)
                                }
                            }
                        else if (total.locus["type",z]=="tri"){
                            if(y==1){
                                tar.ex <- paste(ex.start[ex.start == junc1],ex.end[ex.start == junc1],sep="-")[1]
                                if (info.strand == "+"){
                                    splice.type <- "A3SS"
                                    }
                                else {
                                    splice.type <- "A5SS"
                                    }
                                if (lm.auo.pvalue !="NaN" & lm.auo.pvalue !="NA"){
                                    each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,total.locus["bigrange",z],splice.type,lm.auo.pvalue,perP,tx.gene[1,"GENEID"],"lm")
                                    result <- rbind(result,each.result)
                                    }
                                if (glm.auo.pvalue !="NaN" & glm.auo.pvalue !="NA"){
                                    each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,total.locus["bigrange",z],splice.type,glm.auo.pvalue,perP,tx.gene[1,"GENEID"],"glm")
                                    result <- rbind(result,each.result)
                                    }
                                }
                            else if(y==2){
                                tar.ex <- paste(ex.start[ex.end == junc2],ex.end[ex.end == junc2],sep="-")[1]
                                if (info.strand == "+"){
                                    splice.type <- "A5SS"
                                    }
                                else{
                                    splice.type <- "A3SS"
                                    }
                                if (lm.auo.pvalue !="NaN" & lm.auo.pvalue !="NA"){
                                    each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,total.locus["bigrange",z],splice.type,lm.auo.pvalue,perP,tx.gene[1,"GENEID"],"lm")
                                    result <- rbind(result,each.result)
                                    }
                                if (glm.auo.pvalue !="NaN" & glm.auo.pvalue !="NA"){
                                    each.result <- c(rownames(samp.snp)[x],chrnum,tar.ex,total.locus["bigrange",z],splice.type,glm.auo.pvalue,perP,tx.gene[1,"GENEID"],"glm")
                                    result <- rbind(result,each.result)
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    return (result)
    }
