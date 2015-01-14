sqtlfinder <- function(altInvalue=NULL,overapvalue=NULL,expdata=NULL,snpdata=NULL,method=NULL){
    if (length(overapvalue) == 0){
        return (NULL)
        }
    exonInfo <- altInvalue["exonRange"]
    intronInfo <- altInvalue["intronRange"]
    tx.gene <- altInvalue[["tableBygene"]]
    altintron <- altInvalue[["alterIntron"]]
    rownames(tx.gene) <- tx.gene[,"TXID"]
    predictedSQTL <- NULL
    exon.start <- unlist(start(exonInfo[[1]]))
    exon.end <- unlist(end(exonInfo[[1]]))
    exon.locus <- paste(exon.start,exon.end,sep="-")
    intron.start <- unlist(start(intronInfo[[1]]))
    intron.end <- unlist(end(intronInfo[[1]]))
    intron.locus <- paste(intron.start,intron.end,sep="-")
    names(exon.locus) <- names(exon.end)
    names(intron.locus) <- names(intron.end)
    strandinfo <- as.character(altInvalue[[1]]@strand@values)
    Nchr <- as.character(altInvalue[[1]]@seqnames@values)
    start <- as.character(start(altintron))
    end <- as.character(end(altintron))
    table.start <- table(start)
    table.end <- table(end)
    for (i in 1:length(table.start)){
        total.locus <- NULL
        num.start.locus <- table.start[i]
        start.locus <- names(num.start.locus)
        end.locus <- end[which(start == names(table.start[i]))]
        num.end.locus <- table.end[end.locus]
        if (num.start.locus==1 & length(end.locus)==1){
            if (num.end.locus==1){
                total.locus <- rbind(paste(start.locus,end.locus,sep="-"),"NI")
                NI.range <- as.integer(strsplit(total.locus,"-")[[1]])
                NIjunc1 <- sort(NI.range)[2]+1
                NIjunc2 <- sort(NI.range)[1]-1
                junc1.int <- exon.locus[grep(NIjunc1,exon.locus)]
                junc2.int <- exon.locus[grep(NIjunc2,exon.locus)]
                junc1.sp <- strsplit(junc1.int,"-")
                junc2.sp <- strsplit(junc2.int,"-")
                secondend <- do.call(rbind,junc1.sp)
                firststart <- do.call(rbind,junc2.sp)
                over.start.int <- names(intron.locus[grep(start.locus,intron.locus)])
                over.end.int <- names(intron.locus[grep(end.locus,intron.locus)])
                havename <- unique(c(over.start.int,over.end.int))
                ex.ranges <- unique(names(exon.locus))
                not.ex.name <- ex.ranges[!is.element(ex.ranges,havename)]
                overnotname <- exon.locus[is.element(names(exon.locus),as.character(not.ex.name))]
                se.over.name <- grep(paste(secondend,collapse="|"),overnotname)
                fi.over.name <- grep(paste(firststart,collapse="|"),overnotname)
                NInum <- unique(c(se.over.name,fi.over.name))
                rownames(total.locus) <- c("bigrange","type")
                if (length(NInum) !=0){
                    sqtl.result <- calSignificant(tx.gene,t(unique(t(total.locus))),exon.locus,intron.locus,strandinfo,overapvalue,Nchr,expdata,snpdata,method)
                    predictedSQTL <- rbind(predictedSQTL,sqtl.result)
                    }
                }
            else if (num.end.locus>1){
                multstart <- start[which(end==end.locus)]
                start.locus <- sort(multstart)
                if (length(which(table.start[multstart]!=1)) == 0){
                    total.locus <- lapply(start.locus,function(j){
                        target.ranges <- paste(j,end.locus,sep="-")
                        rbind(target.ranges,"A5SS")
                        })
                    total.locus <- do.call(cbind,total.locus)
                    rownames(total.locus) <- c("bigrange","type")
                    sqtl.result <- calSignificant(tx.gene,t(unique(t(total.locus))),exon.locus,intron.locus,strandinfo,overapvalue,Nchr,expdata,snpdata,method)
                    predictedSQTL <- rbind(predictedSQTL,sqtl.result)
                    }
                }
            }
        else if (num.start.locus>1 & length(end.locus)>1){
            if (!is.element("TRUE",is.element(c(2:10),num.end.locus))){
                end.locus <- sort(end.locus)
                if (length(which(table.end[end.locus]!=1)) == 0){
                    total.locus <- lapply(end.locus,function(j){
                        target.ranges <- paste(start.locus,j,sep="-")
                        rbind(target.ranges,"A3SS")
                        })
                    total.locus <- do.call(cbind,total.locus)
                    rownames(total.locus) <- c("bigrange","type")
                    sqtl.result <- calSignificant(tx.gene,t(unique(t(total.locus))),exon.locus,intron.locus,strandinfo,overapvalue,Nchr,expdata,snpdata,method)
                    predictedSQTL <- rbind(predictedSQTL,sqtl.result)
                    }
                }
            else if (is.element("TRUE",is.element(c(2:10),table.end[end.locus]))){
                left.start.locus <- start[is.element(start,start.locus)]
                right.start.locus <- end[is.element(start,start.locus)]
                left.end.locus <- start[is.element(end,end.locus)]
                right.end.locus <- end[is.element(end,end.locus)]
                start.intron <- paste(left.start.locus,right.start.locus,sep="-")
                end.intron <- paste(left.end.locus,right.end.locus,sep="-")
                total.intron <- unique(append(start.intron,end.intron))
                trans.num <- intron.locus[is.element(intron.locus,total.intron)]
                pre.sort.trans <- sort(trans.num)
                table.trans.num <- table(trans.num)
                sort.trans.num <- table.trans.num[sort(rownames(table.trans.num))]
                sp.trans <- strsplit(rownames(sort.trans.num),"-")
                mat.trans <- do.call(rbind,sp.trans)
                colnames(mat.trans) <- c("start","end")
                start.trans <- mat.trans[which(mat.trans[,1]==start.locus),]
                dis <- as.integer(mat.trans[,"end"])-as.integer(mat.trans[,"start"])
                Nmax <- which(dis==max(dis))
                num1 <- NULL
                num2 <- NULL
                num1 <- sapply(c(1:as.integer(Nmax-1)),function(d){
                    over.trans <- is.element(trans.num,names(sort.trans.num)[d])
                    names(trans.num[over.trans])
                    })
                num2 <- sapply(c(as.integer(Nmax+1):length(sort.trans.num)),function(d){
                    over.trans <- is.element(trans.num,names(sort.trans.num)[d])
                    names(trans.num[over.trans])
                    })
                num1 <- unlist(num1)
                num2 <- unlist(num2)
                kind.num <- unique(append(is.element(num1,num2),is.element(num2,num1)))
                if(length(grep("TRUE",kind.num))==1){
                    low.Nmax <- which(pre.sort.trans==names(sort.trans.num[Nmax]))
                    high.Nmax <- which(pre.sort.trans==names(sort.trans.num[Nmax]))
                    num1.locus <- pre.sort.trans[1:as.integer(min(low.Nmax)-1)]
                    num2.locus <- pre.sort.trans[as.integer(max(high.Nmax)+1):length(pre.sort.trans)]
                    not.over.num1 <- num1.locus[!is.element(num1,num2)]
                    not.over.num2 <- num2.locus[!is.element(num2,num1)]
                    nottotal <- append(table(not.over.num1),table(not.over.num2))
                    tri.name <- names(sort.trans.num)
                    bigrange <- names(sort.trans.num[Nmax])
                    for (d in 1:length(nottotal)){
                        sort.trans.name <- names(sort.trans.num)
                        nottotal.name <- names(nottotal)[d]
                        one.trans <- which(sort.trans.name==nottotal.name)
                        sort.trans.num[one.trans] <- sort.trans.num[nottotal.name]-nottotal[d]
                        }
                    low.Nmax.total <- lapply(c(1:as.integer(Nmax-1)),function(d){
                        target.range <- rbind(sort.trans.name[d],"AS")
                        target.range
                        })
                    high.Nmax.total <- lapply(c(as.integer(Nmax+1):length(sort.trans.name)),function(d){
                        target.range <- rbind(sort.trans.name[d],"AS")
                        target.range
                        })
                    total.locus <- do.call(cbind,c(low.Nmax.total,high.Nmax.total))
                    zero.trans.num <- NULL
                    if (length(which(sort.trans.num==0)) != 0){
                        zero.trans.num <- which(sort.trans.num==0)
                        if (length(which(Nmax > zero.trans.num)) > 0){
                            zero.trans.name <- names(sort.trans.num[zero.trans.num])
                            Nmax <- Nmax-length(unique(zero.trans.name))
                            }
                        }
                    total.locus <- cbind(total.locus,rbind(bigrange,"AS"))
                    total.locus <- cbind(total.locus,rbind(Nmax,"AS"))
                    total.locus <- total.locus[,!is.element(total.locus["bigrange",],unique(zero.trans.num))]
                    rownames(total.locus) <- c("bigrange","type")
                    sqtl.result <- calSignificant(tx.gene,t(unique(t(total.locus))),exon.locus,intron.locus,strandinfo,overapvalue,Nchr,expdata,snpdata,method)
                    predictedSQTL <- rbind(predictedSQTL,sqtl.result)
                    }
                else if(length(grep("FALSE",kind.num))==1){
                    tri.name <- names(sort.trans.num)
                    sp.tri.name <- strsplit(tri.name,"-")
                    tri.mat <- do.call(rbind,sp.tri.name)
                    colnames(tri.mat) <- c("start","end")
                    tri1 <- tri.mat[,"start"]
                    tri2 <- tri.mat[,"end"]
                    tri1.table <- table(tri1)
                    max.tri1 <- max(tri1.table)
                    max.tri1.name <- names(tri1.table[which(tri1.table==max.tri1)])
                    pretotal <- tri.name[grep(paste(max.tri1.name,collapse="|"),tri.name)]
                    tri2.table <- table(tri2)
                    max.tri2 <- max(tri2.table)
                    max.tri2.name <- names(tri2.table[which(tri2.table==max.tri2)])
                    total.locus <- rbind(pretotal[grep(paste(max.tri2.name,collapse="|"),pretotal)],"tri")
                    rownames(total.locus) <- c("bigrange","type")
                    sqtl.result <- calSignificant(tx.gene,total.locus,exon.locus,intron.locus,strandinfo,overapvalue,Nchr,expdata,snpdata,method)
                    predictedSQTL <- rbind(predictedSQTL,sqtl.result)
                    }
                }
            }
        }
    if (method == "boxplot" & length(predictedSQTL) != 0){
        return (predictedSQTL)
        }
    if (method != "boxplot" & length(predictedSQTL) != 0){
        predictedSQTL <- cbind(predictedSQTL,strandinfo)
        colnames(predictedSQTL) <- c("SNP","CHR","targetExon","intron of SNP","type","P.value","per.P.value","gene","method","strand")
        return (unique(predictedSQTL))
        }

    }
