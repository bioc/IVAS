Splicingfinder <- function(GTFdb=NULL,txTable=NULL,calGene=NULL,Ncor=1,out.dir=NULL){
    NItest <- function(exon.info,intron.info,alt.intron.info,tx.gene=NULL){
        NI.exons <- NULL
        alt.intron.info <- alt.intron.info[1,]
        exon.mat <- exon.info
        intron.mat <- intron.info
        NI.range <- alt.intron.info
        NIjunc1 <- sort(NI.range)[2]+1
        NIjunc2 <- sort(NI.range)[1]-1
        NI.exons <- rbind(unique(exon.info[exon.info[,"end"] > NIjunc1 & exon.info[,"start"] < NIjunc2,]))
        if (length(NI.exons) == 0) return (NULL)
        colnames(NI.exons) <- c("start","end")
        NI.intron <- do.call(rbind,lapply(paste(NI.exons[,"start"],NI.exons[,"end"],sep="-"),function(p.NI){
            ni.st.end <- as.integer(strsplit(p.NI,"-")[[1]])
            intron.info[intron.info[,"start"]-1 > ni.st.end[1] & intron.info[,"end"] < ni.st.end[2],]
        }))
        colnames(NI.intron) <- c("start","end")
        new.NIjunc1 <- unique(NI.intron[,2]) + 1
        new.NIjunc2 <- unique(NI.intron[,1]) - 1
        colnames(NI.exons) <- c("start","end")
        targetEX <- paste(NI.exons[,"start"],NI.exons[,"end"],sep="-")
        downEX <- rbind(exon.mat[is.element(exon.mat[,"end"],new.NIjunc2),])
        upEX <- rbind(exon.mat[is.element(exon.mat[,"start"],new.NIjunc1),])
        downEX <- rbind(downEX[is.element(rownames(downEX),intersect(rownames(downEX),rownames(upEX))),])
        upEX <- rbind(upEX[is.element(rownames(upEX),intersect(rownames(downEX),rownames(upEX))),])
        if (nrow(downEX) == 0 | nrow(upEX) == 0)    return (NULL)
        p.downEX <- paste(downEX[,"start"],downEX[,"end"],sep="-")
        p.upEX <- paste(upEX[,"start"],upEX[,"end"],sep="-")
        NI.result <- NULL
        for (i in 1:length(p.upEX)){
            for (j in 1:length(p.downEX)){
                tested.tx <- tx.gene[tx.gene[tx.gene[,"TXID"] == rownames(downEX)[j],"TXSTART"] <= downEX[j,"start"] & tx.gene[tx.gene[,"TXID"] == rownames(downEX)[j],"TXEND"] >= upEX[i,"end"],"TXID"]
                tested.down.num <- rownames(downEX) == rownames(downEX)[j]
                tested.down.ex <- rbind(downEX[tested.down.num,])
                rownames(tested.down.ex) <- rownames(downEX)[tested.down.num]
                tested.down.ex <- rbind(tested.down.ex[is.element(rownames(tested.down.ex),tested.tx),])
                if (nrow(tested.down.ex) == 0)    next
                tested.down.ex <- rbind(tested.down.ex[tested.down.ex[,"end"] < upEX[i,"start"],])
                if (nrow(tested.down.ex) == 0) next
                tested.down.ex <- rbind(tested.down.ex[order(tested.down.ex[,1],decreasing=TRUE)[1],])
                if (upEX[i,"start"] <= tested.down.ex[,"end"]) next
                up.down.EX <- paste(tested.down.ex[,"start"],upEX[i,"end"],sep="-")
                doex.des <- unique(rbind(exon.info[exon.info[,"end"] == tested.down.ex[,"end"],]))
                upex.des <- unique(rbind(exon.info[exon.info[,"start"] == upEX[i,"start"],])) # 넘어 서는거 확인
                doex.des <- paste(sort(paste(doex.des[,"start"],doex.des[,"end"],sep="-")),collapse=",")
                upex.des <- paste(sort(paste(upex.des[,"start"],upex.des[,"end"],sep="-")),collapse=",")
                NIex.des <- paste(sort(unique(c(up.down.EX,targetEX))),collapse=",")
                NI.result <- rbind(NI.result,cbind(rbind(up.down.EX),paste(tested.down.ex,collapse="-"),p.upEX[i],NIex.des,doex.des,upex.des,"IR","possible"))
            }
        }
        if (length(NI.result) == 0) return (NULL)
        colnames(NI.result) <- c("RetainEX","DownEX","UpEX","Retain_des","Do_des","Up_des","Types","status")
        spliced.int <- paste(as.double(do.call(rbind,strsplit(NI.result[,"DownEX"],"-"))[,2])+1,as.double(do.call(rbind,strsplit(NI.result[,"UpEX"],"-"))[,1])-1,sep="-")
        Re.int.test <- is.element(spliced.int,paste(intron.mat[,"start"],intron.mat[,"end"],sep="-"))
        Re.ex.test <- 1:nrow(NI.result) == grep(paste(paste(exon.mat[,"start"],exon.mat[,"end"],sep="-"),collapse="|"),NI.result[,"Retain_des"])
        NI.result[Re.int.test & Re.ex.test,"status"] <- "exist"
        rownames(NI.result) <- 1:nrow(NI.result)
        return (unique(NI.result))
    }
    EStest <- function(exon.info,intron.info,alt.intron.info){
        tx.exon.ranges <- cbind(tapply(exon.info[,"start"],rownames(exon.info),min),tapply(exon.info[,"start"],rownames(exon.info),max))
        colnames(tx.exon.ranges) <- c("start","end")
        alt.intron.info <- alt.intron.info[1,]
        fi.exon.num <- grep(alt.intron.info["start"]-1,exon.info[,"end"])
        se.exon.num <- grep(alt.intron.info["end"]+1,exon.info[,"start"])
        SE.EX <- NULL
        for (i in 1:length(fi.exon.num)){
            fi.exon <- exon.info[fi.exon.num[i],]
            if (length(which(fi.exon["start"] > intron.info[,"start"]-1 & fi.exon["end"] < intron.info[,"end"]-1)) != 0){
                SE.EX <- rbind(SE.EX,fi.exon)
            }
        }
        for (i in 1:length(se.exon.num)){
            se.exon <- exon.info[se.exon.num[i],]
            if (length(which(se.exon["start"] > intron.info[,"start"]-1 & se.exon["end"] < intron.info[,"end"]-1)) != 0){
                SE.EX <- rbind(SE.EX,se.exon)
            }
        }
        if(length(SE.EX) == 0) return (NULL)
        SE.EX <- rbind(SE.EX[SE.EX[,"start"] != SE.EX[,"end"],])
        if(length(SE.EX) == 0) return (NULL)
        SE.EX <- unique(SE.EX)
        test.tx <- unique(unlist(lapply(paste(SE.EX[,"start"],SE.EX[,"end"],sep="-"),function(tx.cal){
            each.tx.cal <- as.double(strsplit(tx.cal,"-")[[1]])
            each.test.tx <- rownames(tx.exon.ranges)[tx.exon.ranges[,"start"] < each.tx.cal[1] & tx.exon.ranges[,"end"] > each.tx.cal[2]]
        })))
        exon.info <- rbind(exon.info[is.element(rownames(exon.info),test.tx),])
        if (nrow(exon.info) == 1) rownames(exon.info) <- test.tx
        intron.info <- rbind(intron.info[is.element(rownames(intron.info),test.tx),])
        if (nrow(intron.info) == 1) rownames(intron.info) <- test.tx
        rownames(SE.EX) <- 1:nrow(SE.EX)
        ES.int.num <- unique(unlist(lapply(paste(SE.EX[,"start"],SE.EX[,"end"],sep="-"),function(tx.cal){
            each.tx.se <- as.double(strsplit(tx.cal,"-")[[1]])
            int.num <- which(intron.info[,"start"]+1 < each.tx.se[1] & intron.info[,"end"]-1 > each.tx.se[2])
            int.num
        })))
        if (length(ES.int.num) == 0) return(NULL)
        ES.int.mat <- cbind(as.double(rownames(intron.info)[ES.int.num]),rbind(intron.info[ES.int.num,]))
        rm.int.num <- NULL
        for (i in 1:nrow(ES.int.mat)){
            num.over.exons <- rownames(exon.info)[exon.info[,"start"] > ES.int.mat[i,"start"]-1 & exon.info[,"end"] < ES.int.mat[i,"end"]+1]
            t.num.over.exons <- table(num.over.exons)
            if (length(which(t.num.over.exons >2)) != 0){
                rm.int.num <- c(rm.int.num,i)
            }
        }
        if(length(rm.int.num) != 0) {ES.int.mat <- rbind(ES.int.mat[-rm.int.num,])}
        if(length(ES.int.mat) == 0) return (NULL)
        colnames(ES.int.mat) <- c("txid","start","end")
        min.int <- min(ES.int.mat[,"start"])
        max.int <- max(ES.int.mat[,"end"])
        final.SE.num <- which(exon.info[,"start"] > min.int-1 & exon.info[,"end"] < max.int+1)
        final.SE.EX <- cbind(as.double(rownames(exon.info)[final.SE.num]),rbind(exon.info[final.SE.num,]))
        colnames(final.SE.EX) <- c("txid","start","end")
        ES.int.num <- lapply(paste(final.SE.EX[,"start"],final.SE.EX[,"end"],sep="-"),function(fse){
            each.fse <- as.double(strsplit(fse,"-")[[1]])
            names(each.fse) <- c("start","end")
            which(intron.info[,"start"]+1 < each.fse["start"] & intron.info[,"end"]-1 > each.fse["end"])
        })
        ES.int.num <- unique(unlist(ES.int.num))
        ES.int.mat <- cbind(as.double(rownames(intron.info)[ES.int.num]),rbind(intron.info[ES.int.num,]))
        colnames(ES.int.mat) <- c("txid","start","end")
        colnames(final.SE.EX) <- c("txid","start","end")
        final.SE.EX <- rbind(final.SE.EX[is.element(final.SE.EX[,"txid"],as.double(test.tx)),])
        t.final.SE.EX <- table(final.SE.EX[,"txid"])
        up2.txid <- as.double(names(t.final.SE.EX)[t.final.SE.EX > 2])
        up1.txid <- as.double(names(t.final.SE.EX)[t.final.SE.EX == 2])
        if (length(ES.int.mat) == 0 | length(final.SE.EX) == 0) return (NULL)
        alt.result <- lapply(paste(ES.int.mat[,1],ES.int.mat[,2],ES.int.mat[,3],sep="-"),function(total.ES.int){
            ES.int <- as.double(strsplit(total.ES.int,"-")[[1]])
            names(ES.int) <- colnames(ES.int.mat)
            alt.ex <- exon.info[rownames(exon.info)==ES.int["txid"] & (exon.info[,"end"] == ES.int["start"]-1 | exon.info[,"start"] == ES.int["end"]+1),]
            alt.ex.result <- cbind(paste(alt.ex[,"start"],alt.ex[,"end"],sep="-"))
        })
        pre.alt.result <- NULL
        if (length(alt.result) != 0){
            for (i in 1:length(alt.result)){
                pre.alt.result <- rbind(pre.alt.result,cbind(alt.result[[i]][1],alt.result[[i]][2]))
            }
        }
        pre.alt.result <- pre.alt.result[!is.element(pre.alt.result[,1],paste(final.SE.EX[,"start"],final.SE.EX[,"end"],sep="-")) & 
                                                                             !is.element(pre.alt.result[,2],paste(final.SE.EX[,"start"],final.SE.EX[,"end"],sep="-")),]
        alt.result <- rbind(unique(pre.alt.result))
        skipped.flank.ex <- alt.result
        up2.es.mat <- NULL
        up2.es.mat.2 <- NULL
        if (length(up1.txid) != 0){
            for (i in 1:length(up1.txid)){
                up1.final.SE.EX <- final.SE.EX[final.SE.EX[,"txid"]==up1.txid[i],]
                min.start <- as.double(min(up1.final.SE.EX[,"start"]))
                max.end <- as.double(max(up1.final.SE.EX[,"end"]))
                le.int <- intron.info[rownames(intron.info) == up1.txid[i] & intron.info[,"end"] == min.start-1,"start"]-1
                ri.int <- intron.info[rownames(intron.info) == up1.txid[i] & intron.info[,"start"] == max.end+1,"end"]+1
                le.ex <- unique(rbind(exon.info[rownames(exon.info) == up1.txid[i] & exon.info[,"end"] == le.int,]))
                ri.ex <- unique(rbind(exon.info[rownames(exon.info) == up1.txid[i] & exon.info[,"start"] == ri.int,]))
                le.ex <- paste(le.ex[,"start"],le.ex[,"end"],sep="-")
                ri.ex <- paste(ri.ex[,"start"],ri.ex[,"end"],sep="-")
                le.ex.des <- paste(sort(unique(paste(exon.info[exon.info[,"end"] == le.int,"start"],exon.info[exon.info[,"end"] == le.int,"end"],sep="-"))),collapse=",")
                ri.ex.des <- paste(sort(unique(paste(exon.info[exon.info[,"start"] == ri.int,"start"],exon.info[exon.info[,"start"] == ri.int,"end"],sep="-"))),collapse=",")
                if (length(le.ex) != 0 & length(ri.ex) != 0){
                    targets.exs <- rbind(c(paste(up1.final.SE.EX[1,-1],collapse="-"),paste(up1.final.SE.EX[2,-1],collapse="-")))
                    up2.es.mat <- rbind(up2.es.mat,cbind(targets.exs,paste(le.ex,collapse=","),paste(ri.ex,collapse=","),targets.exs,le.ex.des,ri.ex.des,"ES"))
                    if (length(alt.result) != 0){
                        up2.es.mat.2 <- cbind(paste(up1.final.SE.EX[1,-1],collapse="-"),paste(up1.final.SE.EX[2,-1],collapse="-"),rbind(alt.result),"ES")
                    }
                }
            }
        }
        inclu.flack.ex <- lapply(paste(final.SE.EX[,"txid"],final.SE.EX[,"start"],final.SE.EX[,"end"],sep="-"),function(total.in.ex){
            in.ex <- as.double(strsplit(total.in.ex,"-")[[1]])
            names(in.ex) <- c("txid","start","end")
            ES.mat.2 <- NULL
            if (length(alt.result) != 0){
                ES.mat.2 <- cbind(paste(in.ex[-1],collapse="-"),"NA",alt.result,"ES")
            }
            ES.mat.2 <- rbind(ES.mat.2,up2.es.mat.2)
            if (length(ES.mat.2) != 0){
                ES.mat.2 <- do.call(rbind,lapply(1:nrow(ES.mat.2),function(each.num){
                    each.ES.mat.2 <- ES.mat.2[each.num,]
                    pre.do.des <- unique(rbind(exon.info[exon.info[,"end"] == unlist(strsplit(each.ES.mat.2[3],"-"))[2],]))
                    pre.up.des <- unique(rbind(exon.info[exon.info[,"start"] == unlist(strsplit(each.ES.mat.2[4],"-"))[1],]))
                    each.ES.mat.2.do.des <- paste(sort(paste(pre.do.des[,"start"],pre.do.des[,"end"],sep="-")),collapse=",")
                    each.ES.mat.2.up.des <- paste(sort(paste(pre.up.des[,"start"],pre.up.des[,"end"],sep="-")),collapse=",")
                    c(each.ES.mat.2[-5],each.ES.mat.2[1:2],each.ES.mat.2.do.des,each.ES.mat.2.up.des,each.ES.mat.2[5])
                }))
            }
            le.intron <- rbind(intron.info[intron.info[,"end"]+1 == in.ex["start"],])
            le.exon.num <- is.element(exon.info[,"end"],le.intron[,"start"]-1)
            le.exon <- cbind(rownames(exon.info)[le.exon.num],rbind(exon.info[le.exon.num,]))
            colnames(le.exon) <- c("txnum","le.start","le.end")
            ri.intron <- rbind(intron.info[intron.info[,"start"]-1 == in.ex["end"],])
            ri.exon.num <- is.element(exon.info[,"start"],ri.intron[,"end"]+1)
            ri.exon <- cbind(rownames(exon.info)[ri.exon.num],rbind(exon.info[ri.exon.num,]))
            colnames(ri.exon) <- c("txnum","ri.start","ri.end")
            inExTx <- in.ex["txid"]
            le.exon <- rbind(le.exon[is.element(le.exon[,"txnum"],inExTx),])
            ri.exon <- rbind(ri.exon[is.element(ri.exon[,"txnum"],inExTx),])
            le.ri.exon <- unique(merge(le.exon,ri.exon,by.x="txnum",by.y="txnum")[,2:5])
            if (nrow(le.ri.exon) != 0){
                ES.result <- do.call(rbind,lapply(1:nrow(le.ri.exon),function(ex.nums){
                    each.le.ri.exon <- rbind(le.ri.exon[ex.nums,])
                    le.des <- unique(rbind(exon.info[exon.info[,"end"] == each.le.ri.exon[,"le.end"],]))
                    ri.des <- unique(rbind(exon.info[exon.info[,"start"] == each.le.ri.exon[,"ri.start"],]))
                    p.le.des <- paste(sort(paste(le.des[,"start"],le.des[,"end"],sep="-")),collapse=",")
                    p.ri.des <- paste(sort(paste(ri.des[,"start"],ri.des[,"end"],sep="-")),collapse=",")
                    cbind(paste(in.ex[-1],collapse="-"),paste(each.le.ri.exon[,"le.start"],each.le.ri.exon[,"le.end"],sep="-"),paste(each.le.ri.exon[,"ri.start"],
                                                                                                                                                                                                                                                     each.le.ri.exon[,"ri.end"],sep="-"),paste(in.ex[-1],collapse="-"),"NA",p.le.des,p.ri.des,"ES")
                }))
                ES.result <- cbind(ES.result[,1],"NA",rbind(ES.result[,-1]))
                ES.result <- rbind(ES.result,up2.es.mat)
                final.ES.result <- NULL
                if (length(ES.mat.2) != 0){
                    for (i in 1:nrow(ES.result)){
                        for (j in 1:nrow(ES.mat.2)){
                            final.ES.result <- rbind(final.ES.result,c(ES.result[i,-c(3:9)],ES.mat.2[j,3],ES.result[i,4],ES.result[i,-c(3:9)],ES.mat.2[j,7],ES.result[i,8],ES.result[i,9]))
                            final.ES.result <- rbind(final.ES.result,c(ES.result[i,-c(4:9)],ES.mat.2[j,4],ES.result[i,5:7],ES.mat.2[j,8],ES.result[i,9]))
                        }
                    }
                }
                ES.result <- unique(rbind(ES.result,final.ES.result,ES.mat.2))
                ES.result
            }
        })
        ES.mat <- do.call(rbind,inclu.flack.ex)
        if (length(ES.mat) == 0) return(NULL)
        ES.mat <- unique(ES.mat)
        sec.mxe.mat <- rbind(ES.mat[ES.mat[,2]!="NA",])
        u.1st.ex <- unique(ES.mat[ES.mat[,2]=="NA",1])
        u.do.ex <- unique(ES.mat[ES.mat[,2]=="NA",3])
        u.up.ex <- unique(ES.mat[ES.mat[,2]=="NA",4])
        final.pre.ES.mat <- NULL
        for (i in 1:length(u.1st.ex)){
            for (j in 1:length(u.do.ex)){
                if (u.1st.ex[i] == u.do.ex[j])    next
                do.des <- unique(ES.mat[ES.mat[,3] == u.do.ex[j],7])
                pre.ES.mat <- do.call(rbind,lapply(u.up.ex,function(each.up.ex){
                    if (u.1st.ex[i] != each.up.ex){
                        up.des <- unique(ES.mat[ES.mat[,4] == each.up.ex,8])
                        cbind(u.1st.ex[i],"NA",u.do.ex[j],each.up.ex,u.1st.ex[i],"NA",do.des,up.des,"ES")
                    }
                }))
                final.pre.ES.mat <- rbind(final.pre.ES.mat,pre.ES.mat)
            }
        }
        ES.mat <- unique(rbind(final.pre.ES.mat,sec.mxe.mat))
        colnames(ES.mat) <- c("TarEX","2ndEX","DownEX","UpEX","Tar_des","2nd_des","Do_des","Up_des","Types")
        final.MXE.cal <- NULL
        pre.MXE.mat <- NULL
        MXE.ES.mat <- NULL
        final.down.ex <- unique(do.call(rbind,strsplit(ES.mat[,"DownEX"],"-")))
        final.up.ex <- unique(do.call(rbind,strsplit(ES.mat[,"UpEX"],"-")))
        colnames(final.down.ex) <- c("start","end")
        colnames(final.up.ex) <- c("start","end")
        u.final.SE.EX <- cbind("1",rbind(unique(final.SE.EX[,c("start","end")])))
        colnames(u.final.SE.EX) <- c("txid","start","end")
        final.MXE.cal <- do.call(rbind,lapply(paste(u.final.SE.EX[,1],u.final.SE.EX[,2],u.final.SE.EX[,3],sep="-"),function(total.es){
            MXE.mat <- NULL
            each.es <- as.double(strsplit(total.es,"-")[[1]])
            names(each.es) <- c("txid","start","end")
            pre.SE.intron <- rbind(ES.mat[is.element(ES.mat[,"TarEX"],paste(each.es["start"],each.es["end"],sep="-")),])
            if (length(pre.SE.intron) != 0){
                SE.intron <- cbind(do.call(rbind,strsplit(pre.SE.intron[,"DownEX"],"-"))[,2],do.call(rbind,strsplit(pre.SE.intron[,"TarEX"],"-"))[,1])
                SE.intron <- rbind(SE.intron,cbind(do.call(rbind,strsplit(pre.SE.intron[,"TarEX"],"-"))[,2],do.call(rbind,strsplit(pre.SE.intron[,"UpEX"],"-"))[,1]))
                u.SE.intron <- rbind(unique(SE.intron))
                MXE.ex.test <- lapply(paste(u.SE.intron[,1],u.SE.intron[,2],sep="-"),function(mxe.intron){
                    mxe.intron <- as.integer(unlist(strsplit(mxe.intron,"-")))
                    do.mxe <- final.down.ex[as.integer(final.down.ex[,"start"]) > mxe.intron[1]-1 & as.integer(final.down.ex[,"end"]) < mxe.intron[2]+1,]
                    up.mxe <- final.up.ex[as.integer(final.up.ex[,"start"]) > mxe.intron[1]-1 & as.integer(final.up.ex[,"end"]) < mxe.intron[2]+1,]
                    rbind(do.mxe,up.mxe)
                })
                MXE.ex.test <- do.call(rbind,MXE.ex.test)
                if (length(MXE.ex.test) != 0){
                    MXE.ex.test <- unique(rbind(MXE.ex.test))
                    MXE.mat <- do.call(rbind,lapply(paste(MXE.ex.test[,1],MXE.ex.test[,2],sep="-"),function(total.other.ex){
                        other.ex <- unlist(strsplit(total.other.ex,"-"))
                        if (as.double(each.es["end"]) < as.double(other.ex[1])){
                            MXE.cal <- cbind(paste(each.es[-1],collapse="-"),paste(other.ex,collapse="-"))
                        }
                        else if (as.double(each.es["start"]) > as.double(other.ex[2])){
                            MXE.cal <- cbind(paste(other.ex,collapse="-"),paste(each.es[-1],collapse="-"))
                        }
                    }))
                }
            }
            MXE.mat
        }))
        if (length(final.MXE.cal) != 0){
            MXE.cal <- final.MXE.cal
            MXE.cal <- unique(MXE.cal)
            if(length(MXE.cal) != 0){
                MXE.ES.mat <- do.call(rbind,lapply(paste(MXE.cal[,1],MXE.cal[,2],sep="-"),function(total.mem){
                    mem <- strsplit(total.mem,"-")[[1]]
                    ES.do <- do.call(rbind,strsplit(ES.mat[,"DownEX"],"-"))
                    ES.up <- do.call(rbind,strsplit(ES.mat[,"UpEX"],"-"))
                    MXE.ES <- unique(rbind(ES.mat[as.double(ES.do[,1]) < min(as.double(mem)) & as.double(ES.up[,2]) > max(as.double(mem)),c("DownEX","UpEX","Do_des","Up_des","Types")]))
                    if (length(MXE.ES) != 0){
                        tar.exs <- matrix(rep(c(paste(mem[1:2],collapse="-"),paste(mem[3:4],collapse="-")),nrow(MXE.ES)),nrow=nrow(MXE.ES),byrow=TRUE)
                        cbind(tar.exs,rbind(MXE.ES[,c("DownEX","UpEX")]),tar.exs,rbind(MXE.ES[,c("Do_des","Up_des","Types")]))
                    }
                }))
                if (length(MXE.ES.mat) != 0){
                    colnames(MXE.ES.mat) <- c("1stEX","2ndEX","DownEX","UpEX","1st_des","2nd_des","Do_des","Up_des","Types")
                    MXE.ES.mat[,"Types"] <- "MXE"
                }
            }
            ES.mat <- rbind(ES.mat,unique(MXE.ES.mat))
        }
        ES.mat <- unique(ES.mat)
        down.ex.e <- as.double(do.call(rbind,strsplit(ES.mat[,"DownEX"],"-"))[,2])
        fi.ex.s.e <- do.call(rbind,strsplit(ES.mat[,"TarEX"],"-"))
        up.ex.e <- as.double(do.call(rbind,strsplit(ES.mat[,"UpEX"],"-"))[,1])
        ES.mat <- rbind(ES.mat[down.ex.e < as.double(fi.ex.s.e[,1]) & as.double(fi.ex.s.e[,2]) < up.ex.e,])
        ES.mat <- cbind(ES.mat,"possible")
        colnames(ES.mat) <- c("1stEX","2ndEX","DownEX","UpEX","1st_des","2nd_des","Do_des","Up_des","Types","status")
        po.test.mat <- paste(do.call(rbind,strsplit(ES.mat[,"DownEX"],"-"))[,2],do.call(rbind,strsplit(ES.mat[,"1stEX"],"-"))[,1],sep="-")
        po.test.mat <- cbind(po.test.mat,paste(do.call(rbind,strsplit(ES.mat[,"1stEX"],"-"))[,2],do.call(rbind,strsplit(ES.mat[,"UpEX"],"-"))[,1],sep="-"))
        po.test.mat <- cbind(po.test.mat,paste(do.call(rbind,strsplit(ES.mat[,"DownEX"],"-"))[,2],do.call(rbind,strsplit(ES.mat[,"UpEX"],"-"))[,1],sep="-"))
        po.test.mat <- cbind(po.test.mat,rbind(ES.mat[,3:4]))
        intron.info.mat <- paste(intron.info[,"start"]-1,intron.info[,"end"]+1,sep="-")
        names(intron.info.mat) <- rownames(intron.info)
        exon.info.mat <- paste(exon.info[,"start"],exon.info[,"end"],sep="-")
        names(exon.info.mat) <- rownames(exon.info)
        po.over.result <- lapply(1:nrow(po.test.mat),function(total.ptm.num){
            ptm <- po.test.mat[total.ptm.num,]
            ptm.test.int <- which(is.element(names(intron.info.mat)[which(intron.info.mat==ptm[1])],names(intron.info.mat)[which(intron.info.mat==ptm[2])])=="TRUE")
            ptm.test.ex <- intersect(names(exon.info.mat)[which(exon.info.mat==ptm[4])],names(exon.info.mat)[which(exon.info.mat==ptm[5])])
            ptm.test.ski.int <- names(intron.info.mat)[which(intron.info.mat==ptm[3])]
            if(length(intersect(ptm.test.ex,ptm.test.ski.int)) != 0 & length(ptm.test.int) != 0){
                total.ptm.num
            }
        })
        po.over.result <- is.element(1:nrow(ES.mat),unlist(po.over.result))
        po.over.result <- po.over.result & (ES.mat[,"2ndEX"] == "NA")
        ES.mat[po.over.result,"status"] <- "exist"
        se.exons.num <- ES.mat[,"2ndEX"]!="NA" & ES.mat[,"Types"]=="ES"# & ES.mat[,"status"]=="exist"
        se.MXE.num <- ES.mat[,"2ndEX"]!="NA" & ES.mat[,"Types"]=="MXE"
        if(length(which(se.exons.num=="TRUE"))!=0){
            es.se.test <- NULL
            es.se.ex <- rbind(ES.mat[se.exons.num,])
            for (i in 1:nrow(es.se.ex)){
                es.se.1 <- which(is.element(names(exon.info.mat)[which(exon.info.mat==es.se.ex[i,1])],names(exon.info.mat)[which(exon.info.mat==es.se.ex[i,2])])=="TRUE")
                es.se.2 <- which(intron.info.mat == paste(unlist(strsplit(es.se.ex[i,3],"-"))[2],unlist(strsplit(es.se.ex[i,4],"-"))[1],sep="-"))
                if (length(es.se.1) != 0 & length(es.se.2) != 0){
                    es.se.test <- c(es.se.test,1==1)
                }
            }
            ES.mat[which(se.exons.num == "TRUE")[es.se.test],"status"] <- "exist"
        }
        if(length(which(se.MXE.num=="TRUE"))!=0){
            mxe.se <- rbind(ES.mat[se.MXE.num,])
            mxe.se.1 <- paste(do.call(rbind,strsplit(mxe.se[,"1stEX"],"-"))[,2],do.call(rbind,strsplit(mxe.se[,"UpEX"],"-"))[,1],sep="-")
            mxe.se.2 <- paste(do.call(rbind,strsplit(mxe.se[,"DownEX"],"-"))[,2],do.call(rbind,strsplit(mxe.se[,"2ndEX"],"-"))[,1],sep="-")
            mxe.se.3 <- paste(do.call(rbind,strsplit(mxe.se[,"DownEX"],"-"))[,2],do.call(rbind,strsplit(mxe.se[,"1stEX"],"-"))[,1],sep="-")
            mxe.se.4 <- paste(do.call(rbind,strsplit(mxe.se[,"2ndEX"],"-"))[,2],do.call(rbind,strsplit(mxe.se[,"UpEX"],"-"))[,1],sep="-")
            mxe.se.test <- is.element(mxe.se.1,intron.info.mat) & is.element(mxe.se.2,intron.info.mat) & is.element(mxe.se.3,intron.info.mat) & is.element(mxe.se.4,intron.info.mat)
            ES.mat[which(se.MXE.num == "TRUE")[mxe.se.test],"status"] <- "exist"
        }
        rownames(ES.mat) <- 1:nrow(ES.mat)
        colnames(ES.mat) <- c("1stEX","2ndEX","DownEX","UpEX","1st_des","2nd_des","Do_des","Up_des","Types","status")
        return (unique(ES.mat))
    }
    ASStest <- function(exon.info,intron.info,alt.intron.info){
        alt.intron.info <- alt.intron.info[1,]
        over.ex <- rbind(exon.info[exon.info[,"end"] > alt.intron.info["end"]+1 & exon.info[,"start"] < alt.intron.info["start"]-1,])
        left.ex.num <- which(exon.info[,"end"] == alt.intron.info["start"]-1)
        right.ex.num <- which(exon.info[,"start"] == alt.intron.info["end"]+1)
        left.ex <- unique(cbind(as.double(rownames(exon.info)[left.ex.num]),rbind(exon.info[left.ex.num,])))
        right.ex <- unique(cbind(as.double(rownames(exon.info)[right.ex.num]),rbind(exon.info[right.ex.num,])))
        colnames(left.ex) <- c("id","start","end")
        colnames(right.ex) <- c("id","start","end")
        left.ex.result <- NULL
        for(i in 1:nrow(left.ex)){
            lx <- intron.info[intron.info[,"start"] < left.ex[i,"start"] & intron.info[,"end"] > left.ex[i,"end"],]
            if(length(lx) != 0){
                left.ex.result <- rbind(left.ex.result,cbind(left.ex[i,"id"],rbind(lx)))
            }
        }
        right.ex.result <- NULL
        for(i in 1:nrow(right.ex)){
            lx <- intron.info[intron.info[,"start"] < right.ex[i,"start"] & intron.info[,"end"] > right.ex[i,"end"],]
            if(length(lx) != 0){
                right.ex.result <- rbind(right.ex.result,cbind(right.ex[i,"id"],rbind(lx)))
            }
        }
        fi.te.num <- exon.info[,"start"] < alt.intron.info["start"]-1 & exon.info[,"end"] > alt.intron.info["start"]-1
        fi.te.ex <- rbind(exon.info[fi.te.num,])
        rownames(fi.te.ex) <- rownames(exon.info)[fi.te.num]
        fi.te.ex.num <- !is.element(paste(fi.te.ex[,"start"],fi.te.ex[,"end"]),paste(over.ex[,"start"],over.ex[,"end"]))
        fi.nm <- rownames(fi.te.ex)[fi.te.ex.num]
        fi.te.ex <- rbind(fi.te.ex[fi.te.ex.num,])
        rownames(fi.te.ex) <- fi.nm
        se.te.num <- intron.info[,"start"] < alt.intron.info["start"] & intron.info[,"end"] > alt.intron.info["start"]
        se.te.int <- rbind(cbind(as.integer(rownames(intron.info)[se.te.num]),rbind(intron.info[se.te.num,])))
        se.te.test <- !is.element(paste(se.te.int[,1],se.te.int[,"start"],se.te.int[,"end"]),paste(left.ex.result[,1],left.ex.result[,"start"],left.ex.result[,"end"]))
        se.te.int <- rbind(se.te.int[se.te.test,])
        th.te.num <- exon.info[,"start"] < alt.intron.info["end"]+1 & exon.info[,"end"] > alt.intron.info["end"]+1
        th.te.ex <- rbind(exon.info[th.te.num,])
        rownames(th.te.ex) <- rownames(exon.info)[th.te.num]
        th.te.ex.num <- !is.element(paste(th.te.ex[,"start"],th.te.ex[,"end"]),paste(over.ex[,"start"],over.ex[,"end"]))
        th.rm <- rownames(th.te.ex)[th.te.ex.num]
        th.te.ex <- rbind(th.te.ex[th.te.ex.num,])
        rownames(th.te.ex) <- th.rm
        fo.te.num <- intron.info[,"start"] < alt.intron.info["end"] & intron.info[,"end"] > alt.intron.info["end"]
        fo.te.int <- rbind(cbind(as.integer(rownames(intron.info)[fo.te.num]),rbind(intron.info[fo.te.num,])))
        fo.te.test <- !is.element(paste(fo.te.int[,1],fo.te.int[,"start"],fo.te.int[,"end"]),paste(right.ex.result[,1],right.ex.result[,"start"],right.ex.result[,"end"]))
        fo.te.int <- rbind(fo.te.int[fo.te.test,])
        ASS.final.result <- NULL
        if(length(fi.te.ex) != 0 | length(se.te.int) != 0){
            tar.tx.num <- which(exon.info[,"end"] == alt.intron.info["start"]-1)
            up.ex.start <- intron.info[is.element(intron.info[,"start"]-1,exon.info[tar.tx.num,"end"]),"end"]+1
            up.tx.exon.num <- which((is.element(rownames(exon.info),names(tar.tx.num)) & is.element(exon.info[,"start"],up.ex.start)) == "TRUE")
            tar.tx <- rbind(exon.info[tar.tx.num,])
            up.tx <- rbind(exon.info[up.tx.exon.num,])
            if (nrow(up.tx) == 1 | nrow(tar.tx) == 1){
                rownames(up.tx) <- rownames(exon.info)[up.tx.exon.num]
                rownames(tar.tx) <- rownames(exon.info)[tar.tx.num]
            }
            each.up.tx <- tapply(up.tx[,"start"],rownames(up.tx),min)
            up.rn <- rownames(up.tx)[is.element(paste(rownames(up.tx),up.tx[,"start"]),paste(names(each.up.tx),each.up.tx))]
            up.tx <- rbind(up.tx[is.element(paste(rownames(up.tx),up.tx[,"start"]),paste(names(each.up.tx),each.up.tx)),])
            rownames(up.tx) <- up.rn
            over.names <- intersect(rownames(tar.tx),rownames(up.tx))
            tar.tx <- rbind(tar.tx[over.names,])
            up.tx <- rbind(up.tx[over.names,])
            tar.ex <- paste(tar.tx[,"start"],tar.tx[,"end"],sep="-")
            tar.up.ex <- paste(up.tx[,"start"],up.tx[,"end"],sep="-")
            ASS.final.result <- NULL
            if (length(fi.te.ex) != 0){
                long.tx <- rbind(fi.te.ex)
                long.int.num <- which(is.element(intron.info[,"start"]-1,long.tx[,"end"])=="TRUE")
                long.up.tx <- rbind(exon.info[is.element(rownames(exon.info),rownames(intron.info)[long.int.num]) & 
                                                                                is.element(exon.info[,"start"],intron.info[long.int.num,"end"]+1),])
                if (nrow(long.up.tx) == 1){
                    rnforup <- rownames(exon.info)[is.element(rownames(exon.info),rownames(intron.info)[long.int.num]) & is.element(exon.info[,"start"],intron.info[long.int.num,"end"]+1)]
                    rownames(long.up.tx) <- rnforup
                }
                each.up.tx <- tapply(long.up.tx[,"start"],rownames(long.up.tx),min)
                long.up.tx <- rbind(long.up.tx[is.element(paste(rownames(long.up.tx),long.up.tx[,"start"]),paste(names(each.up.tx),each.up.tx)),])
                rownames(long.up.tx) <- rownames(long.up.tx)[is.element(paste(rownames(long.up.tx),long.up.tx[,"start"]),paste(names(each.up.tx),each.up.tx))]
                over.names <- intersect(rownames(long.tx),rownames(long.up.tx))
                long.tx <- rbind(long.tx[over.names,])
                long.up.tx <- rbind(long.up.tx[over.names,])
                long.ex <- paste(long.tx[,"start"],long.tx[,"end"],sep="-")
                long.up.ex <- paste(long.up.tx[,"start"],long.up.tx[,"end"],sep="-")
                if (length(long.up.ex) != 0){
                    for (i in 1:length(tar.ex)){
                        for(j in 1:length(long.ex)){
                            s.tar.ex <- unlist(strsplit(tar.ex[i],"-"))
                            s.tar.up.ex <- unlist(strsplit(tar.up.ex[i],"-"))
                            s.long.ex <- unlist(strsplit(long.ex[j],"-"))
                            if (s.tar.ex[2] >= s.long.ex[1] & s.tar.up.ex[1] <= s.long.ex[2]) next
                            s.total.tar.ex <- do.call(rbind,strsplit(tar.ex,"-"))
                            s.total.tar.up.ex <- do.call(rbind,strsplit(tar.up.ex,"-"))
                            
                            Re.rm <- as.double(s.total.tar.up.ex[,1]) <= s.long.ex[2] & as.double(s.total.tar.ex[,2]) >= s.long.ex[1]
                            tarex.nums <- grep(s.tar.ex[2],tar.ex)
                            tarex.des <- paste(sort(unique(tar.ex[is.element(1:length(tar.ex),tarex.nums) & !Re.rm])),collapse=",")
                            tarex.up.des <- paste(sort(unique(tar.up.ex[is.element(1:length(tar.ex),tarex.nums) & !Re.rm])),collapse=",")
                            if (length(which((tarex.nums & !Re.rm) == TRUE)) == 0) next
                            
                            longex.nums <- grep(s.long.ex[2],long.ex)
                            longex.des <- paste(sort(unique(long.ex[longex.nums])),collapse=",")
                            longex.up.des <- paste(sort(unique(long.up.ex[longex.nums])),collapse=",")
                            tar.mat <- cbind(tar.ex[i],long.ex[j],tar.up.ex[i],long.up.ex[j],tarex.des,longex.des,tarex.up.des,longex.up.des,"A5SS")
                            ASS.final.result <- rbind(ASS.final.result,tar.mat)
                        }
                    }
                }
            }
            if (length(se.te.int) != 0){
                se.te.tx <- rownames(intron.info)[is.element(paste(intron.info[,"start"],intron.info[,"end"]),paste(se.te.int[,"start"],se.te.int[,"end"]))]
                long.tx <- rbind(exon.info[is.element(rownames(exon.info),se.te.tx) & is.element(exon.info[,"end"],se.te.int[,"start"]-1),])
                long.up.tx <- rbind(exon.info[is.element(rownames(exon.info),se.te.tx) & is.element(exon.info[,"start"],se.te.int[,"end"]+1),])
                if (nrow(long.tx) == 1 | nrow(long.up.tx) == 1){
                    rownames(long.tx) <- rownames(exon.info)[is.element(rownames(exon.info),se.te.tx) & is.element(exon.info[,"end"],se.te.int[,"start"]-1)]
                    rownames(long.up.tx) <- rownames(exon.info)[is.element(rownames(exon.info),se.te.tx) & is.element(exon.info[,"start"],se.te.int[,"end"]+1)]
                }
                each.up.tx <- tapply(long.up.tx[,"start"],rownames(long.up.tx),min)
                long.up.tx <- rbind(long.up.tx[is.element(paste(rownames(long.up.tx),long.up.tx[,"start"]),paste(names(each.up.tx),each.up.tx)),])
                rownames(long.up.tx) <- rownames(long.up.tx)[is.element(paste(rownames(long.up.tx),long.up.tx[,"start"]),paste(names(each.up.tx),each.up.tx))]
                over.names <- intersect(rownames(long.tx),rownames(long.up.tx))
                long.tx <- rbind(long.tx[over.names,])
                long.up.tx <- rbind(long.up.tx[over.names,])
                long.ex <- paste(long.tx[,"start"],long.tx[,"end"],sep="-")
                long.up.ex <- paste(long.up.tx[,"start"],long.up.tx[,"end"],sep="-")
                if (length(long.up.ex) != 0){
                    for (i in 1:length(tar.ex)){
                        for(j in 1:length(long.ex)){
                            if(as.double(strsplit(long.ex[j],"-")[[1]][2]) > as.double(strsplit(tar.ex[i],"-")[[1]][1])){
                                s.tar.ex <- as.double(unlist(strsplit(tar.ex[i],"-")))
                                s.long.ex <- as.double(unlist(strsplit(long.ex[j],"-")))
                                s.long.up.ex <- unlist(strsplit(long.up.ex[j],"-"))
                                if (s.long.up.ex[1] <= s.tar.ex[2] & s.long.ex[2] >= s.tar.ex[1]) next
                                s.total.long.ex <- do.call(rbind,strsplit(long.ex,"-"))
                                s.total.long.up.ex <- do.call(rbind,strsplit(long.up.ex,"-"))
                                Re.rm <- as.double(s.total.long.ex[,2]) >= s.tar.ex[1] & as.double(s.total.long.up.ex[,1]) <= s.tar.ex[2]
                                
                                tarex.nums <- grep(s.tar.ex[2],tar.ex)
                                tarex.des <- paste(sort(unique(tar.ex[tarex.nums])),collapse=",")
                                tarex.up.des <- paste(sort(unique(tar.up.ex[tarex.nums])),collapse=",")
                                
                                longex.nums <- grep(s.long.ex[2],long.ex)
                                longex.des <- paste(sort(unique(long.ex[is.element(1:length(long.ex),longex.nums) & !Re.rm])),collapse=",")
                                longex.up.des <- paste(sort(unique(long.up.ex[is.element(1:length(long.ex),longex.nums) & !Re.rm])),collapse=",")
                                if (length(which((is.element(1:length(long.ex),longex.nums) & !Re.rm) == TRUE)) == 0) next
                                
                                tar.mat <- cbind(long.ex[j],tar.ex[i],long.up.ex[j],tar.up.ex[i],longex.des,tarex.des,longex.up.des,tarex.up.des,"A5SS")
                                ASS.final.result <- rbind(ASS.final.result,tar.mat)
                            }
                        }
                    }
                }
            }
        }
        if(length(th.te.ex) != 0 | length(fo.te.int) != 0){
            tar.tx.num <- which(exon.info[,"start"] == alt.intron.info["end"]+1)
            down.ex.start <- intron.info[is.element(intron.info[,"end"]+1,exon.info[tar.tx.num,"start"]),"start"]-1
            down.tx.exon.num <- which((is.element(rownames(exon.info),names(tar.tx.num)) & is.element(exon.info[,"end"],down.ex.start)) == "TRUE")
            tar.tx <- rbind(exon.info[tar.tx.num,])
            down.tx <- rbind(exon.info[down.tx.exon.num,])
            if (nrow(down.tx) == 1 | nrow(tar.tx) == 1){
                rownames(down.tx) <- rownames(exon.info)[down.tx.exon.num]
                rownames(tar.tx) <- rownames(exon.info)[tar.tx.num]
            }
            each.down.tx <- tapply(down.tx[,"start"],rownames(down.tx),max)
            down.rn <- rownames(down.tx)[is.element(paste(rownames(down.tx),down.tx[,"start"]),paste(names(each.down.tx),each.down.tx))]
            down.tx <- rbind(down.tx[is.element(paste(rownames(down.tx),down.tx[,"start"]),paste(names(each.down.tx),each.down.tx)),])
            rownames(down.tx) <- down.rn
            over.names <- intersect(rownames(tar.tx),rownames(down.tx))
            tar.tx <- rbind(tar.tx[over.names,])
            down.tx <- rbind(down.tx[over.names,])
            tar.ex <- paste(tar.tx[,"start"],tar.tx[,"end"],sep="-")
            tar.down.ex <- paste(down.tx[,"start"],down.tx[,"end"],sep="-")
            if(length(th.te.ex) != 0){
                long.tx <- rbind(th.te.ex)
                long.int.num <- which(is.element(intron.info[,"end"]+1,long.tx[,"start"])=="TRUE")
                long.down.num <- is.element(rownames(exon.info),rownames(intron.info)[long.int.num]) & is.element(exon.info[,"end"],intron.info[long.int.num,"start"]-1)
                long.down.tx <- rbind(exon.info[long.down.num,])
                rownames(long.down.tx) <- rownames(exon.info)[long.down.num]
                
                each.down.tx <- tapply(long.down.tx[,"start"],rownames(long.down.tx),max)
                long.down.rn <- rownames(long.down.tx)[is.element(paste(rownames(long.down.tx),long.down.tx[,"start"]),paste(names(each.down.tx),each.down.tx))]
                long.down.tx <- rbind(long.down.tx[is.element(paste(rownames(long.down.tx),long.down.tx[,"start"]),paste(names(each.down.tx),each.down.tx)),])
                rownames(long.down.tx) <- long.down.rn
                
                over.names <- intersect(rownames(long.tx),rownames(long.down.tx))
                long.tx <- rbind(long.tx[over.names,])
                long.down.tx <- rbind(long.down.tx[over.names,])
                long.ex <- paste(long.tx[,"start"],long.tx[,"end"],sep="-")
                long.down.ex <- paste(long.down.tx[,"start"],long.down.tx[,"end"],sep="-")
                if (length(long.down.ex) != 0){
                    for (i in 1:length(tar.ex)){
                        for(j in 1:length(long.ex)){
                            s.tar.ex <- as.double(unlist(strsplit(tar.ex[i],"-")))
                            s.long.ex <- as.double(unlist(strsplit(long.ex[j],"-")))
                            s.tar.down.ex <- as.double(unlist(strsplit(tar.down.ex[i],"-")))
                            if (s.tar.down.ex[2] >= s.long.ex[1] & s.tar.ex[1] <= s.long.ex[2]) next
                            s.total.tar.down.ex <- do.call(rbind,strsplit(tar.down.ex,"-"))
                            s.total.tar.ex <- do.call(rbind,strsplit(tar.ex,"-"))
                            Re.rm <- as.double(s.total.tar.down.ex[,2]) >= s.long.ex[1] & as.double(s.total.tar.ex[,1]) <= s.long.ex[2]
                            tarex.nums <- grep(s.tar.ex[1],tar.ex)
                            tarex.des <- paste(sort(unique(tar.ex[is.element(1:length(tar.ex),tarex.nums) & !Re.rm])),collapse=",")
                            tarex.down.des <- paste(sort(unique(tar.down.ex[is.element(1:length(tar.ex),tarex.nums) & !Re.rm])),collapse=",")
                            if (length(which((tarex.nums & !Re.rm) == TRUE)) == 0) next
                            longex.nums <- grep(s.long.ex[1],long.ex)
                            longex.des <- paste(sort(unique(long.ex[longex.nums])),collapse=",")
                            longex.down.des <- paste(sort(unique(long.down.ex[longex.nums])),collapse=",")
                            
                            tar.mat <- cbind(tar.ex[i],long.ex[j],tar.down.ex[i],long.down.ex[j],tarex.des,longex.des,tarex.down.des,longex.down.des,"A3SS")
                            ASS.final.result <- rbind(ASS.final.result,tar.mat)
                        }
                    }
                }
            }
            if(length(fo.te.int) != 0){
                fo.te.tx <- rownames(intron.info)[is.element(paste(intron.info[,"start"],intron.info[,"end"]),paste(fo.te.int[,"start"],fo.te.int[,"end"]))]
                long.tx <- rbind(exon.info[is.element(rownames(exon.info),fo.te.tx) & is.element(exon.info[,"start"],fo.te.int[,"end"]+1),])
                long.down.tx <- rbind(exon.info[is.element(rownames(exon.info),fo.te.tx) & is.element(exon.info[,"end"],fo.te.int[,"start"]-1),])
                if (nrow(long.down.tx) == 1 | nrow(long.tx) == 1){
                    rownames(long.down.tx) <- rownames(exon.info)[is.element(rownames(exon.info),fo.te.tx) & is.element(exon.info[,"end"],fo.te.int[,"start"]-1)]
                    rownames(long.tx) <- rownames(exon.info)[is.element(rownames(exon.info),fo.te.tx) & is.element(exon.info[,"start"],fo.te.int[,"end"]+1)]
                }
                each.down.tx <- tapply(long.down.tx[,"start"],rownames(long.down.tx),max)
                long.down.tx <- rbind(long.down.tx[is.element(paste(rownames(long.down.tx),long.down.tx[,"start"]),paste(names(each.down.tx),each.down.tx)),])
                rownames(long.down.tx) <- rownames(long.down.tx)[is.element(paste(rownames(long.down.tx),long.down.tx[,"start"]),paste(names(each.down.tx),each.down.tx))]
                over.names <- intersect(rownames(long.tx),rownames(long.down.tx))
                long.tx <- rbind(long.tx[over.names,])
                long.down.tx <- rbind(long.down.tx[over.names,])
                long.ex <- paste(long.tx[,"start"],long.tx[,"end"],sep="-")
                long.down.ex <- paste(long.down.tx[,"start"],long.down.tx[,"end"],sep="-")
                if (length(long.down.ex) != 0){
                    for (i in 1:length(tar.ex)){
                        for(j in 1:length(long.ex)){
                            if(as.double(strsplit(long.ex[j],"-")[[1]][1]) < as.double(strsplit(tar.ex[i],"-")[[1]][2])){
                                s.tar.ex <- as.double(unlist(strsplit(tar.ex[i],"-")))
                                tarex.nums <- grep(s.tar.ex[1],tar.ex)
                                tarex.des <- paste(sort(unique(tar.ex[tarex.nums])),collapse=",")
                                tarex.down.des <- paste(sort(unique(tar.down.ex[tarex.nums])),collapse=",")
                                
                                s.long.ex <- as.double(unlist(strsplit(long.ex[j],"-")))
                                s.total.long.ex <- do.call(rbind,strsplit(long.ex,"-"))
                                s.down.long.ex <- as.double(unlist(strsplit(long.down.ex[j],"-")))
                                s.total.down.long.ex <- do.call(rbind,strsplit(long.down.ex,"-"))
                                if (s.down.long.ex[2] >= s.tar.ex[1] & s.long.ex[1] <= s.tar.ex[2]) next
                                Re.rm <- as.double(s.total.down.long.ex[,2]) >= s.tar.ex[1] & as.double(s.total.long.ex[,1]) <= s.tar.ex[2]
                                
                                longex.nums <- grep(s.long.ex[1],long.ex)
                                longex.des <- paste(sort(unique(long.ex[is.element(1:length(long.ex),longex.nums) & !Re.rm])),collapse=",")
                                longex.down.des <- paste(sort(unique(long.down.ex[is.element(1:length(long.ex),longex.nums) & !Re.rm])),collapse=",")
                                if (length(which((is.element(1:length(long.ex),longex.nums) & !Re.rm) == TRUE)) == 0) next
                                tar.mat <- cbind(long.ex[j],tar.ex[i],long.down.ex[j],tar.down.ex[i],longex.des,tarex.des,longex.down.des,tarex.down.des,"A3SS")
                                ASS.final.result <- rbind(ASS.final.result,tar.mat)
                            }
                        }
                    }
                }
            }
        }
        if (length(ASS.final.result) == 0) return (NULL)
        ASS.final.result <- unique(ASS.final.result)
        ASS.final.result <- cbind(ASS.final.result,"possible")
        colnames(ASS.final.result) <- c("ShortEX","LongEX","ShortNeighborEX","LongNeighborEX","Short_des","Long_des",
                                                                        "ShortNeighbor_des","LongNeighbor_des","Types","status")
        A3.num <- which(ASS.final.result[,"Types"] == "A3SS")
        A5.num <- which(ASS.final.result[,"Types"] == "A5SS")
        A3.final.result <- rbind(ASS.final.result[A3.num,])
        A5.final.result <- rbind(ASS.final.result[A5.num,])
        int.test <- paste(intron.info[,"start"]-1,intron.info[,"end"]+1,sep="-")
        if (length(A3.num) != 0){
            A3.long.test.int <- paste(do.call(rbind,strsplit(A3.final.result[,"ShortNeighborEX"],"-"))[,2],do.call(rbind,strsplit(A3.final.result[,"ShortEX"],"-"))[,1],sep="-")
            A3.short.test.int <- paste(do.call(rbind,strsplit(A3.final.result[,"LongNeighborEX"],"-"))[,2],do.call(rbind,strsplit(A3.final.result[,"LongEX"],"-"))[,1],sep="-")
            A3.final.result[is.element(A3.long.test.int,int.test) & is.element(A3.short.test.int,int.test),"status"] <- "exist"
        }
        if (length(A5.num) != 0){
            A5.long.test.int <- paste(do.call(rbind,strsplit(A5.final.result[,"ShortEX"],"-"))[,2],do.call(rbind,strsplit(A5.final.result[,"ShortNeighborEX"],"-"))[,1],sep="-")
            A5.short.test.int <- paste(do.call(rbind,strsplit(A5.final.result[,"LongEX"],"-"))[,2],do.call(rbind,strsplit(A5.final.result[,"LongNeighborEX"],"-"))[,1],sep="-")
            A5.final.result[is.element(A5.long.test.int,int.test) & is.element(A5.short.test.int,int.test),"status"] <- "exist"
        }
        ASS.final.result <- rbind(A5.final.result,A3.final.result)
        ASS.final.result <- rbind(ASS.final.result[!is.na(ASS.final.result[,"ShortNeighborEX"]) & !is.na(ASS.final.result[,"LongNeighborEX"]),])
        return (unique(ASS.final.result))
    }
    CalAlt <- function(altInvalue){
        if (length(altInvalue) == 0)    return (NULL)
        exonInfo <- altInvalue["exonRange"]
        intronInfo <- altInvalue["intronRange"]
        tx.gene <- altInvalue[["tableBygene"]]
        altintron <- altInvalue[["alterIntron"]]
        rownames(tx.gene) <- tx.gene[,"TXID"]
        predictedSQTL <- NULL
        exon.start <- unlist(start(exonInfo[[1]]))
        exon.end <- unlist(end(exonInfo[[1]]))
        exon.mat <- cbind(exon.start,exon.end)
        exon.locus <- paste(exon.start,exon.end,sep="-")
        intron.start <- unlist(start(intronInfo[[1]]))
        intron.end <- unlist(end(intronInfo[[1]]))
        intron.mat <- cbind(intron.start,intron.end)
        intron.locus <- paste(intron.start,intron.end,sep="-")
        alt.mat.int <- cbind(as.integer(start(altintron)),as.integer(end(altintron)))
        names(exon.locus) <- names(exon.end)
        names(intron.locus) <- names(intron.end)
        strandinfo <- as.character(altInvalue[[1]]@strand@values)
        Nchr <- as.character(altInvalue[[1]]@seqnames@values)
        colnames(exon.mat) <- c("start","end")
        colnames(intron.mat) <- c("start","end")
        colnames(alt.mat.int) <- c("start","end")
        uni.tx <- unique(rownames(exon.mat))
        final.exon.mat <- NULL
        for (i in 1:length(uni.tx)){
            each.tx.info <- rbind(exon.mat[rownames(exon.mat)==uni.tx[i],])
            each.tx.mat <- rbind(each.tx.info[order(each.tx.info[,"start"]),])
            if (nrow(each.tx.mat) == 1) rownames(each.tx.mat) <- uni.tx[i]
            final.exon.mat <- rbind(final.exon.mat,each.tx.mat)
        }
        final.intron.mat <- NULL
        for (i in 1:length(uni.tx)){
            each.tx.info <- rbind(intron.mat[rownames(intron.mat)==uni.tx[i],])
            each.tx.mat <- rbind(each.tx.info[order(each.tx.info[,"start"]),])
            if (nrow(each.tx.mat) == 1) rownames(each.tx.mat) <- uni.tx[i]
            final.intron.mat <- rbind(final.intron.mat,each.tx.mat)
        }
        final.result.list <- NULL
        total.list.names <- NULL
        result.n <- 1
        for(i in 1:nrow(alt.mat.int)){
            ES.result <- EStest(final.exon.mat,final.intron.mat,rbind(alt.mat.int[i,]))
            if (newTypes == "anno") ES.result <- rbind(ES.result[ES.result[,"status"] == "exist",])
            else if (newTypes == "new")    ES.result <- rbind(ES.result[ES.result[,"status"] == "possible",])
            if (length(ES.result) != 0){
                final.result.list[[result.n]] <- ES.result
                total.list.names <- c(total.list.names,paste("ES.",length(grep("ES",total.list.names))+1,sep=""))
                result.n <- result.n+1
            }
            ASS.result <- ASStest(final.exon.mat,final.intron.mat,rbind(alt.mat.int[i,]))
            if (newTypes == "anno") ASS.result <- rbind(ASS.result[ASS.result[,"status"] == "exist",])
            else if (newTypes == "new")    ASS.result <- rbind(ASS.result[ASS.result[,"status"] == "possible",])
            if(length(ASS.result) != 0){
                final.result.list[[result.n]] <- ASS.result
                total.list.names <- c(total.list.names,paste("ASS.",length(grep("ASS",total.list.names))+1,sep=""))
                result.n <- result.n+1
            }
            NI.result <- NItest(final.exon.mat,final.intron.mat,rbind(alt.mat.int[i,]),tx.gene=tx.gene)
            if (newTypes == "anno") NI.result <- rbind(NI.result[NI.result[,"status"] == "exist",])
            else if (newTypes == "new")    NI.result <- rbind(NI.result[NI.result[,"status"] == "possible",])
            if(length(NI.result) != 0){
                final.result.list[[result.n]] <- NI.result
                total.list.names <- c(total.list.names,paste("IR.",length(grep("IR",total.list.names))+1,sep=""))
                result.n <- result.n+1
            }
        }
        names(final.result.list) <- total.list.names
        cluster.result <- list(NULL,NULL,NULL)
        names(cluster.result) <- c("ES","ASS","IR")
        if (!is.null(final.result.list)){
            for (each.list.num in 1:length(final.result.list)){
                result.type <- unlist(strsplit(names(final.result.list[each.list.num]),"[.]"))[1]
                cluster.result[[result.type]] <- unique(rbind(cluster.result[[result.type]],final.result.list[[each.list.num]]))
            }
        }
        if (length(cluster.result[["ES"]]) != 0){
            row.names(cluster.result[["ES"]]) <- 1:nrow(rbind(cluster.result[["ES"]]))
            cluster.result[["ES"]] <- cbind(unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),cluster.result[["ES"]])
            colnames(cluster.result[["ES"]])[1:3] <- c("EnsID","Nchr","Strand")
        }
        if (length(cluster.result[["ASS"]]) != 0){
            row.names(cluster.result[["ASS"]]) <- 1:nrow(rbind(cluster.result[["ASS"]]))
            cluster.result[["ASS"]] <- cbind(unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),cluster.result[["ASS"]])
            colnames(cluster.result[["ASS"]])[1:3] <- c("EnsID","Nchr","Strand")
        }
        if (length(cluster.result[["IR"]]) != 0){
            row.names(cluster.result[["IR"]]) <- 1:nrow(rbind(cluster.result[["IR"]]))
            cluster.result[["IR"]] <- cbind(unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),cluster.result[["IR"]])
            colnames(cluster.result[["IR"]])[1:3] <- c("EnsID","Nchr","Strand")
        }
        return(cluster.result)
    }
    GTFdb <- chrseparate(GTFdb,1:22)    
    registerDoParallel(cores=Ncor)
    newTypes <- "both"
    trans.exon.range <- exonsBy(GTFdb,by="tx")
    trans.intron.range <- intronsByTranscript(GTFdb)
    txTable <- select(GTFdb, keys=names(trans.exon.range), columns=c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND"), keytype="TXID")
    txTable <- gsub(" ","",as.matrix(txTable))
    if (length(calGene) != 0){
        txTable <- rbind(txTable[txTable[,"GENEID"] == calGene,])
        }
    Total.chr <- as.character(sort(as.integer(unique(txTable[,"TXCHROM"]))))
    if (length(txTable) != 0){
            txTable[,2] <- do.call(rbind,strsplit(txTable[,2],"[.]"))[,1]
            txTable[,3] <- do.call(rbind,strsplit(txTable[,3],"[.]"))[,1]
            }
    ES.finl.result <- NULL
    ASS.finl.result <- NULL
    IR.finl.result <- NULL
    for (i in 1:length(Total.chr)){
        print (paste("-------------------Processing : chr ",Total.chr[i]," -------------------",sep=""))
        each.chr.db <- chrseparate(GTFdb,Total.chr[i])
        each.chr.names <- names(unlist(trans.exon.range))[as.character(seqnames(unlist(trans.exon.range))) == Total.chr[i]]
        each.chr.exon.range <- trans.exon.range[is.element(names(trans.exon.range),each.chr.names),]
        each.chr.intron.range <- trans.intron.range[is.element(names(trans.intron.range),each.chr.names),]
        transnum <- names(each.chr.exon.range)
        each.chr.txTable <- txTable[txTable[,"TXCHROM"] == Total.chr[i],]
        tx.num <- table(each.chr.txTable[,"GENEID"])
        not.one.gene <- which(tx.num > 1)
        up2gene <- names(not.one.gene)
        called.packages <- c("GenomicRanges","GenomicFeatures")
        j = NULL
        pa.result <- foreach(j=1:length(up2gene),.packages=called.packages) %dopar% {
            strandinfo <- unique(each.chr.txTable[each.chr.txTable[,"GENEID"] == up2gene[j],"TXSTRAND"])
            Altvalue <- findAlternative(up2gene[j],each.chr.txTable,each.chr.exon.range,each.chr.intron.range,Total.chr[i])
            Alt.result <- CalAlt(Altvalue)
            }
        ES.finl.result <- rbind(ES.finl.result,do.call(rbind,lapply(pa.result,function(x) x$"ES")))
        ASS.finl.result <- rbind(ASS.finl.result,do.call(rbind,lapply(pa.result,function(x) x$"ASS")))
        IR.finl.result <- rbind(IR.finl.result,do.call(rbind,lapply(pa.result,function(x) x$"IR")))
        }
    if (length(ES.finl.result) != 0){
        ES.finl.result <- cbind(Index=paste("ES",1:nrow(ES.finl.result),sep=""),rbind(ES.finl.result[,-ncol(ES.finl.result)]))
        rownames(ES.finl.result) <- 1:nrow(ES.finl.result)
        }
    else    ES.finl.result <- matrix("NA")
    if (length(ASS.finl.result) != 0){
        ASS.finl.result <- cbind(Index=paste("ASS",1:nrow(ASS.finl.result),sep=""),rbind(ASS.finl.result[,-ncol(ASS.finl.result)]))
        rownames(ASS.finl.result) <- 1:nrow(ASS.finl.result)
        }
    else    ASS.finl.result <- matrix("NA")
    if (length(IR.finl.result) != 0){
        IR.finl.result <- cbind(Index=paste("IR",1:nrow(IR.finl.result),sep=""),rbind(IR.finl.result[,-ncol(IR.finl.result)]))
        rownames(IR.finl.result) <- 1:nrow(IR.finl.result)
        }
    else    IR.finl.result <- matrix("NA")
    final.total.result <- list(ES.finl.result,ASS.finl.result,IR.finl.result)
    names(final.total.result) <- c("ES","ASS","IR")
    ASdb <- new("ASdb",SplicingModel=final.total.result)
    if (length(out.dir) != 0){
        system(paste("mkdir -p ",out.dir,"/AS_Pattern",sep=""))
        write.table(ES.finl.result,paste(out.dir,"/AS_Pattern/ES_pattern.txt",sep=""),sep='\t',quote=FALSE)
        write.table(ASS.finl.result,paste(out.dir,"/AS_Pattern/ASS_pattern.txt",sep=""),sep='\t',quote=FALSE)
        write.table(IR.finl.result,paste(out.dir,"/AS_Pattern/IR_pattern.txt",sep=""),sep='\t',quote=FALSE)
        }
    return (ASdb)
    }

