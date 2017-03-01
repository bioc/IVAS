RatioFromFPKM <- function(GTFdb=NULL,ASdb=NULL,Total.expdata=NULL,CalIndex=NULL,Ncor=1,out.dir=NULL){
  ES.in.sk.ratio <- NULL
  ASS.in.sk.ratio <- NULL
  IR.in.sk.ratio <- NULL
  calratio <- function(expdata=NULL,in.sk.mat){
    if (length(in.sk.mat) == 0) next
    if (in.sk.mat[1,ncol(in.sk.mat)] == "IR") colnames(in.sk.mat) <- c("Index","exon1","down","up","inclu","skip","Types")
    else  colnames(in.sk.mat) <- c("Index","exon1","exon2","exon1.nei","exon2.nei","inclu","skip","Types")
    final.ratio.mat <- NULL
    for (i in 1:nrow(in.sk.mat)){
      inclu.exp <- rbind(expdata[is.element(rownames(expdata),unlist(strsplit(in.sk.mat[i,"inclu"],"[|]"))),])
      skip.exp <- rbind(expdata[is.element(rownames(expdata),unlist(strsplit(in.sk.mat[i,"skip"],"[|]"))),])
      sum.inclu.exp <- apply(inclu.exp,2,function(x) sum(as.double(x)))
      sum.skip.exp <- apply(skip.exp,2,function(x) sum(as.double(x)))
      in.sk.ratio <- rbind(sum.skip.exp/(sum.inclu.exp+sum.skip.exp))
      in.sk.ratio <- rbind(c(in.sk.mat[i,],in.sk.ratio))
      final.ratio.mat <- rbind(final.ratio.mat,in.sk.ratio)
    }
    if (in.sk.mat[1,ncol(in.sk.mat)] == "IR") colnames(final.ratio.mat) <- c("Index","exons1","down","up","inclusion","skip","Types",colnames(expdata))
    else  colnames(final.ratio.mat) <- c("Index","exons1","exon2","exon1.nei","exon2.nei","inclusion","skip","Types",colnames(expdata))
    rownames(final.ratio.mat) <- 1:nrow(final.ratio.mat)
    return(final.ratio.mat)
  }
  estimateRatio <- function(expdata=NULL,tx.gene=NULL,exon.locus=NULL,test.value=NULL){
    ES.in.sk.ratio <- NULL
    ASS.in.sk.ratio <- NULL
    IR.in.sk.ratio <- NULL
    if(length(test.value[["ES"]]) != 0){
      ES.cluster.result <- rbind(test.value[["ES"]])
      fir.ES.result <- rbind(ES.cluster.result[ES.cluster.result[,"2ndEX"] == "NA",])
      ES.in.sk.mat <- lapply(1:length(fir.ES.result[,"1stEX"]),function(each.nums){
        index.num <- fir.ES.result[each.nums,"Index"]
        each.targets <- fir.ES.result[each.nums,"1stEX"]
        each.ES.result <- rbind(fir.ES.result[each.nums,])
        s.tar.ex <- unlist(strsplit(each.ES.result[,"1st_des"],","))
        s.do.ex <- unlist(strsplit(each.ES.result[,"Do_des"],","))
        s.up.ex <- unlist(strsplit(each.ES.result[,"Up_des"],","))
        test.Do.tx <- names(exon.locus[is.element(exon.locus,s.do.ex) | is.element(exon.locus,s.do.ex)])
        test.up.tx <- names(exon.locus[is.element(exon.locus,s.up.ex) | is.element(exon.locus,s.up.ex)])
        test.tx.ex <- exon.locus[is.element(names(exon.locus),intersect(test.Do.tx,test.up.tx))]
        inclused.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.tar.ex)])
        skipped.tx <-  unique(names(test.tx.ex)[!is.element(names(test.tx.ex),inclused.tx)])
        inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),inclused.tx),"TXNAME"]
        skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),skipped.tx),"TXNAME"]
        test.in.sk.mat <- cbind(Index=index.num,each.targets,"NA",each.ES.result[,"DownEX"],each.ES.result[,"UpEX"],
                                paste(inclused.tx,collapse="|"),paste(skipped.tx,collapse="|"),each.ES.result[,"Types"])
        test.in.sk.mat
      })
      ES.in.sk.mat <- do.call(rbind,ES.in.sk.mat)
      ES.in.sk.ratio <- rbind(calratio(expdata,ES.in.sk.mat))
      
      sec.ES.result <- rbind(ES.cluster.result[ES.cluster.result[,"2ndEX"] != "NA" & ES.cluster.result[,"Types"] == "ES",])
      if (length(sec.ES.result) != 0){
        sec.ES.in.sk.mat <- lapply(1:length(sec.ES.result[,"1stEX"]),function(each.nums){
          index.num <- sec.ES.result[each.nums,"Index"]
          each.targets <- sec.ES.result[each.nums,"1stEX"]
          each.ES.result <- rbind(sec.ES.result[each.nums,])
          s.fi.ex <- unlist(strsplit(each.ES.result[,"1st_des"],","))
          s.se.ex <- unlist(strsplit(each.ES.result[,"2nd_des"],","))
          s.do.ex <- unlist(strsplit(each.ES.result[,"Do_des"],","))
          s.up.ex <- unlist(strsplit(each.ES.result[,"Up_des"],","))
          test.Do.tx <- names(exon.locus[is.element(exon.locus,s.do.ex)])
          test.up.tx <- names(exon.locus[is.element(exon.locus,s.up.ex)])
          test.tx.ex <- exon.locus[is.element(names(exon.locus),intersect(test.Do.tx,test.up.tx))]
          inclused.tx <- intersect(names(test.tx.ex)[is.element(test.tx.ex,s.fi.ex)],names(test.tx.ex)[is.element(test.tx.ex,s.se.ex)])
          skipped.tx <-  unique(names(test.tx.ex)[!is.element(names(test.tx.ex),inclused.tx)])
          inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),inclused.tx),"TXNAME"]
          skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),skipped.tx),"TXNAME"]
          test.in.sk.mat <- cbind(Index=index.num,each.targets,each.ES.result[,"2ndEX"],each.ES.result[,"DownEX"],each.ES.result[,"UpEX"],
                                  paste(inclused.tx,collapse="|"),paste(skipped.tx,collapse="|"),each.ES.result[,"Types"])
          test.in.sk.mat
        })
        sec.ES.in.sk.mat <- do.call(rbind,sec.ES.in.sk.mat)
        ES.in.sk.ratio <- rbind(ES.in.sk.ratio,calratio(expdata,sec.ES.in.sk.mat))
      }
      
      MXE.ES.result <- rbind(ES.cluster.result[ES.cluster.result[,"2ndEX"] != "NA" & ES.cluster.result[,"Types"] == "MXE",])
      if (length(MXE.ES.result) != 0){
        MXE.ES.in.sk.mat <- lapply(1:length(MXE.ES.result[,"1stEX"]),function(each.nums){
          index.num <- MXE.ES.result[each.nums,"Index"]
          each.targets <- MXE.ES.result[each.nums,"1stEX"]
          each.ES.result <- rbind(MXE.ES.result[each.nums,])
          s.fi.ex <- unlist(strsplit(each.ES.result[,"1st_des"],","))
          s.sec.ex <- unlist(strsplit(each.ES.result[,"2nd_des"],","))
          s.do.ex <- unlist(strsplit(each.ES.result[,"Do_des"],","))
          s.up.ex <- unlist(strsplit(each.ES.result[,"Up_des"],","))
          test.Do.tx <- names(exon.locus[is.element(exon.locus,s.do.ex)])
          test.up.tx <- names(exon.locus[is.element(exon.locus,s.up.ex)])
          test.tx.ex <- exon.locus[is.element(names(exon.locus),intersect(test.Do.tx,test.up.tx))]
          fi.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.fi.ex)])
          se.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.sec.ex)])
          inclused.tx <- fi.tx[!is.element(fi.tx,se.tx)]
          skipped.tx <- se.tx[!is.element(se.tx,fi.tx)]
          inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),inclused.tx),"TXNAME"]
          skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),skipped.tx),"TXNAME"]
          test.in.sk.mat <- cbind(Index=index.num,each.targets,rbind(each.ES.result[,c("2ndEX","DownEX","UpEX")]),
                                  paste(inclused.tx,collapse="|"),paste(skipped.tx,collapse="|"),each.ES.result[,"Types"])
          test.in.sk.mat
        })
        MXE.ES.in.sk.mat <- do.call(rbind,MXE.ES.in.sk.mat)
        MXE.ES.in.sk.mat <- rbind(MXE.ES.in.sk.mat[MXE.ES.in.sk.mat[,5] != "" & MXE.ES.in.sk.mat[,6] != "",])
        if (length(MXE.ES.in.sk.mat) != 0){
          ES.in.sk.ratio <- rbind(ES.in.sk.ratio,calratio(expdata,MXE.ES.in.sk.mat))
        }
      }
      ES.in.sk.ratio <- cbind(ES.in.sk.ratio[,"Index"],unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),rbind(ES.in.sk.ratio[,!is.element(colnames(ES.in.sk.ratio),c("Index"))]))
      colnames(ES.in.sk.ratio) <- c("Index","EnsID","Nchr","Strand","1stEX","2ndEX","DownEX","UpEX","Inclusion_TX","Skip_TX","Types",colnames(ES.in.sk.ratio)[-c(1:which(colnames(ES.in.sk.ratio) == "Types"))])
      rownames(ES.in.sk.ratio) <- 1:nrow(ES.in.sk.ratio)
    }
    
    
    if(length(test.value[["ASS"]]) != 0){
      ASS.cluster.result <- rbind(test.value[["ASS"]])
      len.nei <- grep("NeighborEX",colnames(ASS.cluster.result))
      ASS.in.sk.mat <- lapply(1:length(ASS.cluster.result[,"ShortEX"]),function(each.nums){
        index.num <- ASS.cluster.result[each.nums,"Index"]
        each.targets <- ASS.cluster.result[each.nums,"ShortEX"]
        each.ASS.result <- rbind(ASS.cluster.result[each.nums,])
        s.short.ex <- unlist(strsplit(each.ASS.result[,"Short_des"],","))
        s.long.ex <- unlist(strsplit(each.ASS.result[,"Long_des"],","))
        if (length(len.nei) == 1){
          nei.ex <- unlist(strsplit(each.ASS.result[,"Neighbor_des"],","))
          test.tx <- names(exon.locus[is.element(exon.locus,s.short.ex) | is.element(exon.locus,s.long.ex) | is.element(exon.locus,nei.ex)])
          test.tx.ex <- exon.locus[is.element(names(exon.locus),test.tx)]
          shortex <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.short.ex)])
          longex <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.long.ex)])
          inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),longex),"TXNAME"]
          skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),shortex),"TXNAME"]
          test.in.sk.mat <- cbind(Index=index.num,each.targets,each.ASS.result[,"LongEX"],ShortNeighborEX=each.ASS.result[,"NeighborEX"],LongNeighborEX="NA",
                                  paste(skipped.tx,collapse="|"),paste(inclused.tx,collapse="|"),each.ASS.result[,"Types"])
        }
        else if (length(len.nei) == 2){
          short.nei.ex <- unlist(strsplit(each.ASS.result[,"ShortNeighbor_des"],","))
          long.nei.ex <- unlist(strsplit(each.ASS.result[,"LongNeighbor_des"],","))
          test.tx <- c(intersect(names(exon.locus[(is.element(exon.locus,s.short.ex))]),names(exon.locus[(is.element(exon.locus,short.nei.ex))])),
                       intersect(names(exon.locus[(is.element(exon.locus,s.long.ex))]),names(exon.locus[(is.element(exon.locus,long.nei.ex))])))
          test.tx.ex <- exon.locus[is.element(names(exon.locus),test.tx)]
          shortex <- unique(intersect(names(test.tx.ex)[is.element(test.tx.ex,s.short.ex)],names(test.tx.ex)[is.element(test.tx.ex,short.nei.ex)]))
          longex <- unique(intersect(names(test.tx.ex)[is.element(test.tx.ex,s.long.ex)],names(test.tx.ex)[is.element(test.tx.ex,long.nei.ex)]))
          inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),longex),"TXNAME"]
          skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),shortex),"TXNAME"]
          test.in.sk.mat <- cbind(Index=index.num,each.targets,each.ASS.result[,"LongEX"],each.ASS.result[,"ShortNeighborEX"],each.ASS.result[,"LongNeighborEX"],
                                  paste(skipped.tx,collapse="|"),paste(inclused.tx,collapse="|"),each.ASS.result[,"Types"])
        }
        test.in.sk.mat
      })
      ASS.in.sk.mat <- do.call(rbind,ASS.in.sk.mat)
      ASS.in.sk.ratio <- rbind(calratio(expdata,ASS.in.sk.mat))
      ASS.in.sk.ratio <- cbind(ASS.in.sk.ratio[,"Index"],unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),rbind(ASS.in.sk.ratio[,!is.element(colnames(ASS.in.sk.ratio),c("Index"))]))
      rownames(ASS.in.sk.ratio) <- 1:nrow(rbind(ASS.in.sk.ratio))
      colnames(ASS.in.sk.ratio) <- c("Index","EnsID","Nchr","Strand","ShortEX","LongEX","ShortNeighborEX","LongNeighborEX","Short_TX","Long_TX","Types",
                                     colnames(ASS.in.sk.ratio)[-c(1:which(colnames(ASS.in.sk.ratio) == "Types"))])
      if (length(len.nei) == 1){
        ASS.in.sk.ratio <- rbind(ASS.in.sk.ratio[,colnames(ASS.in.sk.ratio) != "LongNeighborEX"])
        colnames(ASS.in.sk.ratio)[which(colnames(ASS.in.sk.ratio) == "ShortNeighborEX")] <- "NeighborEX"
      }
      rownames(ASS.in.sk.ratio) <- 1:nrow(ASS.in.sk.ratio)
    }
    
    if(length(test.value[["IR"]]) != 0){
      IR.cluster.result <- rbind(test.value[["IR"]])
      IR.in.sk.mat <- lapply(1:length(IR.cluster.result[,"RetainEX"]),function(each.nums){
        index.num <- IR.cluster.result[each.nums,"Index"]
        each.targets <- IR.cluster.result[each.nums,"RetainEX"]
        each.IR.result <- rbind(IR.cluster.result[each.nums,])
        s.fi.ex <- unlist(strsplit(each.IR.result[,"Retain_des"],","))
        s.do.ex <- unlist(strsplit(each.IR.result[,"Do_des"],","))
        s.up.ex <- unlist(strsplit(each.IR.result[,"Up_des"],","))
        test.Do.tx <- names(exon.locus[is.element(exon.locus,s.do.ex)])
        test.up.tx <- names(exon.locus[is.element(exon.locus,s.up.ex)])
        test.re.tx <- names(exon.locus[is.element(exon.locus,s.fi.ex)])
        test.tx.ex <- exon.locus[is.element(names(exon.locus),c(intersect(test.Do.tx,test.up.tx),test.re.tx))]
        fi.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.do.ex)])
        se.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.up.ex)])
        re.tx <- unique(names(test.tx.ex)[is.element(test.tx.ex,s.fi.ex)])
        inclused.tx <- unique(intersect(fi.tx,se.tx))
        skipped.tx <- re.tx
        inclused.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),inclused.tx),"TXNAME"]
        skipped.tx <- tx.gene[is.element(as.matrix(tx.gene["TXID"]),skipped.tx),"TXNAME"]
        test.in.sk.mat <- cbind(Index=index.num,each.targets,each.IR.result[,"DownEX"],each.IR.result[,"UpEX"],paste(inclused.tx,collapse="|"),paste(skipped.tx,collapse="|"),"IR")
        test.in.sk.mat
      })
      IR.in.sk.mat <- do.call(rbind,IR.in.sk.mat)
      IR.in.sk.ratio <- rbind(calratio(expdata,IR.in.sk.mat))
      IR.in.sk.ratio <- cbind(IR.in.sk.ratio[,"Index"],unique(tx.gene[,"GENEID"]),unique(tx.gene[,"TXCHROM"]),unique(tx.gene[,"TXSTRAND"]),rbind(IR.in.sk.ratio[,!is.element(colnames(IR.in.sk.ratio),c("Index"))]))
      colnames(IR.in.sk.ratio) <- c("Index","EnsID","Nchr","Strand","RetainEX","DownEX","UpEX","Normal_TX","IR_TX","Types",colnames(IR.in.sk.ratio)[-c(1:which(colnames(IR.in.sk.ratio) == "Types"))])
      rownames(IR.in.sk.ratio) <- 1:nrow(IR.in.sk.ratio)
    }
    final.result.mat <- list(ES.in.sk.ratio,ASS.in.sk.ratio,IR.in.sk.ratio)
    names(final.result.mat) <- c("ES","ASS","IR")
    return (final.result.mat)
  }
  registerDoParallel(cores=Ncor)
  
  
  
  Total.expdata <- as.matrix(Total.expdata)
  pre.ES.alt <- ASdb@"SplicingModel"[["ES"]]
  pre.ASS.alt <- ASdb@"SplicingModel"[["ASS"]]
  pre.IR.alt <- ASdb@"SplicingModel"[["IR"]]
  
  tested.geneid <- unique(c(pre.ES.alt[,is.element(colnames(pre.ES.alt),"EnsID")],pre.ASS.alt[,is.element(colnames(pre.ASS.alt),"EnsID")],pre.IR.alt[,is.element(colnames(pre.IR.alt),"EnsID")]))
  total.chr <- unique(c(pre.ES.alt[,is.element(colnames(pre.ES.alt),"Nchr")],pre.ASS.alt[,is.element(colnames(pre.ASS.alt),"Nchr")],pre.IR.alt[,is.element(colnames(pre.IR.alt),"Nchr")]))
  total.chr <- total.chr[order(as.integer(gsub("chr","",total.chr)))]
  GTFdb <- chrseparate(GTFdb,total.chr)
  total.exon.range <- exonsBy(GTFdb,by="tx")
  total.intron.range <- intronsByTranscript(GTFdb)
  txTable <- try(select(GTFdb, keys=names(total.exon.range), columns=c("TXCHROM","TXNAME","GENEID","TXSTART","TXEND","TXSTRAND"), keytype="TXID"),silent=T)
  
  sub.txTable <- rbind(txTable[is.element(txTable[,"GENEID"],tested.geneid),])
  sub.expdata <- Total.expdata[is.element(rownames(Total.expdata),sub.txTable[,"TXNAME"]),]
  sub.exon.range <- total.exon.range[is.element(names(total.exon.range),as.matrix(sub.txTable[,"TXID"])),]
  sub.intron.range <- total.intron.range[is.element(names(total.intron.range),as.matrix(sub.txTable[,"TXID"])),]
  sampleid <- colnames(sub.expdata)
  called.packages <- c("GenomicRanges","GenomicFeatures")
  AltType <- c("ES","ASS","IR")
  
  final.total.ratio <- list(matrix("NA"),matrix("NA"),matrix("NA"))
  names(final.total.ratio) <- c("ES","ASS","IR")
  j=NULL
  for (i in 1:length(AltType)){
    each.type.result <- slot(ASdb,"SplicingModel")[[AltType[i]]]
    if (length(each.type.result) != 0 & each.type.result[1,1] != "NA" & length(CalIndex) != 0)  each.type.result <- rbind(each.type.result[is.element(each.type.result[,"Index"],CalIndex),])
    if (length(each.type.result) != 0){
      u.genes <- unique(each.type.result[,"EnsID"])
      each.ratio <- NULL
      pa.result <- foreach(j=1:length(u.genes),.packages=called.packages,.combine=rbind) %dopar% {
        each.ES <- rbind(each.type.result[each.type.result[,"EnsID"] == u.genes[j],])
        each.tx.info <- txTable[txTable[,"GENEID"] == unique(each.ES[,"EnsID"]),]
        Exon.info <- findAlternative(unique(each.ES[,"EnsID"]),sub.txTable,sub.exon.range,sub.intron.range,unique(each.ES[,"Nchr"]))$exonRange
        exon.start <- unlist(start(Exon.info))
        exon.end <- unlist(end(Exon.info))
        exon.mat <- cbind(exon.start,exon.end)
        each.exon.locus <- paste(exon.start,exon.end,sep="-")
        names(each.exon.locus) <- names(exon.end)
        expdata <- sub.expdata[intersect(rownames(sub.expdata),each.tx.info[,"TXNAME"]),]
        each.ES <- list(each.ES)
        names(each.ES) <- AltType[i]
        if (length(expdata) != 0){
          each.ratio <- estimateRatio(expdata,each.tx.info,each.exon.locus,each.ES)[[AltType[i]]]
          each.ratio
        }
        else  NULL
      }
      if (length(pa.result) != 0){
        rownames(pa.result) <- c(1:nrow(pa.result))
        final.total.ratio[[i]] <- pa.result
      }
      else final.total.ratio[[i]] <- as.matrix("NA")
    }
    else if (length(each.type.result) == 0){
      final.total.ratio[[i]] <- as.matrix("NA")
    }
  }
  ASdb <- new("ASdb",SplicingModel=ASdb@"SplicingModel",Ratio=final.total.ratio,GroupDiff=ASdb@"GroupDiff",sQTLs=ASdb@"sQTLs",Me.sQTLs=ASdb@"Me.sQTLs",Clinical=ASdb@"Clinical")
  if (length(out.dir) != 0){
    system(paste("mkdir -p ",out.dir,"/AS_Ratio",sep=""))
    write.table(final.total.ratio[["ES"]],paste(out.dir,"/AS_Ratio/ES_Ratio.txt",sep=""),sep='\t',quote=F)
    write.table(final.total.ratio[["ASS"]],paste(out.dir,"/AS_Ratio/ASS_Ratio.txt",sep=""),sep='\t',quote=F)
    write.table(final.total.ratio[["IR"]],paste(out.dir,"/AS_Ratio/IR_Ratio.txt",sep=""),sep='\t',quote=F)
  }
  return(ASdb)
}