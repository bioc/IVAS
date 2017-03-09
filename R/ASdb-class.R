setClass("ASdb",representation(SplicingModel="list",Ratio="list",GroupDiff="list",sQTLs="list",Me.sQTLs="list",Clinical="list"))
ASdb <- new("ASdb")
SplicingModel <- setMethod("show","ASdb",function(object){
    rm.na.f <- function(mat){
        if (length(mat) == 0) return (c(0,0))
        else if (mat[1] == "NA") return (c(0,0))
        else {
            numrow <- nrow(mat)
            numcol <- ncol(mat) - as.integer(which(colnames(mat) == "Types"))
            return (c(numrow,numcol))
        }
    }
    col.row.col <- function(mat){
        if (length(mat) == 0) return ("0 Rows by 0 samples")
        return (paste(paste(mat,c(" Rows"," samples"),sep=""),collapse=" by "))
    }
    rm.na.SplicingModel <- rbind(rm.na.f(object@"SplicingModel"$"ES"),rm.na.f(object@"SplicingModel"$"ASS"),rm.na.f(object@"SplicingModel"$"IR"))
    rownames(rm.na.SplicingModel) <- c("ES","ASS","IR")
    rm.na.ratio <- rbind(rm.na.f(object@"Ratio"$"ES"),rm.na.f(object@"Ratio"$"ASS"),rm.na.f(object@"Ratio"$"IR"))
    rownames(rm.na.ratio) <- c("ES","ASS","IR")
    rm.na.diff <- rbind(rm.na.f(object@"GroupDiff"$"ES"),rm.na.f(object@"GroupDiff"$"ASS"),rm.na.f(object@"GroupDiff"$"IR"))
    rownames(rm.na.diff) <- c("ES","ASS","IR")
    rm.na.sqtls <- rbind(rm.na.f(object@"sQTLs"$"ES"),rm.na.f(object@"sQTLs"$"ASS"),rm.na.f(object@"sQTLs"$"IR"))
    rownames(rm.na.sqtls) <- c("ES","ASS","IR")
    rm.na.Me.sQTLs <- rbind(rm.na.f(object@"Me.sQTLs"$"ES"),rm.na.f(object@"Me.sQTLs"$"ASS"),rm.na.f(object@"Me.sQTLs"$"IR"))
    rownames(rm.na.Me.sQTLs) <- c("ES","ASS","IR")
    rm.na.Clinical <- rbind(rm.na.f(object@"Clinical"$"ES"),rm.na.f(object@"Clinical"$"ASS"),rm.na.f(object@"Clinical"$"IR"))
    rownames(rm.na.Clinical) <- c("ES","ASS","IR")
    
    slot.names <- slotNames(object)
    la.nums <- unlist(lapply(slot.names,function(each.slot)        length(slot(object,each.slot))))
    in.slot.names <- slot.names[la.nums!=0]
    if (is.element("SplicingModel",in.slot.names)){
        pasted.letter <- paste("Splicing Models : ES = ",rm.na.SplicingModel["ES",1]," Rows & ASS = ",rm.na.SplicingModel["ASS",1]," Rows & IR = ",rm.na.SplicingModel["IR",1]," Rows",'\n',sep="")
    }
    if (is.element("Ratio",in.slot.names)){
        pasted.letter <- paste(pasted.letter,paste("Ratio : ES = ",col.row.col(rm.na.ratio["ES",])," & ASS = ",col.row.col(rm.na.ratio["ASS",])," & IR = ",col.row.col(rm.na.ratio["IR",]),'\n',sep=""),sep="")
    }
    if (is.element("GroupDiff",in.slot.names)){
        pasted.letter <- paste(pasted.letter,paste("Differential Ratio : ES = ",rm.na.diff["ES",1]," Rows & ASS = ",rm.na.diff["ASS",1]," Rows & IR = ",rm.na.diff["IR",1]," Rows",'\n',sep=""),sep="")
    }
    if (is.element("sQTLs",in.slot.names)){
        pasted.letter <- paste(pasted.letter,paste("sQTLs : ES = ",rm.na.sqtls["ES",1]," Rows & ASS = ",rm.na.sqtls["ASS",1]," Rows & IR = ",rm.na.sqtls["IR",1]," Rows",'\n',sep=""),sep="")
    }
    if (is.element("Me.sQTLs",in.slot.names)){
        pasted.letter <- paste(pasted.letter,paste("Me-sQTLs : ES = ",rm.na.Me.sQTLs["ES",1]," Rows & ASS = ",rm.na.Me.sQTLs["ASS",1]," Rows & IR = ",rm.na.Me.sQTLs["IR",1]," Rows",'\n',sep=""),sep="")
    }
    if (is.element("Clinical",in.slot.names)){
        pasted.letter <- paste(pasted.letter,paste("Clinical Analysis : ES = ",rm.na.Clinical["ES",1]," Rows & ASS = ",rm.na.Clinical["ASS",1]," Rows & IR = ",rm.na.Clinical["IR",1]," Rows",'\n',sep=""),sep="")
    }
    cat(paste(pasted.letter,"#ASdb object with ",paste(in.slot.names,collapse=" & "),"\n",sep=""))
})

