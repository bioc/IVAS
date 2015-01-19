findAlternative <- function(geneid=NULL,txTable=NULL,totalExrange=NULL,totalInrange=NULL,one.chr=NULL){
    tableBygene <- txTable[which(txTable[,"GENEID"]==geneid),]
    txid <- tableBygene[,"TXID"]
    txstarts <- tableBygene[,"TXSTART"]
    txends <- tableBygene[,"TXEND"]
    each.exon <- totalExrange[as.character(txid)]
    each.intron <- totalInrange[as.character(txid)]
    if (length(txid) > 1 & length(start(unlist(each.intron))) > 0){
        nosigintron <- NULL
        strandinfo <- paste(strand(each.exon)[[1]][1],"",sep="")
        sortintron <- paste(unlist(start(each.intron)),unlist(end(each.intron)),sep="-")
        intronName <- unique(sortintron)
        nosigintron <- lapply(intronName,function(int.na){
            spl.cal <- as.integer(strsplit(int.na,"-")[[1]])
            inclu.tx <- txstarts[(spl.cal[1] < txstarts & spl.cal[2] < txstarts) | (spl.cal[1] > txends & spl.cal[2] > txends)]
            startcount <- length(txid) - length(inclu.tx)
            if (table(sortintron)[int.na] == startcount){
                int.na
                }
            })
        intronName <- unique(unlist(intronName))
        intronName <- intronName[!is.element(intronName,nosigintron)]
        intronrange <- strsplit(intronName,"-")
        In.table <- do.call(rbind,intronrange)
        if (length(intronrange) >0 ){
            colnames(In.table) <- c("start","end")
            irange <- IRanges(start=as.integer(In.table[,"start"]),end=as.integer(In.table[,"end"]))
            alterIntron <- GRanges(seqnames=Rle(one.chr),ranges=irange,strand=strandinfo)
            return.value <- list(alterIntron,tableBygene,each.exon,each.intron)
            names(return.value) <- c("alterIntron","tableBygene","exonRange","intronRange")
            return (return.value)
            }
        }
    return (NULL)
    }

