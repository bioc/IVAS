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
                st.intron <- unlist(start(each.intron))
                en.intron <- unlist(end(each.intron))
                sortintron <- paste(st.intron,en.intron,sep="-")
                intronName <- unique(sortintron)
                for (i in 1:length(intronName)){
                        spl.cal <- as.integer(strsplit(intronName[i],"-")[[1]])
                        startcount <- length(txid) - length(txstarts[(spl.cal[1] < txstarts & spl.cal[2] < txstarts) | (spl.cal[1] > txends & spl.cal[2] > txends)])
                        if (table(sortintron)[intronName[i]] == startcount){
                                nosigintron <- append(nosigintron,intronName[i])
                                }
                        }
                intronName <- intronName[!is.element(intronName,nosigintron)]
                intronrange <- unlist(strsplit(intronName,"-"))
                if (length(intronrange) > 0){
			In.table <- matrix(intronrange,ncol=2,byrow=T)
                	colnames(In.table) <- c("start","end")
                        not.sig <- nrow(In.table)+1
                        for (z in 1:nrow(In.table)){
                                match.num <- NULL
                                match.num <- which(st.intron==In.table[z,1] & en.intron==In.table[z,2])
                                one.num <- NULL
                                last.num <- NULL
                                cal.re <- NULL
                                cal.re2 <- NULL
                                if (is.element(1,match.num) | is.element(length(st.intron),match.num)){
                                        one.num <- which(match.num==1)
                                        last.num <- which(match.num==length(st.intron))
                                        }
                                else {
                                        if (length(which(names(st.intron[match.num]) != names(st.intron[match.num-1]))) != 0){
                                                one.num <- match.num[which(names(st.intron[match.num]) != names(st.intron[match.num-1]))]
                                                }
                                        else if(length(which(names(st.intron[match.num]) != names(st.intron[match.num+1]))) != 0){
                                                last.num <- match.num[which(names(st.intron[match.num]) != names(st.intron[match.num+1]))]
                                                }
                                        }
                                if (length(one.num) > 0){
                                        for (y in 1:length(one.num)){
                                                ex.start <- tableBygene[tableBygene[,1]==names(one.num[y]),"TXSTART"]
                                                cal.re <- which(ex.start > st.intron & ex.start < en.intron)
                                                }
                                        }
                                 else if (length(last.num) > 0){
                                        for (y in 1:length(last.num)){
                                                ex.end <- tableBygene[tableBygene[,1]==names(last.num[y]),"TXEND"]
                                                cal.re2 <- which(ex.end > st.intron & ex.end < en.intron)
                                                }
                                        }
                                if (length(match.num) == 1 & (length(cal.re) + length(cal.re2) > 0)){
                                                not.sig <- c(not.sig,z)
                                                }
                                        }
                        if (length(In.table[-not.sig,1]) > 0){
                                In.table <- In.table[-not.sig,]
				In.table <- matrix(In.table,ncol=2)
				colnames(In.table) <- c("start","end")
                                irange <- IRanges(start=as.integer(In.table[,"start"]),end=as.integer(In.table[,"end"]))
                                alterIntron <- GRanges(seqnames=Rle(one.chr),ranges=irange,strand=strandinfo)
                                return.value <- list(alterIntron,tableBygene,each.exon,each.intron)
                                names(return.value) <- c("alterIntron","tableBygene","exonRange","intronRange")
                                return (return.value)
                                }
                        }
                }
        return (NULL)
        }

