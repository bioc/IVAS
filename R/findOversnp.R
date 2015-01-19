findOversnp <- function(altInvalue=NULL,snprange=NULL){
    Inrange <- altInvalue[["alterIntron"]]
    if (length(Inrange) == 0){
        return (NULL)
        }
    Inoverapmatrix <- NULL
    Startoverapmatrix <- NULL
    Endoverapmatrix <- NULL
    if (length(Inrange) != 0){
        Inoverap <- as.matrix(findOverlaps(snprange,Inrange,select="all"))
        range <- paste(start(Inrange),end(Inrange),sep="-")
        snplist <- unlist(elementMetadata(snprange))
        if (length(Inoverap) != 0){
            over.range <- sapply(Inoverap[,"subjectHits"],function(i){
                range[as.integer(i)]
                })
            over.snp <- sapply(Inoverap[,"queryHits"],function(i){
                snplist[as.integer(i)]
                })
            Inoverapmatrix <- cbind(over.range,over.snp)
            colnames(Inoverapmatrix) <- c("range","snp")
            rownames(Inoverapmatrix) <- c(1:nrow(Inoverapmatrix))
            }
        }
    return (unique(Inoverapmatrix))
    }
