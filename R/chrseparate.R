chrseparate <- function(transdb=NULL,chrname=NULL){
    isActiveSeq(transdb)[seqlevels(transdb)] <- FALSE
    chrsepa <- TRUE
    names(chrsepa) <- chrname
    try.test <- try(isActiveSeq(transdb) <- chrsepa,silent=TRUE)
    if (length(names(try.test)) == 0){
        return (NULL)
        }
    else {
        isActiveSeq(transdb) <- chrsepa
        return (transdb)
        }
    }
