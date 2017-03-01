chrseparate <- function(transdb=NULL,chrname=NULL){
    isActiveSeq(transdb)[seqlevels(transdb)] <- FALSE
    chrsepa <- rep(TRUE,length(chrname))
    names(chrsepa) <- chrname
    chrsepa <- chrsepa[is.element(names(chrsepa),names(isActiveSeq(transdb)))]
    try.test <- try(isActiveSeq(transdb) <- chrsepa,silent=TRUE)
    if (length(names(try.test)) == 0){
        return (NULL)
        }
    else {
        isActiveSeq(transdb) <- chrsepa
        return (transdb)
        }
    }
