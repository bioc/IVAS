chrseparate <- function(GTFdb=NULL,chrname=NULL){
    isActiveSeq(GTFdb)[seqlevels(GTFdb)] <- FALSE
    chrsepa <- rep(TRUE,length(chrname))
    names(chrsepa) <- chrname
    chrsepa <- chrsepa[is.element(names(chrsepa),names(isActiveSeq(GTFdb)))]
    try.test <- try(isActiveSeq(GTFdb) <- chrsepa,silent=TRUE)
    if (length(names(try.test)) == 0){
        return (NULL)
        }
    else {
        isActiveSeq(GTFdb) <- chrsepa
        return (GTFdb)
        }
    }
