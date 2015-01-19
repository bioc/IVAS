saveBplot <- function(sig.sqtl=NULL,expdata=NULL,snpdata=NULL,snplocus=NULL,GTFdata=NULL,outdir=NULL){
    if(length(as.matrix(sig.sqtl)) == 10){
        h.name <- names(sig.sqtl)
        if (length(h.name)==0){h.name <- colnames(sig.sqtl)}
        sig.sqtl <- matrix(sig.sqtl,ncol=10)
        colnames(sig.sqtl)<-h.name
        }
    else{
        h.name <- colnames(sig.sqtl)
        sig.sqtl <- as.matrix(sig.sqtl)
        colnames(sig.sqtl) <- h.name
        sig.sqtl <- gsub(" ","",sig.sqtl)
        }
    Nchr <- unique(sig.sqtl[,"CHR"])
    if (file.exists(outdir) == "FALSE"){
        dir.create(outdir)
        }
    sub.snp.locus <- snplocus[is.element(snplocus[,"SNP"],sig.sqtl[,"SNP"]),]
    sub.snpdata <- snpdata[is.element(rownames(snpdata),sig.sqtl[,"SNP"]),]
    for (n in 1:length(Nchr)){
        each.sig.sqtl <- sig.sqtl[which(sig.sqtl[,"CHR"] == Nchr[n]),]
        each.sig.sqtl <- matrix(each.sig.sqtl,ncol=ncol(sig.sqtl))
        colnames(each.sig.sqtl) <- colnames(sig.sqtl)
        transdb <- chrseparate(GTFdata,Nchr[n])
        exoninfo <- exonsBy(transdb,by="tx")
        introninfo <- intronsByTranscript(transdb)
        txTable <- select(transdb, keys=names(exoninfo), columns=c("TXID","TXNAME","GENEID","TXSTART","TXEND"), keytype="TXID")
        for (s in 1:nrow(each.sig.sqtl)){
            each.sig.gene <- each.sig.sqtl[s,]
            names(each.sig.gene) <- colnames(sig.sqtl)
            Altinfo <- findAlternative(each.sig.gene["gene"],txTable,exoninfo,introninfo,Nchr[n])
            each.snp <- which(sub.snp.locus["SNP"]==each.sig.gene["SNP"])
            each.snplocus <- sub.snp.locus[each.snp,]
            snp.position <- as.integer(as.matrix(each.snplocus["locus"]))
            snp.chr <- as.character(as.matrix(each.snplocus["CHR"]))
            irange <- IRanges(start=snp.position,end=snp.position)
            sig.snp <- GRanges(seqnames=Rle(snp.chr),ranges=irange,metadata=as.character(unlist(each.snplocus["SNP"])))
            overlapsnp <- findOversnp(Altinfo,sig.snp)
            group.matrix <- sqtlfinder(Altinfo,overlapsnp,expdata,sub.snpdata,"boxplot")
            group.matrix <- group.matrix[[each.sig.gene["targetExon"]]]
            ratio <- as.double(group.matrix[,"ratio"])
            names(ratio) <- rownames(group.matrix)
            rsname <- colnames(group.matrix)[2]
            geno <- group.matrix[,rsname]
            geno.names <- names(geno)
            alleles <- names(table(as.matrix(geno)))
            hetero <- sapply(alleles,function(z){
                allele <- strsplit(z,"")[[1]]
                if (allele[1] != allele[2]){
                    z
                    }})
            hetero <- unlist(hetero)
            if (length(hetero)>1){
                geno <- gsub(hetero[1],"hetero",as.character(as.matrix(geno)))
                geno <- matrix(geno,ncol=length(geno.names))
                geno <- gsub(hetero[2],"hetero",as.character(as.matrix(geno)))
                geno <- matrix(geno,ncol=length(geno.names))
                geno <- gsub("hetero",paste(hetero[1],hetero[2],sep="/"),as.character(as.matrix(geno)))
                geno <- matrix(geno,ncol=length(geno.names))
                colnames(geno) <- geno.names
                }
            boxmatrix <- cbind(as.double(ratio),as.character(geno))
            colnames(boxmatrix) <- c("ratio","geno")
            png(filename=paste(outdir,"/",each.sig.gene["gene"],"_",rsname,"_",each.sig.gene["method"],".png",sep=""))
            boxplot(as.double(ratio)~as.character(geno),xlab=(paste("geneID : ",as.character(each.sig.gene["gene"]),sep="")))
            dev.off()
            }
        }
    }
