library(GenomicAlignments)
library(Rsamtools)
library(parallel)

"getcores" <- function(opt_procs){
    # if opt_procs set to 0 then expand to samples by targets
    if( opt_procs == 0 ) opt_procs <- detectCores()
    opt_procs <- min(opt_procs,detectCores())
    return(opt_procs)
}

seqs <- c("NC_004347.2","pBP136_bcA","pMR-1L_hypo","pMR-1S_hypo")

### Given a bam file, read it in, filter results and return a pos_count_matrix
compute_count_matrix <- function(bamfile,mapq = 0){
    #(1) Load the BAM file into a GAlignments object:
    ## Filter on 1) Not primer alignment, 2) duplicate (not sure if working), 3) mapq > 0 [need to test best]
    #param <- ScanBamParam(what=c("seq","mapq","flag"))
    #gal <- readGAlignmentsFromBam(bamfile,index=bamindex, param=param)
    
    param <- ScanBamParam(flag=scanBamFlag(isNotPrimaryRead=NA,isDuplicate=FALSE),what=c("seq","mapq","flag"))
    gal <- readGAlignmentsFromBam(bamfile, param=param)
    
    olen <- length(gal)
    gal <- gal[mcols(gal)$mapq > mapq,] ## Remove Reads not multiply mapped
    qseq <- mcols(gal)$seq  # the query sequences
    print(paste(round((olen-length(gal))/olen*100,1),"% reads removed for too low mapq in file: ",bamfile,sep=""))
    #(2) Use sequenceLayer() to "lay" the query sequences on the reference
    #space. This will remove the parts from the query sequences that
    #correspond to insertions and soft clipping, and it will fill them
    #with - where deletions and/or skipped regions occurred:
    
    
    qseq_on_ref <- sequenceLayer(qseq, cigar(gal),
                                 from="query", to="reference")
    
    qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
    qseq_pos_by_chrom <- splitAsList(start(gal), seqnames(gal))
    
    
    #(3) Compute consensus matrix per chromosome:
    cm_by_chrom <- lapply(names(qseq_pos_by_chrom),
                          function(seqname)
                              consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                                              as.prob=FALSE,
                                              shift=qseq_pos_by_chrom[[seqname]]-1,
                                              width=seqlengths(gal)[[seqname]]))
    cm_by_chrom <- sapply(cm_by_chrom,function(x) x[c("A","C","G","T","N","-"),])
    names(cm_by_chrom) <- levels(seqnames(gal))

    ## 'cm_by_chrom' is a list of consensus matrices. Each matrix
    ## has 18 rows (1 per letter in the DNA alphabet) and 1 column
    ## per chromosome position.
    return(cm_by_chrom)  
}

samples <- read.table("samples.txt",header=T,as.is=T)
bamFolder <- "03-BWA"
mapq = 0

##### Base data set, allelic count data from Bam files
count_data <- mclapply(samples$SAMPLE_ID, function(sample,bf = bamFolder,mq = mapq){
    bamfile=dir(path=file.path(bf,sample),pattern="bam$",full.names = TRUE) 
    return(sapply(bamfile,function(x) compute_count_matrix(x,mapq)))
},mc.cores = getcores(length(samples$SAMPLE_ID)))

names(count_data) <- samples$SAMPLE_ID
count_data <- lapply(count_data,function(x) {names(x) <- seqs;x})

cov_tables <- lapply(seqs,function(x) sapply(lapply(count_data,"[[",x), function(y) summary(colSums(y))))
names(cov_tables) <- seqs

as.numeric.matrix <- function(x) apply(x,2,as.numeric)
monoallelic_sites_tables <- lapply(seqs,function(x) sapply(lapply(count_data,"[[",x), function(y) summary(apply(sweep(y,MARGIN=2,STATS = colSums(y),FUN='/'),2,function(z) any(z==1)))))
names(monoallelic_sites_tables) <- seqs
monoallelic_sites_tables <- lapply(monoallelic_sites_tables,function(x) as.numeric.matrix(x[2:4,]))
monoallelic_sites_tables <- lapply(monoallelic_sites_tables,function(x) {y=round(sweep(x,MARGIN=2,STATS = colSums(x),FUN='/')*100,2);rownames(y) <- c("variable","monomorphic","NA");y})

allelic_counts_tables <- lapply(seqs,function(x) sapply(lapply(count_data,"[[",x), function(y) rowSums(y)))
allelic_counts_tables <- lapply(allelic_counts_tables,function(x) round(sweep(x,MARGIN=2,STATS = colSums(x),FUN='/')*100,2))
names(allelic_counts_tables) <- seqs

Anc <- count_data$Anc
count_data$Anc = NULL

allelic_gain_loss <- function(x,y){
    lapply(seq.int(1,ncol(x)),function(z){
        if (all(c(x[,z]) == 0) | all(c(y[,z]) == 0)) return(NA) ## no coverage in at least one sample
        else {
            return(as.logical(x[,z])-as.logical(y[,z]))
        }
    })
}

gain_loss <- mclapply(count_data, function(sample,ref = Anc,hybrid=TRUE){
    return(sapply(seq.int(1,length(sample)),function(nm) allelic_gain_loss(sample[[nm]],ref[[nm]])))
},mc.cores = getcores(length(count_data)))

gain_loss <- lapply(gain_loss,function(x) {names(x) <- seqs;x})

gain_tables <- lapply(seqs,function(x) sapply(lapply(gain_loss,"[[",x), function(y) table(factor(sapply(y,function(z) sum(z>0)),levels = c(0,1,2,3,4)),useNA="always")))
names(gain_tables) <- seqs

loss_tables <- lapply(seqs,function(x) sapply(lapply(gain_loss,"[[",x), function(y) table(factor(sapply(y,function(z) sum(z<0)),levels = c(0,1,2,3,4)),useNA="always")))
names(loss_tables) <- seqs

## PLOTS NEEDED

sapply(seqs,function(seq){
    png(paste("AllelicGainLoss",seq,".png",sep="-"),width=8,height=8,pointsize = 8, units="in", res=300)    
    par(mfrow=c(2,1))
    par(mar=c(1,4,4,2)+0.1)
    barplot(gain_tables[[seq]][-c(1,6),],space=0.5, axisnames = F,main="Gain Counts")
    mtext("Allelic Gains and Losses by Sample (relative to Anc)",side=3,line=3)
    par(mar=c(6,4,3,2)+0.1)
    x <- barplot(loss_tables[[seq]][-c(1,6),],space=0.5, xaxt="n",main="Loss Counts")
    text(cex=1, x=x+.25, y=par("usr")[3], labels=colnames(gain_tables[[seq]]), adj=c(1,1), xpd=TRUE, srt=60)
    dev.off()
})


fishers_exact <- function(x,y,hybrid=TRUE){
    sapply(seq.int(1,ncol(x)),function(z){
        if (all(c(x[,z]) == 0) | all(c(y[,z]) == 0)) return(NA)
        else if (all(x[,z]/sum(x[,z]) == y[,z]/sum(y[,z]))) return(NA)
        else {
            m <- matrix(c(x[,z],y[,z]),nrow=2,byrow=T)
            return(fisher.test(m,hybrid=hybrid)$p.value)
        }
    })
}

##
p_values <- mclapply(count_data, function(sample,ref = Anc,hybrid=TRUE){
    return(lapply(seq.int(1,length(sample)),function(nm) fishers_exact(sample[[nm]],ref[[nm]],hybrid)))
},mc.cores = getcores(length(count_data)))

p_values <- lapply(p_values,function(x) {names(x) <- seqs;x})

### Fisher's exact test p-values by bp
sapply(seqs,function(seq){
    png(paste(seq,"p-values.png",sep="."),width=10,height=30, pointsize = 8, units="in", res=300)
    par(mfrow=c(8,1))
    tmp <- lapply(p_values,"[[",seq)
    sapply(names(tmp),function(x){
        val <- -log10(tmp[[x]])
        plot(val,col=c("black","green")[as.numeric(val>=5)+1],xlab = "Basepair position",ylab="-log10(p-value)", main= paste(x,"fisher's exact test p-values"))
        abline(h=-log10(1e-5),col="red")
        abline(h=-log10(1e-7),col="blue")
    })
    dev.off()    
})

maf <- mclapply(c(list(Anc=Anc),count_data), function(sample){
    return(lapply(sample,function(nm) apply(nm,2, function(x) max(x)/sum(x))))
},mc.cores = getcores(length(count_data)+1))


sapply(seqs,function(seq){
    png(paste(seq,"maf.png",sep="."),width=10,height=7, pointsize = 8, units="in", res=300)
    boxplot(sapply(maf,"[[",seq),main= paste(seq,"minor allele frequencies"))
    dev.off()
})

#### Major (minor) allele frequences by bp

sapply(seqs,function(seq){
    png(paste(seq,"mafbp.png",sep="."),width=10,height=30, pointsize = 8, units="in", res=300)
    par(mfrow=c(9,1))
    tmp <- lapply(maf,"[[",seq)
    sapply(names(tmp),function(x){
        val <- 1-tmp[[x]]
        plot(val,col=c("black","green")[as.numeric(val>=0.25)+1],xlab = "Basepair position",ylab="1-maf",ylim = c(0,0.6), main= paste(x,"Major Allele Frequency"))
        abline(h=-log10(1e-5),col="red")
        abline(h=-log10(1e-7),col="blue")
    })
    dev.off()    
})


library(qvalue)

qvalue.na <- function(x){
    ind <- which(!is.na(x))
    x[ind] <- qvalue(x[ind])$qvalues
    return(x)
}
tb <- sapply(p_values,sapply,function(x) sum(qvalue.na(x) <= 0.05,na.rm = T))
write.table(tb,"Statistically_diff_multtest.txt",sep="\t",row.names=T,col.names=T,quote=F)

tb2 <- sapply(p_values,sapply,function(x) sum(x <= 1e-5,na.rm = T))
write.table(tb2,"Raw_diff_multtest_1e-5.txt",sep="\t",row.names=T,col.names=T,quote=F)

### Table of differences
pmax <- 1e-3
tb3 <- do.call("rbind", lapply(names(p_values),function(sample){
    do.call("rbind",lapply(seqs, function(genome){
        ind <- which(p_values[[sample]][[genome]] <= pmax )
        if (length(ind) > 0){
            qv <- qvalue.na(p_values[[sample]][[genome]])
            data.frame("Sample"=sample,
                       "Genome"=genome,
                       "bp"=ind,
                       "p-value"= p_values[[sample]][[genome]][ind], 
                       "q-value"=qv[ind],"Alleles_Sample" = t(count_data[[sample]][[genome]][c("A","C","T","G","N","-"),ind]),
                       "Alleles_Anc" = t(Anc[[genome]][c("A","C","T","G","N","-"),ind]))
        }
    }))
})
)

write.table(tb3,"Identified_variants_p1e-3.txt",sep="\t",row.names=F,col.names=T,quote=F)


## CHEM VS BIO COMBINED
efactor <- factor(c("Bio","Bio","Bio","Bio","Chem","Chem","Chem","Chem"))
## "NC_004347.2" "pBP136_bcA"  "pMR-1L_hypo" "pMR-1S_hypo"
Bio = count_data[efactor=="Bio"]
Bio <- lapply(names(Anc),function(x) Reduce("+",lapply(Bio,"[[",x)))
names(Bio) <- names(Anc)

Chem = count_data[efactor=="Chem"]
Chem <- lapply(names(Anc),function(x) Reduce("+",lapply(Chem,"[[",x)))
names(Chem) <- names(Anc)

######################################################################################
joined_cov_tables <- lapply(seqs,function(x) sapply(lapply(list(Anc=Anc,Bio=Bio,Chem=Chem),"[[",x), function(y) summary(colSums(y))))
names(joined_cov_tables) <- seqs

as.numeric.matrix <- function(x) apply(x,2,as.numeric)
joined_monoallelic_sites_tables <- lapply(seqs,function(x) sapply(lapply(list(Anc=Anc,Bio=Bio,Chem=Chem),"[[",x), function(y) summary(apply(sweep(y,MARGIN=2,STATS = colSums(y),FUN='/'),2,function(z) any(z==1)))))
names(joined_monoallelic_sites_tables) <- seqs
joined_monoallelic_sites_tables <- lapply(joined_monoallelic_sites_tables,function(x) as.numeric.matrix(x[2:4,]))
joined_monoallelic_sites_tables <- lapply(joined_monoallelic_sites_tables,function(x) {y=round(sweep(x,MARGIN=2,STATS = colSums(x),FUN='/')*100,2);rownames(y) <- c("variable","monomorphic","NA");y})

joined_allelic_counts_tables <- lapply(seqs,function(x) sapply(lapply(list(Anc=Anc,Bio=Bio,Chem=Chem),"[[",x), function(y) rowSums(y)))
joined_allelic_counts_tables <- lapply(joined_allelic_counts_tables,function(x) round(sweep(x,MARGIN=2,STATS = colSums(x),FUN='/')*100,2))
names(joined_allelic_counts_tables) <- seqs

joined_gain_loss <- mclapply(list(Bio=Bio,Chem=Chem), function(sample,ref = Anc,hybrid=TRUE){
    return(sapply(seq.int(1,length(sample)),function(nm) allelic_gain_loss(sample[[nm]],ref[[nm]])))
},mc.cores = getcores(length(list(Bio=Bio,Chem=Chem))))

joined_gain_loss <- lapply(joined_gain_loss,function(x) {names(x) <- seqs;x})

joined_gain_tables <- lapply(seqs,function(x) sapply(lapply(joined_gain_loss,"[[",x), function(y) table(factor(sapply(y,function(z) sum(z>0)),levels = c(0,1,2,3,4)),useNA="always")))
names(joined_gain_tables) <- seqs

joined_loss_tables <- lapply(seqs,function(x) sapply(lapply(joined_gain_loss,"[[",x), function(y) table(factor(sapply(y,function(z) sum(z<0)),levels = c(0,1,2,3,4)),useNA="always")))
names(joined_loss_tables) <- seqs
######################################################################################

count_data.global = list(Biofilm=Bio,Chemostat=Chem)

p_values.global <- mclapply(count_data.global, function(sample,ref = Anc,hybrid=TRUE){
    return(sapply(seq.int(1,length(sample)),function(nm) fishers_exact(sample[[nm]],ref[[nm]],hybrid)))
},mc.cores = getcores(2))
names(p_values.global) = c("Biofilm","Chemostat")
p_values.global <- lapply(p_values.global,function(x) {names(x) <- seqs;x})

tb <- sapply(p_values.global,sapply,function(x) sum(qvalue.na(x) <= 0.05,na.rm = T))
write.table(tb,"Combined_Statistically_diff_multtest.txt",sep="\t",row.names=T,col.names=T,quote=F)

tb2 <- sapply(p_values.global,sapply,function(x) sum(x <= 1e-5,na.rm = T))
write.table(tb2,"Combined_Raw_diff_multtest_1e-5.txt",sep="\t",row.names=T,col.names=T,quote=F)

### Table of differences
pmax <- 1e-3
tb3 <- do.call("rbind", lapply(names(p_values.global),function(sample){
    do.call("rbind",lapply(seqs, function(genome){
        ind <- which(p_values.global[[sample]][[genome]] <= pmax )
        if (length(ind) > 0){
            qv <- qvalue.na(p_values.global[[sample]][[genome]])
            data.frame("Sample"=sample,
                       "Genome"=genome,
                       "bp"=ind,
                       "p-value"= p_values.global[[sample]][[genome]][ind], 
                       "q-value"=qv[ind],"Alleles_Sample" = t(count_data.global[[sample]][[genome]][c("A","C","T","G","N","-"),ind]),
                       "Alleles_Anc" = t(Anc[[genome]][c("A","C","T","G","N","-"),ind]))
        }
    }))
})
)

write.table(tb3,"Combined_Identified_variants_p1e-3.txt",sep="\t",row.names=F,col.names=T,quote=F)

# qqplots
sapply(seqs,function(x){
    png(paste("ChemVsBioF",x,"QQplots.png",sep="-"),width=8,height=8,pointsize = 8, units="in", res=300)
    tmp <- lapply(p_values.global,"[[",x)
    qqplot(-log10(na.exclude(tmp[[1]])),-log10(na.exclude(tmp[[2]])),xlab="Biofilm (-log10 scale)",ylab="Chemostat (-log10 scale)",       
           main=x)
    #abline(lm(tmp[[2]]~tmp[[1]]),col='red')
    abline(0,1,col='red')
    dev.off()
})

sapply(seqs,function(seq){
    png(paste("Combined",seq,"p-values.png",sep="."),width=20,height=10, pointsize = 8, units="in", res=300)
    par(mfrow=c(2,1))
    tmp <- lapply(p_values.global,"[[",seq)
    sapply(names(tmp),function(x){
        val <- -log10(tmp[[x]])
        plot(val,col=c("black","green")[as.numeric(val>=5)+1],xlab = "Basepair position",ylab="-log10(p-value)", main= paste(x,"fisher's exact test p-values"))
        abline(h=-log10(1e-5),col="red")
        abline(h=-log10(1e-7),col="blue")
    })
    dev.off()    
})

