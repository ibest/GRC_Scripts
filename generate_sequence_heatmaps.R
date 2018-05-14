library(Biostrings)
options(scipen=999)

samples = read.table("../samples.txt",sep="\t",header=T)

dname = sub("./","",dirname(dir(pattern="A.baumannii_ATCC_17978.bam$",recursive=T,full.names=T)))

f <- dir(pattern="A.baumannii_ATCC_17978.bam$",recursive=T,full.names=T)
command <- paste("samtools depth",paste(f,collapse=" "),sep=" ")
depth <- system(command,intern=TRUE)
dp <- strsplit(depth,split="\t")
dpm <- apply(sapply(dp,"[",3:(2+length(f))),1,as.numeric)
pos <- as.numeric(sapply(dp,"[[",2L))
target <- sapply(dp,"[[",1L)

fa <- readDNAStringSet("../References/DOD_ALL_targets.fasta")
newpos <- unlist(lapply(width(fa),seq.int))
newtarget <- rep(sapply(strsplit(names(fa),split=" "),"[[",1L),times=width(fa))
trimmmed_mean <- apply(dpm,2,mean,trim=0.2)
dpm_f <- sweep(dpm,MARGIN=2,STATS=trimmmed_mean,FUN="/")

fdpm <- matrix(0,ncol=ncol(dpm_f),nrow=length(newpos))
fdpm[match(paste(pos,target),paste(newpos,newtarget)),] <- dpm_f

sfdpm <- split(fdpm,newtarget)
sfdpm <- by(fdpm,newtarget,as.matrix)


library(gplots)
library(Heatplus)
library(RColorBrewer)

scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
scaleyellowred[1] = "grey"



plot_htmap <- function(data, name){
    clabel = rep("",nrow(data))
    clabel[seq(1000,nrow(data),by=1000)] = as.character(seq(1000,nrow(data),by=1000))
    sname = paste(dname," ",colSums(round(data) > 0),"bp"," ",round(colMeans(data),2),"rdp",sep="")
    #row.clus <- hclust(dist(t(data)), "aver")
    
    hval = max(7,nrow(samples)*0.25+4)
    wval = max(10.5,nrow(data)/100000*1.5+4)
    png(paste(name,"heatmap.png",sep="_"),width=wval,height=hval, units="in",res=300)
    heatmap.2(t(data),
              #       Rowv = as.dendrogram(row.clus),
              Rowv = FALSE,
              Colv = FALSE,
              dendrogram="none",
              col = scaleyellowred,
              scale="none",
              labCol=clabel,
              labRow=sname,
              cexRow =  0.8,
              margins = c(5, 10),
              key.title="relative depth",
              keysize=0.5,
              #                lmat = rbind( c(4, 3), c(2,1) ), lhei=c(1, 10),lwid=c(1,10),
              lmat = rbind( c(0, 4,0,0,0,0,0,0,0,0), c(2,1,1,1,1,1,1,1,1,1),c(0,3,3,3,3,3,3,3,3,3) ), lhei=c(1, 10, 0.5),lwid=c(0.5,2,1,1,1,1,1,1,1,1),
              trace="none",density.info="none",main=name)
    dev.off()
}

pB10 <- sfdpm$"NC_004840.1"
plot_htmap(pB10,"pB10")

pAB1 <- sfdpm$"NC_009083.1"
plot_htmap(pAB1,"pAB1")

pAB2 <- sfdpm$"NC_009084.1"
plot_htmap(pAB2,"pAB2")

genome <- sfdpm$"NC_009085.1"
plot_htmap(genome,"A.baumannii_ATCC_17978 Genome")
