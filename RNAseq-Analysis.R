suppressPackageStartupMessages(library("parallel"))

samples <- read.table("samples.txt",header=T,as.is=T)

genome_gtf = "/mnt/home/uirig/user_data/genomes/Danio_Rerio/Danio_rerio.Zv9.75.gtf"

# Running TopHat/Bowtie
# bowtie_ref = "/mnt/home/uirig/user_data/genomes/Danio_Rerio/Danio_rerio.Zv9.75"
#
# dir.create("03-Tophat",showWarnings=F)
# tophat_out <- mclapply(samples$SAMPLE_ID, function(sample){
#   try({
#     call <- paste("tophat2 -G",genmoe_gtf, "-p 10 -o ", file.path("03-Tophat",sample) ,bowtie_ref, file.path("02-Cleaned",sample,paste(sample,"_merged_SE.fastq",sep="")))
#     system(call)
#   })
# },mc.cores=6)

# Running BWA
bwa_ref = "/mnt/home/uirig/user_data/genomes/Danio_Rerio/Danio_rerio.Zv9.75.dna.toplevel.fa"

folder = "03-BWA"
dir.create(folder,showWarnings=F)

bwa_out <- mclapply(samples$SAMPLE_ID, function(sample){
    dir.create(file.path(folder,sample))
    try({
        call <- paste("bwa mem -t 10" ,bwa_ref, file.path("02-Cleaned",sample,paste(sample,"_merged_SE.fastq",sep="")), ">", file.path(folder,sample,paste(sample,".sam",sep="")
        ))
        system(call)
    })
},mc.cores=6)

### Generate sorted Bam/Sam files
sam_out <- mclapply(samples$SAMPLE_ID, function(sample,folder="03-BWA"){
    if (!file.exists(file.path(folder,sample,paste(sample,".bam",sep="")))){
        call <- paste("samtools view -bS -o", file.path(folder,sample,paste(sample,".bam",sep="")), file.path(folder,sample,paste(sample,".sam",sep="")))
        system(call)
    }
    # sort by name, convert to SAM for htseq-count
    call <- paste("samtools sort -n",file.path(folder,sample,paste(sample,".bam",sep="")),file.path(folder,sample,paste(sample,"_sn",sep="")))
    system(call)
    call <- paste("samtools view -o",file.path(folder,sample,paste(sample,"_sn.sam",sep="")),file.path(folder,sample,paste(sample,"_sn.bam",sep="")))
    system(call)
    # sort by position and index for IGV
    call <- paste("samtools sort",file.path(folder,sample,paste(sample,".bam",sep="")),file.path(folder,sample,paste(sample,"_s",sep="")))
    system(call)
    call <- paste("samtools view -o",file.path(folder,sample,paste(sample,"_s.sam",sep="")),file.path(folder,sample,paste(sample,"_s.bam",sep="")))
    system(call)
},mc.cores=6)


# Running htseq-count
script = "~/opt/src/HTSeq-0.6.1/build/scripts-2.7/htseq-count"
sam_out <- mclapply(samples$SAMPLE_ID, function(sample,folder="03-BWA"){
    call = paste(script, "-s no -a 10", file.path(folder,sample,paste(sample,"_sn.sam",sep="")) , genome_gtf,">", file.path(folder,sample,paste(sample,"_sn.count",sep="")))
    system(call)
},mc.cores=10)

# (i) Load the edgeR package and use the utility function, readDGE, to read in the COUNT files created from htseq-count:
library("edgeR")
samples <- read.table("samples.txt",header=T,as.is=T)

mfolder = "03-BWA"
folder = "04-HTseqCounts"

countf <- file.path(mfolder,samples$SAMPLE_ID,paste(samples$SAMPLE_ID,"_sn.count",sep=""))
counts = readDGE(countf)$counts

#(ii) Filter weakly expressed and noninformative (e.g., non-aligned) features using a command like:
noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual",
                                "not_aligned","alignment_not_unique")
cpms = cpm(counts)
keep = rowSums(cpms >1) >=3 & !noint
counts = counts[keep,]
# In edgeR, it is recommended to remove features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates (here, n = 3 for the knockdown group).
# (iii) Visualize and inspect the count table as follows:
colnames(counts) = samples$SAMPLE_ID
head( counts[,order(samples$condition,samples$sample)], 5 )
# (iv) Create a DGEList object (edgeR’s container for RNA-seq count data), as follows:
d = DGEList(counts=counts, group=samples$condition)
# (v) Estimate normalization factors using:
d = calcNormFactors(d)
# (vi) Inspect the relationships between samples using a multidimensional scaling (MDS) plot, as shown in Figure 4:
pdf("plots.pdf")
plotMDS(d, labels=samples$SAMPLE_ID,
        col = c("darkgreen","blue")[factor(samples$condition)],cex=0.6, main="MDS")
dev.off()

# (vii) Estimate tagwise dispersion (simple design) using:
d = estimateCommonDisp(d)
d = estimateTagwiseDisp(d)
# (viii) Create a visual representation of the mean-variance relationship using the plotMeanVar (Fig. 5a) and plotBCV (Fig. 5b) functions, as follows:
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE,main="MeanVar")
plotBCV(d,main="BCV")

# (B) edger—complex design
# (i) Follow Step 14A(i–vi).
# (ii) Create a design matrix (see ‘Experimental design’ for further details) to specify the factors that are expected to affect
# expression levels:
design = model.matrix( ~ sample + condition, samples)
design
# (iii) Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR)-adjusted likelihood7,53, as follows:
d2 = estimateGLMTrendedDisp(d, design)
d2 = estimateGLMTagwiseDisp(d2, design)
# (iv) Given the design matrix and dispersion estimates, fit a GLM to each feature:
f = glmFit(d2, design)
# (v) Perform a likelihood ratio test, specifying the difference of interest (here, knockdown versus control, which corresponds to the third column of the above design matrix):
de = glmLRT(f, coef=3)
# (vi) Use the topTags function to present a tabular summary of the differential expression statistics (note that topTags operates on the output of exactTest or glmLRT, but onl
y the latter is shown here):
    tt = topTags(de, n=nrow(d))
head(tt$table)
# (vii) Inspect the depth-adjusted reads per million for some of the top differentially expressed genes:
nc = cpm(d, normalized.lib.sizes=TRUE)
rn = rownames(tt$table)
head(nc[rn,order(samples$condition)],5)
# (viii) Create a graphical summary, such as an M (log-fold change) versus A (log-average expression) plot54, here showing the genes selected as differentially expressed (with
a 5% false discovery rate; Fig. 6):
    deg = rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg,main="Smear")
# (ix) Save the result table as a CSV file (alternative formats are possible) as follows:
write.csv(tt$table, file="toptags_edgeR.csv")
dev.off()
### attach annotation
annof <- "Zv9.ensemble.Mar2014.txt"
anno <- read.table(annof,sep="\t",header=T,comment.char="",quote="",as.is=T)
write.csv(data.frame(tt$table,anno[match(rownames(tt$table),anno$Ensembl.Gene.ID),]),file="toptags_edgeR.annotated.csv")