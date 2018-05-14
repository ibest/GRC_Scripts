
setwd("/mnt/home/uirig/GRC_projects/Forney/Bifidobacteria_genomes")
samples <- read.table("samples.txt",sep="\t",header=T,as.is=T)

newbler_folder <- "03-Assemblies"
anno_folder <- "06-Annotation"

for (i in samples$SAMPLE_ID){
  dir.create(file.path(anno_folder,i))
  call <- paste("prodigal -f gff -a", file.path(anno_folder,i,paste(i,"_Contigs.faa",sep="")),"-i",file.path(newbler_folder,i,"454LargeContigs.fna"), "-o", file.path(anno_folder,i,paste(i,"_Contigs.gff",sep="")),"-s", file.path(anno_folder,i,paste(i,"_scores.txt",sep="")))
  system(call)
}

res <- system(paste("grep -c 'ID='",file.path(anno_folder,"*","*.gff")),intern=T)
splitres <- strsplit(res,split="/|:")

genes <- data.frame(Sample=sapply(splitres,"[[",2L),numGenes=as.numeric(sapply(splitres,"[[",4L)))

# Sample numGenes
# 1     M39     2185
# 2      R1     2157
# 3     R21       63
# 4     R24     2201
# 5     R36     2160
# 6     R37     2245
# 7     R38     2250
# 8      R3     2133
# 9     R56     2192
# 10    R68     2385
# 11    R78     2107

#Rsync to mite for blannotator annotation
call <- paste("rsync -avz -e ssh",file.path(anno_folder,"*","*.faa"), "msettles@bioinfo-mite.ibest.uidaho.edu:~/projects/Forney/Bifidobacteria_genomes/06-Annotation/.")
system(call)

# Annotate genes with BLANNOTATOR:
### on mite  ~/projects/Forney/Bifidobacteria_genomes/06-Annotation
setwd("~/projects/Forney/Bifidobacteria_genomes/06-Annotation")
libs <- dir()
samples <- sapply(strsplit(libs,split="_"),"[[",1L)

for (i in samples){
  dir.create(i)
  file.rename(paste(i,"_Contigs.faa",sep=""),file.path(".",i,paste(i,"_Contigs.faa",sep="")))
  #BLASTDB=/grc/src/BLANNOTATOR/database/
  call <- paste("/grc/src/BLANNOTATOR/RUN_PACKAGE/blannotator.pm --file=",file.path(i,paste(i,"_Contigs.faa",sep="")), " --goa_gene=/grc/src/BLANNOTATOR/database/goa_genes.txt --obo=/grc/src/BLANNOTATOR/database/gene_ontology.1_0.obo --asso=/grc/src/BLANNOTATOR/database/gene_association.goa_uniprot_all_short --database=/grc/src/BLANNOTATOR/database/uniprot_trembl.fasta --run_blast=1 --html=1 --out=",file.path(i,paste(i,"_annotations",sep=""))," --ncpu=14",sep="")
  system(call)
}


#Rsync back to central-storage
call <- paste("rsync -avz -e ssh * msettles@petunia:/mnt/home/uirig/GRC_projects/Forney/Bifidobacteria_genomes/07-Blannotator")
system(call)

