library(VariantAnnotation)
library(Rsamtools)
### Check out replaceLetterAt
snpfile <-"CP1.A.inornata.vcf"
bamfile <- "CP1.A.inornata.bam"
outfile <- "CP1.A.inornata.fasta"

tb <- read.table("04-Bowtie/SummarySample2Targets.txt")
tb$sample <- rownames(tb)

dir.create("05-ConsensusSequences")
apply(tb,1,function(x) {
  sample <- x[4]
  species <- x[3]
  dir.create(file.path("05-ConsensusSequences",sample))
  generateConsensusSequence(file.path("04-Bowtie",sample,paste(sample,species,"vcf",sep=".")),
                            file.path("04-Bowtie",sample,paste(sample,species,"bam",sep=".")),
                            file.path("05-ConsensusSequences",sample,paste(sample,species,"fasta",sep=".")),
                            sample)
})

generateConsensusSequence <- function(snpfile,bamfile,outfile, sample, min_cov_to_call = 2, min_cov_het_call = 10){
  ### read in vcf
  hdr <- scanVcfHeader(snpfile)
  vcf <- readVcf(snpfile,"A.i")
  
  alt.alleles <- (geno(vcf)$GT != "0/1" & geno(vcf)$GT != "1/1")
  ## pull out alleles with are heterozygous other (ie 2 alternates)
  vcf.alt.alleles <- vcf[alt.alleles]
  
  vcf <- vcf[!alt.alleles] 
  ## pull out indels
  indel <- if (class(alt(vcf)) == "DNAStringSetList" ){
              width(unlist(alt(vcf))) != width(ref(vcf))
              } else width(alt(vcf)) != width(ref(vcf))
  
  vcf.indel <- vcf[indel]
  vcf <- vcf[!indel]
  
  #geno(vcf)$GT
  
  #(1) Load the BAM file into a GAlignments object:
  
  param <- ScanBamParam(what=c("seq","qname"))
  gal <- readGAlignmentsFromBam(bamfile, param=param)
  qseq <- mcols(gal)$seq  # the query sequences
  
  #(2) Use sequenceLayer() to "lay" the query sequences on the reference
  #space. This will remove the parts from the query sequences that
  #correspond to insertions and soft clipping, and it will fill them
  #with - where deletions and/or skipped regions occurred:
    
  qseq_on_ref <- sequenceLayer(qseq, cigar(gal),
                                 from="query", to="reference")
  
  #(3) Compute 1 consensus matrix per chromosome:
  
  qseq_on_ref_by_chrom <- splitAsList(qseq_on_ref, seqnames(gal))
  qseq_pos_by_chrom <- splitAsList(start(gal), seqnames(gal))
  
  cm_by_chrom <- lapply(names(qseq_pos_by_chrom),
                        function(seqname)
                          consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                                          as.prob=TRUE,
                                          shift=qseq_pos_by_chrom[[seqname]]-1,
                                          width=seqlengths(gal)[[seqname]]))
  names(cm_by_chrom) <- names(qseq_pos_by_chrom)
  
  cov_by_chrom <- lapply(names(qseq_pos_by_chrom),
                        function(seqname)
                          colSums(consensusMatrix(qseq_on_ref_by_chrom[[seqname]],
                                          as.prob=FALSE,
                                          shift=qseq_pos_by_chrom[[seqname]]-1,
                                          width=seqlengths(gal)[[seqname]])))
  names(cov_by_chrom) <- names(qseq_pos_by_chrom)
  
  
  ## 'cm_by_chrom' is a list of consensus matrices. Each matrix
  ## has 17 rows (1 per letter in the DNA alphabet) and 1 column
  ## per chromosome position.
  
  #(4) Compute the consensus string from each consensus matrix. We'll
  #     put "+" ins the strings wherever there is no coverage for that
  #     position, and "N" where there is coverage but no consensus.
  
  cs_by_chrom <- lapply(names(cm_by_chrom),
    function(seqname) {
      ## Because consensusString() doesn't like consensus
      ## matrices with columns that contain only zeroes (and
      ## you will have columns like for chromosome positions
      ## that don't receive any coverage), we need to "fix"
      ## 'cm' first.
      vcf_idx <-which(seqnames(vcf) == seqname)
      # min coverage of 4
      idx <- cov_by_chrom[[seqname]] < min_cov_to_call
      cm_by_chrom[[seqname]]["+", idx] <- 3
      # min Het coverage, call N otherwise    
      noCallHet <- geno(vcf)$GT[vcf_idx,] == "0/1" & info(vcf)$DP[vcf_idx] < min_cov_het_call
      cm_by_chrom[[seqname]]["N",  start(vcf)[vcf_idx][noCallHet]] <- 2
  
      cm_string <- rownames(cm_by_chrom[[seqname]])[apply(cm_by_chrom[[seqname]],2,which.max)]
    
      # inject het calls
      index_het <- (geno(vcf)$GT[vcf_idx,] == "0/1" & info(vcf)$DP[vcf_idx] >= min_cov_het_call)
      start_het <- start(vcf)[vcf_idx][index_het]
      if(class(alt(vcf[vcf_idx])) == "DNAStringSetList"){
        ref <- as.character(ref(vcf[vcf_idx][index_het]))
        alt <- as.character(unlist(alt(vcf[vcf_idx][index_het])))
        het_code <- names(IUPAC_CODE_MAP)[match(apply(cbind(ref,alt),1, function(x) paste(sort(x),collapse="")),IUPAC_CODE_MAP)]
      } else {
        ref <- as.character(ref(vcf[vcf_idx][index_het]))
        alt <- as.character(alt(vcf[vcf_idx][index_het]))
        het_code <- names(IUPAC_CODE_MAP)[match(apply(cbind(ref,alt),1, function(x) paste(sort(x),collapse="")),IUPAC_CODE_MAP)]
      }
      
      cm_string[start_het] <- het_code
        
      DNAString(paste(cm_string,collapse=""))
      
      ## check that homo other calls are called in the cm_by_chrom
  #    start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]
  #    cm_string[start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]]
  #    cbind(as.character( unlist( alt(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"])) , cm_string[start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]])
  #    DNAString(consensusString(cm, ambiguityMap="N"))
    })
  
  #(5) Write 'cs_by_chrom' to a FASTA file:
  cm_seqs <- DNAStringSet(cs_by_chrom)
  nameSplit <- strsplit(names(cm_by_chrom),split="_:_")
  names(cm_seqs) <- paste(sample,sapply(nameSplit,"[[",1L),sapply(nameSplit,"[[",2L),sep="|")
  writeXStringSet(cm_seqs, outfile)
}


min_cov_to_call = 1
cs_by_chrom <- lapply(names(cm_by_chrom),
  function(seqname) {
    ## Because consensusString() doesn't like consensus
    ## matrices with columns that contain only zeroes (and
    ## you will have columns like for chromosome positions
    ## that don't receive any coverage), we need to "fix"
    ## 'cm' first.
    idx <- cov_by_chrom[[seqname]] < min_cov_to_call
    cm_by_chrom[[seqname]]["N", idx] <- 3
#    cm_by_chrom[[seqname]]["+", idx] <- 3
    
    cm_string <- rownames(cm_by_chrom[[seqname]])[apply(cm_by_chrom[[seqname]],2,which.max)]
    
    DNAString(paste(cm_string,collapse=""))
    
    ## check that homo other calls are called in the cm_by_chrom
    #    start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]
    #    cm_string[start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]]
    #    cbind(as.character( unlist( alt(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"])) , cm_string[start(vcf)[vcf_idx][geno(vcf)$GT[vcf_idx,] == "1/1"]])
    #    DNAString(consensusString(cm, ambiguityMap="N"))
})

