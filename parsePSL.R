################ Function to read and manipulate BLAT PSL and PSLx files
### 
require(bitops)

### PSLx header names
colnamesPSLX <- c(
  "matches",           # Number of bases that match that aren't repeats
  "misMatches",        # Number of bases that don't match
  "repMatches",        # Number of bases that match but are part of repeats
  "nCount",            # Number of 'N' bases
  "qNumInsert",        # Number of inserts in query
  "qBaseInsert",       # Number of bases inserted in query
  "tNumInsert",        # Number of inserts in target
  "tBaseInsert",       # Number of bases inserted in target
  "strand",            # + or - for query strand, optionally followed by + or – for target strand
  "qName",             # Query sequence name
  "qSize",             # Query sequence size
  "qStart",            # Alignment start position in query
  "qEnd",              # Alignment end position in query
  "tName",             # Target sequence name
  "tSize",             # Target sequence size
  "tStart",            # Alignment start position in target
  "tEnd",              # Alignment end position in target
  "blockCount",        # Number of blocks in alignment. A block contains no gaps.
  "blockSizes",        # Size of each block in a comma separated list
  "qStarts",           # Start of each block in query in a comma separated list
  "tStarts",           # Start of each block in target in a comma separated list
  "qseq",					 # Query Sequence
  "tseq"					 # Target Sequence
)
### PSL is a subset of PSLx
colnamesPSL <- colnamesPSLX[-c(22,23)]

## column types
colCX <- rep("numeric",23)
colCX[c(9,10,14,19,20,21,22,23)] <- "character"
colC <- colCX[-c(22,23)]


#######################################################################
## Read in PSL and PSLx files
readPSL <- function(filename,pslHeader=TRUE,removeSeqs = TRUE,skip=5, quote="", ...){    
    pslx <- as.logical(regexpr("pslx$",filename)>0)
    skip=ifelse(pslHeader,skip,0)
    if(pslx){ 
    	colN=colnamesPSLX;colC=colCX
    }else{colN=colnamesPSL;colC=colC}
    tab <- read.table(filename, quote=quote, skip=skip, header=FALSE,col.names=colN, colClasses=colC, ...)
    if(removeSeqs & pslx)  tab <- tab[,-c(22,23)]
    tab
}

#######################################################################
## determine if the psl is in protein space
pslIsProtein <- function(psl)
  ##is psl a protein psl (are it's blockSizes and scores in protein space) 
{
  psl <- psl[1,]
  lastBlock = psl$blockCount;
  tStarts <- as.numeric(unlist(strsplit(as.character(psl$tStarts), ',')))
  blockSizes <- as.numeric(unlist(strsplit(as.character(psl$blockSizes), ',')))
  return  (((psl$strand == '+' ) &&
              (psl$tEnd == tStarts[lastBlock] + 
                 3*blockSizes[lastBlock])) ||
             ((psl$strand == '-') &&
                (psl$tStart == (psl$tEnd-(tStarts[lastBlock] + 
                                            3*blockSizes[lastBlock])))));
}


#######################################################################
## Calculate coverage
pslCoverage <- function(psl) return(((psl$matches + psl$misMatches + psl$repMatches + psl$nCount)/psl$qSize) * 100)
  
#######################################################################
## Calculate identity
pslIdentity <- function(psl) return(100.0 - pslCalcMilliBad(psl) * 0.1)

# Here is the source for pslCalcMilliBad:
#The complexity in milliBad arises primarily from how it handles inserts. 
#Ignoring the inserts, the calculation is simply mismatches expressed as parts per thousand. 
#However, the algorithm factors in insertion penalties as well,
#which are relatively weak compared to say blasts but still present. 
#When huge inserts are allowed (which is necessary to accommodate introns), 
#it is typically necessary to resort to logarithms like this calculation does.
pslCalcMilliBad <- function( psl, isMrna=FALSE)
{
  # Calculate badness in parts per thousand. 
  sizeMul <- ifelse(pslIsProtein(psl), 3, 1);

    qAliSize <- sizeMul * (psl$qEnd - psl$qStart);
    tAliSize <- psl$tEnd - psl$tStart
    aliSize = pmin(qAliSize, tAliSize);
    aliSize <- as.numeric(!(aliSize<= 0))

    sizeDif = qAliSize - tAliSize;
    sizeDif[sizeDif < 0] <- ifelse(isMrna,
            0,
            -sizeDif[sizeDif < 0])

    insertFactor = psl$qNumInsert;
    if (!isMrna)
        insertFactor = insertFactor + psl$tNumInsert;

    total = (sizeMul * (psl$match + psl$repMatch + psl$misMatch));
    if (any(total == 0)) stop("total is equal to zero")
    milliBad <- (1000 * (psl$misMatch*sizeMul + insertFactor + 
            round(3*log(1+sizeDif)))) / total
             
    return(aliSize*milliBad)
}

#######################################################################
## Calculate Score
pslScore <- function( psl )
{
    sizeMul <- ifelse(pslIsProtein(psl), 3, 1);
    return ( sizeMul * (psl$match + bitShiftR(psl$repMatch,1)) -
             sizeMul * psl$misMatch - psl$qNumInsert - psl$tNumInsert )
}

