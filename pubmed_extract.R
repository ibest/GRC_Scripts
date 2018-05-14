
library("RCurl")
library("XML")
library("plyr")
library("ggplot2")
library("directlabels")


########################
# Download PubMed Data #
########################
PubMedTrend <- function(query, yrStart=1950, yrMax=2009) {
  
  ### Some error checking ###
  if (is.numeric(yrStart) == FALSE || is.numeric(yrMax) == FALSE) stop("One of the year values is not numeric")
  if (yrStart < 1800) stop(paste("Sure you want to look for hits from the 18th century (yrStart = " ,yrStart, ")?\n", sep=""))
  this.year <- Sys.time()
  this.year <- as.integer(format(this.year, "%Y"))
  if (yrMax > this.year) stop(paste("Are you from the future? Please check your year interval; yrMax =",yrMax,"\n"))
  if (yrMax < yrStart) stop("yrMax is smaller than yrMin!")
  query_new <- c(query,"total_counts"="")
  ### Start main search function ###
  getCount <- function(query.term) {
    # convert spaces to '+'
    query.gsub <- gsub(" ", "+", query.term)
    # convert some characters to brower friendly text (better to be safe than sorry)
    query.gsub <- gsub('"', "%22", query.gsub)
    query.gsub <- gsub("\\[", "%5B", query.gsub)
    query.gsub <- gsub("\\]", "%5D", query.gsub)
    # add progressbar
    pb <- txtProgressBar(min = yrStart, max = yrMax, style = 3)
    # create empty data frame
    df <- data.frame(NULL)
    cat("Searching for: ", query.term,"\n")
    
    # Start retrieval loop
    for(i in yrStart:yrMax) {
      # tell progressbar how it's going
      setTxtProgressBar(pb, i)
      # add publication date [dp] to query
      query.parsed <- paste(query.gsub, "+AND+",i, "%5Bppdat%5D", sep="")
      # Get XML with number of hits for query.parsed
      success = FALSE
      while(!success){
        pub.esearch <- try(getURL(paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&term=", 
                                  query.parsed, sep = "")))
        if(inherits(pub.esearch,"try-error")){
          print("NCBI error")
          Sys.sleep(1.0)
        }else{
          
          success=TRUE
        }
      }
      # Parse XML
      pub.esearch <- xmlTreeParse(pub.esearch, asText = TRUE)
      # Get number of hits from XML
      pub.count <- as.numeric(xmlValue(pub.esearch[["doc"]][["eSearchResult"]][["Count"]]))
      
      # Don't add anything if count is 0
      if (pub.count != 0) df <- rbind(df, data.frame("year" = i, "count" = pub.count))
      # Wait 0.5 sec
      Sys.sleep(1.0)
    }
    # close progressbar
    close(pb)
    return(df)
  } 
  # Run getCount() for all query terms
  df <- ldply(query_new, getCount)
  
  ### Calculate relative frequencies ###
  # load file with pubmed total citations from 1947-2009
#  load(file="total_table")
  # match year
  total.table <- df[which(df$.id == "total_counts"),]
  df <- df[which(df$.id %in% names(query)),]
  match <- match(df$year, total.table$year)
  # add total count
  df$total_count <- total.table$count[match]
  # compute relative count times 10 000, i.e. show number of matches per 1 million PubMed citations
  df$relative <- (df$count / df$total_count) * 10000
  cat("\nAll done!")
  return(df)
}

#######################
### Show total hits ###
#######################
PubTotalHits <- function(args=FALSE) {
  # Get column total for query 'x'
  GetCount <- function(x) {
    df <- data.frame("search_name" = x, "total_hits" = colSums(df[df$.id == x,][3]))
  }
  # Index all query names
  query.index <- unique(df$.id)
  # Use GetCount() for every term in 'query.index' and return as data.frame
  df <- ldply(query.index, GetCount)
  # if argument is 'query' add full query instead of query name. 
  # if there is no argument specified both name and query will be shown
  if (args == "query" || args == FALSE) {
    # remove names
    names(query) <- NULL
    # add queries to df
    df <- cbind(df, "query" = query)
    # reorder columns
    df <- df[,c(1,3,2)]
    # remove 'names' if we only want queries
    if (args == "query") df <- df[-1]
  }
  return(df)
}

query <- c("seq"= "rna seq", microarray="expression microarray")
#query <- c("cbt"= "cognitive behav* psychotherap*[tiab] OR cognitive behav* therap*[tiab]", 
#           "pdt" = "psychodynamic therap*[tiab] OR psychodynamic psychotherap*[tiab]",
#           "psychoanalytic" = "psychoanalytic therap*[tiab] OR psychoanalytic psychoterap*[tiab]", 
#           "ssri" = "selective serotonin reuptake inhibitor*[tiab]",
#           "mindfulness" = "mindfulness[tiab]")
pubmedDF <- PubMedTrend(query, yrStart=2000, yrMax=2013)





### AREA PLOT ###
ggplot(df, aes(year, relative, group=.id, fill=.id)) + 
  geom_area() +
  opts(title=paste("Area Plot of PubMed Publications per Year\nfor", paste(names(query), collapse = ", "))) +
  xlab("year") +
  ylab("Publications per 1 million PubMed articles") +
  scale_fill_brewer()

### LINE PLOTS ###

# RAW
ggplot(df, aes(year, relative, group=.id, color=.id)) + 
  geom_line(show_guide=F) + 
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = paste("Pubmed hits for", paste(names(query), collapse = ", ")))

pdf("Pubmed_hits.pdf",height=6,width=8,pointsize=12)
p <- ggplot(df, aes(year, relative, group=.id, color=.id)) + 
  geom_line(alpha = I(7/10), show_guide=F,size=2) +
#  stat_smooth(size=2, span=0.3, se=F, show_guide=F) + 
  xlab("Publication year") +
  ylab("Publications per 1 million PubMed articles") +
  opts(title = "Pubmed hits for ('rna seq' and 'expression microarray')") +
  xlim(2000,2018)
direct.label(p, list("last.bumpup", hjust = -0.1, vjust = 0))
dev.off()

sra <- read.table("../../Documents/Presentations/ISAFG-2013/sra_stat.csv.txt",sep=",",header=T,as.is=T)
sra$Tb <- sra$bases/1e12
sra$Tboe <- sra$open_access_bases/1e12
sra$Date <- as.Date(sra$date,"%m/%d/%Y")

bases <- data.frame(sra[,c("Date","Tb")],group="bases")
oabases <- data.frame(sra[,c("Date","Tboe")],group="open access bases")
names(oabases) <- names(bases)                    
sraTB <- rbind(bases,oabases)

pdf("SRA_hits.pdf",height=6,width=8,pointsize=12)
# SMOOTHED
p <- ggplot(sraTB, aes(Date, Tb, group=group, color=group)) + 
  geom_line(alpha = I(7/10), show_guide=F,size=2) +
  xlab("Submission year") +
  ylab("Terabytes") +
  opts(title = "Terabytes of sequence uploaded to the SRA per year") +
  xlim(sraTB[1,"Date"],sraTB[1631,"Date"]+2*365)
direct.label(p, list("last.bumpup", hjust = -0.1, vjust = 0))

dev.off()







######
#source('pubmed_trend.r')
pubmed_trend <- function(search.str = 'Sex+Characteristics[mh] AND Pain[mh]', year.span=1970:2011) {
  require(XML)
  require(RCurl)
  
  results <- NULL
  tmpf <- "./tempfile.xml"
  ## clean before
  system(paste("rm", tmpf))
  
  for(i in year.span){
    queryString <- paste(search.str, ' AND ', i, '[dp]', sep="")
    print(paste('queryString:', queryString))
    sysString <- paste('./pubmed_trend.pl "', queryString,'"', sep="")
    system(sysString)
    
    xml <- xmlTreeParse(tmpf, useInternalNodes=TRUE)
    pubTerm <- as.numeric(xmlValue(getNodeSet(xml, "//Count")[[1]]))
    print(paste("#______num pub for",i,":",pubTerm))
    rm(xml)
    results <- append(results, pubTerm)
    ## avoid being kicked out!
    Sys.sleep(1)
  }
  names(results) <- year.span
  ## clean after
  system(paste("rm", tmpf))
  
  return(results)
}

sex.pub <- pubmed_trend(search.str = 'Sex+Characteristics[mh] AND Pain[mh]', year.span=1970:2011)
analgesic.pub <- pubmed_trend(search.str = 'Sex+Characteristics[mh] AND Analgesics[mh]', year.span=1970:2011)

#source('plot_bar.r')
plot_bar <- function(x=sex.pub, linecol="royalblue", cols, addArg=TRUE) {
  bp <- barplot(x, col=cols, add=addArg)
  fit <- stats::lowess(x, f=1/3)
  lines(x=bp, fit$y, col=linecol, lwd=3)  
}
library("RColorBrewer")

pdf(file='sex_pain.pdf', height=8, width=8)
par(las=1)
colorfunction = colorRampPalette(brewer.pal(9, "Reds"))
mycolors = colorfunction(length(sex.pub))
plot_bar(x=sex.pub, linecol="#525252", cols=mycolors, addArg=FALSE)

colorfunction = colorRampPalette(brewer.pal(9, "Blues"))
mycolors = colorfunction(length(analgesic.pub))
plot_bar(x=analgesic.pub, linecol='black', cols=mycolors, addArg=TRUE)
title('Number of publication per year')
legend('topleft',
       legend=c('Sex and Pain', 'Sex and Analgesics'),
       fill=c("red", "blue"),
       bty="n",
       cex=1.1
)
dev.off()