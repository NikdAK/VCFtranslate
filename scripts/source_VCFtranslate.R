# Load libraries
print("Loading packages")
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(systemPipeR))
suppressPackageStartupMessages(library(splitstackshape))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(here))

# some package functions are overwritten but we need exactly those ones (which are not selected by :: in code) 
getSequence <- biomaRt::getSequence
translate <- Biostrings::translate
between <- data.table::between

registerDoParallel(cores = 2)

# Init
organism <- "human"
ensemblVersion <- 92

# set versions
release.version <- "VCFtranslate-1.0"

#set date
analysis.date <- format(Sys.time(), "%d-%m-%Y-%H:%M")

# combine to header element which will be included in results file
versionControlTag <- paste0(release.version,",","ensembl-", ensemblVersion, ",", "analysis-",analysis.date)

# check if output directory exists, else create
resultsDir <- "results_out/"
ifelse(!dir.exists(file.path(resultsDir)), dir.create(file.path(resultsDir)), FALSE)

# Do pre settings
if (organism == "human"){
  db <- useEnsembl("ensembl",dataset="hsapiens_gene_ensembl", version=ensemblVersion)
} else {
  stop("organism not implemented and tested yet")
}

print("Loading functions")

# MAIN FUNCTION
VCFtranslate <- function(vcf.DT, vcf.file, seqType.input){
  
  if(!dir.exists(resultsDir)){stop("Folder results_out not found in workdir")}
  
  print(paste0("Using ensembl version ", ensemblVersion))

  # Init lookups
  print("Initializing lookup tables from biomart")
  tscrIDs <- vcf.DT[, unique(FEATUREID.noversion)]
  
  # for large datasets we need to split this into smaller batches (n=500), since biomart will break due to a timeout 
  tscrIDs.chunks <- split(tscrIDs, ceiling(seq_along(tscrIDs) / 500))
  
  # get all exons for all transcripts
  print("(1/4) Downloading transcript information")
  exonDT <<- data.table(getBM(attributes=c('ensembl_transcript_id','rank', 'ensembl_exon_id', 'exon_chrom_start', 'exon_chrom_end', 'strand',
                                          'transcript_start', 'transcript_end', 'transcript_biotype', 'gene_biotype'), 
                             filters = 'ensembl_transcript_id', values=tscrIDs, mart=db))
  setkey(exonDT, ensembl_transcript_id)

  # get all sequences for all transcripts (coding, cdna)
  print("(2/4) Downloading transcript sequences")
  # sequenceDT <<- data.table(getSequence(id=tscrIDs, type="ensembl_transcript_id", seqType=seqType.input, mart = db))
  sequenceDT <<- rbindlist(lapply(1:length(tscrIDs.chunks), function(x){
    seq.tmp <- data.table(getSequence(id=tscrIDs.chunks[[x]], type="ensembl_transcript_id", seqType=seqType.input, mart = db))
    seq.tmp
  }))
  setkey(sequenceDT, ensembl_transcript_id)
  
  
  
  # get all sequences for splice site mutated transcripts (exons + introns)
  print("(3/4) Downloading intron sequences for transcripts with splice site mutations")
  spliceTscrIDs <<- vcf.DT[grep("splice", EFFECT), unique(FEATUREID.noversion)]
  if(length(spliceTscrIDs) > 0){
    intronexonSeqs <<- data.table(getSequence(id=spliceTscrIDs, type="ensembl_transcript_id", seqType="transcript_exon_intron", mart = db))
    setkey(intronexonSeqs, ensembl_transcript_id)
  }
  
  # get UTRs 
  print("(4/4) Downloading UTR sequences")
  UTR3seq <<- rbindlist(lapply(1:length(tscrIDs.chunks), function(x){
    utr3.tmp <- getSequence(id=tscrIDs.chunks[[x]], type="ensembl_transcript_id", seqType="3utr", mart = db)
    utr3.tmp
  }))
  
  UTR5seq <<- rbindlist(lapply(1:length(tscrIDs.chunks), function(x){
    utr5.tmp <- getSequence(id=tscrIDs.chunks[[x]], type="ensembl_transcript_id", seqType="5utr", mart = db)
    utr5.tmp
  }))
  
  UTRseqs <<- data.table(merge(UTR3seq, UTR5seq, all=T))
  setkey(UTRseqs, ensembl_transcript_id)

  # Init files
  
  vcf.filename <- gsub("(.*)\\..*","\\1", basename(vcf.file))
  outfile.name <- paste0("results_out/",vcf.filename,"_results.tempfasta")
  errorfile.name <- paste0("results_out/",vcf.filename,"_error.txt")
  
  print("Processing in parallel...")
  
  catchreturns <- foreach(rowN = 1:nrow(vcf.DT), .packages = c('data.table'), .combine = rbind ) %dopar% {
    
    if(debugMode){
      print(rowN)
    }

    # Init variables ----
    mutDT <- vcf.DT[rowN,]
    header.NOTEtag <- NA
    
    tscrID <- mutDT[, FEATUREID]
    tscrID <- gsub("\\..*","", tscrID)
    
    # error + extension flag for all runs
    notCorrect <- FALSE
    isExtended <- FALSE
    
    # Init variables
    mutInfo <- mutDT[, HGVS_C]
    # regex to get elements from mutInfo
    regex <- ".\\.(\\d*_\\d*|\\d*)([a-z]*)([A-Z]*)"
    
    # special case splice region variant (is fine for del and sub, because the REF is ignored BUT wrong for ins, because ref can be intronic)
    isSpliceRegion <- grepl("splice_region", mutDT[, EFFECT])
    isInDel <- grepl("ins", mutInfo) & grepl("del", mutInfo)
    isInser <- grepl("ins", mutInfo)
    isDelet <- grepl("del", mutInfo)
    isSubst <- grepl(">", mutInfo)
    isDup <- grepl("dup",  mutInfo)
    isSpliceAD <- grepl("acceptor|donor", mutDT$EFFECT)
    isSpliceA <- grepl("acceptor", mutDT$EFFECT)
    isSpliceD <- grepl("donor", mutDT$EFFECT)

    # change regex for substitution (>)
    if(isSubst){
      regex <- ".\\.(\\d*_\\d*|\\d*)([A-Z]*)>([A-Z]*)"
      mutSubRef <- gsub(regex,"\\2",mutInfo)
    }

    # catch cryptic variant calling outputs
    # examples:
    # "c.-3_77delACCATGGAGGCGGCGGTCGCTGCTCCGCGTCCCCGGCTGCTCCTCCTCGTGCTGGCGGCGGCGGCGGCGGCGGCGGCGGC"
    # "n.370_*1delAAA"
    if (grepl(".\\.(\\+|\\-)[0-9]+|\\*", mutInfo)){
      notCorrect <- TRUE
      if(debugMode){warn(paste0("Invalid HGVS_C: ",mutInfo))}
      return(NA)
    }
    
    # again check for correct position string (for the downstream isSpliceAD & x functions)
    if(isSpliceAD){
      cdsPos <- as.numeric(gsub("..(\\d*).*","\\1",mutInfo))
      if (is.na(cdsPos)){
        notCorrect <- TRUE
        if(debugMode){warn("Invalid HGVS_C")}
        # what the hell are cases like this? c.-22-1G>A
        return(NA)
      } 
    }
    
    # SPLICE AD FOR INDEL (insertion with deletion at the same time)
    if(isSpliceAD & isInDel){

      caseRegex <- "(.\\.)(.*)_(.*)(del.*)"
      caseVal <- gsub(caseRegex,"\\3",mutInfo)
      
      # case1 = c.240-3_240-1delCAGinsTTT
      # case2 = c.738-1_738delGGinsAA
      # case3 = c.7434_7434+2delGGTinsAAA
      
      # case1
      if(grepl("-", caseVal)){
        calcRegex <- "(.\\.)(\\d*)(-|\\+)(\\d*)_(\\d*)(-|\\+)(\\d)(.*)"
        startString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- as.numeric(gsub(calcRegex,"\\2",mutInfo))
        startCalc <-  as.numeric(gsub(calcRegex,"\\4",mutInfo))
        endRaw <- as.numeric(gsub(calcRegex,"\\5",mutInfo))
        endCalc <- as.numeric(gsub(calcRegex,"\\7",mutInfo))
        rest <- gsub(calcRegex,"\\8",mutInfo)
        
        startNew <- startRaw-startCalc
        endNew <- endRaw-endCalc
        
      } else if (grepl("\\+", caseVal)){
        #case3
        calcRegex <- "(.\\.)(\\d*)_(\\d*)(-|\\+)(\\d*)(.*)"
        startString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- as.numeric(gsub(calcRegex,"\\2",mutInfo))
        startCalc <- 0
        endCalc <-  as.numeric(gsub(calcRegex,"\\5",mutInfo))
        endRaw <- as.numeric(gsub(calcRegex,"\\3",mutInfo))
        rest <- gsub(calcRegex,"\\6",mutInfo)
        
        startNew <- startRaw+startCalc
        endNew <- endRaw+endCalc
        
      } else {
        # case2
        calcRegex <- "(.\\.)(\\d*)(-|\\+)(\\d*)_(\\d*)(.*)"
        startString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- as.numeric(gsub(calcRegex,"\\2",mutInfo))
        startCalc <-  as.numeric(gsub(calcRegex,"\\4",mutInfo))
        endRaw <- as.numeric(gsub(calcRegex,"\\5",mutInfo))
        endCalc <- 0
        rest <- gsub(calcRegex,"\\6",mutInfo)
        
        startNew <- startRaw-startCalc
        endNew <- endRaw-endCalc
      }

      mutInfo <- paste0(startString,startNew,"_", endNew,rest)
      cdsPos <- startNew
      # regex is reset in isIns 
    }
    
    # splice AD for substitution
    if(isSpliceAD & isSubst){

      calcRegex <- "(.\\.)(\\d*)((-|\\+)\\d*)(.*)" # sub
      startString <- gsub(calcRegex,"\\1",mutInfo)
      startRaw <- as.numeric(gsub(calcRegex,"\\2",mutInfo))
      startCalc <-  gsub(calcRegex,"\\3",mutInfo)
      rest <- gsub(calcRegex,"\\5",mutInfo)
      startNew <- startRaw+eval(parse(text=startCalc)) # calculate new start
      
      # rewrite mutInfo
      mutInfo <- paste0(startString,startNew,rest)
      
      # trset variables
      regex <- ".\\.(\\d*_\\d*|\\d*)([A-Z]*)>([A-Z]*)"
      
      mutSubRef <- gsub(regex,"\\2",mutInfo)
      cdsPos <- startNew
      
      # example c.1-2A>C (ENST00000295702.8, splice acceptor mutation, but actually 5'UTR)
      if (cdsPos < 0){
        notCorrect <- TRUE
        if(debugMode){warn(paste0("Invalid HGVS_C: ",mutInfo))}
        return(NA)
      }
    }
    
    # splice AD for insertion
    if(isSpliceAD & isInser & !isInDel){
      regex <- ".\\.(.*)_.*(ins)(.*)"
      # we do not need to change cdsPos, because the IncludeMutRetainedIntron() function
      # greps the closes exon or the hit one, so both coordinates are be fine
    
      # but we need a new mutStart and mutEnd for some cases 
      # we will rewrite the mutInfo 
      calcRegex <- "(.\\.)(.*)_(.*)(ins.*)"
      
      firstString <- gsub(calcRegex,"\\1",mutInfo)
      startRaw <- gsub(calcRegex,"\\2",mutInfo)
      endRaw <- gsub(calcRegex,"\\3",mutInfo)
      rest <- gsub(calcRegex,"\\4",mutInfo)
      
      startNew <- eval(parse(text=startRaw)) # calculate from string
      endNew <- eval(parse(text=endRaw)) # calculate from string

      # can there be cases without _ ?
      if(!grepl("_", mutInfo)){
        stop("debug this A4")
      }
      
      mutInfo <- paste0(firstString,startNew,"_",endNew,rest)
    }
    
    # splice AD for deletion
    if(isSpliceAD & isDelet & !isInDel){
      # we do not need to change cdsPos, because the IncludeMutRetainedIntron() function
      # greps the closes exon or the hit one, so both coordinates are be fine
      # but we need a new mutStart and mutEnd (this regex is different than the ins above)
      # we will rewrite the mutInfo 

      # two cases
      #n.425delA # case 1
      #c.3+5454_3+5455delTT # case 2
      
      # case 1
      if(grepl("_", mutInfo)){
        regex <- ".\\.(.*_.*)(del)(.*)"
        calcRegex <- "(.\\.)(.*)_(.*)(del.*)"

        firstString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- gsub(calcRegex,"\\2",mutInfo)
        endRaw <- gsub(calcRegex,"\\3",mutInfo)
        rest <- gsub(calcRegex,"\\4",mutInfo)
        
        startNew <- eval(parse(text=startRaw)) # calculate from string
        endNew <- eval(parse(text=endRaw)) # calculate from string
        
        mutInfo <- paste0(firstString,startNew,"_",endNew,rest)
      } else {
        # case 2
        regex <- ".\\.(.*)(del)(.*)"
        calcRegex <- "(.\\.)(.*)(del.*)"
        
        firstString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- gsub(calcRegex,"\\2",mutInfo)
        rest <- gsub(calcRegex,"\\3",mutInfo)
        
        startNew <- eval(parse(text=startRaw)) # calculate from string
        mutInfo <- paste0(firstString,startNew,"_",rest)
      }
    }

    # dup is the almost the same as ins but the insertions is done on/after mutEnd
    if(isSpliceAD & isDup){

      # but we need a new mutStart and mutEnd (this regex is different than the ins above) and we will rewrite the mutInfo 
      
      # there are 2 cases: with and without _
      # mutInfo <- c.938-4_938-3dupTT # case1
      # mutInfo <- c.1741-3dupT # case2
      
      # case 1
      if(grepl("_", mutInfo)){
        regex <- ".\\.(.*)_.*(dup)(.*)"
        calcRegex <- "(.\\.)(.*)_(.*)(dup.*)"
        
        firstString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- gsub(calcRegex,"\\2",mutInfo)
        endRaw <- gsub(calcRegex,"\\3",mutInfo)
        rest <- gsub(calcRegex,"\\4",mutInfo)
        
        startNew <- eval(parse(text=startRaw)) # calculate from string
        endNew <- eval(parse(text=endRaw)) # calculate from string
        
        mutInfo <- paste0(firstString,startNew,"_",endNew,rest)
      } else {
        # case 2
        regex <- ".\\.(.*)(dup)(.*)"
        calcRegex <- "(.\\.)(.*)(dup.*)"
        
        firstString <- gsub(calcRegex,"\\1",mutInfo)
        startRaw <- gsub(calcRegex,"\\2",mutInfo)
        rest <- gsub(calcRegex,"\\3",mutInfo)

        startNew <- eval(parse(text=startRaw)) # calculate from string
        mutInfo <- paste0(firstString,startNew,"_",rest)
      }
    }
    
    #special case c.1450_1454delGATGCinsTTGT
    if(isInDel){
      #save the ins information for later
      insInfo <- gsub(".*(ins[A-Z]*).*","\\1",mutInfo)
      
      # set only the del and change after del was performed
      mutInfo <- gsub("ins[A-Z]*","",mutInfo)
      
      #reset regex (due to splice acceptor donor changes with indel as mutation)
      regex <- ".\\.(\\d*_\\d*|\\d*)([a-z]*)([A-Z]*)"
    }
    
    # another special case: insertion on first bp of an exon but REF is on intron so it is +1 (but the mutation is NOT splice)
    # normally c.2361_2362insAGAT (without +- on mutEnd)
    # new: c.11160_11160+1insAGAGCAT
    # so is ins, not splice and has +or- 
    
    # additionally, the check for correctness will fail since the REF is intronic 
    # so we will skip it and add a note to the header
    specialInsCase.intronRef <- FALSE
    
    if(isInser && !isSpliceAD &&  grepl("\\+|-", mutInfo)){
      calcRegex <- "(.\\.)(.*)_(.*)(ins.*)"
      
      firstString <- gsub(calcRegex,"\\1",mutInfo)
      startRaw <- gsub(calcRegex,"\\2",mutInfo)
      endRaw <- gsub(calcRegex,"\\3",mutInfo)
      rest <- gsub(calcRegex,"\\4",mutInfo)
      
      startNew <- eval(parse(text=startRaw)) # calculate from string
      endNew <- eval(parse(text=endRaw)) # calculate from string
      
      mutInfo <- paste0(firstString,startNew,"_",endNew,rest)
    
      specialInsCase.intronRef <- TRUE
    }
    
    #temp debug
    if(isInser && !grepl("_", mutInfo)){
      stop("debugD1")
    }
    
    tempSplit <- gsub(regex,"\\1",mutInfo) # get the mutation window 
    tempSplit <- as.numeric(strsplit(tempSplit, split="_")[[1]]) # we only need start
    # for cases without _, it will simply take the start
    
    # temp debug
    if(any(is.na(tempSplit))){
      stop("debug this A3")
    }
    
    mutStart.firstSave <- tempSplit[1] # we need this for debug prints
    mutStart <- tempSplit[1]
    
    mutEnd <- tempSplit[length(tempSplit)] 
    #can only be of length 1 or 2, but 2 can be missing so this gives then 1 for both
    # this will be equal to mutStart for ins
    
    mutType <- gsub(regex,"\\2",mutInfo) # get the mutation type (del, ins, dup)
    mutString <- gsub(regex,"\\3",mutInfo) # get mutation sequence (to delete or to insert)
    
    # HGVS_P given for coding, but not for non-coding
    if (mutDT[, HGVS_P] != "NA"){
      mutAApos <- as.numeric(mutDT[, gsub(".\\.[a-zA-Z]*(\\d*).*","\\1",HGVS_P)]) # get mutation sequence (to delete or to insert)
    } else {
      mutAApos <- ceiling(mutStart/3) # bugfix: round is now ceiling
    }
    
    # set mutType correctly (change from mutSubRef (1 base) to sub, its more tidy)
    if(isSubst){
      mutType <- "sub"
    }
    
    # Get sequence ----
    
    seq <- as.character(sequenceDT[ensembl_transcript_id == tscrID, 1])
    
    isIntronInclusion <- FALSE
    # include intron for splice donor acceptor
    if(isSpliceAD){
      intronFun.out <- IncludeMutRetainedIntron(mutDT,seqType.input)
      
      isIntronInclusion <- TRUE
      
      if(is.null(intronFun.out[[1]])){
        if(debugMode){warn(paste0("Version problem (exon boundaries changed), ignoring: ", tscrID))}
        return(NA)
        }
      
      seq <- intronFun.out[[1]]
      intronLength <- intronFun.out[[2]]
      included.intron <- intronFun.out[[3]]
      
      # include intron to relative start end position in sequence
      # only for acceptor? (because there mutpos is at the END of intron)
      # donor should be fine
      
      if(isSpliceA && isSpliceD && isSubst){
        # there are some cryptic cases of donor + acceptor mutation and we do not know if we include intron or not
        # this solves most of the cases, the others will be discarded in mutation check
        # all cases are printed to warning_cases.txt
      } else if(isSpliceA){
        # normal case
        mutStart <- mutStart+intronLength
        mutEnd <- mutEnd+intronLength
      }
      
    }
    
    seq.old <- seq # resave original sequences after intron inclusion

    # Check if transcript found ----
    #if(isEmpty(seqL[,1])){
    if(length(seq)==0){
      if(debugMode){warn(paste0("Transcript not found in Ensembl! Is the transcript ID input correct for: ", organism), "?")}
      return(NA)
    }
    
    # Check for CDS 5' and/or 3' incomplete ----
    noStart <- subseq(seq.old,1,3) != "ATG" #first 3 bases not ATG
    
    # check for stop codon at end
    tripplets <- CutTripplets(seq,1) # cut sequence into tripplets
    lastCodon <- tripplets[length(tripplets)] # get last codon
    noStop <- lastCodon != "TAA" & lastCodon != "TAG" & lastCodon != "TGA" # check if stop codon is in frame
    
    # Only check last without frame if donor acceptor case
    if(isSpliceAD){
      stopC <- subseq(seq.old,nchar(seq.old)-2,nchar(seq.old))
      noStop <- stopC != "TAA" & stopC != "TAG" & stopC != "TGA"
    }
    
    if(seqType.input == "coding" & (noStart | noStop ) ){
      
      if(debugMode){warn(paste0("Incomplete CDS in biomart annotation (no start/stop), ignoring: ", tscrID))}
      
      incompleteCDS.filename <- paste0("results_out/",vcf.filename,"_incompleteCDS.txt")
      write.table(mutDT,file=incompleteCDS.filename,append=TRUE, col.names = F, row.names = F, quote = F)
      return(NA)
    }
    
    # Check mutation sequence correctness ----
    # First check if the mutString is actually present at that position (transcription version check)
    if (mutType == "del"){
      if(substr(seq.old, mutStart, mutEnd) != mutString){
        
        if(debugMode){warn(paste0("Correct position for deletion NOT found in given sequence, ignoring: ", tscrID))}
        notCorrect <- TRUE
      }
    }
    
    if (mutType == "ins" & !isSpliceRegion){
      
      # in this case we need the strand information to adjust REF
      
      strand <- exonDT[ensembl_transcript_id == tscrID, unique(strand)]

      REF.given <- mutDT[, REF] # for +strand 
      REF.seq <- substr(seq.old, mutStart, mutStart) # for +strand
      
      if(strand == "-1"){
        # set other REF
        REF.given <- chartr("ATGC","TACG",REF.given)
        REF.seq <- substr(seq.old, mutEnd, mutEnd) # last position as ref because -strand flip
      }
      
      if (REF.given != REF.seq && !specialInsCase.intronRef){
        if(debugMode){warn(paste0("Correct position for insertion NOT found in given sequence, ignoring: ", tscrID))}
        notCorrect <- TRUE 
      }
    }
    
    if (mutType == "dup"){
      # sequence at insertion site should be given
      if(substr(seq.old, mutStart, mutStart+nchar(mutString)-1) != mutString){
        
        if(debugMode){warn(paste0("Correct position for deletion NOT found in given sequence, ignoring: ", tscrID))}

        notCorrect <- TRUE 
      }
    }
    
    if (mutType == "sub"){
      if(substr(seq.old, mutStart, mutEnd) != mutSubRef){
        
        if(debugMode){warn(paste0("Correct position for substitution NOT found in given sequence, ignoring: ", tscrID))}

        notCorrect <- TRUE 
      }
    }
    
    if(notCorrect){
      write.table(mutDT,file=errorfile.name,append=TRUE, col.names = F, row.names = F, quote = F)
      return(NA)
    }
    
    # Perform mutation calculation ------
    # INSERTION
    # insert mutString exactly onto mutStart, so keep the base at mutStart position in left block and +1 in the right block
    
    # DELETION
    if (mutType == "del"){
      # to delete simply shrink the sequence around the position
      seq <- paste0(substr(seq.old, 1, mutStart-1), substr(seq.old, mutEnd+1, nchar(seq.old)))
    }
    
    # Catch special case of insertion after deletion annotated in 1 mutation
    if(isInDel){
      # insInfo is the new mutInfo
      mutType <- "ins"
      mutString <- gsub("ins([A-Z]*)","\\1",insInfo) # get mutation sequence (to delete or to insert)
      mutStart <- mutStart-1
      if(debugMode){warn("Special case of ins+del! Please double check the mutation result")}
    }
    
    #INSERTION
    if (mutType == "ins"){
      # use seq and not seq.old in this one, because of del+ins case
      seq <- paste0(substr(seq, 1, mutStart), mutString, substr(seq, mutStart+1, nchar(seq)))
      
      # mutString is NOT the position of the mutation, but of the reference 1 before it!
      mutStart <- mutStart+1
    }
    
    # DUPLICATION
    # dup is the almost the same as ins but the insertions is done on/after mutEnd
    if (mutType == "dup"){
      seq <- paste0(substr(seq.old, 1, mutEnd), mutString, substr(seq.old, mutEnd+1, nchar(seq.old)))
    }
    
    # SUBSTITUTION
    if (mutType == "sub"){
      # to substitutie simply bind the mutString into the position
      substr(seq, mutStart, mutEnd) <- mutString
    }
    
    
    # Find ORFs ------
    
    # Check if mutation in coding removed start and/or stop codons, then read from/into UTR (extending seq by UTR sequence)
    noStart <- subseq(seq,1,3) != "ATG" #first 3 bases not ATG
    
    #tripplets <- CutTripplets(seq,1) # cut sequence into tripplets
    
    lastCodon <- tripplets[length(tripplets)] # get last codon
    noStop <- lastCodon != "TAA" & lastCodon != "TAG" & lastCodon != "TGA" # check if stop codon is in frame
    
    # Only check last without frame if donor acceptor case
    if(isSpliceAD){
      stopC <- subseq(seq.old,nchar(seq.old)-2,nchar(seq.old))
      noStop <- stopC != "TAA" & stopC != "TAG" & stopC != "TGA"
    }

    seq.saveLen <- nchar(seq) # save length of mutated sequence to compare against ORF extension
    if(seqType.input == "coding" & (noStart & noStop)){
      # problematic case: frameshift due to to deletion in start codon AND no stop codon (unlikely but possible?)
      # add 5' and 3' UTR sequence to seq and try to find ORFs
      
      UTR3seq <- as.character(UTRseqs[ensembl_transcript_id == tscrID, `3utr`])
      UTR5seq <- as.character(UTRseqs[ensembl_transcript_id == tscrID, `5utr`])

      seq <- paste0(UTR5seq,seq,UTR3seq)
      
      isExtended <- TRUE
      print("start and stop codons lost! trying to extend ORF into UTRs")
    } else if(seqType.input == "coding" & noStart){
      # add 5' UTR sequence to seq (and later try to find ORFs)
      
      UTR5seq <- as.character(UTRseqs[ensembl_transcript_id == tscrID, `5utr`])
      seq <- paste0(UTR5seq,seq)
      
      isExtended <- TRUE
      if(debugMode){print("start codon lost! trying to extend ORF into 5'UTR")}
    } else if(seqType.input == "coding" & noStop){
      # add 3' UTR sequence to seq (and later try to find ORFs)
      
      UTR3seq <- as.character(UTRseqs[ensembl_transcript_id == tscrID, `3utr`])
      seq <- paste0(seq,UTR3seq)
      
      isExtended <- TRUE
      if(debugMode){print("stop codon lost! trying to extend ORF into 3'UTR")}
    }
    
    # check if UTRs were not found (we do this after the processing, since those cases are very rare and we do not want to lookup for all)
    if(exists("UTR3seq")){
      if(length(UTR3seq) == 0){
        notCorrect <- TRUE
        if(debugMode){warn("3UTR for transcript not found in biomart")}
        rm(UTR3seq) # just in case I forget and it stays in global envir
        return(NA)
      }
    }
    if(exists("UTR5seq")){
      if(length(UTR5seq) == 0){
        notCorrect <- TRUE
        if(debugMode){warn("3UTR for transcript not found in biomart")}
        rm(UTR5seq) # just in case I forget and it stays in global envir
        return(NA)
      }
    }
    
    # get sequence of ORFs and translate
    postprocessed.out <- Postprocess.Mutation(seq, seqType.input,isExtended,mutStart,isSpliceAD)

    AA.list <- postprocessed.out[[2]]
    DNAseq.list <- postprocessed.out[[1]]
    ORF.DT <- DefineORFs(seq,mutStart,seqType.input,isExtended,isSpliceAD)

    # Remove the stop codon * (if not empty)
    if(!is.null(AA.list)){
      AA.list <- subseq(AA.list,1,nchar(AA.list)-1)
    }
    
    # print if seq was extende into UTR
    if(!is.null(AA.list) & isExtended){
      extend.length <- nchar(AA.list)*3
      if(debugMode){print(paste0("ORF extended ",abs(extend.length-seq.saveLen)," bp into UTR"))}
    }
   
    
    # Large loop (get peptide, check synonymous, print out) ----- 
    
    # First, check if any ORF was found
    if(is.null(ORF.DT)){ # if no ORF found (e.g. start/stop loss) then print into ORF_lost file
      
      transcriptLost.filename <- paste0("results_out/",vcf.filename,"_ORF_lost.txt")
      
      ORFlost <- paste(mutDT[,CHROM], mutDT[,POS], mutDT[,REF], mutDT[,ALT], mutDT[,GENE],
                          mutDT[,EFFECT],mutDT[,IMPACT], mutDT[,FEATUREID], mutDT[,HGVS_C], mutDT[,HGVS_P],mutDT[,CDS_POS],mutDT[,AA_POS],
                          seqType.input, vcf.filename,sep="|")
      ORFlost <- paste0(">",ORFlost)
      ORFlost <- paste0(ORFlost,"|ORFlost") # ORflost is ORF column, there is no flag/note column
      
      write(ORFlost,file=transcriptLost.filename,append=TRUE)
      
      if(debugMode){warn(paste0("No ORF found, skipping: ", tscrID))}
      return(NA)
    }
    
    # For non-codings (type=cdna) there can be multiple ORFs, which have to be debug printed individually (for synonymous check)
    # for coding, there will be only 1 element in DNAseq.list and we skip it
    nORFS <- length(DNAseq.list)
    
    # save the original variables for the ORF loop
    xsave.seq <- copy(seq)
    xsave.seq.old <- copy(seq.old)
    xsave.mutStart <- copy(mutStart)
    xsave.mutEnd <- copy(mutEnd)
    xsave.mutAApos <- copy(mutAApos)
    xsave.mutStart.firstSave <- mutStart.firstSave
    xsave.mutDT <- copy(mutDT)

    # for every ORF  
    for(ORFid in 1:nORFS){
      # for coding, a splice acceptor (and probably donor) mutation will lead to intron inclusion
      # in that intron there can be a premature stop codon, so the ORF will be shorter and NOT INCLUDE the mutation site
      # which breaks all outprint functions, so we use an additional outprint function
      if(isSpliceAD){
        if(mutStart > ORF.DT[ORFid, endPos]){
          # premature stop! we can not print debug, because mutation is not in ORF 
          
          seq <- DNAseq.list[ORFid] # get the ORF DNA sequence

          # get AA sequence
          AAseq.new <- AA.list[ORFid] #translate(DNAString(x=seq), if.fuzzy.codon = "solve")
          
          # get gene and transcript biotype
          Gbiotype <- exonDT[ensembl_transcript_id == tscrID, unique(gene_biotype)]
          Tbiotype <- exonDT[ensembl_transcript_id == tscrID, unique(transcript_biotype)]
          
          # set header
          headerLine <- paste(paste(mutDT[, .(CHROM, POS, REF, ALT, GENE, EFFECT, IMPACT, FEATUREID, HGVS_C, HGVS_P,
                                              CDS_POS, AA_POS, G5, AC, AN, AF,
                                              TumorVF, TumorAD, TumorRD, NormalAD, NormalRD, patientID, SOURCE, RESCUED)], collapse="|"), 
                              seqType.input, Gbiotype, Tbiotype, versionControlTag, sep="|")
          headerLine <- paste0(">",headerLine)
          
          # set ORF column
          ORF.window <- paste0("[",ORF.DT[ORFid,paste0(startPos,"-",endPos)], "]")
          header.ORFtag <- paste0("ORF",ORFid,ORF.window)

          # set the note column (was NA or rescued)
          # we do not want NA before the following tags (e.g. prematureStop), but we want "rescued-premature stop"
          header.NOTEtag <- c(header.NOTEtag, "prematureIntronStop")
          header.NOTEtag <- header.NOTEtag[!is.na(header.NOTEtag)] # we remove the NA here
          header.NOTEtag <- paste0(header.NOTEtag, collapse = ",")
          
          # combine final headerLine
          headerLine <- paste(headerLine, header.ORFtag, header.NOTEtag, sep="|")
          outputLine <- paste0("&START&",headerLine, "&NEWLINE&",as.character(AAseq.new),"&END&")
          
          # write to file
          conn <- file( paste0(outfile.name ,"_pid", Sys.getpid()) , open = "a" )
          write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
          close( conn )

          # write to report file
          
          fullreport.filename <- paste0("results_out/",vcf.filename,"_fullReport.txt")
          sometext <- "Complex case of premature stop codon due to intron inclusion in splice site mutation (full sequences in _report_PrematureIntronStop.txt"
          outputLine <- paste0("&START&",headerLine, "&NEWLINE&",sometext,"&END&", "\n")
          
          conn <- file( paste0(fullreport.filename ,"_pid", Sys.getpid()) , open = "a" )
          write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
          close( conn )
          
          # also write DNA and AA into a new file
          PrematureIntronStop.filename <- paste0("results_out/",vcf.filename,"_report_PrematureIntronStop.txt")
          
          outputLine <- paste(headerLine,
                              paste0("<Included Intron: ",included.intron, " (length=",intronLength,")"),
                              paste0("<Original MutPos: ",mutDT[,HGVS_C], "; New MutPos (Original+Intron): ", mutStart),
                              "<original DNA sequence (full length transcript)",
                              as.character(seq.old),
                              "<mutated DNA sequence (short ORF, until premature stop)",
                              as.character(seq),
                              "<mutated peptide AA sequence (short ORF, until premature stop)",
                              as.character(AAseq.new),
                              sep="&NEWLINE&")
          outputLine <- paste0("&START&",outputLine,"&END&","\n")
          
          conn <- file( paste0(PrematureIntronStop.filename ,"_pid", Sys.getpid()) , open = "a" )
          write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
          close( conn )

          if(debugMode){print(paste0("Premature STOP after intron inclusion for splice site mutation in: ", tscrID, " (results in file without print)"))}
          # after we printed it we "skip"
          # we dont need a synonymous check
          return(NA)
        }
      }

      # for every non-codings (cdna): set the data
      # only for non splice acceptor/donor mutations
      if(seqType.input == "cdna"){
        seq <- DNAseq.list[ORFid] # get the ORF DNA sequence
        
        mutStart <- mutStart - ORF.DT[ORFid, startPos]+1
        mutEnd <- mutEnd - ORF.DT[ORFid, startPos]+1
        mutAApos <- ceiling(mutStart/3)
        
        # this is a dirty indel catch downstream, we bypass this for now
        mutStart.firstSave <- mutStart
        
        # problem large mutations: n.690_691insAAAAAATAAATAAATAATAAATCCCGTAAT
        if(nchar(seq.old) < ORF.DT[ORFid, endPos]){
          seqdif <- ORF.DT[ORFid, endPos]-nchar(seq.old) # ORF is x larger than seq.old
          seq.old <- subseq(seq.old, ORF.DT[ORFid, startPos], ORF.DT[ORFid, endPos]-seqdif)
        } else {
          # classic get ORF window old seq
          seq.old <- subseq(seq.old, ORF.DT[ORFid, startPos], ORF.DT[ORFid, endPos])
        }
      }
      
      
      # Full debug outprints ----
      if (TRUE == TRUE){
        # Get original sequence with positions
        x <- 10
        
        # if mutation insertion is very large, increase x
        if (nchar(mutString) > 10){
          x <- nchar(mutString)
        }
        
        # Catch AA at end so that x would be out of range
        len <- nchar(seq.old)
        if (mutStart+x > len){
          x <- len-mutStart
        } else if (mutStart-x < 0){
          x <- mutStart-1
        }
        
        #reset mutStart due to indel case (what?)
        if(isInDel){
          mutStart <- mutStart.firstSave
        }
        
        # Get original sequence
        substr.old <- substr(seq.old,mutStart-x,mutStart+x)
        substr.old <- data.table(do.call(rbind, strsplit(substr.old, "")))
        names(substr.old) <- as.character((mutStart-x):(mutStart+x))[1:ncol(substr.old)] # catch shorter then 21 because of start/end
        # Get mutated sequence with positions
        substr <- substr(seq,mutStart-x,mutStart+x)
        substr <- data.table(do.call(rbind, strsplit(substr, "")))
        names(substr) <- as.character((mutStart-x):(mutStart+x))[1:ncol(substr)]
        
        # Combine and show mutation window with *
        seqCompare <- rbind(substr.old,substr, fill=T)
        
        # on very long mutations, we need a manual curation
        # e.g. n.28_79delTCGCTGTCGCCGCCGCCGCCGCCCGCCGCCGCCGCCGCCGCCGCCGCCGCTG
        
        # get the number of the last position in list
        tempLen <- length(seqCompare)
        windowPrintEnd <- as.numeric(names(seqCompare[, ..tempLen]))
        
        # set correct mut end for indels
        if(isInser){
          mutEnd <- mutStart + nchar(mutString)-1
        }
        
        # compare if the mutation would be out of bound
        mutationPrintEnd <- mutStart+nchar(mutString)
        
        if(windowPrintEnd < mutationPrintEnd){
          #newEnd <- length(seqCompare)
          setnames(seqCompare, paste0(mutStart:windowPrintEnd), paste0(mutStart:windowPrintEnd, "*"))
        } else {
          setnames(seqCompare, paste0(mutStart:mutEnd), paste0(mutStart:mutEnd, "*"))
        }
        
        # get original AA sequence
        seqString.old <- DNAString(seq.old)
        options(warn=-1)
        AAseq.old <- translate(seqString.old, if.fuzzy.codon = "solve")
        
        # get new AA sequence
        seqString <- DNAString(x=seq)
        AAseq.new <- translate(seqString, if.fuzzy.codon = "solve")
        options(warn=0)
        
        x <- 10
        # for mutation in codon 1, change AA to 1
        if(mutStart <= 3){
          mutAApos <- 1
        } else if(mutStart >= nchar(seq)-3){
          x <- 0
        }
        # Catch AA at end so that x would be out of range
        len <- min(length(AAseq.new),length(AAseq.old))
        
        # there are some cases of splice sites, where len < mutAApos
        if (mutAApos > len){
          mutAApos <- len
        }
        
        if (mutAApos+x >= len){
          #x <- len-mutAApos
          x <- 0
        } else if (mutAApos-x <= 0){
          x <- 0
        }
        
        # show AA sequence
        substrAA.old <- as.character(substr(AAseq.old, mutAApos-x, mutAApos+x))
        substrAA.old <- data.table(do.call(rbind, strsplit(substrAA.old, "")))
        names(substrAA.old) <- as.character((mutAApos-x):(mutAApos+x))
        setnames(substrAA.old, as.character(mutAApos), paste0(mutAApos,"*"))
        
        substrAA.new <- as.character(substr(AAseq.new, mutAApos-x, mutAApos+x))
        substrAA.new <- data.table(do.call(rbind, strsplit(substrAA.new, "")))
        names(substrAA.new) <- as.character((mutAApos-x):(mutAApos+x))
        setnames(substrAA.new, as.character(mutAApos), paste0(mutAApos,"*"))
        
        AACompare <- rbind(substrAA.old,substrAA.new, fill=T)
        
        # now update the missing information for noncodings (cdna)
        if(seqType.input == "cdna"){
          oldA <- AACompare[1, grep("\\*", colnames(AACompare)), with=F]
          newA <- AACompare[2, grep("\\*", colnames(AACompare)), with=F]
          mutDT[,HGVS_P := paste0(oldA,">",newA)]
          mutDT[,CDS_POS := mutStart]
          mutDT[,AA_POS := mutAApos]
          mutDT[,CDS_LEN := nchar(seq)]
          mutDT[,AA_LEN := nchar(AAseq.new)] #or nchar(seq)/3
        }
        # also update coding is the information is missing
        mutDT[is.na(CDS_POS),CDS_POS := mutStart]
        mutDT[is.na(AA_POS),AA_POS := mutAApos]
        mutDT[is.na(CDS_LEN),CDS_LEN := nchar(seq)]
        mutDT[is.na(AA_LEN),AA_LEN := nchar(AAseq.new)]
        
        # get gene and transcript biotype
        Gbiotype <- exonDT[ensembl_transcript_id == tscrID, unique(gene_biotype)]
        Tbiotype <- exonDT[ensembl_transcript_id == tscrID, unique(transcript_biotype)]
        
        # set header
        headerLine <- paste(paste(mutDT[, .(CHROM, POS, REF, ALT, GENE, EFFECT, IMPACT, FEATUREID, HGVS_C, HGVS_P,
                                            CDS_POS, AA_POS, G5, AC, AN, AF, 
                                            TumorVF, TumorAD, TumorRD, NormalAD, NormalRD, patientID, SOURCE, RESCUED)], collapse="|"), 
                            seqType.input, Gbiotype, Tbiotype, versionControlTag, sep="|")
        headerLine <- paste0(">",headerLine)
        
        # set ORF column
        ORF.window <- paste0("[",ORF.DT[ORFid,paste0(startPos,"-",endPos)], "]")
        header.ORFtag <- paste0("ORF",ORFid,ORF.window)

        # add additional note if intron was included (these cases are very inlikely, since most of the times there is a stop in the included intron)
        if(isIntronInclusion){
          addHeader1 <- paste0("IntronInclusion[intron:",included.intron,"]")
          header.NOTEtag <- c(header.NOTEtag, addHeader1)
        }
        
        # add additional note if special case ins (REF on intron because last position of exon)
        # this cases are not checked for correctness
        if(specialInsCase.intronRef){
          header.NOTEtag <- c(header.NOTEtag, "InsertionRefInIntronNotValidated")
        }

        # remove NA if there are additional tags
        if(length(header.NOTEtag[!is.na(header.NOTEtag)]) > 0){
          header.NOTEtag <- header.NOTEtag[!is.na(header.NOTEtag)] # we remove the NA here
        }
        # combine all NOTE tags
        header.NOTEtag <- paste0(header.NOTEtag, collapse = ",")
        
        # combine final headerLine
        headerLine <- paste(headerLine, header.ORFtag, header.NOTEtag, sep="|")
        
        # Print to console
        if (debugMode == TRUE){
          mutDT.withORFid <- cbind(mutDT, ORF=paste0("ORF",ORFid))
          print(mutDT.withORFid)
          print(seqCompare)
          print(AACompare)
        }
        
        if (debugPrintMode == TRUE){
          
          fullreport.filename <- paste0("results_out/",vcf.filename,"_fullReport.txt")
          
          outputLine <- paste(headerLine,
                              paste(colnames(seqCompare), collapse = "\t"),
                              paste(seqCompare[1,], collapse = "\t"),
                              paste(seqCompare[2,], collapse = "\t"),
                              paste(colnames(AACompare), collapse = "\t"),
                              paste(AACompare[1,], collapse = "\t"),
                              paste(AACompare[2,], collapse = "\t"),
                              sep="&NEWLINE&")
          outputLine <- paste0("&START&",outputLine,"&END&","\n")
          
          
          conn <- file( paste0(fullreport.filename ,"_pid", Sys.getpid()) , open = "a" )
          write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
          close( conn )

          # also write full seq (DNA and AA) to a file
          fullsequences.filename <- paste0("results_out/",vcf.filename,"_fullSequences.txt")
          
          outputLine <- paste(headerLine,
                              "<original DNA sequence (full ORF)",
                              seq.old,
                              "<mutated DNA sequence (full ORF)",
                              seq,
                              "<original peptide sequence (full ORF)",
                              as.character(AAseq.old),
                              "<mutated peptide sequence (full ORF)",
                              as.character(AAseq.new),
                              sep="&NEWLINE&")
          outputLine <- paste0("&START&",outputLine,"&END&","\n")
          
          conn <- file( paste0(fullsequences.filename ,"_pid", Sys.getpid()) , open = "a" )
          write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
          close( conn )
        } 
      }
      
      # Check for synonymous mutation (only in substitution), set flag for outprints ----
      synonymousFlag <- FALSE #init synonymous mutation flag
      AA.pair <- AACompare[, c(colnames(AACompare) == paste0(mutAApos,"*")), with=F] # get * column
  
      # check for synonymous amino acid on mutAApos* for substitutions or del+ins (which are not on splice acceptor or donor)
      if(AA.pair[1] == AA.pair[2] & (mutType=="sub" | isInDel) & !isSpliceAD){
        synonymousFlag <- TRUE
      }
      
      # Outprints -----

      # get amino acid sequence output (for the current ORFid, which is always 1 for every coding)
      AAseq.new.outprint <- AA.list[ORFid]
      
      # just a sanity check (should never happen)
      if(AAseq.new.outprint != AA.list[ORFid]){
        stop("ERROR1: debug this")
        
        # where are stop codons?
        gregexpr(pattern ="\\*",AAseq.new)
        ts <- as.character(subseq(AAseq.new, 1,875-1)) #-1 because stop codon
        identical(AA.list, ts)
        
      }
      
      # start outprints
      if (synonymousFlag) { # check if mutation is synonymous
        
        synonymous.filename <- paste0("results_out/",vcf.filename,"_synonymous.txt")
        
        outputLine <- paste0("&START&",headerLine, "&NEWLINE&",as.character(AAseq.new.outprint),"&END&", "\n")
        
        conn <- file( paste0(synonymous.filename ,"_pid", Sys.getpid()) , open = "a" )
        write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
        close( conn )

        if(debugMode){warn(paste0("synonymous mutation found for: ", tscrID))}
        
      } else { # else we have a relevant mutation which we print
        # write to results file
        
        outputLine <- paste0("&START&",headerLine, "&NEWLINE&",as.character(AAseq.new.outprint),"&END&", "\n")
        
        conn <- file( paste0(outfile.name ,"_pid", Sys.getpid()) , open = "a" )
        write.table( data.table(outputLine) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
        close( conn )
      }
      
      # End of large loop -----
      
      # reset the original variables for the ORF loop
        seq <- copy(xsave.seq)
        seq.old <- copy(xsave.seq.old)
        mutStart <- copy(xsave.mutStart)
        mutEnd <- copy(xsave.mutEnd)
        mutAApos <- copy(xsave.mutAApos)
        mutStart.firstSave <- xsave.mutStart.firstSave
        mutDT <- copy(xsave.mutDT)
    }
    
    # Free up memory ------
    gc(verbose = F)
    
  } 
}

# MAIN FUNCTION WRAPPER
RunVCFtranslate <- function(input.file){
  
  print(paste0("INPUT: ", input.file))
  start_time <- Sys.time()
  
  input.list <- Stratify.VCF(input.file)
  #input.list[[3]] # data for protein-coding variants
  #input.list[[3]] # data for non-protein-coding variants
  #input.list[[3]] # summary of selection and filtering
  
  vcf.DT <- input.list[[1]] # coding
  if(nrow(vcf.DT) > 0){
    print("Translating VCF for coding transcripts")
    VCFtranslate(vcf.DT, input.file, seqType.input="coding")
    print ("Finished translating coding transcripts")
  } else {
    print("All VCF lines for coding were filtered out")
  }
  
  vcf.DT <- input.list[[2]] # non protein coding
  if(nrow(vcf.DT) > 0){
    print("Translating VCF for non-coding transcripts")
    VCFtranslate(vcf.DT, input.file, seqType.input="cdna")
    print ("Finished translating non-coding transcripts")
  } else {
    print("All VCF lines for non-coding were filtered out")
  }
  
  # combine all temp files 
  CombineAndCleanup(input.file)
  # generate final files
  PostprocessOutput(input.file)
  
  end_time <- Sys.time()
  print(end_time - start_time)
}

# SUPPORT FUNCTIONS
debugMode <<- F # print mutation in console
debugPrintMode <<- F #  print mutation to fullReport.txt file

# fast cutting of a character sequence into tripplets + naming
FastCutTripplets <- function(tempSeq){
  sst <- strsplit(tempSeq, "")[[1]]
  out <- paste0(sst[c(TRUE, FALSE,FALSE)], sst[c(FALSE, TRUE,FALSE)], sst[c(FALSE, FALSE, TRUE)])
  names(out) <- 1:length(out) # set tripplet position as name
  
  # apparently this does also repeat the last element even if it is not a tripplet, so we remove that difference
  lenDif <- (length(out)*3)-nchar(tempSeq)
  if(lenDif != 0){
    lastTrip <- out[length(out)] # get last tripplet
    out[length(out)] <- substr(lastTrip, 1, 3-lenDif) # remove the difference 
  }

  # they should always be identical, but lets double check it
  #identical(out, out.old)
  Tripplet.lastBP <- out[length(out)]
  
  # get the last x bp of tempSeq corresponding to the last tripplet of tripplet list
  True.LastBP <- substr(tempSeq, nchar(tempSeq)-nchar(Tripplet.lastBP)+1, nchar(tempSeq))
  if(Tripplet.lastBP != True.LastBP){
    stop("New function for tripplet cutting is buggy!")
  }
  
  return(out)
}

# slightly modified version of FastCutTripplets
CutTripplets <- function(seq,StartCodonPos){
  tempSeq <- subseq(seq, StartCodonPos,nchar(seq))
  
  out <- FastCutTripplets(tempSeq)
  return(out)
}

# Find the ORFs (slower, old version)
DefineORFs <- function(seq,mutPos,seqType.input,isExtended, isSpliceAD){

  # List of possible start codons (without frame!)
  ATG.list <- as.data.table(vmatchPattern("ATG", seq)) 
  start.list <- ATG.list[,start]
  
  ORF.DT <- data.table() # Init return table
  
  #Quick fix: for coding, drop all except start codon 1 
  #but only IF it was not UTR extended (isExtended in loop) due to start/stop loss
  if(seqType.input == "coding" & !isExtended){
    start.list <- start.list[1]
  }
  
  # only ORF startPos which are smaller (before) the mutPos can result in a valid ORF AROUND the mutPos
  start.list <- start.list[start.list <= mutPos]
  
  # for every possible start codon, get all possible stop codons within the same frame (tripplets)
  for (StartCodonPos in start.list){
    
    # cut sequence infront of ATG, keep everything after it
    tempSeq <- subseq(seq, StartCodonPos,nchar(seq))
    
    # Get Tripplets as list (fast new method)
    out <- FastCutTripplets(tempSeq)
    
    # get list of stop positions
    endCodon.list <- out[out=="TGA" | out =="TAA" | out == "TAG"]
    
    # check if there is a stop codon, skip to next if not (invalid ORF)
    if(isEmpty(endCodon.list)){
      next
    }
    
    # only the smallest value is valid since elsewise there would be a stop codon within an ORF
    stopCodon <- min(as.numeric(names(endCodon.list)))
    
    # Find sequence position from tripplet position
    StartCodonSeq <- StartCodonPos #start codon position in original sequence
    StopCodonSeq <- StartCodonPos-1+stopCodon*3 #stop codon position in original sequence

    #get ORF coordinates (start codon given by large loop)
    ORF.tripplets <- out[1:stopCodon]
    ORF.seq <- paste(ORF.tripplets, collapse="") # this is the actual sequence, but we need the positions to work on the mutated seq later 
    
    temp <- data.table(startPos=StartCodonSeq, endPos=StopCodonSeq, ORF.width=nchar(ORF.seq))
    ORF.DT <- rbind(ORF.DT,temp)
  }
  
  if(nrow(ORF.DT) == 0){
    return(NULL)
  }
  
  ORF.DT <- ORF.DT[ORF.width >= 18] # at least 6 amino acids 
  
  # only keep ORFs around mutation (start checked above, but now endPos)
  # except it is a splice donor / acceptor mutation! because we can have a premature STOP within the retained intron
  # this does work for protein codings, but we have to do some filtering for the non-coding (else we get many fragment ORFs)
  
  if(!isSpliceAD){
    # normal case, no splice region mutation (only keep ORFs around mutation)
    ORF.DT <- ORF.DT[startPos <= mutPos & endPos >= mutPos]
  } else if (seqType.input=="coding"){
    # do nothing, keep the found ORF 
    # which is probably a shorter version due to premature stop, but not necessarily if there is no stop in intron
    # we catch this in main function
  } else {
    # for non-coding transcripts with spliceAD mutation and intron inclusion, we need to filter out irrelevant ORFs not around prev exon
    # prevExonPos.intronInclusion is the cdsPos, which is the relative position of the previous exon (the one before the included intron)
    
    ORF.DT <- ORF.DT[startPos <= prevExonPos.intronInclusion & endPos >= prevExonPos.intronInclusion]
    # we do not need to catch donor mutations, because the mutation should be prevExon +-2, so will be within ORF
  }
  
  if(nrow(ORF.DT) == 0){
    return(NULL)
  }
  
  # same stop codons = same frame --> select longest ORF (all others are elements of this long one)
  outORF.DT <- ORF.DT[, min(startPos), by="endPos"]
  setcolorder(outORF.DT, c(2,1))
  setnames(outORF.DT, "V1", "startPos")
  
  outORF.DT
  return(outORF.DT)
}

# Support function for IncludeMutRetainedIntron, flipping coordinates for -strand
FlipStartStop.IfReverseStrand <- function(DT){
  if(DT[1, strand == -1]){
    DT[, tempEND := exon_chrom_end]
    DT[, exon_chrom_end := exon_chrom_start]
    DT[, exon_chrom_start := tempEND]
    DT[, tempEND := NULL]
    
    setnames(DT, "transcript_start", "TEMP_start")
    setnames(DT, "transcript_end", "transcript_start")
    setnames(DT, "TEMP_start", "transcript_end")
    
  }
  setcolorder(DT, c('ensembl_transcript_id','ensembl_exon_id','exon_chrom_start','exon_chrom_end','strand','rank', 
                    'transcript_start','transcript_end'))
  #print(DT)
  return(DT)
}

# Large function to include intron for acceptor donor splice site mutations
IncludeMutRetainedIntron <- function(mutDT,seqType.input){
  # Function block ----
  tscrID <- mutDT[, gsub("(.*)\\.\\d.*","\\1",FEATUREID)]
  
  # we need to reset cdsPos
  cdsPos <- as.numeric(mutDT[, gsub("..(\\d*).*","\\1",HGVS_C)])
  
  # for snps (X>Y) it will always be acceptor pos-1 and donor pos+1 (also for ins and del?)
  
  # now we need -1 for the following cases
  # c.13243-2A>C acceptor
  # c.178-1G>A acceptor (both strands)
  # c.350-5_353delCACAGATGC acceptor (forward)
  # c.296-8_296delGCCCTAAGG acceptor (reverse)
  # c.4-4_4delCCAGG acceptor (forward)
  # c.938-4_938-3dupTT acceptor
  
  # acceptor need +1, donor NOT
  if(grepl("-", mutDT$HGVS_C)){
    cdsPos <- cdsPos-1
  }
  
  # but keep cdsPos unchanged for:
  # c.4149+1G>T donor (reverse)
  # c.144+1G>C donor
  # c.468+1G>T donor (forward)
  # c.1271_1271+1insAA donor
  # c.2944+1_2944+2delGT donor
  
  # are there cases like this?
  # ENST00000442948.3 c.313-1_313delAT (intron is deleted, actually another wild type but whatever)
  # ENST00000538324.2 c.259-1_259insG
  if(grepl("-", mutDT$HGVS_C) && grepl("donor", mutDT$EFFECT)){
    #warning("debug C1")
    write.table("acceptor + donor mutation (small intron?): double check this",file="results_out/warning_cases.txt",append=TRUE, col.names = F, row.names = F, quote = F)
    write.table(vcf.file,file="results_out/warning_cases.txt",append=TRUE, col.names = F, row.names = F, quote = F)
    write.table(mutDT,file="results_out/warning_cases.txt",append=TRUE, col.names = F, row.names = F, quote = F)
    write.table("",file="results_out/warning_cases.txt",append=TRUE, col.names = F, row.names = F, quote = F)
  }
  
  
  # get tscrDT from subselection of large table
  tscrDT <- exonDT[ensembl_transcript_id == tscrID]
  
  Tstrand <- tscrDT[1, strand]
  
  # start end needs to be flipped for reverse strand (annoying biomart)
  FlipStartStop.IfReverseStrand(tscrDT)
  
  # sort by exon rank
  setorder(tscrDT, "rank")
  
  # get exon seq (coding, cdna)
  CodingSeq <- as.character(sequenceDT[ensembl_transcript_id == tscrID, 1])

  # get sequence intron and exon
  seq <- as.character(intronexonSeqs[ensembl_transcript_id == tscrID, 1])
  
  Tstart <- tscrDT[1,transcript_start] 
  Tend <-  tscrDT[1,transcript_end] 
  
  # get mutGenom to check intron/exon hit
  mutGenomPos <- mutDT$POS
  
  isAcceptor <- mutDT[,grepl("acceptor", EFFECT)]   # ACCEPTOR = start of exon (1-2 bp in intron)
  isDonor <-  mutDT[,grepl("donor", EFFECT)]   # DONOR = end of exon (1-2 bp)
  isBoth <- isAcceptor && isDonor # for very small introns we get both tags
  
  # there are examples with missing acceptor donor tag (those are EXCLUDED in stratifyVCF)
  # splice_region_variant&non_coding_transcript_exon_variant (K72C1G.Mutect2.indel.txt, ENST00000520657.1)
  if(!isDonor && !isAcceptor){
    stop(paste0("no tag found for: ", tscrID))
  }
  
  # get dt with all introns
  intronDT <- data.table()
  
  # build intronDT
  if (Tstrand == "1"){
    # row -1 so last exon is skipped
    for(r in 1:(nrow(tscrDT)-1)){
      intron.start <- tscrDT[r, exon_chrom_end+1]
      intron.end <- tscrDT[r+1, exon_chrom_start-1]
      
      # sometimes ensembl skips exon 1 (probably if not used)
      prevExon <- tscrDT[r,rank]
      nextExon <- tscrDT[r+1,rank]
      id <- paste0(prevExon,"-",nextExon)
      
      intronDT <- rbind(intronDT, data.table(intron_id=id, intron_start=intron.start, intron_end=intron.end, 
                                             previous_exon=prevExon, next_exon=nextExon))
    }
  } else if (Tstrand == "-1"){
    # - strand
    
    # row -1 so last exon is skipped
    for(r in 1:(nrow(tscrDT)-1)){
      intron.start <- tscrDT[r, exon_chrom_end-1]
      intron.end <- tscrDT[r+1, exon_chrom_start+1]
      
      # sometimes ensembl skips exon 1 (probably if not used)
      prevExon <- tscrDT[r,rank]
      nextExon <- tscrDT[r+1,rank]
      id <- paste0(prevExon,"-",nextExon)
      
      intronDT <- rbind(intronDT, data.table(intron_id=id, intron_start=intron.start, intron_end=intron.end, 
                                             previous_exon=prevExon, next_exon=nextExon))
    }
  }
  
  # update mutGenomPos if we have both tags (acceptor and donor) due to small intron
  # so we get the intronHit and everything gets easier
  # but this is very experimental, these cases will be reported in warning_cases.txt
  if(isBoth){
    mutGenomPos <- mutGenomPos+1
  }
  
  # check if mutGenomPos (POS) is in intron or exon
  if (Tstrand == "1"){
    intronHit <- intronDT[between(mutGenomPos, intronDT$intron_start, intronDT$intron_end)]
    exonHit <- tscrDT[between(mutGenomPos, tscrDT$exon_chrom_start, tscrDT$exon_chrom_end)]
  } else if (Tstrand == "-1"){
    intronHit <- intronDT[between(mutGenomPos, intronDT$intron_end, intronDT$intron_start)]
    exonHit <- tscrDT[between(mutGenomPos, tscrDT$exon_chrom_end, tscrDT$exon_chrom_start)]
  }
  
  # there should only be only 1 hit, but who knows?
  if(nrow(exonHit) > 1 | nrow(intronHit) > 1) {
    stop("More than 1 position found in function IncludeMutRetainedIntron()")
  }
  
  isIntron <- F
  isExon <- F
  if(nrow(intronHit) == 1){
    isIntron <- T
    exonN <- intronHit$previous_exon
  } else if(nrow(exonHit) == 1){   # it may be exon hit (large in/del over exon boundaries + intron)
    isExon <- T
    #this is fine for donor
    exonN <- exonHit$rank
    
    # for exon hit we need to go 1 back for acceptors
    if(isAcceptor){
      exonN <- (exonHit$rank)-1
    }
  }
  # now we need to got the previous exon if intron
  # if we find an intron hit (between exons above) we automatically have the exonN for all cases (donor/acceptor, forward/reverse strand)
  # exonN IS ALWAYS the previous exon (exon before the matched intron)
  
  # check if exon was found
  if(isEmpty(exonN)){
    warn("Exon for splice site not found.")
    return(NULL)
  }
  
  # for exonN, get the intron to the next exon
  exonStart <- tscrDT[rank == exonN, exon_chrom_start]
  exonEnd <- tscrDT[rank == exonN, exon_chrom_end]
  
  # get exon intron sequences (main function)
  if(Tstrand == "-1"){
    # this is the closes exon 
    seq.exonStart <- Tstart-exonStart+1 #Tend - E4start
    seq.exonEnd <- Tend-exonEnd-1 #Tend - E4end
    exonSeq <- subseq(seq,seq.exonStart,seq.exonEnd) # this is the previous exon
    
    # start of next exon
    next.exonStart <- tscrDT[rank == exonN+1, exon_chrom_start]
    
    # this is intron between exon 4 and exon 5
    seq.intronStart <- Tend-exonEnd
    seq.intronEnd <- Tend-next.exonStart-2
    intronSeq <- subseq(seq,seq.intronStart,seq.intronEnd)
  }
  if(Tstrand == "1"){
    
    # this is the PREVIOUS exon
    seq.exonStart <- abs(Tstart-exonStart)+1 # relative start of sequence
    seq.exonEnd <- abs(Tstart-exonEnd)+1
    exonSeq <- subseq(seq,seq.exonStart,seq.exonEnd)
    
    # no matter if it is donor or acceptor, we have exonN (always previous) and the next one
    
    next.exonStart <- tscrDT[rank == exonN+1, exon_chrom_start] # start of next exon
    
    # this is intron between exonN and next exon
    seq.intronStart <- seq.exonEnd+1
    seq.intronEnd <- abs(Tstart-next.exonStart)
    intronSeq <- subseq(seq,seq.intronStart,seq.intronEnd)
  }
  
  # should not be needed, but checks for *
  if(is.na(cdsPos)) { # should not be needed, because main function skips invaid with next (this is mainly for debugging)
    seq <- NULL
  } else if(nchar(CodingSeq) < cdsPos){
    seq <- NULL
  } else {
    # set new sequence (combined exon intron exon)
    seq <- paste0(subseq(CodingSeq,1,cdsPos),intronSeq,subseq(CodingSeq,cdsPos+1,nchar(CodingSeq)))
  }
  
  # assign to global (for DefineUTR function to check if the ORF is around the previous exon)
  assign("prevExonPos.intronInclusion", cdsPos, envir=.GlobalEnv)
  
  # DEBUG outprints ----
  debugOutprint <- F
  
  Included.Intron <- paste0(exonN,"-",exonN+1)
  Intron.Length <- nchar(intronSeq) 
  if(debugMode){print(paste0("Including intron ", Included.Intron))}
  
  if(debugOutprint == T){
    
    # for cds and the ends we need a smaller print window
    printWindow <- 10
    if(printWindow > cdsPos){
      printWindow <- 1
    }
    
    print(mutDT)
    print(paste0("exonN (previous) is: ",exonN, " (length ", nchar(exonSeq), ")"))
    print(paste0(subseq(exonSeq, 1,printWindow), "...",subseq(exonSeq, nchar(exonSeq)-printWindow,nchar(exonSeq)))) # this is with UTR, but our seq is not
    
    print("First part before Intron")
    tmpSeq <- subseq(CodingSeq,1,cdsPos)
    # for very small exons (< 10 bp) we need a smaller window
    if(printWindow > nchar(tmpSeq)){
      print(tmpSeq)
    } else {
      print(paste0(subseq(tmpSeq, 1,printWindow), "...",subseq(tmpSeq, nchar(tmpSeq)-printWindow,nchar(tmpSeq))))
    }
    
    Included.Intron <- paste0(exonN,"-",exonN+1)
    print(paste0("INTRON ", Included.Intron, " (length ", Intron.Length, ")"))
    
    # for very small introns (< 10 bp) we need a smaller window
    if(printWindow > Intron.Length){
      print(intronSeq)
    } else {
      print(paste0(subseq(intronSeq, 1,printWindow), "...",subseq(intronSeq, nchar(intronSeq)-printWindow,nchar(intronSeq))))
    }
    
    print(paste0("NEXT EXON: ", exonN+1, " until end of CDS"))
    tempseqx <- subseq(CodingSeq,cdsPos+1,nchar(CodingSeq))
    
    # for very small exons (< 10 bp) we need a smaller window
    if(printWindow > nchar(tempseqx)){
      print(tempseqx)
    } else {
      print(paste0(subseq(tempseqx, 1,printWindow), "...",subseq(tempseqx, nchar(tempseqx)-printWindow,nchar(tempseqx))))
    }
    
    
    if(tscrDT[rank == min(rank), rank] == exonN || tscrDT[rank == max(rank), rank] == exonN+1){
      print("UTR not included in exon for the output seq, this is just to debug the exon finding!")
    }
    
  }
  
  #return without print ----
  intron.length <- nchar(intronSeq)
  invisible(list(seq,intron.length,Included.Intron))
}

# Support the building of vcf.DT, get strand and biotype information
GetBiotype <- function(DT){
  IDs <- gsub("\\..*","", DT[, FEATUREID])
  
  tscrDT <- getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','transcript_biotype'), 
                  filters = 'ensembl_transcript_id', values=IDs, mart=db)
  
  tscrDT <- data.table(tscrDT)
  tscrDT <- unique(tscrDT) # this is important! getBM somehow sometimes returns duplicate rows, which are both found in loop and cause an error
  
  # insert NAs if ensembl entry not found
  if(nrow(tscrDT)==0){
    tscrDT <- rbind(tscrDT, list("NA","NA","NA"))
    #sapply(tscrDT, class)
  }
  
  
  DT[, ensembl_transcript_id := gsub("\\..*","", FEATUREID)]
  tscrDT <- merge(tscrDT, DT)
  
  tscrDT[, unique(transcript_biotype)]
  
  return(tscrDT)
}

# Function to return warnings AS THEY HAPPEN in loop
warn <- function(x) {
  suppressWarnings(w <- warning(x))
  message(w)
}

# Find ORFs and return AA sequnce
Postprocess.Mutation <- function(seq, type,isExtended,mutStart,isSpliceAD){
  ORF.DT <- DefineORFs(seq,mutStart,type,isExtended,isSpliceAD)
  
  if(is.null(ORF.DT)){
    return(NULL)
  }
  
  # for coding, only select longest ORF
  if(type=="coding"){
    ORF.DT[, width := endPos-startPos]
    ORF.DT<- ORF.DT[which.max(width)]
  }
  
  ORFseqAA.list <- c()
  ORFseq.list <- c()
  for(i in 1:nrow(ORF.DT)){
    ORF.start <- ORF.DT[i,startPos]
    ORF.end <- ORF.DT[i,endPos]
    ORFseq <- subseq(seq,ORF.start,ORF.end)
    ORFseqString <- DNAString(ORFseq)
    ORFseqAA <- as.character(translate(ORFseqString, if.fuzzy.codon = "solve"))
    
    ORFseqAA.list <- c(ORFseqAA.list,ORFseqAA)
    ORFseq.list <- c(ORFseq.list,ORFseq)
  }
  returnlist <- list(ORFseq.list,ORFseqAA.list)
  return(returnlist)
}

Stratify.VCF <- function(vcf.file){
  vcf.raw <- fread(vcf.file)
  
  names(vcf.raw) <- gsub("ANN[*].","",names(vcf.raw), fixed=T)
  vcf.raw[, FEATUREID.noversion := gsub("\\..*","", FEATUREID)]
  
  # rename some columns: GEN[Tumor].RD --> TumorRD
  colnames(vcf.raw) <- gsub("GEN\\[(.*)\\].(.*)","\\1\\2",colnames(vcf.raw))
  
  # we need at least NA columns for
  if(!"CDS_LEN" %in% colnames(vcf.raw)){
    vcf.raw[, CDS_LEN := NA]
  }
  if(!"AA_LEN" %in% colnames(vcf.raw)){
    vcf.raw[, AA_LEN := NA]
  }
  if(!"CDS_POS" %in% colnames(vcf.raw)){
    vcf.raw[, CDS_POS := NA]
  }
  if(!"AA_POS" %in% colnames(vcf.raw)){
    vcf.raw[, AA_POS := NA]
  }
  # set NA columns for SNP information (G5 database below 5% is true, allele counts, total counts AN, allele frequency AF)
  if(!"G5" %in% colnames(vcf.raw)){
    vcf.raw[, G5 := NA]
  }
  if(!"AC" %in% colnames(vcf.raw)){
    vcf.raw[, AC := NA]
  }
  if(!"AN" %in% colnames(vcf.raw)){
    vcf.raw[, AN := NA]
  }
  if(!"AF" %in% colnames(vcf.raw)){
    vcf.raw[, AF := NA]
  }
  if(!"TumorVF" %in% colnames(vcf.raw)){
    vcf.raw[, TumorVF := NA]
  }
  if(!"TumorRD" %in% colnames(vcf.raw)){
    vcf.raw[, TumorRD := NA]
  }
  if(!"TumorAD" %in% colnames(vcf.raw)){
    vcf.raw[, TumorAD := NA]
  } else {
    # FILTER
    vcf.raw <- vcf.raw[TumorAD > 0]
  }
  if(!"NormalRD" %in% colnames(vcf.raw)){
    vcf.raw[, NormalRD := NA]
  }
  if(!"NormalAD" %in% colnames(vcf.raw)){
    vcf.raw[, NormalAD := NA]
  } else {
    # the mutation should not be present in reference, so filter on GEN[Normal].AD
    vcf.raw <- vcf.raw[NormalAD <= 1]
  }
  if(!"RESCUED" %in% colnames(vcf.raw)){
    vcf.raw[, RESCUED := NA]
  }
  
  # only select relevant columns
  vcf.raw <- vcf.raw[, .(CHROM, POS, REF, ALT, GENE, EFFECT, IMPACT, FEATUREID, HGVS_C, HGVS_P, CDS_POS, CDS_LEN, AA_POS, AA_LEN, FEATUREID.noversion,
                         G5, AC, AN, AF, TumorVF, TumorRD, TumorAD, NormalAD, NormalRD, patientID, SOURCE, RESCUED)]
  
  # all available effects before filtering
  stratDT <- data.table(EFFECT=vcf.raw[, unique(EFFECT)])
  
  # remove empty effects
  vcf.raw <- vcf.raw[EFFECT != ""]
  
  # set missing HGVS_P to NA
  vcf.raw[HGVS_P == "", HGVS_P := "NA"]
  
  # rescue acceptor donor splice
  saveDonAc <- vcf.raw[grep("donor|acceptor", EFFECT)]
  
  # exclude irrelevant EFFECTS (by grep, so all combination and subtypes are found)
  blacklist <- c("downstream_gene_variant","intergenic", "upstream_gene_variant", "intron", "UTR", "intragenic",
                 " ","TF_binding_site_variant", "sequence_feature","synonymous_variant", "fusion", "miRNA", "snoRNA", "structural_interaction_variant")
  blacklist <- paste(blacklist,collapse="|")
  hits <- grep(blacklist, vcf.raw$EFFECT)
  vcf.raw[!hits, unique(EFFECT)]
  excludedElements <- vcf.raw[hits, unique(EFFECT)]
  vcf.raw <- vcf.raw[!hits]
  
  # manually remove splice, full word
  vcf.raw <- vcf.raw[EFFECT != "splice_region_variant"]
  # same for stop retained
  vcf.raw <- vcf.raw[EFFECT != "stop_retained_variant"]
  
  # also remove every combination of "splice_region_variant", because they are full intronic 
  # splice_region_variant&non_coding_transcript_exon_variant
  # disruptive_inframe_deletion&splice_region_variant
  # frameshift_variant&splice_region_variant
  vcf.raw <- vcf.raw[!grep("splice_region_variant", EFFECT)]
  
  excludedElements <- c(excludedElements, "splice_region_variant", "stop_retained_variant")
  
  # check donor acceptor splice mutation for retained introns
  if(nrow(saveDonAc) != 0){
    saveDonAc.bio <- GetBiotype(saveDonAc)
    saveDonAc[, ensembl_transcript_id := NULL]
    
    pc.list <- saveDonAc.bio[transcript_biotype == "protein_coding", FEATUREID]
    saveDonAc.pc <- saveDonAc[FEATUREID %in% pc.list]
    
    npc.list <- saveDonAc.bio[transcript_biotype != "protein_coding", FEATUREID]
    saveDonAc.npc <- saveDonAc[FEATUREID %in% npc.list]
  }
  
  
  # set coding variants
  # structural_interaction_variant ??
  coding.keywords <- c("frameshift_variant","missense_variant",
                       "stop", "start", "coding_sequence_variant", "inframe", "amino_acid", "splice_acceptor", "splice_donor")
  coding.keywords <- paste(coding.keywords,collapse="|")
  c.hits <- grep(coding.keywords, vcf.raw$EFFECT)
  vcf.coding <- vcf.raw[c.hits]
  vcf.coding[CDS_POS == "-1"]
  
  #rebind the splice donor acceptor
  if(nrow(saveDonAc) != 0){
    vcf.coding <- rbind(vcf.coding,saveDonAc.pc)
  }
  
  
  # set non-protein coding variants
  NOTcoding.keywords <- c("non_coding")
  NOTcoding.keywords <- paste(NOTcoding.keywords,collapse="|")
  nc.hits <- grep(NOTcoding.keywords, vcf.raw$EFFECT)
  vcf.Ncoding <- vcf.raw[nc.hits]
  
  #rebind the splice donor acceptor
  if(nrow(saveDonAc) != 0){
    vcf.Ncoding <- rbind(vcf.Ncoding,saveDonAc.npc)
  }
  
  # Check if keywords were missed
  foundHits <- c(grep(NOTcoding.keywords, vcf.raw$EFFECT),grep(coding.keywords, vcf.raw$EFFECT))
  if(nrow(vcf.raw[!foundHits])!=0){
    print("MISSED effect keywords:")
    print(vcf.raw[!foundHits])
  }
  
  # Show stratification
  stratDT[EFFECT=="", biotype := "excluded"]
  stratDT[EFFECT %in% excludedElements, biotype := "excluded"]
  stratDT[EFFECT %in% excludedElements, biotype := "excluded"]
  stratDT[EFFECT %in% vcf.raw[c.hits, unique(EFFECT)], biotype := "coding"]
  stratDT[EFFECT %in% vcf.raw[nc.hits, unique(EFFECT)], biotype := "non-P-coding"]
  stratDT[is.na(biotype), biotype := "included"]
  setorder(stratDT, EFFECT)
  
  outList <- list(vcf.coding,vcf.Ncoding,stratDT)
}

# merge all files from the differend PIDs
CombineAndCleanup <- function(vcf.file){
  
  vcf.filename <- gsub("(.*)\\..*","\\1", basename(vcf.file))
  
  # we will move all temp files into a temp dir in case we need them
  tempdirname <- paste0("results_out/tempFiles/",vcf.filename, "_", format(Sys.time(), "%Y-%m-%d_%H:%M"))
  dir.create(tempdirname, recursive = T)
  
  # merge all created tempfiles
  all.patterns <- c("results.tempfasta_pid","synonymous.txt_pid", "fullSequences.txt_pid", "fullReport.txt_pid", "report_PrematureIntronStop.txt_pid")
  
  # only select existing file patterns
  all.outfiles <- list.files("results_out/", pattern = "pid")
  all.outfiles <- unique(gsub(".*_(.*_pid).*", "\\1", all.outfiles))
  all.outfiles <- paste(all.outfiles, collapse = "|")
  all.patterns <- all.patterns[grep(all.outfiles, all.patterns)]
  
  for(mypattern in all.patterns){
    all.files <- list.files("results_out/", pattern=mypattern, full.names = T)
    
    # just in case also grep vcf filename
    all.files <- all.files[grep(vcf.filename, all.files)]
    
    mylist <- lapply(all.files, fread, header=F, sep=NULL, blank.lines.skip=T)
    mydata <- rbindlist( mylist )
    
    outfilename <- gsub("_pid.*","",all.files[1])
    fwrite(mydata, outfilename, sep="\t", col.names = F, quote = F, row.names = F)
    
    # move all files
    for(myfile in all.files){
      newfilename <- paste0(tempdirname, "/", basename(myfile))
      file.rename(myfile, newfilename)
    }
  }
  # we actually do not need the tmp files anymore (+tempfastq)
  unlink("results_out/tempFiles/", recursive = T)
}

# remove sanity check &START& &END&, replace &NEWLINE& by a true newline and generate summary file and some statistics
PostprocessOutput <- function(vcf.file){
  vcf.filename <- gsub("(.*)\\..*","\\1", basename(vcf.file))
  tempfasta.filename <- paste0("results_out/",vcf.filename,"_results.tempfasta")
  truefasta.filename <- paste0("results_out/",vcf.filename,"_results.fasta")
  summary.filename <- paste0("results_out/",vcf.filename,"_results_summary.tsv")
  summarystatistics.filename <- paste0("results_out/",vcf.filename,"_results_summary_statistics.txt")
  
  resdt <- fread(tempfasta.filename, sep="\n", header = F, blank.lines.skip=T)
  
  # check for correct format
  # every row should start with &START& and end with &END&
  if(nrow(resdt[!grep("^&START&", V1)]) > 0 | nrow(resdt[!grep("&END&$", V1)]) > 0){
    
    resdt[!grep("^&START&", V1)]
    resdt[!grep("&END&$", V1)]
    
    stop("PROBLEM 1")
  }
  
  resdt[, rowID := 1:.N]
  
  # remove START and END tags
  resdt <- resdt[, gsub("&START&", "", V1), by="rowID"]
  resdt <- resdt[, gsub("&END&", "", V1), by="rowID"]
  
  # now there should be no more START END left
  if(nrow(resdt[grep("&START&", V1)]) > 0 | nrow(resdt[grep("&END&$", V1)]) > 0){
    stop("PROBLEM 2")
  }
  
  resdt[, rowID := NULL]
  
  # split into newline at &NEWLINE&
  fastaOut <- cSplit(resdt, "V1", "&NEWLINE&", direction = "long", type.convert = F)
  
  # print results file
  fwrite(fastaOut, truefasta.filename, col.names = F, quote = F)

  # also print summary file
  summaryDT <- fastaOut[grep(">", V1)]
  
  summaryDT <- cSplit(summaryDT, "V1", "|", type.convert = F)
  
  VCFtranslate.headernames <- c("CHROM", "POS", "REF", "ALT", "GENE", "EFFECT", "IMPACT", "FEATUREID", "HGVS_C", "HGVS_P",
                                "CDS_POS", "AA_POS", "G5", "AC", "AN", "AF", "TumorVF", "TumorAD", "TumorRD", "NormalAD", "NormalRD",
                                "patientID", "SOURCE", "RESCUED", "seqType","geneBiotype", "transcriptBiotype", "versionControl", "ORF", "NOTE")
  
  names(summaryDT) <- VCFtranslate.headernames

  # change factors to character
  changeCols <- colnames(summaryDT)[which(as.vector(summaryDT[,lapply(.SD, class)]) == "factor")]
  summaryDT[,(changeCols):= lapply(.SD, as.character), .SDcols = changeCols]
  
  summaryDT[, CHROM := gsub(">", "", CHROM)]
  
  # add mutation type
  summaryDT[grepl(">", HGVS_C), mutationType := "substitution"]
  summaryDT[grepl("dup", HGVS_C), mutationType := "insertion"]
  summaryDT[grepl("dup", HGVS_C), mutationSubType := "duplication"]
  summaryDT[grepl("ins", HGVS_C) & !grepl("del", HGVS_C), mutationType := "insertion"]
  summaryDT[grepl("del", HGVS_C) & !grepl("ins", HGVS_C), mutationType := "deletion"]
  summaryDT[grepl("del", HGVS_C) & grepl("ins", HGVS_C), mutationType := "indel"] 
  
  # check if indel is multi-substitution (DNV, TNV, ...)
  indels <- unique(summaryDT[mutationType == "indel", "HGVS_C"])
  indels[, del := gsub(".*del([A-Z]*).*", "\\1", HGVS_C)]
  indels[, ins := gsub(".*ins([A-Z]*).*", "\\1", HGVS_C)]
  indels[, delN := nchar(del)]
  indels[, insN := nchar(ins)]
  indels[delN == insN, mutType := "multi-substitution"]
  indels[, mutSubType := paste0(delN, "NV")]
  summaryDT <- merge(summaryDT, indels[, .(HGVS_C, mutType, mutSubType)], by="HGVS_C", all.x = T)
  summaryDT[HGVS_C %in% indels$HGVS_C, `:=` (mutationType = mutType, mutationSubType = mutSubType), by="POS"]
  summaryDT <- summaryDT[, -c("mutType", "mutSubType")]
  
  setcolorder(summaryDT, VCFtranslate.headernames)
  
  # print summary file
  fwrite(summaryDT, summary.filename, sep="\t", col.names = T)
  
  # some statistics
  conn <- file( summarystatistics.filename , open = "a" )
  
  posDT <- unique(summaryDT[, .(CHROM, POS, EFFECT, mutationType, mutationSubType, SOURCE)])
  
  uniquemutations <- nrow(unique(posDT[, .(CHROM, POS)]))
  write(paste0("number of unique mutations = ", uniquemutations),
        conn, append=F)
  write("note: these are unique genomic mutation sites, so the same mutation for multiple transcripts = 1 mutation (RNA+exome combined)", conn, append=T)
  
  write("\n", conn, append=T)
  write("with these mutation called by:", conn, append=T)
  sourceDT <- posDT[, paste0(SOURCE, collapse=";"), by=.(CHROM, POS)]
  sourceDT[, calledby := V1]
  sourceDT[grep(";", V1), calledby := "both"]
  sourceTAB <- sort(table(sourceDT[, "calledby"]), decreasing = T)
  write.table(sourceTAB, conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
  
  effectDT <- unique(posDT[, .(CHROM, POS, EFFECT)]) # remove the duplicared entries due to caller etc
  effectDT <- effectDT[, paste0(EFFECT, collapse="+++"), by=.(CHROM,POS)] # collapse multiple effects for 1 mutation
  effectTAB <- sort(table(effectDT$V1), decreasing = T)
  
  write("\n", conn, append=T)
  write("with the mutation classes distributed like:", conn, append=T)
  write.table( effectTAB , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
  write("note: mutation sites with overlapping genes/transcripts can have multiple classes, which are combined with +++", conn, append=T)
  
  muttypeDT <- unique(posDT[, .(CHROM, POS, mutationType)]) # remove the duplicared entries due to caller etc
  muttypeDT <- muttypeDT[, paste0(mutationType, collapse="+++"), by=.(CHROM,POS)] # collapse multiple effects for 1 mutation
  muttypeTAB <- sort(table(muttypeDT$V1), decreasing = T)
  
  write("\n", conn, append=T)
  write("with the mutation types distributed like:", conn, append=T)
  write.table( muttypeTAB , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
  write("note: mutation sites with overlapping genes/transcripts can have multiple types, which are combined with +++", conn, append=T)
  write("note: 1 genomic position can have multiple different mutations, e.g. G>A and GG>AA with different amino acid changes", conn, append=T)

  genesDT <- unique(summaryDT[, .(GENE, geneBiotype)])
  transcrDT <- unique(summaryDT[, .(FEATUREID, transcriptBiotype)])
  
  write("\n", conn, append=T)
  write(paste0("number of affected genes = ", length(genesDT$GENE), "\n"), conn, append=T)
  
  write("note: one mutation can affect multiple transcripts/genes (genomic overlaps)", conn, append=T)
  write("note: one transcript/gene can have multiple mutations", conn, append=T)
  
  write("with the biotype distributed like:", conn, append=T)
  write.table( table(genesDT$geneBiotype) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
  
  write("\n", conn, append=T)
  write(paste0("number of affected transcripts = ", length(transcrDT$FEATUREID), "\n"), conn, append=T)
  write("with the biotype distributed like:", conn, append=T)
  write.table( table(transcrDT$transcriptBiotype) , conn , append = TRUE , col.names = FALSE, quote = F, row.names = F)
  
  close( conn )
  
  file.remove(tempfasta.filename)
}










