library(shiny)
library(shinyjs)
library(shinythemes)
library(shinybusy)
library(sangeranalyseR)
library(seqinr)
options(repos = BiocManager::repositories())


ui = fluidPage(
  useShinyjs(),
  theme = shinytheme("cerulean"),
  navbarPage("SequenceR", tabPanel("Home",sidebarPanel(
    numericInput(inputId = "q_ts", label = "Quality threshold", value = 20, min = 1, step = 0.5),
    numericInput(inputId = "gap_p", label = "Gap penalty", value = 18, min = 0),
    fileInput(inputId = "forward_sequence", label = "Forward sequence(s) (.ab1)", multiple = TRUE), 
    fileInput(inputId = "reverse_sequence", label = "Reverse sequence(s) (.ab1)", multiple = TRUE),
    fileInput(inputId = "reference", label = "Reference sequence (.fasta)"),
    actionButton("submit", "Submit", class="btn btn-primary"),
    actionButton("reset", "Clear")
    
  ), mainPanel(add_busy_spinner(spin = "double-bounce", position = 'bottom-right',color = "red"),tabsetPanel(type = "tabs",
                           tabPanel("Alignment", htmlOutput('alignment'),tags$style(type="text/css",
                                                                                    "span { white-space:nowrap; }",
                                                                                    "pre { display: inline-block; }"),
                                    
                                   
                                    verbatimTextOutput("amb_pos"),
                                    downloadButton(outputId = "contig", "Download  contig", class = "btn-block")),
                           
                           tabPanel("Summary", textOutput("coverage"), tableOutput("summary_used"), tableOutput("summary_unused"))
  ))),
  tabPanel("About", 
           div(includeMarkdown("dula.Rmd")),
           h4("Author: Anis Mansouri"),
           h4("E-Mail: anis.mansouri.dz@gmail.com.")))
  
)#,



my_BrowseSeqs= function(myXStringSet,
                        htmlFile=paste(tempdir(),"/myXStringSet.html",sep=""),
                        openURL=interactive(),
                        colorPatterns=TRUE,
                        highlight=NA,
                        patterns=c("-", alphabet(myXStringSet, baseOnly=TRUE)),
                        colors=substring(rainbow(length(patterns), v=0.8, start=0.9, end=0.7), 1, 7),
                        colWidth=Inf,
                        ...) {
  
  # error checking
  if (!is(myXStringSet, "XStringSet"))
    stop("myXStringSet must be an XStringSet.")
  if (length(myXStringSet)==0)
    stop("No sequence information to display.")
  type <- switch(class(patterns),
                 `DNAStringSet` = 1L,
                 `RNAStringSet` = 2L,
                 `AAStringSet` = 3L,
                 `list` = -1L,
                 0L)
  if (type > 0L) {
    patterns <- as.character(patterns)
    if (type==1L) { # DNAStringSet
      patterns <- gsub("M",
                       "[ACM]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("R",
                       "[AGR]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("W",
                       "[ATW]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("S",
                       "[CGS]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("Y",
                       "[CTY]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("K",
                       "[GTK]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("V",
                       "[ACGMRSV]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("H",
                       "[ACTMWYH]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("D",
                       "[AGTRWKD]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("B",
                       "[CGTSYKB]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("N",
                       "[ACGTMRWSYKVHDBN]",
                       patterns,
                       fixed=TRUE)
    } else if (type==2L) { # RNAStringSet
      patterns <- gsub("M",
                       "[ACM]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("R",
                       "[AGR]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("W",
                       "[AUW]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("S",
                       "[CGS]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("Y",
                       "[CUY]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("K",
                       "[GUK]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("V",
                       "[ACGMRSV]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("H",
                       "[ACUMWYH]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("D",
                       "[AGURWKD]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("B",
                       "[CGUSYKB]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("N",
                       "[ACGUMRWSYKVHDBN]",
                       patterns,
                       fixed=TRUE)
    } else { # AAStringSet
      patterns <- gsub("B",
                       "[NDB]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("Z",
                       "[QEZ]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("J",
                       "[ILJ]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("X",
                       "[ARNDCQEGHILKMFPSTWYVUOBJZX]",
                       patterns,
                       fixed=TRUE)
      patterns <- gsub("*",
                       "\\*",
                       patterns,
                       fixed=TRUE)
    }
  } else if (type==0) { # character vector
    if (!is.null(patterns) && !is.character(patterns))
      stop("patterns must be a character vector.")
    if (any(grepl("=|\"|<|>|[1-9]|[a-z]", patterns)))
      stop("patterns cannot contain numbers, lower case characters, or the characters (=, <, >, \").")
  }
  if (type < 0) {
    if (length(myXStringSet) != length(patterns))
      stop("patterns is not the same length as myXStringSet.")
  } else {
    w <- which(patterns %in% c("?", "*", "+", "."))
    if (length(w) > 0)
      patterns[w] <- paste("\\", patterns[w], sep="")
    if (length(colors) != length(patterns) || !is.character(colors))
      stop("colors must be a character vector of the same length as patterns.")
  }
  # check that the file exist
  if (is.character(htmlFile)) {
    htmlfile <- file(htmlFile, "w")
    on.exit(close(htmlfile))
  } else {
    stop("'htmlFile' must be a character string or connection.")
  }
  if (!is.logical(openURL) || is.na(openURL))
    stop("openURL must be TRUE or FALSE.")
  if (!is.logical(colorPatterns) && !is.numeric(colorPatterns))
    stop("colorPatterns must be a logical or numeric.")
  if (is.numeric(colorPatterns)) {
    if ((length(colorPatterns) %% 2) == 1 || length(colorPatterns) == 0)
      stop("colorPatterns must specify all start and endpoints.")
    if (any((colorPatterns[seq(2, length(colorPatterns), 2)] - colorPatterns[seq(1, length(colorPatterns), 2)]) < 0))
      stop("colorPatterns specifies a negative range.")
    if (any(colorPatterns <= 0))
      stop("colorPatterns must be a positive numeric.")
    if (length(colorPatterns) > 2)
      if (any((colorPatterns[seq(3, length(colorPatterns), 2)] - colorPatterns[seq(2, length(colorPatterns) - 2, 2)]) <= 0))
        stop("Ranges specified in colorPatterns must be non-overlapping.")
    if (max(colorPatterns) > max(width(myXStringSet)))
      stop("Ranges specified in colorPatterns are out-of-bounds.")
  }
  if (is.numeric(colorPatterns) & !is.infinite(colWidth))
    stop("colWidth must be Inf if colorPatterns is numeric.")
  if (is.null(names(myXStringSet)))
    names(myXStringSet) <- 1:length(myXStringSet)
  names(myXStringSet) <- gsub("\t", " ", names(myXStringSet), fixed=TRUE)
  if (!is.na(highlight)) {
    if (highlight < 0 || highlight > length(myXStringSet))
      stop("highlight must be 0 or the index of a sequence in myXStringSet.")
    if (highlight != floor(highlight))
      stop("highlight be be a whole number.")
  }
  # add a consensus sequence to myXStringSet
  
  l <- length(myXStringSet)
  
  
  # convert the XStringSet to character strings
  html <- as.character(myXStringSet)
  
  if (type < 0) { # record positions of letter colors
    s <- strsplit(html, "", fixed=TRUE)
    n <- nchar(html)
    v <- vector("list", length(html)) # ignore the consensus
    for (i in seq_len(length(html))) {
      if (!is.matrix(patterns[[i]]))
        stop("All elements of the list patterns contain a single matrix.")
      if (nrow(patterns[[i]]) != 3)
        stop("All elements of the list patterns must be a matrix with 3 rows.")
      if (any(patterns[[i]] < 0) || any(patterns[[i]] > 1))
        stop("All elements of patterns[[", i, "]] are not between 0 and 1.")
      v[[i]] <- character(n[i])
      if (any(patterns[[i]] > 1) || any(patterns[[i]] < 0))
        stop("All values of patterns[[", i, "]] must be between 0 and 1.")
      w <- which(s[[i]] %in% c(LETTERS, letters, "*"))
      if (length(w) != ncol(patterns[[i]]))
        stop("The number of columns in patterns[[", i, "]] is different than the number of letters in myXStringSet[", i, "]")
      if (length(w) > 0)
        v[[i]][w] <- apply(patterns[[i]],
                           2,
                           function(x) {
                             rgb(x[1], x[2], x[3])
                           })
      if (is.numeric(colorPatterns)) {
        start <- 0L
        cPs <- c(colorPatterns, length(v[[i]]) + 1)
        for (j in seq_len(length(cPs))) {
          if (j %% 2) {
            end <- cPs[j]
            if (end > length(v[[i]]))
              end <- length(v[[i]]) + 1L
            if (start < end - 1)
              v[[i]][(start + 1L):(end - 1L)] <- ""
          } else {
            start <- cPs[j]
            if (start > length(v[[i]]) + 1)
              break
          }
        }
      }
    }
  }
  
  if (!is.na(highlight)) { # highlight a sequence
    if (highlight==0) {
      highlight <- length(html)
      index <- 1:(length(html) - 1L)
    } else {
      index <- (1:(length(html) - 1L))[-highlight]
    }
    html <- sapply(html, strsplit, split="", fixed=TRUE)
    for (i in index) {
      L <- min(length(html[[highlight]]), length(html[[i]]))
      w <- which(html[[i]][1:L]==html[[highlight]][1:L])
      if (length(w) > 0)
        html[[i]][w] <- "\u00B7"
    }
    html <- sapply(html, paste, collapse="")
  }
  
  # pad shorter sequences with spaces
  maxW <- max(width(myXStringSet))
  if (maxW==0)
    stop("No sequence information to display.")
  if (colWidth > maxW)
    colWidth <- maxW
  for (i in seq_len(l)) {
    html[i] <- paste(html[i],
                     paste(rep(" ", maxW - nchar(html[i])), collapse=""),
                     sep="")
  }
  
  # create a legend that gives position
  if (maxW < 20) {
    if (maxW < 10) {
      counter <- maxW
    } else {
      counter <- 10
    }
  } else {
    counter <- 20
  }
  offset <- (counter - 1) - nchar(maxW)
  if (offset < 0)
    offset <- 0
  legend <- paste(paste(rep(" ", offset), collapse=""), format(seq(counter, maxW, by=counter)), collapse="")
  counter <- ceiling(counter/2)
  tickmarks <- paste(paste(rep("'", counter - 1), collapse=""), rep("|", floor(maxW/counter)), collapse="", sep="")
  tickmarks <- paste(tickmarks, paste(rep("'", maxW - counter*floor(maxW/counter)), collapse=""), sep="")
  
  # split the html into multiple pages
  starts <- seq(1, maxW, colWidth)
  stops <- starts + colWidth - 1L
  stops[length(stops)] <- maxW
  temp <- character((length(html) + 5L)*length(starts))
  count <- 1L
  for (i in 1:length(starts)) {
    temp[count:(count + 1L)] <- substring(c(legend, tickmarks), starts[i], stops[i])
    count <- count + 2L
    temp[c(count:(count + length(html) - 2L), count + length(html))] <- substring(html, starts[i], stops[i])
    count <- count + length(html) + 3L
  }
  html <- temp
  
  # add the cumulative nucleotide lengths on the right
  if (length(starts) > 1) {
    lengths <- numeric(length(starts)*l)
    myLengths <- character(length(starts)*(l + 5L))
    for (i in 1:length(starts)) {
      s <- substring(myXStringSet,
                     starts[i],
                     stops[i])
      lengths[((i - 1L)*l + 1L):(i*l)] <- nchar(s) - rowSums(letterFrequency(BStringSet(s),
                                                                             c("-", " ", ".", "+")))
      if (i > 1)
        lengths[((i - 1L)*l + 1L):(i*l)] <- lengths[((i - 1L)*l + 1L):(i*l)] + lengths[((i - 2L)*l + 1L):((i - 1L)*l)]
      myLengths[((i - 1L)*(l + 5L) + 3L):(i*(l + 5L) - 4L)] <- lengths[((i - 1L)*l + 1L):(i*l - 1L)]
      myLengths[(i*(l + 5L) - 2L)] <- lengths[i*l]
    }
  } else {
    lengths <- width(myXStringSet) - rowSums(letterFrequency(myXStringSet, "-"))
    myLengths <- c("", "", lengths[-seq(l, length(lengths), l)], "", lengths[seq(l, length(lengths), l)], "", "")
  }
  
  if (is.numeric(colorPatterns) || colorPatterns) {
    if (type < 0) {
      s <- strsplit(html, "", fixed=TRUE)
      for (i in seq_along(v)) {
        count <- 2L
        for (j in seq_along(starts)) {
          w <- which(v[[i]][starts[j]:stops[j]] != "")
          if (length(w) > 0)
            s[[i + count]][w] <- paste("<span style=\"color:#FFF;background:",
                                       v[[i]][w + starts[j] - 1],
                                       "\">",
                                       s[[i + count]][w],
                                       "</span>",
                                       sep="")
          count <- count + length(myXStringSet) + 5L
        }
      }
      html <- sapply(s, paste, collapse="")
    } else {
      patterns <- paste("(",
                        patterns,
                        ifelse(nchar(patterns)==1, "+)", ")"),
                        sep="")
      classes <- paste("<span class=\"_",
                       seq_along(patterns),
                       "\">\\1</span>",
                       sep="")
      if (is.numeric(colorPatterns)) {
        htm <- substring(html, 0, colorPatterns[1] - 1)
        for (i in seq(1, length(colorPatterns), 2)) {
          htmi <- substring(html, colorPatterns[i], colorPatterns[i + 1])
          for (j in seq_along(colors)) {
            htmi <- gsub(patterns[j],
                         classes[j],
                         htmi)
          }
          end <- ifelse(i==(length(colorPatterns) - 1),
                        max(nchar(html)),
                        colorPatterns[i + 2] - 1)
          
          htm <- paste(htm, htmi, substring(html, colorPatterns[i + 1] + 1, end), sep="")
        }
        html <- htm
      } else {
        for (j in seq_along(colors)) {
          html <- gsub(patterns[j],
                       classes[j],
                       html)
        }
      }
    }
    
    html <- paste(html,
                  myLengths,
                  "", # post-spacing
                  sep="    ")
    
    # add the legend and consensus sequence
    html <- paste("", # pre-spacing
                  format(c("", "", names(myXStringSet)[1:(l-1)], "", "Contig", "", ""), justify="right"),
                  html,
                  sep="    ")
    
    styles <- character()
    for (i in seq_along(colors)) {
      styles <- paste(styles,
                      "span._", i,
                      " {background:", colors[i],
                      "; color:#FFF;} ",
                      sep="")
    }
    styles <- paste("<style type=text/css> ",
                    styles,
                    "</style>",
                    sep="")
    html <- c('<html><head><meta http-equiv="Content-Type" content="text/html; charset=utf-8"></head>',styles,"<pre>", html, "</pre></html>")
  } else {
    html <- paste(html,
                  myLengths,
                  "", # post-spacing
                  sep="    ")
    
    # add the legend and consensus sequence
    html <- paste("", # pre-spacing
                  format(c("", "", names(myXStringSet)[1:(l-1)], "","Contig",  "", ""), justify="right"),	
                  html,
                  sep="    ")
    
    html <- c("<html>","<pre>", html, "</pre></html>")
  }
  
  # replace unicode 'middle dot' with html entity
  html <- gsub("\u00B7", "&#183;", html, fixed=TRUE)
  
  writeLines(html, htmlfile)
  
  
  
  invisible(htmlFile)
}

ends_trim=function(Forward_files, Reverse_files, trim_cutoff5, trim_cutoff3, peak_cutoff = 0.33) {
  used_files=c(Forward_files, Reverse_files)
  read_set=DNAStringSet()
  for (f in 1:length(used_files)) {print(used_files[f])
    if (f<= length(Forward_files)) {
      #load the read
      F_read=makeBaseCalls(sangerseqR::readsangerseq(used_files[f]), ratio = peak_cutoff)
      iupac_F_read=mostConsensusString(AlignSeqs(c(DNAStringSet(F_read@primarySeq),DNAStringSet(F_read@secondarySeq)),verbose=FALSE))
      
      #Covert the read
      raw_F_read=str_split(iupac_F_read$contig, '')
      raw_F_read_scores=ifelse(raw_F_read[[1]]=="A" | raw_F_read[[1]]== "C" | raw_F_read[[1]]=="G" |      raw_F_read[[1]]=="T", 1, -1)
      
      #Find the starting trimming point
      start_point=NULL
      s_criterion=0
      for (i in 1:length(raw_F_read_scores)) {
        s_criterion=s_criterion+raw_F_read_scores[i]
        if (raw_F_read_scores[i]==-1) {s_criterion=0}
        if (s_criterion==trim_cutoff5) {start_point=c(start_point,i-(trim_cutoff5-1))}
      }
      if (is.null(start_point)) {start_point=30} else {
      diffs=diff(start_point)
      if (length(start_point)==1 | length(start_point)==2) {s_p=which(start_point==start_point[1])} else {
        if (diffs[1]>50) {s_p=which(start_point==start_point[1])} else {s_p=which(diffs==sort(diffs, decreasing = TRUE)[1])}}}
      #Find terminal trimmin point
      end_point=NULL
      e_criterion=0
      for (i in 1:length(raw_F_read_scores)) {
        e_criterion=e_criterion+rev(raw_F_read_scores)[i]
        if (rev(raw_F_read_scores)[i]==-1) {e_criterion=0}
        if (e_criterion==trim_cutoff3) {end_point=c(end_point,i-(trim_cutoff3-1))}}
      if (is.null(end_point)) {end_point=30}
      
      #Add the trimmed F read to the read set
      print(c(start_point[s_p], length(raw_F_read_scores)-end_point[1]))
      trimmed_F_read=iupac_F_read$contig[start_point[s_p]:(length(raw_F_read_scores)-end_point[1])]
      read_set=c(read_set, DNAStringSet(trimmed_F_read))}
    
    ##############################################
    
    else {
      
      #load the read
      R_read=makeBaseCalls(sangerseqR::readsangerseq(used_files[f]), ratio = peak_cutoff)
      reverse_p=reverseComplement(R_read@primarySeq)
      reverse_s=reverseComplement(R_read@secondarySeq)
      iupac_R_read=mostConsensusString(AlignSeqs(c(DNAStringSet(reverse_p),DNAStringSet(reverse_s)),verbose=FALSE))
      
      #Covert the read
      raw_R_read=str_split(iupac_R_read$contig, '')
      raw_R_read_scores=ifelse(raw_R_read[[1]]=="A" | raw_R_read[[1]]== "C" | raw_R_read[[1]]=="G" |      raw_R_read[[1]]=="T", 1, -1)
      
      #Find the starting trimming point
      start_point=NULL
      s_criterion=0
      for (i in 1:length(raw_R_read_scores)) {
        s_criterion=s_criterion+raw_R_read_scores[i]
        if (raw_R_read_scores[i]==-1) {s_criterion=0}
        if (s_criterion==trim_cutoff5) {start_point=c(start_point,i-(trim_cutoff5-1))}
      }
      if (is.null(start_point)) {start_point=30} else {
      diffs=diff(start_point)
      if (length(start_point)==1 | length(start_point)==2) {s_p=which(start_point==start_point[1])} else {
        if (diffs[1]>50) {s_p=which(start_point==start_point[1])} else {s_p=which(diffs==sort(diffs, decreasing = TRUE)[1])}}}
      #Find terminal trimmin point
      end_point=NULL
      e_criterion=0
      for (i in 1:length(raw_R_read_scores)) {
        e_criterion=e_criterion+rev(raw_R_read_scores)[i]
        if (rev(raw_R_read_scores)[i]==-1) {e_criterion=0}
        if (e_criterion==trim_cutoff3) {end_point=c(end_point,i-(trim_cutoff3-1))}}
      if (is.null(end_point)) {end_point=30}
      
      #Add the trimmed R read to the read set
      print(c(start_point[s_p], length(raw_R_read_scores)-end_point[1]))
      trimmed_R_read=iupac_R_read$contig[start_point[s_p]:(length(raw_R_read_scores)-end_point[1])]
      read_set=c(read_set, DNAStringSet(trimmed_R_read))}    
    names(read_set)[f]=used_files[f]
    
  }
  
  return(list("readset" = read_set))
}


make.readset.dula <- function(fwd.fnames, rev.fnames, trim = TRUE, trim.cutoff = 0.0001, max.secondary.peaks = NULL, secondary.peak.ratio = 0.33, min.length = 1, processors = NULL){
  
  processors = 1
  
  fwd.dat = parallel::mclapply(fwd.fnames, loadread.dula, trim, trim.cutoff, revcomp = FALSE, max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, min.length = min.length, processors = 1, mc.cores = processors)
  rev.dat = parallel::mclapply(rev.fnames, loadread.dula, trim, trim.cutoff, revcomp = TRUE,  max.secondary.peaks = max.secondary.peaks, secondary.peak.ratio = secondary.peak.ratio, min.length = min.length, processors = 1, mc.cores = processors)
  
  fwd.reads = lapply(fwd.dat, function(x) x[["read"]])
  rev.reads = lapply(rev.dat, function(x) x[["read"]])
  names(fwd.reads) = fwd.fnames
  names(rev.reads) = rev.fnames
  
  # remove the NULL reads
  all.reads = c(fwd.reads, rev.reads)
  all.reads = Filter(Negate(is.null), all.reads)
  readset = DNAStringSet(all.reads)
  
  # build the summary data frame just as in summarise.abi.folder
  #fwd.summaries = lapply(fwd.dat, function(x) x[["summary"]])
  #rev.summaries = lapply(rev.dat, function(x) x[["summary"]])
  #all.summaries = c(fwd.summaries, rev.summaries)
  #all.summaries = do.call(rbind, all.summaries)
  #abi.fnames = unlist(c(fwd.fnames, rev.fnames))
  #folder.names = basename(dirname(abi.fnames))
  #file.names = basename(abi.fnames)
  #all.summaries = cbind.data.frame("file.path" = as.character(abi.fnames), "folder.name" = as.character(folder.names), "file.name" = file.names, all.summaries, stringsAsFactors = FALSE)
  
  
  used.reads = names(readset)
  #all.summaries$read.included.in.readset = all.summaries$file.path %in% used.reads
  
  
  return(list("readset" = readset)) #, "read.summaries" = all.summaries)
  
}

loadread.dula <- function(fname, trim, trim.cutoff, revcomp, max.secondary.peaks, secondary.peak.ratio, min.length, processors){
  
  read.abi = sangerseqR::read.abif(fname)
  
  s = summarise.abi.file(read.abi, trim.cutoff, secondary.peak.ratio, processors = processors)
  
  summary = s$summary
  
  # here we store the secondary peaks by storing them in a single sequence
  # as ambiguity codes. Note, this is version of a read that we use.
  # So in this package, a read has an ambiguit whereever there's a 
  # secondary peak
  d = c(DNAStringSet(s$read@primarySeq), DNAStringSet(s$read@secondarySeq))
  read = ConsensusSequence(d)[[1]]
  
  if(trim == TRUE){
    trim.start = summary["trim.start"]
    trim.finish = summary["trim.finish"]
    sp = summary["trimmed.secondary.peaks"]
    
  }else if(trim == FALSE){
    trim.start = 1
    trim.finish = length(read)
    sp = summary["raw.secondary.peaks"]
  }
  
  # FILTER read based on user specified limits
  read = read[trim.start:trim.finish]
  
  if(!is.null(max.secondary.peaks)){
    if(sp > max.secondary.peaks){
      read = NULL
    }
  }
  
  if(length(read) < min.length){
    read = NULL
  }    
  
  if(!is.null(read)) {
    if(revcomp == TRUE){
      read = reverseComplement(read)
    }
  }
  return(list('read' = read, summary = summary))
  
}


my_IUPAC_CODE_MAP <- c(
  A="A",
  C="C",
  G="G",
  T="T",
  M="AC",
  R="AG",
  W="AT",
  S="CG",
  Y="CT",
  K="GT",
  V="ACG",
  H="ACT",
  D="AGT",
  B="CGT",
  N="ACGT",
  A="AM",
  C="CM",
  A="AR",
  G="GR",
  A="AW",
  T="TW",
  C="CS",
  G="GS",
  C="CY",
  T="TY",
  G="GK",
  T="TK",
  T="KT",
  T="ST",
  #
  R="RV",
  A="AH",
  A="AY",
  Y="YB",
  G="GB",
  S="SB",
  A="AR",
  C="CR",
  G="RS",
  G="KR",
  G="GM",
  T="AKT",
  G="DG",
  C="BC"
)

my_IUPAC_CODE_MAP_2 <- c(
  A="A",
  C="C",
  G="G",
  T="T",
  M="AC",
  R="AG",
  W="AT",
  S="CG",
  Y="CT",
  K="GT",
  V="ACG",
  H="ACT",
  D="AGT",
  B="CGT",
  N="ACGT",
  M="AM",
  M="CM",
  R="AR",
  R="GR",
  W="AW",
  W="TW",
  S="CS",
  S="GS",
  S="SG",
  Y="CY",
  Y="YC",
  Y="TY",
  K="GK",
  K="TK",
  K="KT",
  S="ST",
  #
  R="RV",
  H="AH",
  Y="AY",
  Y="YB",
  B="GB",
  S="SB",
  R="AR",
  R="CR",
  G="RS",
  G="KR",
  M="GM",
  T="AKT",
  D="DG",
  B="BC"
)


mostConsensusString <- function(seqs){
  conm <- consensusMatrix(seqs)/length(seqs)
  conm <- conm[-c(16,17,18),]
  #conm <- conm[,-(which(colSums(conm)==0))]
  row_maxes <- apply(conm, 2, max)
  new_seq <- rep('+', ncol(conm))
  for (i in 1:ncol(conm)){
    lets <- sort(row.names(conm)[conm[,i] == row_maxes[i]])
    lets <- paste(lets, sep="", collapse="")
    if (nchar(lets) >= 2) {lets <- gsub("[M,R,W,S,Y,K,V,H,D,B,N]", "", lets)}
    if (lets == "") {lets="N"}
    if (lets == "M" || lets == "R" || lets == "W" || lets == "S" || lets == "Y" || lets == "K" || lets == "V" || lets == "H" || lets == "D" || lets == "B" || lets == "N") {
      new_seq[i] <- lets
    } else {
      amb_char <- names(my_IUPAC_CODE_MAP)[my_IUPAC_CODE_MAP == lets]
      new_seq[i] <- amb_char}
    
  }
  return(c(list("contig"=DNAString(paste(new_seq, sep="", collapse=""))),list("consensus_matrix"= conm)))
}

mostConsensusString_2 <- function(seqs){
  conm <- consensusMatrix(seqs)/length(seqs)
  conm <- conm[-c(16,17,18),]
  #conm <- conm[,-(which(colSums(conm)==0))]
  row_maxes <- apply(conm, 2, max)
  new_seq <- rep('+', ncol(conm))
  for (i in 1:ncol(conm)){
    lets <- sort(row.names(conm)[conm[,i] == row_maxes[i]])
    lets <- paste(lets, sep="", collapse="")
    if (nchar(lets) > 2) {lets <- gsub("[M,R,W,S,Y,K,V,H,D,B,N]", "", lets)}
    if (lets == "") {lets="N"}
    if (lets == "M" || lets == "R" || lets == "W" || lets == "S" || lets == "Y" || lets == "K" || lets == "V" || lets == "H" || lets == "D" || lets == "B" || lets == "N") {
      new_seq[i] <- lets
    } else {
      amb_char <- names(my_IUPAC_CODE_MAP_2)[my_IUPAC_CODE_MAP_2 == lets]
      new_seq[i] <- amb_char}
    
  }
  return(c(list("contig"=DNAString(paste(new_seq, sep="", collapse=""))),list("consensus_matrix"= conm)))
}


server = function(input, output, session) {
  
  all_data=reactive({
    
    options("digits"=2)
    # Check
    if(is.null(input$reference)) return(NULL)
    
    # Upload reference sequence
    reference_sequence=toupper(as.character(seqinr::read.fasta(input$reference$datapath,seqonly = TRUE)))
    names(reference_sequence)='Reference'
    
    # Select reads
    F_files=input$forward_sequence$datapath
    F_names=input$forward_sequence$name
    R_files=input$reverse_sequence$datapath
    R_names=input$reverse_sequence$name
    q_threshold=input$q_ts
    used_F_files=NULL
    used_F_names=NULL
    unused_F_files=NULL
    unused_F_names=NULL
    means_F=NULL
    means_unused_F=NULL
    for (i in 1:length(F_files)){
      if (mean(sangerseqR::read.abif(F_files[i])@data[["PCON.2"]])>q_threshold) {used_F_files=c(used_F_files, F_files[i])
      used_F_names=c(used_F_names, F_names[i])
      means_F=c(means_F, mean(sangerseqR::read.abif(F_files[i])@data[["PCON.2"]]))} else {unused_F_files=c(unused_F_files, F_files[i])
      unused_F_names=c(unused_F_names, F_names[i])
      means_unused_F=c(means_unused_F, mean(sangerseqR::read.abif(F_files[i])@data[["PCON.2"]]))}
      
    }
    used_R_files=NULL
    used_R_names=NULL
    unused_R_files=NULL
    unused_R_names=NULL
    means_R=NULL
    means_unused_R=NULL
    for (i in 1:length(R_files)){
      if (mean(sangerseqR::read.abif(R_files[i])@data[["PCON.2"]])>q_threshold) {used_R_files=c(used_R_files, R_files[i])
      used_R_names=c(used_R_names, R_names[i])
      means_R=c(means_R, mean(sangerseqR::read.abif(R_files[i])@data[["PCON.2"]]))} else {unused_R_files=c(unused_R_files, R_files[i])
      unused_R_names=c(unused_R_names, R_names[i])
      means_unused_R=c(means_unused_R, mean(sangerseqR::read.abif(R_files[i])@data[["PCON.2"]]))}
    }
    good_files=c(used_F_names,used_R_names)
    good_means=c(means_F,means_R)
    good_df=cbind(good_files, good_means)
    if (is.null(good_df)==FALSE) {
    colnames(good_df)=c("Used file","Average Phred score")} 
    
    bad_files=c(unused_F_names,unused_R_names)
    bad_means=c(means_unused_F,means_unused_R)
    bad_df=cbind(bad_files, bad_means)
    if (is.null(bad_df)==FALSE) {
    colnames(bad_df)=c("Unused file","Average score")}
    tc_3=15
    pc=0.33
    if ((sum(good_means)/length(good_means))<26) {tc_3=13
    pc=0.4}
    
    # Trimming

    trimmed_readset=ends_trim(Forward_files = used_F_files, Reverse_files = used_R_files, trim_cutoff5 = 20, trim_cutoff3 = tc_3, peak_cutoff = pc)
    
    read_set=c(trimmed_readset$readset,DNAStringSet(reference_sequence))
    aln_read_set=AlignSeqs(read_set, verbose = FALSE, gapOpening = -(input$gap_p))
    contig_c1=mostConsensusString(aln_read_set[-length(aln_read_set)])
    contig_c2=mostConsensusString_2(aln_read_set[-length(aln_read_set)])
    consensus_matrix=contig_c1$consensus_matrix
    ambigous_positions=NULL
    for (i in 1:ncol(consensus_matrix)) {
      sum_dna=sum(consensus_matrix[c(1:4),i])
      sum_iupac=sum(consensus_matrix[c(5:10),i])
      if (sum_dna==sum_iupac & (sum_dna+sum_iupac)>0) {ambigous_positions=c(ambigous_positions, i)}
    }
    
    
    not_N_c1= str_locate_all(contig_c1$contig, pattern = "[^N]")[[1]][,1]
    not_N_c2= str_locate_all(contig_c2$contig, pattern = "[^N]")[[1]][,1]
    contig_c1=contig_c1$contig[min(not_N_c1):max(not_N_c1)]
    contig_c2=contig_c2$contig[min(not_N_c2):max(not_N_c2)]
    contig_c1=DNAStringSet(contig_c1)
    contig_c2=DNAStringSet(contig_c2)
    names(contig_c1)="Contig"
    names(contig_c2)="Contig"
    trimed_contig_c1=trim.to.reference(contig_c1, DNAString(reference_sequence))
    trimed_contig_c2=trim.to.reference(contig_c2, DNAString(reference_sequence))
    aln_read_set_c1=AlignSeqs(c(read_set, contig_c1), verbose = FALSE, gapOpening = -(input$gap_p))
    aln_read_set_c2=AlignSeqs(c(read_set, contig_c2), verbose = FALSE, gapOpening = -(input$gap_p))
    
    F_R_names=c(used_F_names, used_R_names)
    
    for (i in 1:length(F_R_names)) {names(aln_read_set_c1)[i]=F_R_names[i]}
    
    return(list("aln_c1"=aln_read_set_c1, "aln_c2"=aln_read_set_c2, "cont_c1"=contig_c1, "cont_c2"=contig_c2,"trimed_c1"=trimed_contig_c1, "trimed_c2"=trimed_contig_c2, "ambig_pos"=ambigous_positions, 
           "good_df"=good_df, "bad_df"=bad_df))
    
    })
  
  observeEvent(input$reset, {
    shinyjs::reset("forward_sequence")
    shinyjs::reset("reverse_sequence")
    shinyjs::reset("reference")
  })
  
  
  
    
    
  output$alignment= renderUI({
      if (input$submit>0) {isolate(all_data())} else {return(h4("Server is ready to process the reads"))}
    
      column(width = 12,
             includeHTML(my_BrowseSeqs(all_data()$aln_c1,
                                       colors = c('#57e127','blue','#434749','red','#dfc5db'),patterns =c('A','C','G','T','-'),
                                       openURL = FALSE)),
             style = paste("height:100%; ",
                           "overflow-y: hidden;",
                           "overflow-x: scroll;"))
      })
  

   
  
    output$amb_pos= renderText({
      if (input$submit>0) {isolate(all_data())
        shinyjs::show("contig")
        as.character(all_data()$ambig_pos)
        
        }
      
      })
    
    observe({
      shinyjs::hide("contig")})
  
    
    output$contig <- downloadHandler(
      filename = function() {
        paste("contig", ".fasta", sep = "")
      },
      content = function(file) {
        write.fasta(sequences = list(all_data()$trimed_c1, all_data()$trimed_c2), file.out = file, names = c("contig ", "contig (alternative form)"))
      }
    )
  
      
    
    output$summary_used= renderTable({
      if (input$submit>0) {isolate(all_data())}
      all_data()$good_df
    })
    output$summary_unused= renderTable({
      if (input$submit>0) {isolate(all_data())}
      all_data()$bad_df
    })
    
    
      
      
    
   
     
    
  
}
shinyApp(ui, server)
