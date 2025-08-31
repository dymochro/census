# This function retrieves the copy numbers of proteins of interest in mESCs


getScaffoldCopyNumbers <- function(proteins=NA, setAtMin=TRUE, nucleosome_content) {
  # load iBAQ values from Atlasi et al, 2020, Nat Commun 11: 1617, PMID 32238817, https://doi.org/10.1038/s41467-020-15449-9.
  forward_proteome <- rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15449-9/MediaObjects/41467_2020_15449_MOESM7_ESM.xlsx", which="Forward_proteome")
  reverse_proteome <- rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15449-9/MediaObjects/41467_2020_15449_MOESM7_ESM.xlsx", which="Reverse_proteome")
  
  # get relative abundance based on iBAQ values
  fwd <- forward_proteome[,c("Gene...1", "iBAQ")]
  rev <- reverse_proteome[,c("Gene names...1", "iBAQ")]
  fwd$iBAQ <- fwd$iBAQ/sum(fwd$iBAQ)
  rev$iBAQ <- rev$iBAQ/sum(rev$iBAQ)
  
  # merge forward/reverse proteomes into one data frame
  colnames(fwd) <- c("Gene", "iBAQ, fwd")
  colnames(rev) <- c("Gene", "iBAQ, rev")
  combi <- merge(fwd, rev, all=T)
  
  # remove rows without gene name
  combi <- combi[!is.na(combi[,"Gene"]),]
  
  # fill missing values of either fwd or rev
  quantile_fwd <- min(combi[,c("iBAQ, fwd")],na.rm=TRUE)
  quantile_rev <- min(combi[,c("iBAQ, rev")],na.rm=TRUE)
  combi$`iBAQ, fwd` <- unname(sapply(combi[,c("iBAQ, fwd")],function(x) {if(is.na(x)){return(quantile_fwd)} else{return(x)}}))
  combi$`iBAQ, rev` <- unname(sapply(combi[,c("iBAQ, rev")],function(x) {if(is.na(x)){return(quantile_rev)} else{return(x)}}))

  # get coefficient to convert iBAQ values to copy numbers
  conversion <- getIBAQtoCopyNumberCoeff(nucleosome_content)
  slope <- conversion[[1]]
  slope_error <- conversion[[2]]
  
  # average forward/reverse values and convert to copy numbers
  combi$`ibaq, avg` <- apply(combi[,c("iBAQ, fwd", "iBAQ, rev")], 1, function(x) mean(x))
  combi$`ibaq, ratio_avg` = (combi$`ibaq, avg`/sum(combi$`ibaq, avg`))
  combi$`CopyNumber, avg` = combi$`ibaq, ratio_avg`*slope

  # retrieve copy numbers
  cn <- sapply(proteins, function(x) {
    # match gene names, accounting for the fact that some names in the table are separated by ";"
    indices <- sapply(combi[,"Gene"], function(y) sum(tolower(unlist(strsplit(y,";")))==tolower(x))>0)
      
    # sum copy numbers if there are multiple matches
    sum(combi[indices, "CopyNumber, avg"])
  })

  # calculate errors of copy numbers
  combi$`ibaq, ratio_avg_error` <- apply(combi[,c("iBAQ, fwd", "iBAQ, rev")], 1, function(x) sem(x))

  cn_error <- sapply(proteins, function(x) {
    # match gene names, accounting for the fact that some names in the table are separated by ";"
    indices <- sapply(combi[,"Gene"], function(y) sum(tolower(unlist(strsplit(y,";")))==tolower(x))>0)
      
    # calculate error for each entry
    err <- sqrt((combi[indices, "ibaq, ratio_avg"]*slope_error)^2 + (slope*combi[indices, "ibaq, ratio_avg_error"])^2)
      
    # calculate combined error
    return(sqrt(sum(err^2)))
  }) 

    
  if(setAtMin) {
    # set a minimum value for undetected proteins
    min_abundance <- round(min(combi$`CopyNumber, avg`))
    print(paste0("Undetected proteins set at minimum abundance (",min_abundance,"):"))
    print(proteins[which(cn==0)])

    # attribute minimum value to undetected proteins, assuming a relative error of 10%
    cn <- sapply(cn, function(x) {if(x<1) {return(min_abundance)} else {return(x)}})
    cn_error <- sapply(cn_error, function(x) {if(x<1) {return(0.1*min_abundance )} else {return(x)}})
  }

  # return copy numbers and their errors
  res <- cbind(unname(round(cn)), unname(round(cn_error)))
  return(res)
}


## ###################
## accessory functions
## ###################

# standard error of the mean
sem <- function(x) sd(x, na.rm=T)/sqrt(length(na.omit(x)))

# conversion between IBAQ and absolute copy number
getIBAQtoCopyNumberCoeff <- function(nucleosome_content) {
  # load iBAQ values from Atlasi et al, 2020, Nat Commun 11: 1617, PMID 32238817, https://doi.org/10.1038/s41467-020-15449-9.
  forward_proteome <- rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15449-9/MediaObjects/41467_2020_15449_MOESM7_ESM.xlsx", which="Forward_proteome")
  reverse_proteome <- rio::import("https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-15449-9/MediaObjects/41467_2020_15449_MOESM7_ESM.xlsx", which="Reverse_proteome")
  fwd <- forward_proteome[,c("Gene...1", "iBAQ")]
  rev <- reverse_proteome[,c("Gene names...1", "iBAQ")]
  
  # merge forward/reverse proteomes into one data frame
  colnames(fwd) <- c("Gene", "iBAQ, fwd")
  colnames(rev) <- c("Gene", "iBAQ, rev")
  combi <- merge(fwd, rev, all=T)
  
  # remove rows without gene name
  combi <- combi[!is.na(combi[,"Gene"]),]
  
  # average forward/reverse values and normalize
  combi$`iBAQ, avg` <- apply(combi[,c("iBAQ, fwd", "iBAQ, rev")], 1, function(x) mean(x, na.rm=T))
  combi$`iBAQ, ratio` <- combi$`iBAQ, avg`/sum(combi$`iBAQ, avg`)

  # use absolute copy numbers for CTCF, Rad21 amd H4 to find relationship between iBAQ and copy number
  number_h4 <- 2*nucleosome_content[which(nucleosome_content$Organelle=="Nucleus"),]$TotalNucleosome_number # each nucleosome contains two copies of H4
  
  ## Cattoglio et al, 2019, eLife 8: e40164, PMID 31205001, https://doi.org/10.7554/eLife.40164
  number_ctcf <- 217200 
	number_rad21 <- 109400
  
	## Huseyin & Klose, 2021, Nat Commun 12: 887, PMID 33563969, https://doi.org/10.1038/s41467-021-21130-6
	number_rnf2 <- 62423 # Rnf2 = Ring1B
	number_rybp <- 30000
	number_pcgf6 <- 14000
	
  # fit a linear function
	# y ~ a*x + b (copy number ~ a*ibaq/ibaq_tot + b)
	y <- c(number_ctcf, number_h4, number_rad21, number_rnf2, number_rybp, number_pcgf6) 
  x <- c(combi[combi[,"Gene"]=="Ctcf","iBAQ, ratio"], combi[grepl("Hist1h4a", combi[,"Gene"]),"iBAQ, ratio"], combi[combi[,"Gene"]=="Rad21","iBAQ, ratio"], combi[combi[,"Gene"]=="Rnf2","iBAQ, ratio"], combi[combi[,"Gene"]=="Rybp","iBAQ, ratio"], combi[combi[,"Gene"]=="Pcgf6","iBAQ, ratio"])

	fit_res <- nls(y ~ a*x, start=c(a=10^6), lower=c(a=10), upper=c(a=10^10), alg="port", control=list(maxiter=100, tol=1e-05, minFactor=1e-07, printEval=FALSE, warnOnly=TRUE))
  coefficients <- summary(fit_res)$coefficients
  
  # plot(x, y, log="xy", xlim=c(1e-5,1e-1), ylim=c(2e4,2e8))
  # lines(x, predict(fit_res), col="blue")
  return(coefficients)
}