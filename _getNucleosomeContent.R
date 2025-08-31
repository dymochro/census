## determine nucleosome content of organelles of interest

getNucleosomeContent <- function(path, info) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  
  ## genome size
  # http://genomewiki.ucsc.edu/index.php/Mm10_Genome_size_statistics
  genome_size <- 2730871774 # bp
  
  
  ## nucleosome repeat length
  # Teif et al, 2012, Nat Struct Mol Biol 19: 1185-1192.
  # PMID 23085715, https://doi.org/10.1038/nsmb.2419, Fig. 7b
  nrl <- 186.1 # bp
  nrl_error <- 0.4
  
  
  ## total number of nucleosomes (diploid genome)
  number_nucleosomes <- 2*genome_size/nrl
  number_nucleosomes_error <- nrl_error*2*genome_size/nrl^2
  
  ## relative nucleosome content in nucleoli
  # Ballmer et al, 2022, Nucleic Acids Res 51: 117–143.
  # PMID 36533441, https://doi.org/10.1093/nar/gkac1159, Fig. 4
  # quantify in Cbx1 channel, segment in DAPI channel
  
  # load image data
  img1 <- flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341600&s=157&r=7&c=1")),90))
  img1 <- abind(img1, flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341600&s=157&r=7&c=2")),90)), along=1)
  
  img_bg1 <- img1[47:58,15:149,1]
  img_fbl1 <- img1[326:460,15:149,1]
  img_dapi1 <- img1[47:181,15:149,1]-mean(img_bg1)
  img_dapi1[img_dapi1<0] <- 0
  dapi_enrichment_nucleolus1 <- getEnrichment(img_fbl1, img_dapi1, 0.02, 0.2, erode_disc=13, dilate_disc_foci=21, img4nuc=TRUE, outlier_factor=0, min_size_nucleus=250)
  
  img2 <- flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341609&s=157&r=1&c=1")),90))
  img2 <- abind(img2, flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341609&s=157&r=1&c=2")),90)), along=1)
  
  img_bg2 <- img2[53:63,51:176,1]
  img_fbl2 <- img2[310:434,51:176,1]
  img_dapi2 <- img2[53:177,51:176,1]-mean(img_bg2)
  img_dapi2[img_dapi2<0] <- 0
  dapi_enrichment_nucleolus2 <- getEnrichment(img_fbl2, img_dapi2, 0.02, 0.2, erode_disc=13, dilate_disc_foci=15, img4nuc=TRUE, outlier_factor=0, min_size_nucleus=250)
  
  nucleosome_enrichment_nucleolus <- mean(c(dapi_enrichment_nucleolus1[[1]],dapi_enrichment_nucleolus2[[1]]))
  nucleosome_enrichment_nucleolus_error <- sd(c(dapi_enrichment_nucleolus1[[1]],dapi_enrichment_nucleolus2[[1]]))
  
  ## relative nucleosome content in heterochromatin foci
  # Yu, Sun, Tan, Pan et al, 2021, Nat Commun 12: 6365.
  # PMID 34753899, https://doi.org/10.1038/s41467-021-26576-2, Fig. 3c (DAPI)
  target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-021-26576-2/MediaObjects/41467_2021_26576_Fig3_HTML.png"
  img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target, .opts = list(ssl.verifypeer = FALSE))),90))
  img_dapi <- img[1171:1394,50:267,3]
  dapi_enrichment_heterochromatin <- getEnrichment(img_dapi, img_dapi, 1, 1.5, outlier_factor=1.2)
  
  # Muller-Ott et al, 2014, Mol Syst Biol 10: 746.
  # PMID 25134515, https://doi.org/10.15252/msb.20145377, Fig. 1
  dapi_to_nucleosomes_het <- 2.7/1.8 # from the comparison between H2A-RFP versus DAPI enrichment in heterochromatin foci
  dapi_to_nucleosomes_het_error <- sqrt((0.1/1.8)^2 + (2.7/1.8^2*0.3)^2)
    
  nucleosome_enrichment_heterochromatin <- dapi_enrichment_heterochromatin[[1]]/dapi_to_nucleosomes_het
  nucleosome_enrichment_heterochromatin_error <- sqrt((dapi_enrichment_heterochromatin[[2]]/dapi_to_nucleosomes_het)^2 + (dapi_enrichment_heterochromatin[[1]]/dapi_to_nucleosomes_het^2*dapi_to_nucleosomes_het_error)^2)
  
  
  
  ## relative nucleosome content in transcriptional condensates
  # Sabari et al, 2018, Science 361: eaar3958.
  # PMID 29930091, https://doi.org/10.1126/science.aar3958, Fig. 1A (Med1)
  target <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6092193/bin/nihms-983922-f0001.jpg"
  target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/963d/6092193/89fc5bd52e09/nihms-983922-f0001.jpg"
  img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
  img_dapi <- img[140:230, 158:257, 3] # select blue channel
  img_dapi[65:91,85:100] <- 0
  img_med1 <- img[140:230, 158:257, 2] # select green channel
  img_med1[65:91,85:100] <- 0
  dapi_enrichment_txn <- getEnrichment(img_med1, img_dapi, 0.8, 0.18, outlier_factor=0.1, img4nuc=TRUE)

  nucleosome_enrichment_txn <- dapi_enrichment_txn[[1]]
  nucleosome_enrichment_txn_error <- dapi_enrichment_txn[[2]]
  
  
  
  ## relative nucleosome content in Polycomb bodies
  # Tardat, Albert et al, 2015, Mol Cell 58: 157-171.
  # PMID 25801166, https://doi.org/10.1016/j.molcel.2015.02.013, Fig. 1C (Cbx2)
  target <- "https://ars.els-cdn.com/content/image/1-s2.0-S1097276515001276-gr1_lrg.jpg"
  img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
  img_dapi <- img [61:277, 1366:1581, 2] # select green channel (source image is grayscale)
  img_cbx2 <- img [61:277, 1598:1813, 2] # select green channel (source image is grayscale)
  dapi_enrichment_polycomb <- getEnrichment(img_cbx2, img_dapi, 1, 1.5, outlier_factor=1.5, img4nuc=TRUE)
  
  nucleosome_enrichment_polycomb <- dapi_enrichment_polycomb[[1]]
  nucleosome_enrichment_polycomb_error <- dapi_enrichment_polycomb[[2]]
  
  
  
  ## volume of the nucleoplasm (nucleus without nucleoli and heterochromatin foci) 
  vnp <- info$TotalVolume[which(info$Organelle=="Nucleoplasm")]-info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]
  vnp_error <- sqrt(info$TotalVolume_Err[which(info$Organelle=="Nucleoplasm")]^2+info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]^2)
    
  ## weighted volume average
  wv <- vnp + nucleosome_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")] + nucleosome_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]
  wv_error <- sqrt(vnp_error^2 + (nucleosome_enrichment_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (nucleosome_enrichment_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus")])^2 +(nucleosome_enrichment_heterochromatin_error*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")])^2 + (nucleosome_enrichment_heterochromatin*info$TotalVolume_Err[which(info$Organelle=="Heterochromatin focus")])^2)
  
  ## number of nucleosomes in heterochromatin
  number_nucleosomes_heterochromatin <- nucleosome_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv
  number_nucleosomes_heterochromatin_error <- sqrt((nucleosome_enrichment_heterochromatin_error*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_heterochromatin*info$Volume_Err[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv^2)^2)
  total_number_nucleosomes_heterochromatin <- nucleosome_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv
  total_number_nucleosomes_heterochromatin_error <- sqrt((nucleosome_enrichment_heterochromatin_error*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_heterochromatin*info$TotalVolume_Err[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_nucleosomes/wv^2)^2)
  
  ## number of nucleosomes in nucleoli
  number_nucleosomes_nucleolus <- nucleosome_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv
  number_nucleosomes_nucleolus_error <- sqrt((nucleosome_enrichment_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv^2)^2)
  total_number_nucleosomes_nucleolus <- nucleosome_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv
  total_number_nucleosomes_nucleolus_error <- sqrt((nucleosome_enrichment_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_nucleosomes/wv^2)^2)
  
  ## number of nucleosomes in GC/DFC
  number_nucleosomes_nucleolus_gc <- number_nucleosomes_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  number_nucleosomes_nucleolus_gc_error <- sqrt((number_nucleosomes_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_nucleosomes_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, GC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  number_nucleosomes_nucleolus_dfc <- number_nucleosomes_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  number_nucleosomes_nucleolus_dfc_error <- sqrt((number_nucleosomes_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_nucleosomes_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  total_number_nucleosomes_nucleolus_gc <- number_nucleosomes_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  total_number_nucleosomes_nucleolus_gc_error <- sqrt((number_nucleosomes_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_nucleosomes_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  total_number_nucleosomes_nucleolus_dfc <- number_nucleosomes_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  total_number_nucleosomes_nucleolus_dfc_error <- sqrt((number_nucleosomes_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_nucleosomes_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_nucleosomes_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  
  ## number of nucleosomes in small transcriptional condensates
  number_nucleosomes_txn_small <- nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv
  number_nucleosomes_txn_small_error <- sqrt((nucleosome_enrichment_txn_error*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$Volume_Err[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv^2)^2)
  total_number_nucleosomes_txn_small <- nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv
  total_number_nucleosomes_txn_small_error <- sqrt((nucleosome_enrichment_txn_error*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$TotalVolume_Err[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_nucleosomes/wv^2)^2)
  
  ## number of nucleosomes in large transcriptional condensates
  number_nucleosomes_txn_large <- nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv
  number_nucleosomes_txn_large_error <- sqrt((nucleosome_enrichment_txn_error*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$Volume_Err[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv^2)^2)
  total_number_nucleosomes_txn_large <- nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv
  total_number_nucleosomes_txn_large_error <- sqrt((nucleosome_enrichment_txn_error*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$TotalVolume_Err[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_nucleosomes/wv^2)^2)
  
  ## number of nucleosomes in Polycomb bodies
  number_nucleosomes_polycomb <- nucleosome_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv
  number_nucleosomes_polycomb_error <- sqrt((nucleosome_enrichment_polycomb_error*info$Volume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_polycomb*info$Volume_Err[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv^2)^2)
  total_number_nucleosomes_polycomb <- nucleosome_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv
  total_number_nucleosomes_polycomb_error <- sqrt((nucleosome_enrichment_polycomb_error*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_polycomb*info$TotalVolume_Err[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv)^2 + (nucleosome_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_nucleosomes_error/wv)^2 + (wv_error*nucleosome_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_nucleosomes/wv^2)^2)
  
  ## assemble and return results
  result <- data.frame(
    Organelle = c("Nucleolus", "Nucleolus, GC", "Nucleolus, DFC", "Txn condensate, small", "Txn condensate, large", "Heterochromatin focus", "Polycomb body", "Nucleus"),
    Nucleosome_enrichment = c(nucleosome_enrichment_nucleolus, nucleosome_enrichment_nucleolus, nucleosome_enrichment_nucleolus, nucleosome_enrichment_txn, nucleosome_enrichment_txn, nucleosome_enrichment_heterochromatin, nucleosome_enrichment_polycomb, 1),
    Nucleosome_enrichment_Err = c(nucleosome_enrichment_nucleolus_error, nucleosome_enrichment_nucleolus_error, nucleosome_enrichment_nucleolus_error, nucleosome_enrichment_txn_error, nucleosome_enrichment_txn_error, nucleosome_enrichment_heterochromatin_error, nucleosome_enrichment_polycomb_error, 0),
    Nucleosome_number = c(number_nucleosomes_nucleolus, number_nucleosomes_nucleolus_gc, number_nucleosomes_nucleolus_dfc, number_nucleosomes_txn_small, number_nucleosomes_txn_large, number_nucleosomes_heterochromatin, number_nucleosomes_polycomb, number_nucleosomes),
    Nucleosome_number_Err = c(number_nucleosomes_nucleolus_error, number_nucleosomes_nucleolus_gc_error, number_nucleosomes_nucleolus_dfc_error, number_nucleosomes_txn_small_error, number_nucleosomes_txn_large_error, number_nucleosomes_heterochromatin_error, number_nucleosomes_polycomb_error, number_nucleosomes_error),
    TotalNucleosome_number = c(total_number_nucleosomes_nucleolus, total_number_nucleosomes_nucleolus_gc, total_number_nucleosomes_nucleolus_dfc, total_number_nucleosomes_txn_small, total_number_nucleosomes_txn_large, total_number_nucleosomes_heterochromatin, total_number_nucleosomes_polycomb, number_nucleosomes),
    TotalNucleosome_number_Err = c(total_number_nucleosomes_nucleolus_error, total_number_nucleosomes_nucleolus_gc_error, total_number_nucleosomes_nucleolus_dfc_error, total_number_nucleosomes_txn_small_error, total_number_nucleosomes_txn_large_error, total_number_nucleosomes_heterochromatin_error, total_number_nucleosomes_polycomb_error, number_nucleosomes_error),
    stringsAsFactors = FALSE)

  return(result)
}
