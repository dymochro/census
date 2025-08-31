## determine RNA content of organelles of interest

getRNAContent <- function(path, info) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  
  ## total amount of cellular RNA
  # Calabrese et al, 2007, Proc Natl Acad Sci USA 104: 18097-18102.
  # PMID 17989215, https://doi.org/10.1073/pnas.0709193104, Methods section (first paragraph)
  total_rna <- 20 # pg
  
  ## fraction of nuclear RNA
  # Piwnicka et al, 1983, Cytometry 3: 269-275.
  # PMID 6185286, https://doi.org/10.1002/cyto.990030407, Abstract (in different cell types)
  nuclear_rna_fraction <- 0.15
  nuclear_rna_fraction_error <- 0.05
  
  ## total amount of nuclear RNA
  nuclear_rna <- total_rna*nuclear_rna_fraction # pg
  nuclear_rna_error <- total_rna*nuclear_rna_fraction_error # pg
  
  ## total number of RNAs (in equivalents of ribonucleotides)
  number_rnas <- nuclear_rna/340*6.022*10^11
  number_rnas_error <- nuclear_rna_error/340*6.022*10^11
  
  ## average length of RNA
  mrna_length <- 2790 # ribonucleotides
  rrna_length <- 14000 # ribonucleotides
  
  
  ## relative RNA content in heterochromatin foci and nucleoli
  # Beagrie, Thieme, Annunziatella et al, 2023, Nat Methods 20: 1037-1047.
  # PMID 37336949, https://doi.org/10.1038/s41592-023-01903-1, Extended Data Fig. 4e
  target <- "https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41592-023-01903-1/MediaObjects/41592_2023_1903_Fig9_ESM.jpg"
  img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target, .opts = list(ssl.verifypeer = FALSE))),90))
  img_crop <- img[1332:1998,1009:1675,]
  img_crop[abs(img_crop[,,2]-img_crop[,,1])<0.3 & (img_crop[,,1]>0.4*max(img_crop[,,1]))] <- 0
  img_sytorna <- img_crop[,,2]
  img_dapi <- img_crop[,,3]-img_crop[,,1]
  img_dapi <- img_dapi - min(img_dapi)
  sytorna_enrichment_heterochromatin <- getEnrichment(img_dapi, img_sytorna, 0.85, 1.3, outlier_factor=1, min_size_condensate=5, min_size_nucleus=600)
  
  rna_enrichment_heterochromatin <- sytorna_enrichment_heterochromatin[[1]]
  rna_enrichment_heterochromatin_error <- sytorna_enrichment_heterochromatin[[2]]
  
  
  img_nucleolus <- img_sytorna
  img_nucleolus[EBImage::erode(sytorna_enrichment_heterochromatin[[3]][[1]],EBImage::makeBrush(21, shape='disc'))==0] <- 0
  sytorna_enrichment_nucleolus <- getEnrichment(img_nucleolus, img_sytorna, 0.85, 1.6, outlier_factor=1.5, min_size_condensate=5, min_size_nucleus=600)
  
  rna_enrichment_nucleolus <- sytorna_enrichment_nucleolus[[1]]
  rna_enrichment_nucleolus_error <- sytorna_enrichment_nucleolus[[2]]
  
  
  ## relative RNA content in transcriptional condensates
  # generate a mask for the nucleoplasm without nucleoli
  img_sytorna_mask <- img_sytorna
  img_sytorna_mask[sytorna_enrichment_heterochromatin[[3]][[1]]<1] <- 0
  nuclei <- segmentNuclei(img_nucleolus, 0.85*EBImage::otsu(img_nucleolus), 3, 3, 600)
  for(i in 1:nuclei[[2]]) { # remove nucleoli
    nucleoli <- segmentCondensates(ref=img_nucleolus, img=img_sytorna, mask_nucleus=nuclei[[1]], noi=as.numeric(names(table(nuclei[[1]]))[i+1]), seg_factor=1.6, dilate_disc=3, erode_disc=3, min_size=5, outlier_factor=1.5)
    img_sytorna_mask[EBImage::dilate(nucleoli[[2]],EBImage::makeBrush(11, shape='disc'))>0] <- 0
  }
  # quantify maximum values in the nucleoplasm outside of nucleoli
  sytorna_enrichment_max_in_nucleoplasm <- getEnrichment(img_sytorna_mask, img_sytorna, 0.85, 1.6, outlier_factor=1, min_size_condensate=2, min_size_nucleus=600)
  rna_enrichment_txn <- sytorna_enrichment_max_in_nucleoplasm[[1]]
  rna_enrichment_txn_error <- sytorna_enrichment_max_in_nucleoplasm[[2]]
  
  
  
  ## relative RNA content in Polycomb bodies
  rna_enrichment_polycomb <- 1
  rna_enrichment_polycomb_error <- 0
  
  
  
  ## volume of the nucleoplasm (nucleus without nucleoli and heterochromatin foci) 
  vnp <- info$TotalVolume[which(info$Organelle=="Nucleoplasm")]-info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]
  vnp_error <- sqrt(info$TotalVolume_Err[which(info$Organelle=="Nucleoplasm")]^2+info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]^2)
    
  ## weighted volume average
  wv <- vnp + rna_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")] + rna_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]
  wv_error <- sqrt(vnp_error^2 + (rna_enrichment_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (rna_enrichment_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus")])^2 +(rna_enrichment_heterochromatin_error*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")])^2 + (rna_enrichment_heterochromatin*info$TotalVolume_Err[which(info$Organelle=="Heterochromatin focus")])^2)
  
  ## number of RNAs in heterochromatin
  number_rnas_heterochromatin <- rna_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv
  number_rnas_heterochromatin_error <- sqrt((rna_enrichment_heterochromatin_error*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv)^2 + (rna_enrichment_heterochromatin*info$Volume_Err[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv)^2 + (rna_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_heterochromatin*info$Volume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv^2)^2)
  total_number_rnas_heterochromatin <- rna_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv
  total_number_rnas_heterochromatin_error <- sqrt((rna_enrichment_heterochromatin_error*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv)^2 + (rna_enrichment_heterochromatin*info$TotalVolume_Err[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv)^2 + (rna_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_heterochromatin*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]*number_rnas/wv^2)^2)
  
  ## number of RNAs in nucleoli
  number_rnas_nucleolus <- rna_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_rnas/wv
  number_rnas_nucleolus_error <- sqrt((rna_enrichment_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus")]*number_rnas/wv)^2 + (rna_enrichment_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus")]*number_rnas/wv)^2 + (rna_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_nucleolus*info$Volume[which(info$Organelle=="Nucleolus")]*number_rnas/wv^2)^2)
  total_number_rnas_nucleolus <- rna_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_rnas/wv
  total_number_rnas_nucleolus_error <- sqrt((rna_enrichment_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_rnas/wv)^2 + (rna_enrichment_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_rnas/wv)^2 + (rna_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus")]*number_rnas/wv^2)^2)
  
  ## number of RNAs in GC/DFC
  number_rnas_nucleolus_gc <- number_rnas_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  number_rnas_nucleolus_gc_error <- sqrt((number_rnas_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_rnas_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_rnas_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, GC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  number_rnas_nucleolus_dfc <- number_rnas_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  number_rnas_nucleolus_dfc_error <- sqrt((number_rnas_nucleolus_error*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_rnas_nucleolus*info$Volume_Err[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_rnas_nucleolus*info$Volume[which(info$Organelle=="Nucleolus, DFC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  total_number_rnas_nucleolus_gc <- number_rnas_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  total_number_rnas_nucleolus_gc_error <- sqrt((number_rnas_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_rnas_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus, GC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_rnas_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  total_number_rnas_nucleolus_dfc <- number_rnas_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")]
  total_number_rnas_nucleolus_dfc_error <- sqrt((number_rnas_nucleolus_error*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (number_rnas_nucleolus*info$TotalVolume_Err[which(info$Organelle=="Nucleolus, DFC")]/info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 + (info$TotalVolume_Err[which(info$Organelle=="Nucleolus")]*number_rnas_nucleolus*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]/(info$TotalVolume[which(info$Organelle=="Nucleolus")])^2 )^2 )
  
  ## number of RNAs in small transcriptional condensates
  number_rnas_txn_small <- rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv
  number_rnas_txn_small_error <- sqrt((rna_enrichment_txn_error*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$Volume_Err[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv^2)^2)
  total_number_rnas_txn_small <- rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv
  total_number_rnas_txn_small_error <- sqrt((rna_enrichment_txn_error*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$TotalVolume_Err[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, small")]*number_rnas/wv^2)^2)
  
  ## number of RNAs in large transcriptional condensates
  number_rnas_txn_large <- rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv
  number_rnas_txn_large_error <- sqrt((rna_enrichment_txn_error*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$Volume_Err[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_txn*info$Volume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv^2)^2)
  total_number_rnas_txn_large <- rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv
  total_number_rnas_txn_large_error <- sqrt((rna_enrichment_txn_error*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$TotalVolume_Err[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv)^2 + (rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_txn*info$TotalVolume[which(info$Organelle=="Txn condensate, large")]*number_rnas/wv^2)^2)
  
  ## number of RNAs in Polycomb bodies
  number_rnas_polycomb <- rna_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_rnas/wv
  number_rnas_polycomb_error <- sqrt((rna_enrichment_polycomb_error*info$Volume[which(info$Organelle=="Polycomb body")]*number_rnas/wv)^2 + (rna_enrichment_polycomb*info$Volume_Err[which(info$Organelle=="Polycomb body")]*number_rnas/wv)^2 + (rna_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_polycomb*info$Volume[which(info$Organelle=="Polycomb body")]*number_rnas/wv^2)^2)
  total_number_rnas_polycomb <- rna_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_rnas/wv
  total_number_rnas_polycomb_error <- sqrt((rna_enrichment_polycomb_error*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_rnas/wv)^2 + (rna_enrichment_polycomb*info$TotalVolume_Err[which(info$Organelle=="Polycomb body")]*number_rnas/wv)^2 + (rna_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_rnas_error/wv)^2 + (wv_error*rna_enrichment_polycomb*info$TotalVolume[which(info$Organelle=="Polycomb body")]*number_rnas/wv^2)^2)
  
  ## total number of RNA molecules in the nucleus, using a weighted average
  wnumber <- number_rnas_nucleolus/rrna_length + (number_rnas-number_rnas_nucleolus)/mrna_length
  wnumber_error <- sqrt((number_rnas_nucleolus_error/rrna_length)^2 + (number_rnas_error/mrna_length)^2 + (number_rnas_nucleolus_error/mrna_length)^2)
  
  
  ## assemble and return results
  result <- data.frame(
    Organelle = c("Nucleolus", "Nucleolus, GC", "Nucleolus, DFC", "Txn condensate, small", "Txn condensate, large", "Heterochromatin focus", "Polycomb body", "Nucleus"),
    RNuc_enrichment = c(rna_enrichment_nucleolus, rna_enrichment_nucleolus, rna_enrichment_nucleolus, rna_enrichment_txn, rna_enrichment_txn, rna_enrichment_heterochromatin, rna_enrichment_polycomb, 1),
    RNuc_enrichment_Err = c(rna_enrichment_nucleolus_error, rna_enrichment_nucleolus_error, rna_enrichment_nucleolus_error, rna_enrichment_txn_error, rna_enrichment_txn_error, rna_enrichment_heterochromatin_error, rna_enrichment_polycomb_error, 0),
    RNuc_number = c(number_rnas_nucleolus, number_rnas_nucleolus_gc, number_rnas_nucleolus_dfc, number_rnas_txn_small, number_rnas_txn_large, number_rnas_heterochromatin, number_rnas_polycomb, number_rnas),
    RNuc_number_Err = c(number_rnas_nucleolus_error, number_rnas_nucleolus_gc_error, number_rnas_nucleolus_dfc_error, number_rnas_txn_small_error, number_rnas_txn_large_error, number_rnas_heterochromatin_error, number_rnas_polycomb_error, number_rnas_error),
    TotalRNuc_number = c(total_number_rnas_nucleolus, total_number_rnas_nucleolus_gc, total_number_rnas_nucleolus_dfc, total_number_rnas_txn_small, total_number_rnas_txn_large, total_number_rnas_heterochromatin, total_number_rnas_polycomb, number_rnas),
    TotalRNuc_number_Err = c(total_number_rnas_nucleolus_error, total_number_rnas_nucleolus_gc_error, total_number_rnas_nucleolus_dfc_error, total_number_rnas_txn_small_error, total_number_rnas_txn_large_error, total_number_rnas_heterochromatin_error, total_number_rnas_polycomb_error, number_rnas_error),
    RNA_number = c(number_rnas_nucleolus/rrna_length, number_rnas_nucleolus_gc/rrna_length, number_rnas_nucleolus_dfc/rrna_length, number_rnas_txn_small/mrna_length, number_rnas_txn_large/mrna_length, number_rnas_heterochromatin/mrna_length, number_rnas_polycomb/mrna_length, wnumber),
    RNA_number_Err = c(number_rnas_nucleolus_error/rrna_length, number_rnas_nucleolus_gc_error/rrna_length, number_rnas_nucleolus_dfc_error/rrna_length, number_rnas_txn_small_error/mrna_length, number_rnas_txn_large_error/mrna_length, number_rnas_heterochromatin_error/mrna_length, number_rnas_polycomb_error/mrna_length, wnumber_error),
    TotalRNA_number = c(total_number_rnas_nucleolus/rrna_length, total_number_rnas_nucleolus_gc/rrna_length, total_number_rnas_nucleolus_dfc/rrna_length, total_number_rnas_txn_small/mrna_length, total_number_rnas_txn_large/mrna_length, total_number_rnas_heterochromatin/mrna_length, total_number_rnas_polycomb/mrna_length, wnumber),
    TotalRNA_number_Err = c(total_number_rnas_nucleolus_error/rrna_length, total_number_rnas_nucleolus_gc_error/rrna_length, total_number_rnas_nucleolus_dfc_error/rrna_length, total_number_rnas_txn_small_error/mrna_length, total_number_rnas_txn_large_error/mrna_length, total_number_rnas_heterochromatin_error/mrna_length, total_number_rnas_polycomb_error/mrna_length, wnumber_error),
    stringsAsFactors = FALSE)

  return(result)
}
