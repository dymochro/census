## determine enrichment based on a microscopy image

getEnrichment <- function(ref, img, seg_factor_nucleus, seg_factor_condensates, dilate_disc=3, erode_disc=3, dilate_disc_foci=0, erode_disc_foci=0, img4nuc=FALSE, outlier_factor=1, min_size_nucleus=50, min_size_condensate=1) {
  # segment nucleus
  if(!img4nuc) {nuclei <- segmentNuclei(ref, threshold=EBImage::otsu(ref)*seg_factor_nucleus, dilate_disc, erode_disc, min_size=min_size_nucleus)} else
               {nuclei <- segmentNuclei(img, threshold=EBImage::otsu(ref)*seg_factor_nucleus, dilate_disc, erode_disc, min_size=min_size_nucleus)}
  enrichment <- rep(0, nuclei[[2]])
  display(nuclei[[1]])
  
  # quantify median enrichment in condensates (for each segmented nucleus)
  if(dilate_disc_foci==0) {dilate_disc_foci <- dilate_disc}
  if(erode_disc_foci==0) {erode_disc_foci <- erode_disc}
  quanti <- quantify(ref=ref, img=img, nuc=nuclei, seg_factor=seg_factor_condensates, dilate_disc_foci, erode_disc_foci, min_size=min_size_condensate, outlier_factor=outlier_factor)
  #EBImage::display(quanti[[3]][[1]])

  # return enrichment
  return(list(mean(quanti[[1]], na.rm=T), mean(quanti[[2]], na.rm=T), nuclei, quanti[[3]]))
}


## #########################################
## accessory functions used for segmentation
## #########################################

segmentNuclei <- function(image, threshold, dilate_disc, erode_disc, min_size) {
  # make threshold image
  mask_nucleus <- image > threshold
  
  # erode and dilate the segmented area
  mask_nucleus <- EBImage::dilate(mask_nucleus, EBImage::makeBrush(dilate_disc, shape='disc'))
  mask_nucleus <- EBImage::erode(mask_nucleus, EBImage::makeBrush(erode_disc, shape='disc'))
  mask_nucleus <- EBImage::fillHull(mask_nucleus) # fill holes in the mask
  mask_nucleus <- EBImage::bwlabel(mask_nucleus) # find connected sets of pixels
  
  # remove small regions
  regions_nuc <- table(mask_nucleus)
  small_regions <- which(regions_nuc<min_size)-1
  mask_nucleus[mask_nucleus %in% small_regions] <- 0
  
  # count nuclei
  regions_nuc <- table(mask_nucleus)
  nuclei_number <- length(regions_nuc)-1
  
  # return results
  return(list(mask_nucleus, nuclei_number))
}

segmentCondensates <- function(ref, img, mask_nucleus, noi, seg_factor, dilate_disc, erode_disc, min_size, outlier_factor) {
  # ref: Reference image, used to segment condensates
  # img: Image showing protein of interest
  # noi: nucleus of interest
  # seg_factor: scaling factor for segmentation threshold (for condensates)
  # outlier_factor: scaling factor for segmentation threshold (for nucleoplasm)
  
  # select nucleus of interest in reference image
  ref_noi <- ref
  ref_noi[mask_nucleus != noi] <- NA
  
  # select nucleus of interest in image with condensates (will be used to remove out-of-focus condensates from the nucleoplasm)
  img_noi <- img
  img_noi[mask_nucleus != noi] <- NA
  
  # make threshold image
  threshold_condensates <- EBImage::otsu(ref_noi)*seg_factor
  mask_condensate <- ref_noi > threshold_condensates
  
  # erode and dilate the segmented area
  mask_condensate <- EBImage::erode(mask_condensate, EBImage::makeBrush(erode_disc, shape='disc'))
  mask_condensate <- EBImage::dilate(mask_condensate, EBImage::makeBrush(dilate_disc, shape='disc'))
  mask_condensate <- EBImage::fillHull(mask_condensate) # fill holes in the mask
  mask_condensate <- EBImage::bwlabel(mask_condensate) # find connected sets of pixels
  
  # remove small regions
  regions_cond <- table(mask_condensate)
  small_regions <- which(regions_cond<min_size)-1
  mask_condensate[mask_condensate %in% small_regions] <- 0
  
  # count condensates
  regions_cc <- table(mask_condensate)
  condensate_number <- as.numeric(names(regions_cc))[-1]
  
  # select the nucleoplasm
  mask_nucleoplasm <- img_noi < threshold_condensates/outlier_factor
  mask_nucleoplasm <- EBImage::erode(mask_nucleoplasm,EBImage::makeBrush(2.5, shape='disc') )
  mask_nucleoplasm[mask_condensate > 0] <- 0
  mask_nucleoplasm[is.na(ref_noi)] <- 0
  mask_nucleoplasm <- EBImage::bwlabel(mask_nucleoplasm)
  
  # select largest region as nucleoplasm
  regions_np <- table(mask_nucleoplasm)
  mask_nucleoplasm[(which(regions_np==max(regions_np[2:length(regions_np)]))-1)] <- 1
  
  # return results
  return(list(condensate_number, mask_condensate, mask_nucleoplasm))
}

quantify <- function(ref, img, nuc, seg_factor, dilate_disc, erode_disc, min_size, outlier_factor=4, bckg = 0) {
  # initialize result variable
  result <- rep(NaN, nuc[[2]])
  error <- rep(NA, nuc[[2]])
  image_with_masks <- list()
  ref_with_masks <- list()
  
  # loop through nuclei, segment condensates based on image 'ref' and quantify enrichment based on image 'img'
  for(nuc_i in 1:nuc[[2]]) {
    # select nucleus of interest
    noi <- as.numeric(names(table(nuc[[1]]))[nuc_i+1])
    # segment condensates in nucleus of interest
    
    cond <- segmentCondensates(ref=ref, img=img, mask_nucleus=nuc[[1]], noi=noi, seg_factor=seg_factor, dilate_disc=dilate_disc, erode_disc=erode_disc, min_size=min_size, outlier_factor=outlier_factor)
    
    # make image with masks
    iwm <- EBImage::toRGB(EBImage::normalize(img))
    iwm <- EBImage::paintObjects(cond[[3]], iwm, opac=c(0.2,0.2), col=c("green","green"), thick=T)
    iwm <- EBImage::paintObjects(nuc[[1]]==noi, iwm, opac=c(0.5,0), col=c("magenta","magenta"), thick=T)
    iwm <- EBImage::paintObjects(cond[[2]], iwm, opac=c(0.5,0), col=c("yellow","yellow"), thick=T)
    image_with_masks[[nuc_i]] <- iwm
    
    # make ref image with masks
    iwm2 <- EBImage::toRGB(EBImage::normalize(ref))
    iwm2 <- EBImage::paintObjects(cond[[2]], iwm2, opac=c(0.5,0), col=c("yellow","yellow"), thick=T)
    ref_with_masks[[nuc_i]] <- iwm2
    
    # define nuclear panel
    mask_nucleus <- (nuc[[1]] == noi)
    mask_condensate <- cond[[2]]
    mask_nucleoplasm <- cond[[3]]
    condensate_number <- cond[[1]]
    # print(condensate_number)
    
    if(length(condensate_number)>0){
      # quantify intensity/area of the nucleus
      int_nucleus <- mean(img[mask_nucleus>0], na.rm = TRUE) # intensity of the nucleus
      area_nucleus <- EBImage::computeFeatures.shape(mask_nucleus) # area of the nucleus (in pixels)
      
      # quantify intensity/area of the nucleoplasm
      int_nucleoplasm <- mean(img[mask_nucleoplasm !=0], na.rm = TRUE) # intensity of the nucleoplasm
      area_nucleoplasm <- EBImage::computeFeatures.shape(mask_nucleoplasm)[1] # area of the nucleoplasm (in pixels)
      # quantify intensity/area of the condensates
      int_condensates <- sapply(condensate_number, function(x) {mean(img[mask_condensate==x], na.rm = TRUE) }) # intensity of individual condensates
      area_condensates <- EBImage::computeFeatures.shape(mask_condensate)[,1] # size of individual condensates (in pixels)
      # sd_all_condensates <- sd(img[mask_condensate>0], na.rm = TRUE) # sd over pixels in all condensates
      # print(sd_all_condensates)
      
      # sd_enrichment = sd_all_condensates/int_nucleoplasm
      # print(sd_enrichment)
      # calculate mean/sd of enrichment
      #print(condensate_number)
      enrichment_median <- median((int_condensates-bckg)/(int_nucleoplasm-bckg))
      #print(enrichment_median)
      enrichment_mean <- mean((int_condensates-bckg)/(int_nucleoplasm-bckg))
      #print(enrichment_mean)
      enrichment_sem <- (sd((int_condensates-bckg))/(int_nucleoplasm-bckg))/sqrt(length(int_condensates))
      #print(enrichment_sem)
      if(is.na(enrichment_sem)){
        enrichment_sem=0
      }
      # store enrichment
      result[nuc_i] <- enrichment_median
      error[nuc_i] <- enrichment_sem
    }
    
  }
  if(mean(error,na.rm=T) == 0) {
    error = rep(sd(result, na.rm=T)/sqrt(length(result)), length(result))
  }
  
  # return result
  return(list(enrichment = result, error = error, img = image_with_masks, ref=ref_with_masks))
}
