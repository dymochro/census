# These functions retrieve the enrichment level for proteins of interest in their condensates in mESCs
# The results are based on multiple papers listed below

# Segmentation results can be displayed using 'display(quanti[[4]][[i]])', with i = index of segmented cell nucleus

getPolycombEnrichment <- function(proteins) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
    
  # initialize result variable
  enrichment <- rep(NaN, length(proteins))
  enrichment_error <- rep(NaN, length(proteins))
    
  # loop through proteins and determine enrichment
  for(i in 1:length(proteins)) {
    
    if(tolower(proteins[i])=="cbx2") {
      # Tatavosian et al, 2019, J Biol Chem 294: 1451-1463.
      # PMID 30514760, https://doi.org/10.1074/jbc.RA118.006620, Fig. 1
      # quantify in Cbx2 channel, segment in Cbx2 channel
              
      # load image data
      target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/c7f3/6364756/3958b6adade1/zbc0051999880001.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      img_grayscale <- apply(img, c(1,2), mean)
      cbx2 <- img_grayscale[4:104, 329:427]

      # quantify enrichment
      quanti <- getEnrichment(ref=cbx2, img=cbx2, 1, 1, dilate_disc=11, erode_disc=9, min_size_nucleus=1000, dilate_disc_foci=1.5, erode_disc_foci=1, min_size_condensate=3, outlier_factor=3)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
      
    if(tolower(proteins[i])=="phc1") {
      # Tatavosian et al, 2019, J Biol Chem 294: 1451-1463.
      # PMID 30514760, https://doi.org/10.1074/jbc.RA118.006620, Fig. 1
      # quantify in Phc1 channel, segment in Cbx2 channel
            
      # load image data
      target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/c7f3/6364756/3958b6adade1/zbc0051999880001.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      img_grayscale <- apply(img, c(1,2), mean)
            
      # crop relevant panels
      cbx2 <- img_grayscale[4:104, 329:427]
      phc1 <- img_grayscale[121:221, 329:427]

      # quantify enrichment
      quanti <- getEnrichment(ref=cbx2, img=phc1, 0.9, 0.6, img4nuc=TRUE, dilate_disc=11, erode_disc=5, min_size_nucleus=800, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=1, outlier_factor=1)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="ube2i") {
      # Tessier et al, 2022, Nat Commun 13: 5726.
      # PMID 36175410, https://doi.org/10.1038/s41467-022-33147-6, Suppl. Fig. 2e (bottom left panel)
      # quantify and segment in Ube2i/Ubc9 channel
      
      # load image data
      target <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-33147-6/MediaObjects/41467_2022_33147_MOESM1_ESM.pdf"
      img <- EBImage::flop(EBImage::rotate(pdftools::pdf_render_page(RCurl::getURLContent(target),page=3,dpi=300,numeric=T),90))
      img_grayscale <- apply(img[,,1:3], c(1,2), mean)
      
      # crop relevant panels
      ube2i <- img_grayscale[213:412,1762:1948]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=ube2i, img=ube2i, 1, 1.5, img4nuc=TRUE, dilate_disc=11, erode_disc=5, min_size_nucleus=800, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=1, outlier_factor=1.5)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="sumo3") {
      # Liu et al, 2020, EMBO J 39: e103697
      # PMID 32395866, https://doi.org/10.15252/embj.2019103697, Suppl. Fig. EV3A
      # quantify and segment in Sumo2/3 channel
      
      # load image data
      target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/73a5/7327501/bf86cc8f40ed/EMBJ-39-e103697-g007.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      img_grayscale <- img[,,1]
      
      # crop relevant panels
      sumo3 <- img_grayscale[155:219,88:151]
      pml <- img_grayscale[88:152,88:151]
      
      # exclude PML bodies (which do not colocalize with Polycomb bodies)
      sumo3[pml>0.05] <- 0
      
      # quantify enrichment
      quanti <- getEnrichment(ref=sumo3, img=sumo3, 1.5, 1.7, img4nuc=FALSE, dilate_disc=3, erode_disc=3, min_size_nucleus=200, dilate_disc_foci=1, erode_disc_foci=1, min_size_condensate=3, outlier_factor=0.5)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="cbx4") {
      # the majority of Cbx4 is considered to be in Polycomb bodies
      enrichment[i] <- Inf
      enrichment_error[i] <- 0
    }
  }
  
  # return result
  return(list(Enrichment = enrichment, Enrich_Err = enrichment_error))
}

getHeterochromatinEnrichment <- function(proteins) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  
  # initialize result variable
  enrichment <- rep(NaN, length(proteins))
  enrichment_error <- rep(NaN, length(proteins))

    # loop through proteins and determine enrichment
    for(i in 1:length(proteins)) {
      
      if(tolower(proteins[i])=="cbx5") {
        # He et al, 2015, Cell Stem Cell 17: 273-286.
        # PMID 26340527, https://doi.org/10.1016/j.stem.2015.07.022, Fig. 6
        # quantify in Cbx5 channel, segment in DAPI channel
        
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S1934590915003215-gr6_lrg.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        img_grayscale <- apply(img, c(1,2), mean)
        
        # crop relevant panels
        cbx5 <- img_grayscale[900:1247, 2200:2565]
        dapi <- img_grayscale[1279:1636, 1796:2161]
        
        # remove captions
        cbx5[0:133,0:70] <- NA
        dapi[0:133,0:70] <- NA
        
        # register DAPI channel
        dapi <- EBImage::Image(RNiftyReg::niftyreg(dapi,cbx5)$image)
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=cbx5, 0.25, 1.35, dilate_disc=5, erode_disc=5, min_size_nucleus=250, min_size_condensate=50)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
      
      if(tolower(proteins[i])=="hmga2") {
        # Yu et al, 2014, Cell Rep 6: 684-697.
        # PMID 24508460, https://doi.org/10.1016/j.celrep.2014.01.014, Fig. 1
        # quantify in Hmga2 channel, segment in DAPI channel
        
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S221112471400031X-gr1_lrg.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        hmga2 <- img[1209:1549, 841:1195, 2]
        dapi <- img[1209:1549, 841:1195, 3]
        
        # remove captions
        hmga2 <- hmga2[,1:300]
        dapi <- dapi[,1:300]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=hmga2, 0.5, 1.5, dilate_disc=5.5, erode_disc=7.5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=20, outlier_factor=1.5)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
        if(is.na(enrichment_error[i])) {enrichment_error[i] <- 0}
      }

      if(tolower(proteins[i])=="cbx1") {
        # Ballmer et al, 2022, Nucleic Acids Res 51: 117–143.
        # PMID 36533441, https://doi.org/10.1093/nar/gkac1159, Fig. 6B
        # quantify in Cbx1 channel, segment in DAPI channel
        
        # load image data
        img_row1 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341618&s=157&r=3&c=1")),90))
        img_row1 <- abind(img_row1, flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341618&s=157&r=3&c=2")),90)), along=1)
        img_row2 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341618&s=157&r=4&c=1")),90))
        img_row2 <- abind(img_row2, flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=341618&s=157&r=4&c=2")),90)), along=1)
        img <- abind(img_row1, img_row2, along=2)
        
        # crop relevant panels
        cbx1 <- img[188:337, 197:348, 1]
        dapi <- img[33:182, 197:348, 3]
        
        # remove scale bars
        cbx1[99:145, 131:141] <- 0
        dapi[99:145, 131:141] <- 0
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=cbx1, 0.6, 1.5, dilate_disc=5, erode_disc=7.5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=20, outlier_factor=2)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }

      if(tolower(proteins[i])=="trim28") {
        # Ding et al, 2018, J Biol Chem 293: 2711-2724.
        # PMID 29284678, https://doi.org/10.1074/jbc.RA117.000959, Fig. 5C 
        # quantify in Kap1/Trim28 channel, segment in DAPI channel
        
        # load image data
        target <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5827453/bin/zbc0081881910005.jpg"  
        target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/80ac/5827453/3167f5ed882c/zbc0081881910005.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        trim28 <- img[320:463, 388:515,1]
        dapi <- img[170:313, 388:515,1]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=trim28, 0.5, 1.5, dilate_disc=5, erode_disc=3.5, dilate_disc_foci=1, erode_disc_foci=1, min_size_nucleus=100, min_size_condensate=1, outlier_factor=0.8)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }

      if(tolower(proteins[i])=="mecp2") {
        # Li et al 2020, Nature 586: 440-444.
        # PMID 32698189, https://doi.org/10.1038/s41586-020-2574-4, Fig. 1
        # quantify in Mecp2 channel, segment in DAPI channel
        
        # load image data
        target <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7735819/bin/nihms-1646018-f0001.jpg"
        target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/8e58/7735819/ade1e1d259d1/nihms-1646018-f0001.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        mecp2 <- img[32:127, 20:114, 2] - img[32:127, 20:114, 3]
        dapi <- img[134:229, 20:114, 3] - img[134:229, 20:114, 2]
        mecp2[mecp2<0] <- 0
        dapi[dapi<0] <- 0
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=mecp2, 0.4, 0.9, dilate_disc=5, erode_disc=7.5, dilate_disc_foci=2.5, erode_disc_foci=3, min_size_nucleus=900, min_size_condensate=10, outlier_factor=2)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
      
      if(tolower(proteins[i])=="suv39h1") {
        # Endogenous Suv39h1 not visible, Tosolini et al, 2018, https://doi.org/10.1038/s41598-018-23822-4 (Fig. 2B)
        # Enrichment thus determined based on GFP-Suv39h1
        
        # Mallm et al 2015, Cell Reports 11: 1-12.
        # PMID 26051938, https://doi.org/10.1016/j.celrep.2015.05.015, Fig. 3A
        # quantify and segment in Suv39h1 channel
        
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S2211124715005276-gr3.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        suv39h1 <- img[33:97, 157:221, 1]
         
        # quantify enrichment
        quanti <- getEnrichment(ref=suv39h1, img=suv39h1, 0.7, 1.2, dilate_disc=5, erode_disc=7.5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=10, outlier_factor=1.5)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
      
      if(tolower(proteins[i])=="kmt5c") {
        # Hahn et al 2013, Genes Dev 27: 859-872.
        # PMID 23599346, https://doi.org/10.1101/gad.210377.112, Fig. S2
        # quantify in Kmt5c channel, segment in DAPI channel
        
        # load image data
        target <- "https://genesdev.cshlp.org/content/suppl/2013/04/11/gad.210377.112.DC1/Supplemental_figure_S2.pdf"
        img <- EBImage::flop(EBImage::rotate(pdftools::pdf_render_page(RCurl::getURLContent(target),page=1,dpi=300,numeric=T),90))
        
        # crop relevant panels
        kmt5c <- img[549:884, 2413:2749, 1]
        dapi <- img[209:544, 2413:2749, 1]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=dapi, img=kmt5c, 1, 1.8, dilate_disc=5, erode_disc=7.5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=100, outlier_factor=3)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }

      if(tolower(proteins[i])=="ube2i") {
        # Tessier et al, 2022, Nat Commun 13: 5726.
        # PMID 36175410, https://doi.org/10.1038/s41467-022-33147-6, Suppl. Fig. 2e (bottom left panel)
        # quantify and segment in Ube2i/Ubc9 channel (no heterochromatin reference, max. enrichment taken)
        
        # load image data
        target <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-33147-6/MediaObjects/41467_2022_33147_MOESM1_ESM.pdf"
        img <- EBImage::flop(EBImage::rotate(pdftools::pdf_render_page(RCurl::getURLContent(target),page=3,dpi=300,numeric=T),90))
        img_grayscale <- apply(img[,,1:3], c(1,2), mean)
        
        # crop relevant panels
        ube2i <- img_grayscale[213:412,1762:1948]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=ube2i, img=ube2i, 1, 1.5, img4nuc=TRUE, dilate_disc=11, erode_disc=5, min_size_nucleus=800, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=1, outlier_factor=1.5)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
    }
    
  # return result
  return(list(Enrichment = enrichment, Enrich_Err = enrichment_error))
}

getTxnEnrichment <- function(proteins) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  library(magick)
  
  # initialize result variable
  enrichment <- rep(NaN, length(proteins))
  enrichment_error <- rep(NaN, length(proteins))
  
  # loop through proteins and determine enrichment
  for(i in 1:length(proteins)) {

      if(tolower(proteins[i])=="cdk9") {
        # Ghamari et al, 2013, Genes Dev 27: 767–777.
        # PMID 23592796, https://doi.org/10.1101/gad.216200.113, Movie S1
        # quantify in Cdk9 channel, segment in Cdk9 channel
        
        # load image data
        target <- "https://genesdev.cshlp.org/content/suppl/2013/04/16/27.7.767.DC1/Supplemental_Movie_S1.avi"
        mov <- RCurl::getURLContent(target, binary = TRUE)
        tf <- tempfile()
        readr::write_file(mov, tf)
        img_series <- EBImage::readImage(av::av_video_images(tf))
        img <- apply(img_series, c(1,2), mean)
        
        # crop relevant panel
        cdk9 <- img[1:255, 27:280]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=cdk9, img=cdk9, 0.9, 1.4, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=10, outlier_factor=1)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
    
      if(tolower(proteins[i])=="ccnt1") {
        # Jang et al, 2005, Mol Cell 19: 523-34
        # PMID 16109376, https://doi.org/10.1016/j.molcel.2005.06.027, Fig. 4 (3T3 cells)
        # quantify in Ccnt1 channel, segment in Brd4 channel
        
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S1097276505014322-gr4_lrg.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      
        # crop relevant panels
        ccnt1 <- img[63:453,93:606,2]
        brd4 <- img[521:911,93:606,1]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=brd4, img=ccnt1, 0.1, 2, img4nuc=T, dilate_disc=7.5, erode_disc=7.5, dilate_disc_foci=3, erode_disc_foci=3, min_size_nucleus=900, min_size_condensate=100, outlier_factor=2)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
    }
    
      if(tolower(proteins[i])=="pou5f1"){
        # Verneri et al, 2020, Sci Rep 10: 5195.
        # PMID 32251342, https://doi.org/10.1038/s41598-020-62235-0, Fig. 2
        # quantify in Oct4/Pou5f1 channel, segment in Oct4/Pou5f1 channel

        # load image data
        target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41598-020-62235-0/MediaObjects/41598_2020_62235_Fig2_HTML.png"
        img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
        
        # crop relevant panel
        oct4 <- img[185:430,95:330,2]
        
        # quantify enrichment
        quanti <- getEnrichment(ref=oct4, img=oct4, 1.1, 1.4, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=5, outlier_factor=1)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
    
      if(tolower(proteins[i])=="brd3") {
        # Banani et al, 2022, Dev Cell 57: 1776-1788.e8
        # PMID 35809564, https://doi.org/10.1016/j.devcel.2022.06.010, Fig. 3C
        # quantify in Brd3 channel, segment in Brd3 channel
        
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S1534580722004506-gr3_lrg.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panel
        brd3 <- img[1516:1962,961:1407,2] - img[1516:1962,961:1407,3]
        brd3[brd3<0] <- 0
        
        # quantify enrichment
        quanti <- getEnrichment(ref=brd3, img=brd3, 0.3, 2, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=100, outlier_factor=1)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
  
      if(tolower(proteins[i])=="brd4") {
        # Sabari et al, 2018, Science 361:eaar3958.
        # PMID 29930091, https://doi.org/10.1126/science.aar3, Fig. 1
        # quantify in Brd4 channel, segment in Brd4 channel
        
        # load image data
        target <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6092193/bin/nihms-983922-f0001.jpg"
        target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/963d/6092193/89fc5bd52e09/nihms-983922-f0001.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        brd4 <- img[138:241, 50:146, 2]-img[138:241, 50:146, 1]
        brd4[brd4<0] <- 0
        
        # quantify enrichment
        quanti <- getEnrichment(ref=brd4, img=brd4, 0.2, 1, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=5, outlier_factor=1.5)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
      
      if(tolower(proteins[i])=="med1") {
        # Sabari et al, 2018, Science 361:eaar3958.
        # PMID 29930091, https://doi.org/10.1126/science.aar3, Fig. 1
        # quantify in Med1 channel, segment in Med1 channel
        
        # load image data
        target <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6092193/bin/nihms-983922-f0001.jpg"
        target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/963d/6092193/89fc5bd52e09/nihms-983922-f0001.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
        # crop relevant panels
        med1 <- img[138:241, 156:255, 2]-img[138:241, 156:255, 1]
        med1[med1<0] <- 0
        dapi <- img[138:241, 156:255, 3]-img[138:241, 156:255, 1] # for segmentation of the nucleus
        dapi[dapi<0] <- 0
        
        # quantify enrichment
        quanti <- getEnrichment(ref=0.8*med1+0.2*dapi, img=med1, 0.2, 0.7, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=5, outlier_factor=1.5)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }

      if(tolower(proteins[i])=="smarca4") {
        # Efroni et al, 2008, Cell Stem Cell 2: 437-447.
        # PMID 18462694, https://doi.org/10.1016/j.stem.2008.03.021, Fig. 6D
        # quantify in Brg1/Smarca4 channel, segment in Brg1/Smarca4 channel
      
        # load image data
        target <- "https://ars.els-cdn.com/content/image/1-s2.0-S1934590908001616-gr6_lrg.jpg"
        img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      
        # crop relevant panel
        smarca4 <- img[885:1524,1667:2173,2]
        
        # remove scale bar
        smarca4[466:630,434:497] <- 0
      
        # quantify enrichment
        quanti <- getEnrichment(ref=smarca4, img=smarca4, 1, 1.5, dilate_disc=5, erode_disc=5, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=900, min_size_condensate=2, outlier_factor=1.2)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }
    
      if(tolower(proteins[i])=="polr2a") {
        # Guo et al, 2019, Nature 572:543-548.
        # PMID 31391587, https://doi.org/10.1038/s41586-019-1464-0, Fig. 3
        # quantify in Polr2a channel, segment in Med1 channel
        
        # load image data
        img_row1 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=551140&s=91&r=1&c=1")),90))
        img_row2 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=551140&s=91&r=2&c=1")),90))
        for(c in 2:4) {
          img_row1 <- abind(img_row1, flop(EBImage::rotate(jpeg::readJPEG(getURLContent(paste0("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=551140&s=91&r=1&c=",c))),90)), along=1)
          img_row2 <- abind(img_row2, flop(EBImage::rotate(jpeg::readJPEG(getURLContent(paste0("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=551140&s=91&r=2&c=",c))),90)), along=1)
        }
        img <- abind(img_row1, img_row2, along=2)
        img_grayscale <- abs(img[,,1]-img[,,2])
        
        # crop relevant panels
        polr2a <- img_grayscale[480:807,148:474]
        polr2a[polr2a<0] <- 0
        med1 <- img_grayscale[132:459,148:474]
        med1[med1<0] <- 0

        # quantify enrichment
        quanti <- getEnrichment(ref=0.8*polr2a+0.2*img[480:807,148:474,2], img=polr2a, 1, 2, dilate_disc=21, erode_disc=11, dilate_disc_foci=1.5, erode_disc_foci=1.5, min_size_nucleus=9000, min_size_condensate=30, outlier_factor=2)
        enrichment[i] <- quanti[[1]]
        enrichment_error[i] <- quanti[[2]]
      }   
    
      if(tolower(proteins[i])=="mllt3" | tolower(proteins[i])=="nsd2") {
        # the majority of these (lowly abundant) proteins are considered to be present in txn condensates
        enrichment[i] <- Inf
        enrichment_error[i] <- 0
      }
    }
    
  # return result
  return(list(Enrichment = enrichment, Enrich_Err = enrichment_error))
}

getNucleolusDFCEnrichment <- function(proteins) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  
  # initialize result variable
  enrichment <- rep(NaN, length(proteins))
  enrichment_error <- rep(NaN, length(proteins))
  
  # loop through proteins and determine enrichment
  for(i in 1:length(proteins)) {
    
    if(tolower(proteins[i])=="fbl") {
      # Yu, Sun, Tan, Pan et al, 2021, Nat Commun 12: 6365
      # PMID 34753899, https://doi.org/10.1038/s41467-021-26576-2, Fig. 3d
      # quantify and segment in Fbl channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-021-26576-2/MediaObjects/41467_2021_26576_Fig3_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
        
      # crop relevant panels
      fbl <- img[360:630, 565:810,1]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=fbl, img=fbl, 0.3, 1.6, dilate_disc=3, erode_disc=5, min_size_nucleus=200, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=5, outlier_factor=2)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="polr1b") {
      # Jiang et al, 2020, Genome Biol 21: 158
      # PMID 32616013, https://doi.org/10.1186/s13059-020-02067-3, Fig. 1C
      # quantify and segment in Polr1b channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-02067-3/MediaObjects/13059_2020_2067_Fig1_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      polr1b <- img[301:412, 493:594,2]
      dapi <- img[184:295, 493:594,2]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=0.5*polr1b+0.5*dapi, img=polr1b, 0.3, 1.1, dilate_disc=3, erode_disc=5, min_size_nucleus=200, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=5, outlier_factor=1)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }

    if(tolower(proteins[i])=="lin28a") {
      # Sun, Yu, Zhao et al, 2021, Protein Cell 13: 490–512
      # PMID 34331666, https://doi.org/10.1007/s13238-021-00864-5, Fig. 1C
      # quantify and segment in Lin28a channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs13238-021-00864-5/MediaObjects/13238_2021_864_Fig1_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
        
      # crop relevant panels
      lin28a <- img[71:339, 649:892,2]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=lin28a, img=lin28a, 0.4, 1, dilate_disc=3, erode_disc=5, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=3, min_size_condensate=10, outlier_factor=2)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }  
    
    if(tolower(proteins[i])=="eloa") {
      # Ardehali et al, 2021, J Biol Chem 296: 100202
      # PMID 33334895, https://doi.org/10.1074/jbc.RA120.015877, Fig. 5C (3T3 cells)
      # quantify and segment in Eloa channel
      
      # load image data
      target <- "https://ars.els-cdn.com/content/image/1-s2.0-S0021925820001982-gr5.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
        
      # crop relevant panels
      eloa <- img[19:177, 401:530,1]
      
      # remove text
      eloa[1:60,115:130] <- 0
      
      # quantify enrichment
      quanti <- getEnrichment(ref=eloa, img=eloa, 0.55, 1, dilate_disc=5, erode_disc=1.5, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=6.5, min_size_condensate=10, outlier_factor=2)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }  
  }
  
  # return result
  return(list(Enrichment = enrichment, Enrich_Err = enrichment_error))
}

getNucleolusGCEnrichment <- function(proteins) {
  ## load library and function for image segmentation
  source(paste0(path,"_getEnrichment.R"))
  library(EBImage)
  
  # initialize result variable
  enrichment <- rep(NaN, length(proteins))
  enrichment_error <- rep(NaN, length(proteins))
  
  # loop through proteins and determine enrichment
  for(i in 1:length(proteins)) {
    
    if(tolower(proteins[i])=="ncl") {
      # Yu, Sun, Tan, Pan et al, 2021, Nat Commun 12: 6365
      # PMID 34753899, https://doi.org/10.1038/s41467-021-26576-2, Fig. 3b
      # quantify and segment in Ncl channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-021-26576-2/MediaObjects/41467_2021_26576_Fig3_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      ncl <- img[651:867, 52:267, 3]
      dapi <- img[430:646, 52:267, 3]
      
      # remove captions
      ncl[0:217,0:29] <- 0 # remove text
      dapi[0:217,0:29] <- 0 # remove same area from the image that is used to make the nuclear mask
      
      # quantify enrichment
      quanti <- getEnrichment(ref=0.5*ncl+0.5*dapi, img=ncl, 0.4, 1.2, dilate_disc=3, erode_disc=3, min_size_nucleus=800, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=400, outlier_factor=2)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="npm1") {
      # Yu, Sun, Tan, Pan et al, 2021, Nat Commun 12: 6365
      # PMID 34753899, https://doi.org/10.1038/s41467-021-26576-2, Fig. 3b
      # quantify and segment in Npm1 channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-021-26576-2/MediaObjects/41467_2021_26576_Fig3_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      npm1 <- img[1403:1622, 51:266, 3]
      dapi <- img[1173:1392, 51:266, 3]
      
      # remove captions
      npm1[1:220,1:26] <- 0 # remove text
      dapi[1:220,1:26] <- 0 # remove same area from the image that is used to make the nuclear mask
      
      # quantify enrichment
      quanti <- getEnrichment(ref=0.5*npm1+0.5*dapi, img=npm1, 1, 1.2, dilate_disc=3, erode_disc=3, min_size_nucleus=800, dilate_disc_foci=3, erode_disc_foci=3, min_size_condensate=600, outlier_factor=1.6)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="surf6") {
      # Moraleva et al, 2017, Cell Cycle 16: 1979–1991
      # PMID 28873013, https://doi.org/10.1080/15384101.2017.1371880, Fig. 1G (3T3)
      # quantify and segment in Surf6 channel
      
      # load image data
      target <- "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/4a0f/5638364/4e55307e773d/kccy-16-20-1371880-g001.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      surf6 <- img[373:570,367:514,1]-img[373:570,367:514,2]
      surf6[surf6<0] <- 0
      dapi <- img[371:568,527:674,3]-img[371:568,527:674,1]
      dapi[dapi<0] <- 0
      
      # quantify enrichment
      quanti <- getEnrichment(ref=0.5*surf6+0.5*dapi, img=surf6, 0.6, 1.4, dilate_disc=3, erode_disc=3, min_size_nucleus=50, dilate_disc_foci=1, erode_disc_foci=1, min_size_condensate=10, outlier_factor=1.6)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="lin28a") {
      # Sun, Yu, Zhao et al, 2021, Protein Cell 13: 490–512
      # PMID 34331666, https://doi.org/10.1007/s13238-021-00864-5, Fig. 1C
      # quantify and segment in Lin28a channel
      
      # load image data
      target <- "https://media.springernature.com/full/springer-static/image/art%3A10.1007%2Fs13238-021-00864-5/MediaObjects/13238_2021_864_Fig1_HTML.png"
      img <- EBImage::flop(EBImage::rotate(png::readPNG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      lin28a <- img[71:339, 649:892,2]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=lin28a, img=lin28a, 0.4, 0.8, dilate_disc=3, erode_disc=5, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=3, min_size_condensate=10, outlier_factor=1.5)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }  
    
    if(tolower(proteins[i])=="eloa") {
      # Ardehali et al, 2021, J Biol Chem 296: 100202
      # PMID 33334895, https://doi.org/10.1074/jbc.RA120.015877, Fig. 5C (3T3 cells)
      # quantify and segment in Eloa channel
      
      # load image data
      target <- "https://ars.els-cdn.com/content/image/1-s2.0-S0021925820001982-gr5.jpg"
      img <- EBImage::flop(EBImage::rotate(jpeg::readJPEG(RCurl::getURLContent(target)),90))
      
      # crop relevant panels
      eloa <- img[19:177, 401:530,1]
      
      # remove text
      eloa[1:60,115:130] <- 0
      
      # quantify enrichment
      quanti <- getEnrichment(ref=eloa, img=eloa, 0.55, 0.65, dilate_disc=5, erode_disc=1.5, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=6.5, min_size_condensate=100, outlier_factor=1)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="rpl5") {
      # Yu et al, 2006, Mol Cell Biol 26: 3798–3809
      # PMID 16648475, https://doi.org/10.1128/MCB.26.10.3798-3809.2006, Fig. 7B (3T3 cells)
      # quantify and segment in Rpl5 channel
      
      # load image data
      img_row1 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=12&c=1")),90))
      img_row2 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=13&c=1")),90))
      img_row3 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=14&c=1")),90))
      for(c in 2:4) {
        img_row1 <- abind(img_row1, flop(EBImage::rotate(jpeg::readJPEG(getURLContent(paste0("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=12&c=",c))),90)), along=1)
        img_row2 <- abind(img_row2, flop(EBImage::rotate(jpeg::readJPEG(getURLContent(paste0("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=13&c=",c))),90)), along=1)
        img_row3 <- abind(img_row3, flop(EBImage::rotate(jpeg::readJPEG(getURLContent(paste0("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=465465&s=2&r=14&c=",c))),90)), along=1)
      }
      img <- abind(img_row1, img_row2, img_row3, along=2)
      img_grayscale <- abs(img[,,1]-img[,,2])
      
      # crop relevant panels
      rpl5 <- img[(152+480):865, 33:765, 2]
      nuc <- img[(152+480):865, 33:765, 1]
      nuc[1:5,110:460] <- 1
      
      # remove dashed line
      mask <- img[(152+480):865, 33:765, 1]
      rpl5[mask>0.2] <- 0
      
      # quantify enrichment
      quanti <- getEnrichment(ref=0.2*nuc+0.8*rpl5, img=rpl5, 1, 1.6, dilate_disc=21, erode_disc=19, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=6.5, min_size_condensate=600, outlier_factor=1)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="rpl23a") {
      # Plafker & Macara, 2002, Mol Cell Biol 22: 1266–1275.
      # PMID 11809816, https://doi.org/10.1128/MCB.22.4.1266-1275.2002, Fig. 2C (BHK cells)
      # quantify and segment in Rpl23a channel
      
      # load image data
      img_row1 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=112768&s=10&r=2&c=3")),90))
      img_row2 <- flop(EBImage::rotate(jpeg::readJPEG(getURLContent("https://www.ncbi.nlm.nih.gov/corecgi/tileshop/tileshop.fcgi?p=PMC3&id=112768&s=10&r=3&c=3")),90))
      img <- abind(img_row1, img_row2, along=2)
      
      # crop relevant panels
      rpl23a <- img[30:256, 207:450, 2]
      
      # quantify enrichment
      quanti <- getEnrichment(ref=rpl23a, img=rpl23a, 0.4, 1.2, dilate_disc=3, erode_disc=3, min_size_nucleus=900, dilate_disc_foci=1, erode_disc_foci=6.5, min_size_condensate=200, outlier_factor=2)
      enrichment[i] <- quanti[[1]]
      enrichment_error[i] <- quanti[[2]]
    }
    
    if(tolower(proteins[i])=="gnl2") {
      # the majority of this protein is considered to be in the nucleolus
      enrichment[i] <- Inf
      enrichment_error[i] <- 0
    }
  }
  
  # return result
  return(list(Enrichment = enrichment, Enrich_Err = enrichment_error))
}