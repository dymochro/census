## information about membrane-less organelles of interest in mESCs (number, diameter, volume)

# functions to convert between radii/diameters and volumes of spheres
getVolumeSphere <- function(r) {
	v <- (4/3) * pi * (r^3)
	return(v)
}

getVolumeErrorSphere <- function(r, re) {
  dv_dr <- 4 * pi * (r^2)
  return(dv_dr*re)
}


getMLOInfo <- function() {
  ## cell nucleus
  # Faro-Trindade and Cook, 2006, Mol Biol Cell 17: 2910-20.
  # PMID 16624866, https://doi.org/10.1091/mbc.e05-11-1024, Table 1

  nucleus_number <- 1
  nucleus_number_error <- 0

  nucleus_volume <- 840 + 480 # µm^3
  nucleus_volume_error <- 240 + 134 # µm^3

  nucleus_total_volume <- nucleus_number*nucleus_volume # µm^3
  nucleus_total_volume_error <- sqrt((nucleus_number_error*nucleus_volume)^2 + (nucleus_number*nucleus_volume_error)^2) # µm^3

  
  ## nucleoplasm
  # Faro-Trindade and Cook, 2006, Mol Biol Cell 17: 2910-20.
  # PMID 16624866, https://doi.org/10.1091/mbc.e05-11-1024, Table 1
  
  nucleoplasm_number <- 1
  nucleoplasm_number_error <- 0
  
  nucleoplasm_volume <- 840 # µm^3
  nucleoplasm_volume_error <- 240 # µm^3
  
  nucleoplasm_total_volume <- nucleoplasm_number*nucleoplasm_volume # µm^3
  nucleoplasm_total_volume_error <- sqrt((nucleoplasm_number_error*nucleoplasm_volume)^2 + (nucleoplasm_number*nucleoplasm_volume_error)^2) # µm^3
  

  ## nucleolus, total
  # Gupta and Santoro, 2020, Stem Cell Reports 15: 1206-1219.
  # PMID 32976768, https://doi.org/10.1016/j.stemcr.2020.08.012

  nucleolus_number <- 1
  nucleolus_number_error <- 0

  # Faro-Trindade and Cook, 2006, Mol Biol Cell 17: 2910-20.
  # PMID 16624866, https://doi.org/10.1091/mbc.e05-11-1024, Table 1

  nucleolus_volume <- 480 # µm^3
  nucleolus_volume_error <- 134 # µm^3

  nucleolus_total_volume <- nucleolus_number*nucleolus_volume # µm^3
  nucleolus_total_volume_error <- sqrt((nucleolus_number_error*nucleolus_volume)^2 + (nucleolus_number*nucleolus_volume_error)^2) # µm^3


  ## nucleolus, GC and DFC (assuming that relative sizes between GC and DFC are similar in mESCs and HeLa cells)
  # Caragine et al, 2019, eLife 8: e47533.
  # PMID 31769409, https://doi.org/10.7554/eLife.47533, section "Physical nature of nucleolar subcompartments" following Fig. 3

  hela_gc_fraction <- 0.9
  hela_dfc_fraction <- 0.1

  # Shan et al, 2023, Nature 615: 526-534.
  # PMID 36890225, https://doi.org/10.1038/s41586-023-05767-5, Extended Data Fig. 9g

  dfc_number <- 43
  dfc_number_error <- 3

  gc_number <- 1
  gc_number_error <- 0

  gc_volume <- nucleolus_volume*hela_gc_fraction # µm^3
  gc_volume_error <- nucleolus_volume_error*hela_gc_fraction # µm^3
  dfc_volume <- nucleolus_volume*hela_dfc_fraction/dfc_number # µm^3
  dfc_volume_error <- sqrt((nucleolus_volume_error*hela_dfc_fraction/dfc_number)^2 + (nucleolus_volume*hela_dfc_fraction/dfc_number^2*dfc_number_error)^2) # µm^3

  gc_total_volume <- gc_number*gc_volume # µm^3
  gc_total_volume_error <- sqrt((gc_number_error*gc_volume)^2 + (gc_number*gc_volume_error)^2) # µm^3

  dfc_total_volume <- dfc_number*dfc_volume # µm^3
  dfc_total_volume_error <- sqrt((dfc_number_error*dfc_volume)^2 + (dfc_number*dfc_volume_error)^2) # µm^3



  ## transcriptional condensates
  # Sabari et al, 2018, Science 361: eaar3958
  # PMID 29930091, https://doi.org/10.1126/science.aar3958, Table S1 (MED1)

  txn_number <- 983
  txn_number_error <- 102

  # Cho et al, 2018, Science 361: 412-415.
  # PMID 29930094, https://doi.org/10.1126/science.aar4199, text around Fig. 1B/C (Mediator)

  txn_small_diameter <- 0.2 # µm, 100 nm radius
  txn_large_diameter <- 0.6 # µm, 300 nm radius
  txn_large_number <- 14
  txn_small_number <- txn_number - txn_large_number

  # assume 20% error margins (no errors indicated in Cho et al, 2018) 
  txn_small_number_error <- 0.2*txn_small_number
  txn_large_number_error <- 0.2*txn_large_number
  txn_small_diameter_error <- 0.2*txn_small_diameter # µm
  txn_large_diameter_error <- 0.2*txn_large_diameter # µm

  txn_small_volume <- getVolumeSphere(txn_small_diameter/2) # µm^3
  txn_small_volume_error <- getVolumeErrorSphere(txn_small_diameter/2, txn_small_diameter_error/2) # µm^3
  txn_large_volume <- getVolumeSphere(txn_large_diameter/2) # µm^3
  txn_large_volume_error <- getVolumeErrorSphere(txn_large_diameter/2, txn_large_diameter_error/2) # µm^3

  txn_small_total_volume <- txn_small_number*txn_small_volume # µm^3
  txn_small_total_volume_error <- sqrt((txn_small_number_error*txn_small_volume)^2 + (txn_small_number*txn_small_volume_error)^2) # µm^3
  txn_large_total_volume <- txn_large_number*txn_large_volume # µm^3
  txn_large_total_volume_error <- sqrt((txn_large_number_error*txn_large_volume)^2 + (txn_large_number*txn_large_volume_error)^2) # µm^3



  ## heterochromatin foci
  # Rausch, Weber et al, 2020, Nucleic Acids Res 48: 12751-12777.
  # PMID 33264404, https://doi.org/10.1093/nar/gkaa1124, Suppl. Table 9

  heterochromatin_number <- 10 # median for undiff. mESCs
  heterochromatin_number_error <- 1.24 # sem for undiff. mESCs

  heterochromatin_volume <- 4.8 # µm^3, median for undiff. mESCs
  heterochromatin_volume_error <- 0.63 # µm^3, sem for undiff. mESCs

  heterochromatin_total_volume <- heterochromatin_number*heterochromatin_volume # µm^3
  heterochromatin_total_volume_error <- sqrt((heterochromatin_number_error*heterochromatin_volume)^2 + (heterochromatin_number*heterochromatin_volume_error)^2) # µm^3



  ## Polycomb bodies
  # Tatavosian et al, 2019, J Biol Chem 294: 1451-1463.
  # PMID 30514760, https://doi.org/10.1074/jbc.RA118.006620, Fig. 4f/g

  polycomb_number <- 13
  polycomb_number_error <- 4

  polycomb_area <- 0.133 # µm^2
  polycomb_area_error <- 0.033 # µm^2
  polycomb_volume <- getVolumeSphere(sqrt(polycomb_area/pi)) # µm^3
  polycomb_volume_error <- getVolumeErrorSphere(sqrt(polycomb_area/pi), polycomb_area_error/2/sqrt(polycomb_area*pi)) # µm^3

  polycomb_total_volume <- polycomb_number*polycomb_volume # µm^3
  polycomb_total_volume_error <- sqrt((polycomb_number_error*polycomb_volume)^2 + (polycomb_number*polycomb_volume_error)^2) # µm^3



  ## assemble and return information about the nucleus and the MLOs of interest
  info <- data.frame(
    Organelle = c("Nucleolus", "Nucleolus, GC", "Nucleolus, DFC", "Txn condensate, small", "Txn condensate, large", "Heterochromatin focus", "Polycomb body", "Nucleoplasm", "Nucleus"),
    Volume = c(nucleolus_volume, gc_volume, dfc_volume, txn_small_volume, txn_large_volume, heterochromatin_volume, polycomb_volume, nucleoplasm_volume, nucleus_volume),
    Volume_Err = c(nucleolus_volume_error, gc_volume_error, dfc_volume_error, txn_small_volume_error, txn_large_volume_error, heterochromatin_volume_error, polycomb_volume_error, nucleoplasm_volume_error, nucleus_volume_error),
    Number = c(nucleolus_number, gc_number, dfc_number, txn_small_number, txn_large_number, heterochromatin_number, polycomb_number, nucleoplasm_number, nucleus_number),
    Number_Err = c(nucleolus_number_error, gc_number_error, dfc_number_error, txn_small_number_error, txn_large_number_error, heterochromatin_number_error, polycomb_number_error, nucleoplasm_number_error, nucleus_number_error),
    TotalVolume = c(nucleolus_total_volume, gc_total_volume, dfc_total_volume, txn_small_total_volume, txn_large_total_volume, heterochromatin_total_volume, polycomb_total_volume, nucleoplasm_total_volume, nucleus_total_volume),
    TotalVolume_Err = c(nucleolus_total_volume_error, gc_total_volume_error, dfc_total_volume_error, txn_small_total_volume_error, txn_large_total_volume_error, heterochromatin_total_volume_error, polycomb_total_volume_error, nucleoplasm_total_volume_error, nucleus_total_volume_error),
    stringsAsFactors = FALSE)

    # units: all volumes given in µm^3
  
  return(info)
}
