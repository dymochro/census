# This function calculates the copy numbers of proteins in MLOs of interest

getScaffoldCopyNumbersInMLOs <- function(organelle, info, data) {
  # define variables
  enrich <- data$Enrichment_MLO
  ntot <- data$Number_per_cell
  vnuc <- info$Volume[which(info$Organelle=="Nucleus")]
    
  if(organelle=="Txn condensate") {
    vorg <- info$Volume[which(info$Organelle=="Txn condensate, small")]
    vtot <- info$TotalVolume[which(info$Organelle=="Txn condensate, small")]+info$TotalVolume[which(info$Organelle=="Txn condensate, large")]
  }
  else {
    vorg <- info$Volume[which(info$Organelle==organelle)]
    vtot <- info$TotalVolume[which(info$Organelle==organelle)]
  }
    
  # calculate copy numbers in a small transcriptional condensate
  copy_number <- enrich*ntot*vorg/((enrich-1)*vtot+vnuc)
    
  # error, enrichment
  err_enrich <- (ntot*vorg/((enrich-1)*vtot+vnuc) - enrich*ntot*vorg*vtot/((enrich-1)*vtot+vnuc)^2) * data$Enrichment_MLO_Err
    
  # error, copy number per cell
  err_ntot <- enrich*vorg/((enrich-1)*vtot+vnuc) * data$Number_per_cell_Err
      
  # error, MLO volume
  if(organelle=="Txn condensate") {err_vol <- enrich*ntot/((enrich-1)*vtot+vnuc) * info$Volume_Err[which(info$Organelle=="Txn condensate, small")]}
  else {err_vol <- enrich*ntot/((enrich-1)*vtot+vnuc) * info$Volume_Err[which(info$Organelle==organelle)]}
      
  # error, MLO total volume
  if(organelle=="Txn condensate") {err_totvol <- -(enrich-1)*enrich*ntot*vorg/((enrich-1)*vtot+vnuc)^2 * (info$TotalVolume_Err[which(info$Organelle=="Txn condensate, small")]+info$TotalVolume_Err[which(info$Organelle=="Txn condensate, large")])} 
  else {err_totvol <- -(enrich-1)*enrich*ntot*vorg/((enrich-1)*vtot+vnuc)^2 * info$TotalVolume_Err[which(info$Organelle==organelle)]}
    
  # error, nuclear volume
  err_nucvol <- -enrich*ntot*vorg/((enrich-1)*vtot+vnuc)^2 * info$Volume_Err[which(info$Organelle=="Nucleus")]
    
  # final error
  err_copy_number <- sqrt(err_enrich^2 + err_ntot^2 + err_vol^2 + err_totvol^2 + err_nucvol^2)
    
    
  # calculate copy numbers for proteins with infinite enrichment values
  indices_inf <- which(is.infinite(enrich))
  if(organelle=="Txn condensate") {
    copy_number[indices_inf] <- ntot[indices_inf]*vorg/vtot
    err_copy_number[indices_inf] <- sqrt((data$Number_per_cell_Err[indices_inf]*vorg/vtot)^2 + (ntot[indices_inf]*info$Volume_Err[which(info$Organelle=="Txn condensate, small")]/vtot)^2 + (-ntot[indices_inf]*vorg/vtot^2*(info$TotalVolume_Err[which(info$Organelle=="Txn condensate, small")]+info$TotalVolume_Err[which(info$Organelle=="Txn condensate, large")]))^2)
  }
  else {
    copy_number[indices_inf] <- ntot[indices_inf]/info$Number[which(info$Organelle==organelle)]
    err_copy_number[indices_inf] <- sqrt( (data$Number_per_cell_Err[indices_inf]/info$Number[which(info$Organelle==organelle)])^2 + (-ntot[indices_inf]/info$Number[which(info$Organelle==organelle)]^2*info$Number_Err[which(info$Organelle==organelle)])^2)
  }
  
  # return result
  return(c(copy_number, err_copy_number))
}
  