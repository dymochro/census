## determine scaffold protein content of organelles of interest

getScaffoldProteinContent <- function(info, gc, dfc, txn, het, pcg) {
  ## calculate volume fractions occupied by scaffold proteins
  gc_spvf <- sum(gc$Number_per_MLO*getAAVolumes(gc$protein_sequence))/info$Volume[which(info$Organelle=="Nucleolus, GC")]
  dfc_spvf <- sum(dfc$Number_per_MLO*getAAVolumes(dfc$protein_sequence))/info$Volume[which(info$Organelle=="Nucleolus, DFC")]
  txn_spvf <- sum(txn$Number_per_MLO*getAAVolumes(txn$protein_sequence))/info$Volume[which(info$Organelle=="Txn condensate, small")]
  het_spvf <- sum(het$Number_per_MLO*getAAVolumes(het$protein_sequence))/info$Volume[which(info$Organelle=="Heterochromatin focus")]
  pcg_spvf <- sum(pcg$Number_per_MLO*getAAVolumes(pcg$protein_sequence))/info$Volume[which(info$Organelle=="Polycomb body")]
  nuc_spvf <- (gc_spvf*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]+dfc_spvf*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]+
               txn_spvf*(info$TotalVolume[which(info$Organelle=="Txn condensate, small")]+info$TotalVolume[which(info$Organelle=="Txn condensate, large")])+
               het_spvf*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]+pcg_spvf*info$TotalVolume[which(info$Organelle=="Polycomb body")]) / 
               info$TotalVolume[which(info$Organelle=="Nucleus")]

  gc_mf <- sum(gc$Number_per_MLO*getAAMasses(gc$protein_sequence))/info$Volume[which(info$Organelle=="Nucleolus, GC")]
  dfc_mf <- sum(dfc$Number_per_MLO*getAAMasses(dfc$protein_sequence))/info$Volume[which(info$Organelle=="Nucleolus, DFC")]
  txn_mf <- sum(txn$Number_per_MLO*getAAMasses(txn$protein_sequence))/info$Volume[which(info$Organelle=="Txn condensate, small")]
  het_mf <- sum(het$Number_per_MLO*getAAMasses(het$protein_sequence))/info$Volume[which(info$Organelle=="Heterochromatin focus")]
  pcg_mf <- sum(pcg$Number_per_MLO*getAAMasses(pcg$protein_sequence))/info$Volume[which(info$Organelle=="Polycomb body")]
  nuc_mf <- (gc_mf*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]+dfc_mf*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")]+
                 txn_mf*(info$TotalVolume[which(info$Organelle=="Txn condensate, small")]+info$TotalVolume[which(info$Organelle=="Txn condensate, large")])+
                 het_mf*info$TotalVolume[which(info$Organelle=="Heterochromatin focus")]+pcg_mf*info$TotalVolume[which(info$Organelle=="Polycomb body")]) / 
    info$TotalVolume[which(info$Organelle=="Nucleus")]
  
  ## assemble and return results
  result <- data.frame(
    Organelle = c("Nucleolus", "Nucleolus, GC", "Nucleolus, DFC", "Txn condensate, small", "Txn condensate, large", "Heterochromatin focus", "Polycomb body", "Nucleus"),
    ScaffoldProteinVolumeFraction = c((gc_spvf*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]+dfc_spvf*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")])/info$Volume[which(info$Organelle=="Nucleolus")], gc_spvf, dfc_spvf, txn_spvf, txn_spvf, het_spvf, pcg_spvf, nuc_spvf),
    ScaffoldProteinMassConcentration = 1e12*c((gc_mf*info$TotalVolume[which(info$Organelle=="Nucleolus, GC")]+dfc_mf*info$TotalVolume[which(info$Organelle=="Nucleolus, DFC")])/info$Volume[which(info$Organelle=="Nucleolus")], gc_mf, dfc_mf, txn_mf, txn_mf, het_mf, pcg_mf, nuc_mf),
    stringsAsFactors = FALSE)

  return(result)
}

getAAVolumes <- function(scaffolds) {
  # load library
  library(seqinr)
  
  # get partial specific volumes of each amino acid
  data(aaindex)
  partial_specific_volumes_per_aa <- aaindex$BULH740102$I # cm^3/g
  names(partial_specific_volumes_per_aa) <- seqinr::a(names(partial_specific_volumes_per_aa))
  
  # get molecular mass of each amino acid
  mol_mass_per_aa <- unlist(lapply(names(partial_specific_volumes_per_aa), function(x) seqinr::pmw(x))) # Dalton
  mol_mass_per_aa <- mol_mass_per_aa*(1.6605e-24) # g
  
  # get volume of each amino acid
  volumes_per_aa <- partial_specific_volumes_per_aa*mol_mass_per_aa # cm^3
  volumes_per_aa <- volumes_per_aa*1e12 # µm^3
  
  # get volumes of amino acids of scaffold proteins
  volumes <- numeric()
  for(i in 1:length(scaffolds)) {
    volumes[i] <- unlist(lapply(strsplit(scaffolds[i], split=""), function(x) sum(volumes_per_aa[x])))
  }
  
  # return volumes (given in µm^3)
  return(volumes)
}

getAAMasses <- function(scaffolds) {
  # load library
  library(seqinr)
  
  # get partial specific volumes of each amino acid
  data(aaindex)
  partial_specific_volumes_per_aa <- aaindex$BULH740102$I # cm^3/g
  names(partial_specific_volumes_per_aa) <- seqinr::a(names(partial_specific_volumes_per_aa))
  
  # get molecular mass of each amino acid
  mol_mass_per_aa <- unlist(lapply(names(partial_specific_volumes_per_aa), function(x) seqinr::pmw(x))) # Dalton
  mol_mass_per_aa <- mol_mass_per_aa*(1.6605e-24) # g
  names(mol_mass_per_aa) <- names(partial_specific_volumes_per_aa)
  
  # get volumes of amino acids of scaffold proteins
  masses <- numeric()
  for(i in 1:length(scaffolds)) {
    masses[i] <- unlist(lapply(strsplit(scaffolds[i], split=""), function(x) sum(mol_mass_per_aa[x])))
  }
  
  # return masses (given in mg)
  return(masses*1000)
}
