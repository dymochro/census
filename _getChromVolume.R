## retrieve the volume occupied by amino acids and nucleotides of a nucleosome (with linker histone)

getChromVolume <- function() {
  ## calculate volume fractions occupied by histone proteins
  histone_seq <- c(
    H1 = "MSETAPVAQAASTATEKPAAAKKTKKPAKAAAPRKKPAGPSVSELIVQAVSSSKERSGVSLAALKKSLAAAGYDVEKNNSRIKLGLKSLVNKGTLVQTKGTGAAGSFKLNKKAESKAITTKVSVKAKASGAAKKPKKTAGAAAKKTVKTPKKPKKPAVSKKTSKSPKKPKVVKAKKVAKSPAKAKAVKPKASKAKVTKPKTPAKPKKAAPKKK", #https://www.uniprot.org/uniprot/P43275
    H3 = "MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEACEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA", #https://www.uniprot.org/uniprot/P68433
    H4 = "MSGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", #https://www.uniprot.org/uniprot/P62806
    H2A = "MSGPTKRGGKARAKVKSRSSRAGLQFPVGRVHRLLRQGNYAQRIGAGAPVYLAAVLEYLTAEVLELAGNAARDNKKTRITPRHLQLAIRNDEELNKLLGRVTIAQGGVLPNIQAVLLPKKTESHKSQTK", #https://www.uniprot.org/uniprot/Q8CGP4
    H2B = "MPEPSKSAPAPKKGSKKAISKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIASEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK" #https://www.uniprot.org/uniprot/Q64475
  )
  
  # volumes of histone proteins
  vol <- getAAVolumes(histone_seq)
  
  # volume of protein portion of a nucleosome, including the histone octamer and one linker histone
  vprot <- 2*sum(vol[2:5]) + vol[1]
  
  # volume of 186 bp of DNA, nucleosome repeat length in mESCs
  # Teif et al, 2012, Nat Struct Mol Biol 19: 1185-1192.
  # PMID 23085715, https://doi.org/10.1038/nsmb.2419, Fig. 7b
  vdna_cm3 <- 0.55 * 327*(1.6605e-24) * 2*186 # cm^3
  vdna <- vdna_cm3 * 1e12 # µm^3
  
  # molecular weight of deoxyribonucleotide monophosphates: 327 Da
  # partial specific volume of DNA: 0.55 cm^3/g
  
  return(vprot + vdna)
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
