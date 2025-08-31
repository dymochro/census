# Retrieve sizes of scaffold proteins

getScaffoldSizes <- function(scaffolds, vers="v1", minlen=15, dssp_path = "/opt/anaconda3/envs/intel_env/bin/mkdssp") {
  # load library
  library(bio3d)
  
  # calculate size based on AlphaFold
  size_afold <- rep(0, length(scaffolds$uniprotswissprot))
  size_relaxed <- rep(0, length(scaffolds$uniprotswissprot))
  size_expanded <- rep(0, length(scaffolds$uniprotswissprot))
  size_denatured <- rep(0, length(scaffolds$uniprotswissprot))
  
  for(i in 1:length(scaffolds$uniprotswissprot)) {
    # read pdb structure
    pdb <- bio3d::read.pdb(paste0("https://alphafold.ebi.ac.uk/files/AF-",scaffolds$uniprotswissprot[i],"-F1-model_",vers,".pdb")) 
    
    # make list with the mass of each atom
    mass <- rep(12, length(pdb$xyz)/3)
    mass[substr(pdb$atom[,"elety"], 1, 1) == "N"] <- 14
    mass[substr(pdb$atom[,"elety"], 1, 1) == "H"] <- 1
    mass[substr(pdb$atom[,"elety"], 1, 1) == "O"] <- 16 
    mass[substr(pdb$atom[,"elety"], 1, 1) == "S"] <- 32
    
    # calculate the diameter of gyration in nm (rgyr gives the radius of gyration in Angstrom)
    size_afold[i] <- 2*bio3d::rgyr(pdb, mass)/10
    
    
    ## consider more extended conformations
    
    # get secondary structure annotation
    sse_pdb <- bio3d::dssp(pdb, exefile = dssp_path)$sse
    
    # get disordered patches that are longer than 'minlen'
    index_disorder <- sse_pdb == " " # dssp will return " " for disordered regions
    patches <- bio3d::rle2(index_disorder)  # determine patches with the same structure
    ends <- patches$ind[patches$values==T & patches$lengths<minlen] # end positions of disordered patches shorter than 'minlen'
    starts <- ends - patches$lengths[patches$values==T & patches$lengths<minlen] # start positions of disordered patches shorter than 'minlen'
    index_disorder[unlist(mapply(":", starts, ends))] <- F # remove disordered patches shorter than 'minlen'  
    
    # get the lengths of ordered and disordered patches
    patches <- bio3d::rle2(index_disorder)
    lengths_disordered_patches <- patches$lengths[patches$values==T]
    lengths_ordered_patches <- patches$lengths[patches$values==F]
    ends_ordered_patches <- patches$ind[patches$values==F]
    starts_ordered_patches <- ends_ordered_patches - lengths_ordered_patches
    
    # calculate the sizes of disordered patches
    Rh_disordered <- 0.249*lengths_disordered_patches^0.509 # formula for Rh of IDPs from Marsh & Forman-Kay: https://pubmed.ncbi.nlm.nih.gov/20483348/
    Rg_disordered <- 1.58*Rh_disordered # ratio from Clisby & Duenweg: https://arxiv.org/abs/2001.03138
    
    Rg_ordered <- rep(0, length(lengths_ordered_patches))
    for(j in 1:length(lengths_ordered_patches)) {
      temp_structure <- bio3d::trim(pdb, atom.select(pdb, resno=starts_ordered_patches[j]:ends_ordered_patches[j]))
      mass <- rep(12, length(temp_structure$xyz)/3)
      mass[substr(temp_structure$atom[,"elety"], 1, 1) == "N"] <- 14
      mass[substr(temp_structure$atom[,"elety"], 1, 1) == "H"] <- 1
      mass[substr(temp_structure$atom[,"elety"], 1, 1) == "O"] <- 16 
      mass[substr(temp_structure$atom[,"elety"], 1, 1) == "S"] <- 32
      
      # calculate the diameter of gyration in nm (rgyr gives the radius of gyration in Angstrom)
      Rg_ordered[j] <- bio3d::rgyr(temp_structure, mass)/10
    }
    
    # calculate the diameter of gyration for relaxed and expanded conformations
    diameter_random_walk <- 2*sqrt(sum(Rg_ordered^2)+sum(Rg_disordered^2)) # random walk with steps corresponding to the sizes of ordered and disordered patches
    diameter_extended <- 2*(sum(Rg_ordered)+sum(Rg_disordered)) # linear arrangement of ordered and disordered patches
    
    size_relaxed[i] <- diameter_random_walk
    size_expanded[i] <- mean(c(diameter_extended,diameter_random_walk))
    
    # calculate the diameter of gyration for the denatured protein chain
    # formula for Rh of denatured proteins from Marsh & Forman-Kay: https://pubmed.ncbi.nlm.nih.gov/20483348/
    # conversion from Rh to Rg from Clisby & Duenweg: https://arxiv.org/abs/2001.03138
    size_denatured[i] <- 2*1.58*0.233*nchar(scaffolds$protein_sequence[i])^0.549
    
    # ensure that the relaxed and expanded sizes are not smaller than the size obtained from AlphaFold
    size_relaxed[i] <- max(size_relaxed[i],size_afold[i])
    size_expanded[i] <- max(size_expanded[i],size_afold[i])
    
    # ensure that relaxed and expanded sizes do not exceed the size of the denatured protein
    size_relaxed[i] <- min(size_relaxed[i],size_denatured[i])
    size_expanded[i] <- min(size_expanded[i],size_denatured[i])
  }
  
  # return sizes
  return(c(size_afold, size_relaxed, size_expanded))
}


