## Retrieve information about condensates reconstituted in vitro (as reference)
## uses the function getAAVolumes that should be loaded

getReferenceInfo <- function(refs) {
  # prepare results object
  res <- refs
  
  # concentrations
  # Npm1
  # Ferrolino et al, 2018, Nat Commun 9: 5064.
  # PMID 30498217, https://doi.org/10.1038/s41467-018-07530-1, Supplementary Table 1 (5% PEG)
  
  # Cbx5
  # Her et al, 2022, Nucleic Acids Res 50: 12702–12722.
  # PMID 36537242, https://doi.org/10.1093/nar/gkac1194, Supplementary Figure 3d (HP1)
  
  res$Concentration <- c(4.6, 10) # mM
  res$Concentration_Err <- c(0.2, 1) # mM
  
  res$MassConcentration <- c(4.6*32575/1000, 10*22181/1000) # mg/mL
  res$MassConcentration_Err <- c(0.2*32575/1000, 1*22181/1000) # mg/mL
  
  # volume fractions
  volumes <- getAAVolumes(res$protein_sequence) # µm^3
  res$ScaffoldProteinVolumeFraction <- volumes*res$Concentration*10^(-3)*(6.022*10^23)*10^(-15)
  res$ScaffoldProteinVolumeFraction_Err <- volumes*res$Concentration_Err*10^(-3)*(6.022*10^23)*10^(-15)
  
  # intermolecular distances
  phi_afold <- 4/3*pi*(res$Size_AlphaFold/2)^3/(1e9*volumes/res$ScaffoldProteinVolumeFraction)
  res$Distance_AlphaFold <- exp(8*phi_afold)*res$Size_AlphaFold/2/3/phi_afold^(1/3)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)
  res$Distance_AlphaFold <- res$Distance_AlphaFold*(phi_afold<0.65)
  
  phi_relaxed <- 4/3*pi*(res$Size_relaxed/2)^3/(1e9*volumes/res$ScaffoldProteinVolumeFraction)
  res$Distance_relaxed <- exp(8*phi_relaxed)*res$Size_relaxed/2/3/phi_relaxed^(1/3)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)
  res$Distance_relaxed <- res$Distance_relaxed*(phi_relaxed<0.65)
  
  phi_expanded <- 4/3*pi*(res$Size_expanded/2)^3/(1e9*volumes/res$ScaffoldProteinVolumeFraction)
  res$Distance_expanded <- exp(8*phi_expanded)*res$Size_expanded/2/3/phi_expanded^(1/3)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)
  res$Distance_expanded <- res$Distance_expanded*(phi_expanded<0.65)
  
  err_phi_afold <- res$ScaffoldProteinVolumeFraction_Err * 4/3*pi*(res$Size_AlphaFold/2)^3/(1e9*volumes)
  res$Distance_AlphaFold_Err <- err_phi_afold*(res$Size_AlphaFold/2/9/phi_afold^(4/3)*(6*phi_afold^(1/3)-exp(8*phi_afold)*(24*phi_afold-1)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)))
  res$Distance_AlphaFold_Err <- res$Distance_AlphaFold_Err
  
  err_phi_relaxed <- res$ScaffoldProteinVolumeFraction_Err * 4/3*pi*(res$Size_relaxed/2)^3/(1e9*volumes)
  res$Distance_relaxed_Err <- err_phi_relaxed*(res$Size_relaxed/2/9/phi_relaxed^(4/3)*(6*phi_relaxed^(1/3)-exp(8*phi_relaxed)*(24*phi_relaxed-1)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)))
  res$Distance_relaxed_Err <- res$Distance_relaxed_Err
  
  err_phi_expanded <- res$ScaffoldProteinVolumeFraction_Err * 4/3*pi*(res$Size_expanded/2)^3/(1e9*volumes)
  res$Distance_expanded_Err <- err_phi_expanded*(res$Size_expanded/2/9/phi_expanded^(4/3)*(6*phi_expanded^(1/3)-exp(8*phi_expanded)*(24*phi_expanded-1)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)))
  res$Distance_expanded_Err <- res$Distance_expanded_Err
  
  # return result
  return(res)
}

