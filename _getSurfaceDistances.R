# Retrieve intermolecular distances between surfaces of neighboring scaffold proteins
# Considering one species or multiple species together

getMultipleSpeciesSurfaceDistances <- function(info, scaffolds, organelle, rna, nuc) {
  # calculate distances between surfaces of neighboring proteins, considering all protein species
  ravg_afold <- (sum(scaffolds$Number_per_MLO*(scaffolds$Size_AlphaFold/2)^3)/sum(scaffolds$Number_per_MLO))^(1/3)
  phi_afold <- 4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]/sum(scaffolds$Number_per_MLO))
  distances_allp_afold <- exp(8*phi_afold)*ravg_afold/3/phi_afold^(1/3)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)
  distances_allp_afold <- distances_allp_afold*(phi_afold<0.65)
  
  err_phi_afold <- sqrt(sum(scaffolds$Number_per_MLO_Err^2) * (4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allp_afold <- sqrt(err_phi_afold^2*(ravg_afold/9/phi_afold^(4/3)*(6*phi_afold^(1/3)-exp(8*phi_afold)*(24*phi_afold-1)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)))^2)
  err_distances_allp_afold <- err_distances_allp_afold
  
  ravg_relaxed <- (sum(scaffolds$Number_per_MLO*(scaffolds$Size_relaxed/2)^3)/sum(scaffolds$Number_per_MLO))^(1/3)
  phi_relaxed <- 4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]/sum(scaffolds$Number_per_MLO))
  distances_allp_relaxed <- exp(8*phi_relaxed)*ravg_relaxed/3/phi_relaxed^(1/3)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)
  distances_allp_relaxed <- distances_allp_relaxed*(phi_relaxed<0.65)
  
  err_phi_relaxed <- sqrt(sum(scaffolds$Number_per_MLO_Err^2) * (4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allp_relaxed <- sqrt(err_phi_relaxed^2*(ravg_relaxed/9/phi_relaxed^(4/3)*(6*phi_relaxed^(1/3)-exp(8*phi_relaxed)*(24*phi_relaxed-1)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)))^2)
  err_distances_allp_relaxed <- err_distances_allp_relaxed
  
  ravg_expanded <- (sum(scaffolds$Number_per_MLO*(scaffolds$Size_expanded/2)^3)/sum(scaffolds$Number_per_MLO))^(1/3)
  phi_expanded <- 4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]/sum(scaffolds$Number_per_MLO))
  distances_allp_expanded <- exp(8*phi_expanded)*ravg_expanded/3/phi_expanded^(1/3)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)
  distances_allp_expanded <- distances_allp_expanded*(phi_expanded<0.65)
  
  err_phi_expanded <- sqrt(sum(scaffolds$Number_per_MLO_Err^2) * (4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allp_expanded <- sqrt(err_phi_expanded^2*(ravg_expanded/9/phi_expanded^(4/3)*(6*phi_expanded^(1/3)-exp(8*phi_expanded)*(24*phi_expanded-1)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)))^2)
  err_distances_allp_expanded <- err_distances_allp_expanded
  
  # prepare data frame with resulting distances
  res <- as.data.frame(c(distances_allp_afold, distances_allp_relaxed, distances_allp_expanded))
  names(res) <- "Distance, all proteins"
  
  # add errors to data frame with results
  res$`Distance_Err, all proteins` <- c(err_distances_allp_afold, err_distances_allp_relaxed, err_distances_allp_expanded)
  
   
  # calculate distances between surfaces of neighboring proteins, considering all protein species and RNA
  rna_size <- 2*getRNASize(rna[which(rna$Organelle==organelle),]$RNuc_number/rna[which(rna$Organelle==organelle),]$RNA_number)
  rna_size_expanded <- 2*getRNASizeExpanded(rna[which(rna$Organelle==organelle),]$RNuc_number/rna[which(rna$Organelle==organelle),]$RNA_number)
  rna_size_relaxed <- mean(c(rna_size,rna_size_expanded))
  
  ravg_afold <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_AlphaFold/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))^(1/3)
  phi_afold <- 4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))
  distances_allpr_afold <- exp(8*phi_afold)*ravg_afold/3/phi_afold^(1/3)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)
  distances_allpr_afold <- distances_allpr_afold*(phi_afold<0.65)
  
  err_phi_afold <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2) * (4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allpr_afold <- sqrt(err_phi_afold^2*(ravg_afold/9/phi_afold^(4/3)*(6*phi_afold^(1/3)-exp(8*phi_afold)*(24*phi_afold-1)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)))^2)
  err_distances_allpr_afold <- err_distances_allpr_afold
  
  ravg_relaxed <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_relaxed/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size_relaxed/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))^(1/3)
  phi_relaxed <- 4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))
  distances_allpr_relaxed <- exp(8*phi_relaxed)*ravg_relaxed/3/phi_relaxed^(1/3)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)
  distances_allpr_relaxed <- distances_allpr_relaxed*(phi_relaxed<0.65)
  
  err_phi_relaxed <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2) * (4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allpr_relaxed <- sqrt(err_phi_relaxed^2*(ravg_relaxed/9/phi_relaxed^(4/3)*(6*phi_relaxed^(1/3)-exp(8*phi_relaxed)*(24*phi_relaxed-1)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)))^2)
  err_distances_allpr_relaxed <- err_distances_allpr_relaxed
  
  ravg_expanded <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_expanded/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size_expanded/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))^(1/3)
  phi_expanded <- 4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number))
  distances_allpr_expanded <- exp(8*phi_expanded)*ravg_expanded/3/phi_expanded^(1/3)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)
  distances_allpr_expanded <- distances_allpr_expanded*(phi_expanded<0.65)
  
  err_phi_expanded <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2) * (4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allpr_expanded <- sqrt(err_phi_expanded^2*(ravg_expanded/9/phi_expanded^(4/3)*(6*phi_expanded^(1/3)-exp(8*phi_expanded)*(24*phi_expanded-1)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)))^2)
  err_distances_allpr_expanded <- err_distances_allpr_expanded
  
  # add distances to data frame with results
  res$`Distance, all proteins and RNAs` <- c(distances_allpr_afold, distances_allpr_relaxed, distances_allpr_expanded)
  # add errors to data frame with results
  res$`Distance_Err, all proteins and RNAs` <- c(err_distances_allpr_afold, err_distances_allpr_relaxed, err_distances_allpr_expanded)
  
  
  # calculate distances between surfaces of neighboring proteins, considering all protein species and RNA and nucleosomes
  nuc_size <- 11 # nucleosome diameter in nm
  
  ravg_afold <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_AlphaFold/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size/2)^3+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number*(nuc_size/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))^(1/3)
  phi_afold <- 4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))
  distances_allprn_afold <- exp(8*phi_afold)*ravg_afold/3/phi_afold^(1/3)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)
  distances_allprn_afold <- distances_allprn_afold*(phi_afold<0.65)
  
  err_phi_afold <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number^2) * (4/3*pi*ravg_afold^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allprn_afold <- sqrt(err_phi_afold^2*(ravg_afold/9/phi_afold^(4/3)*(6*phi_afold^(1/3)-exp(8*phi_afold)*(24*phi_afold-1)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)))^2)
  err_distances_allprn_afold <- err_distances_allprn_afold
  
  ravg_relaxed <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_relaxed/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size_relaxed/2)^3+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number*(nuc_size/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))^(1/3)
  phi_relaxed <- 4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))
  distances_allprn_relaxed <- exp(8*phi_relaxed)*ravg_relaxed/3/phi_relaxed^(1/3)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)
  distances_allprn_relaxed <- distances_allprn_relaxed*(phi_relaxed<0.65)
  
  err_phi_relaxed <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number^2) * (4/3*pi*ravg_relaxed^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allprn_relaxed <- sqrt(err_phi_relaxed^2*(ravg_relaxed/9/phi_relaxed^(4/3)*(6*phi_relaxed^(1/3)-exp(8*phi_relaxed)*(24*phi_relaxed-1)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)))^2)
  err_distances_allprn_relaxed <- err_distances_allprn_relaxed
  
  ravg_expanded <- ((sum(scaffolds$Number_per_MLO*(scaffolds$Size_expanded/2)^3)+rna[which(rna$Organelle==organelle),]$RNA_number*(rna_size_expanded/2)^3+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number*(nuc_size/2)^3)/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))^(1/3)
  phi_expanded <- 4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]/(sum(scaffolds$Number_per_MLO)+rna[which(rna$Organelle==organelle),]$RNA_number+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number))
  distances_allprn_expanded <- exp(8*phi_expanded)*ravg_expanded/3/phi_expanded^(1/3)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)
  distances_allprn_expanded <- distances_allprn_expanded*(phi_expanded<0.65)
  
  err_phi_expanded <- sqrt((sum(scaffolds$Number_per_MLO_Err^2)+rna[which(rna$Organelle==organelle),]$RNA_number^2+nuc[which(nuc$Organelle==organelle),]$Nucleosome_number^2) * (4/3*pi*ravg_expanded^3/(1e9*info$Volume[which(info$Organelle==organelle)]))^2)
  err_distances_allprn_expanded <- sqrt(err_phi_expanded^2*(ravg_expanded/9/phi_expanded^(4/3)*(6*phi_expanded^(1/3)-exp(8*phi_expanded)*(24*phi_expanded-1)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)))^2)
  err_distances_allprn_expanded <- err_distances_allprn_expanded
  
  # add distances to data frame with results
  res$`Distance, all proteins and RNAs and nucleosomes` <- c(distances_allprn_afold, distances_allprn_relaxed, distances_allprn_expanded)
  # add errors to data frame with results
  res$`Distance_Err, all proteins and RNAs and nucleosomes` <- c(err_distances_allprn_afold, err_distances_allprn_relaxed, err_distances_allprn_expanded)
  
  # return results
  return(res)
}

getOneSpeciesSurfaceDistances <- function(info, scaffolds, organelle) {
  # calculate distances between surfaces of neighboring proteins, considering one protein species at a time
  phi_afold <- 4/3*pi*(scaffolds$Size_AlphaFold/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)]/scaffolds$Number_per_MLO)
  distances_afold <- exp(8*phi_afold)*scaffolds$Size_AlphaFold/2/3/phi_afold^(1/3)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)
  distances_afold <- distances_afold*(phi_afold<0.65)
  
  phi_relaxed <- 4/3*pi*(scaffolds$Size_relaxed/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)]/scaffolds$Number_per_MLO)
  distances_relaxed <- exp(8*phi_relaxed)*scaffolds$Size_relaxed/2/3/phi_relaxed^(1/3)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)
  distances_relaxed <- distances_relaxed*(phi_relaxed<0.65)
  
  phi_expanded <- 4/3*pi*(scaffolds$Size_expanded/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)]/scaffolds$Number_per_MLO)
  distances_expanded <- exp(8*phi_expanded)*scaffolds$Size_expanded/2/3/phi_expanded^(1/3)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)
  distances_expanded <- distances_expanded*(phi_expanded<0.65)
  
  err_phi_afold <- scaffolds$Number_per_MLO_Err * 4/3*pi*(scaffolds$Size_AlphaFold/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)])
  err_distances_afold <- err_phi_afold*(scaffolds$Size_AlphaFold/2/9/phi_afold^(4/3)*(6*phi_afold^(1/3)-exp(8*phi_afold)*(24*phi_afold-1)*gamma(1/3)*pgamma(8*phi_afold,1/3,1,lower=FALSE)))
  err_distances_afold <- err_distances_afold
  
  err_phi_relaxed <- scaffolds$Number_per_MLO_Err * 4/3*pi*(scaffolds$Size_relaxed/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)])
  err_distances_relaxed <- err_phi_relaxed*(scaffolds$Size_relaxed/2/9/phi_relaxed^(4/3)*(6*phi_relaxed^(1/3)-exp(8*phi_relaxed)*(24*phi_relaxed-1)*gamma(1/3)*pgamma(8*phi_relaxed,1/3,1,lower=FALSE)))
  err_distances_relaxed <- err_distances_relaxed
  
  err_phi_expanded <- scaffolds$Number_per_MLO_Err * 4/3*pi*(scaffolds$Size_expanded/2)^3/(1e9*info$Volume[which(info$Organelle==organelle)])
  err_distances_expanded <- err_phi_expanded*(scaffolds$Size_expanded/2/9/phi_expanded^(4/3)*(6*phi_expanded^(1/3)-exp(8*phi_expanded)*(24*phi_expanded-1)*gamma(1/3)*pgamma(8*phi_expanded,1/3,1,lower=FALSE)))
  err_distances_expanded <- err_distances_expanded
  
  # return distances
  return(c(distances_afold, distances_relaxed, distances_expanded, err_distances_afold, err_distances_relaxed, err_distances_expanded))
}


# radius of gyration of RNA (in nm)
# Hyeon et al, 2006, J Chem Phys 125: 194905, PMID 17129165, https://doi.org/10.1063/1.2364190
getRNASize <- function(N) {
  # N is the number of nucleotides
  return((5.5*N^(1/3))/10)
}

# radius of gyration of RNA as worm-like chain (in nm)
getRNASizeExpanded <- function(N, lc_monomer=0.696) {
  # N is the number of nucleotides
  # lc_monomer is the contour lenght per monomer, which is 0.696 nm for nucleotides: Chi et al, 2013, Physica A, 392: 1072-1079, https://doi.org/10.1016/j.physa.2012.09.022
  lc <- lc_monomer*N # contour length of the entire RNA molecules
  lp <- (1.5*N^0.33)/10 # persistence length of RNA: Hyeon et al, 2006, J Chem Phys 125: 194905, PMID 17129165, https://doi.org/10.1063/1.2364190
  Rg2 <- (1/3)*lp*lc-(lp^2)+2*((lp^3)/lc)*(1-(lp/lc)*(1-exp(-lc/lp))) # radius of gyration
  return(sqrt(Rg2))
}

