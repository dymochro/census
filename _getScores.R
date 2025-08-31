## determine scores for an MLO of interest
## scores reflect how well the properties of an MLO are in line with different condensation mechanisms

getScores <- function(distances) {
  # score intermolecular distances
  if(distances$`Distance, all proteins`[2]<=2.2) {dist_p_score <- 1} else {dist_p_score <- exp(-(distances$`Distance, all proteins`[2]-2.2)/2.2)}
  if(distances$`Distance, all proteins and RNAs`[2]<=2.2) {dist_pr_score <- 1} else {dist_pr_score <- exp(-(distances$`Distance, all proteins and RNAs`[2]-2.2)/2.2)}
  if(distances$`Distance, all proteins and RNAs and nucleosomes`[2]<=2.2) {dist_prn_score <- 1} else {dist_prn_score <- exp(-(distances$`Distance, all proteins and RNAs and nucleosomes`[2]-2.2)/2.2)}
  
  # return result
  return(c(dist_p_score,dist_pr_score,dist_prn_score)) 
}