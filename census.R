### Molecular census of membrane-less organelles ###

# load functions
path <- "/Users/census/census_code/"
path <- "/Users/fabian/Documents/Projects/Cheryn/_draft/_submission/GenomeBiology/_revisions2/census_Rcode/"
source(paste0(path,"_getScores.R"))
source(paste0(path,"_getMLOInfo.R"))
source(paste0(path,"_getScaffolds.R"))
source(paste0(path,"_getRNAContent.R"))
source(paste0(path,"_getChromVolume.R"))
source(paste0(path,"_getScaffoldSizes.R"))
source(paste0(path,"_getReferenceInfo.R"))
source(paste0(path,"_getSurfaceDistances.R"))
source(paste0(path,"_getNucleosomeContent.R"))
source(paste0(path,"_getScaffoldCopyNumbers.R"))
source(paste0(path,"_getScaffoldProteinContent.R"))
source(paste0(path,"_getScaffoldEnrichmentsInMLOs.R"))
source(paste0(path,"_getScaffoldCopyNumbersInMLOs.R"))

# general information about membrane-less organelles of interest (number, diameter, volume) from literature
MLO_info <- getMLOInfo()

# nucleosome content of organelles of interest
nucleosome_content <- getNucleosomeContent(path, MLO_info)

# RNA content of organelles of interest
rna_content <- getRNAContent(path, MLO_info)

# proteome of organelles of interest
nucleolus_gc_scaffolds <- getScaffolds("Nucleolus, GC")
nucleolus_dfc_scaffolds <- getScaffolds("Nucleolus, DFC")
txn_scaffolds <- getScaffolds("Txn condensate")
heterochromatin_scaffolds <- getScaffolds("Heterochromatin focus")
polycomb_scaffolds <- getScaffolds("Polycomb body")

# protein copy numbers per cell
nucleolus_gc_scaffolds[c("Number_per_cell","Number_per_cell_Err")] <- getScaffoldCopyNumbers(nucleolus_gc_scaffolds$external_gene_name, nucleosome_content=nucleosome_content)
nucleolus_dfc_scaffolds[c("Number_per_cell","Number_per_cell_Err")] <- getScaffoldCopyNumbers(nucleolus_dfc_scaffolds$external_gene_name, nucleosome_content=nucleosome_content)
txn_scaffolds[c("Number_per_cell","Number_per_cell_Err")] <- getScaffoldCopyNumbers(txn_scaffolds$external_gene_name, nucleosome_content=nucleosome_content)
heterochromatin_scaffolds[c("Number_per_cell","Number_per_cell_Err")] <- getScaffoldCopyNumbers(heterochromatin_scaffolds$external_gene_name, nucleosome_content=nucleosome_content)
polycomb_scaffolds[c("Number_per_cell","Number_per_cell_Err")] <- getScaffoldCopyNumbers(polycomb_scaffolds$external_gene_name, nucleosome_content=nucleosome_content)

# protein enrichments in organelles of interest
nucleolus_gc_scaffolds[c("Enrichment_MLO","Enrichment_MLO_Err")] <- getNucleolusGCEnrichment(nucleolus_gc_scaffolds$external_gene_name)
nucleolus_dfc_scaffolds[c("Enrichment_MLO","Enrichment_MLO_Err")] <- getNucleolusDFCEnrichment(nucleolus_dfc_scaffolds$external_gene_name)
txn_scaffolds[c("Enrichment_MLO","Enrichment_MLO_Err")] <- getTxnEnrichment(txn_scaffolds$external_gene_name)
heterochromatin_scaffolds[c("Enrichment_MLO","Enrichment_MLO_Err")] <- getHeterochromatinEnrichment(heterochromatin_scaffolds$external_gene_name)
polycomb_scaffolds[c("Enrichment_MLO","Enrichment_MLO_Err")] <- getPolycombEnrichment(polycomb_scaffolds$external_gene_name)

# protein copy numbers in organelles of interest
nucleolus_gc_scaffolds[c("Number_per_MLO","Number_per_MLO_Err")] <- getScaffoldCopyNumbersInMLOs("Nucleolus, GC", MLO_info, nucleolus_gc_scaffolds)
nucleolus_dfc_scaffolds[c("Number_per_MLO","Number_per_MLO_Err")] <- getScaffoldCopyNumbersInMLOs("Nucleolus, DFC", MLO_info, nucleolus_dfc_scaffolds)
txn_scaffolds[c("Number_per_MLO","Number_per_MLO_Err")] <- getScaffoldCopyNumbersInMLOs("Txn condensate", MLO_info, txn_scaffolds)
heterochromatin_scaffolds[c("Number_per_MLO","Number_per_MLO_Err")] <- getScaffoldCopyNumbersInMLOs("Heterochromatin focus", MLO_info, heterochromatin_scaffolds)
polycomb_scaffolds[c("Number_per_MLO","Number_per_MLO_Err")] <- getScaffoldCopyNumbersInMLOs("Polycomb body", MLO_info, polycomb_scaffolds)

# scaffold protein content of organelles of interest
scaffold_protein_content <- getScaffoldProteinContent(MLO_info, nucleolus_gc_scaffolds, nucleolus_dfc_scaffolds, txn_scaffolds, heterochromatin_scaffolds, polycomb_scaffolds)

# determine protein sizes
nucleolus_gc_scaffolds[c("Size_AlphaFold","Size_relaxed","Size_expanded")] <- getScaffoldSizes(nucleolus_gc_scaffolds)
nucleolus_dfc_scaffolds[c("Size_AlphaFold","Size_relaxed","Size_expanded")] <- getScaffoldSizes(nucleolus_dfc_scaffolds)
txn_scaffolds[c("Size_AlphaFold","Size_relaxed","Size_expanded")] <- getScaffoldSizes(txn_scaffolds)
heterochromatin_scaffolds[c("Size_AlphaFold","Size_relaxed","Size_expanded")] <- getScaffoldSizes(heterochromatin_scaffolds)
polycomb_scaffolds[c("Size_AlphaFold","Size_relaxed","Size_expanded")] <- getScaffoldSizes(polycomb_scaffolds)

# calculate intermolecular distances in organelles of interest (considering one species at a time)
nucleolus_gc_scaffolds[c("Distance_AlphaFold","Distance_relaxed","Distance_expanded","Distance_AlphaFold_Err","Distance_relaxed_Err","Distance_expanded_Err")] <- getOneSpeciesSurfaceDistances(MLO_info, nucleolus_gc_scaffolds, "Nucleolus, GC")
nucleolus_dfc_scaffolds[c("Distance_AlphaFold","Distance_relaxed","Distance_expanded","Distance_AlphaFold_Err","Distance_relaxed_Err","Distance_expanded_Err")] <- getOneSpeciesSurfaceDistances(MLO_info, nucleolus_dfc_scaffolds, "Nucleolus, DFC")
txn_scaffolds[c("Distance_AlphaFold","Distance_relaxed","Distance_expanded","Distance_AlphaFold_Err","Distance_relaxed_Err","Distance_expanded_Err")] <- getOneSpeciesSurfaceDistances(MLO_info, txn_scaffolds, "Txn condensate, small")
heterochromatin_scaffolds[c("Distance_AlphaFold","Distance_relaxed","Distance_expanded","Distance_AlphaFold_Err","Distance_relaxed_Err","Distance_expanded_Err")] <- getOneSpeciesSurfaceDistances(MLO_info, heterochromatin_scaffolds, "Heterochromatin focus")
polycomb_scaffolds[c("Distance_AlphaFold","Distance_relaxed","Distance_expanded","Distance_AlphaFold_Err","Distance_relaxed_Err","Distance_expanded_Err")] <- getOneSpeciesSurfaceDistances(MLO_info, polycomb_scaffolds, "Polycomb body")

# adjust distances for Med1 and Polr2b in transcriptional condensates because they are part of well-known complexes
# Mediator-PIC complex: Diameter of 40 nm (PMID 24550107)
# RNA Pol II complex: Diameter of 15 nm (PMIDs 30575770, 40436841, 31155237)
# The number of Polr2a molecules is reduced by one as the Mediator-PIC complex, which is considered separately, already includes one Polr2a
txn_scaffolds_complexes <- txn_scaffolds
txn_scaffolds_complexes$Size_AlphaFold[which(txn_scaffolds_complexes$external_gene_name=="Med1")] <- 40
txn_scaffolds_complexes$Size_relaxed[which(txn_scaffolds_complexes$external_gene_name=="Med1")] <- 40*txn_scaffolds$Size_relaxed[which(txn_scaffolds$external_gene_name=="Med1")]/txn_scaffolds$Size_AlphaFold[which(txn_scaffolds$external_gene_name=="Med1")]
txn_scaffolds_complexes$Size_expanded[which(txn_scaffolds_complexes$external_gene_name=="Med1")] <- 40*txn_scaffolds$Size_expanded[which(txn_scaffolds$external_gene_name=="Med1")]/txn_scaffolds$Size_AlphaFold[which(txn_scaffolds$external_gene_name=="Med1")]
txn_scaffolds_complexes$Number_per_MLO[which(txn_scaffolds_complexes$external_gene_name=="Polr2a")] <- txn_scaffolds_complexes$Number_per_MLO[which(txn_scaffolds_complexes$external_gene_name=="Polr2a")]-1
txn_scaffolds_complexes$Size_AlphaFold[which(txn_scaffolds_complexes$external_gene_name=="Polr2a")] <- 15
txn_scaffolds_complexes$Size_relaxed[which(txn_scaffolds_complexes$external_gene_name=="Polr2a")] <- 15*txn_scaffolds$Size_relaxed[which(txn_scaffolds$external_gene_name=="Polr2a")]/txn_scaffolds$Size_AlphaFold[which(txn_scaffolds$external_gene_name=="Polr2a")]
txn_scaffolds_complexes$Size_expanded[which(txn_scaffolds_complexes$external_gene_name=="Polr2a")] <- 15*txn_scaffolds$Size_expanded[which(txn_scaffolds$external_gene_name=="Polr2a")]/txn_scaffolds$Size_AlphaFold[which(txn_scaffolds$external_gene_name=="Polr2a")]

# calculate intermolecular distances in organelles of interest (considering all species together)
nucleolus_gc_distances_multi <- getMultipleSpeciesSurfaceDistances(MLO_info, nucleolus_gc_scaffolds, "Nucleolus, GC", rna_content, nucleosome_content) 
nucleolus_dfc_distances_multi <- getMultipleSpeciesSurfaceDistances(MLO_info, nucleolus_dfc_scaffolds, "Nucleolus, DFC", rna_content, nucleosome_content) 
txn_distances_multi <- getMultipleSpeciesSurfaceDistances(MLO_info, txn_scaffolds_complexes, "Txn condensate, small", rna_content, nucleosome_content) 
heterochromatin_distances_multi <- getMultipleSpeciesSurfaceDistances(MLO_info, heterochromatin_scaffolds, "Heterochromatin focus", rna_content, nucleosome_content) 
polycomb_distances_multi <- getMultipleSpeciesSurfaceDistances(MLO_info, polycomb_scaffolds, "Polycomb body", rna_content, nucleosome_content) 

# compute scores
nucleolus_gc_scores <- getScores(nucleolus_gc_distances_multi)
nucleolus_dfc_scores <- getScores(nucleolus_dfc_distances_multi)
txn_scores <- getScores(txn_distances_multi)
heterochromatin_scores <- getScores(heterochromatin_distances_multi)
polycomb_scores <- getScores(polycomb_distances_multi)

# gather data about condensates reconstituted in vitro
refs <- subset(nucleolus_gc_scaffolds[nucleolus_gc_scaffolds$external_gene_name=="Npm1",], select=c("uniprotswissprot","external_gene_name","protein_sequence","Size_AlphaFold","Size_relaxed","Size_expanded"))
refs <- rbind.data.frame(refs,subset(heterochromatin_scaffolds[heterochromatin_scaffolds$external_gene_name=="Cbx5",], select=c("uniprotswissprot","external_gene_name","protein_sequence","Size_AlphaFold","Size_relaxed","Size_expanded")))
reference_scaffolds <- getReferenceInfo(refs)

  
############################
### display some results ###
############################

# display general information about MLOs (number, volume)
print.data.frame(subset(MLO_info, select=c(Organelle,Volume,Number,TotalVolume)), row.names=FALSE)

# display nucleosomal content of MLOs
print.data.frame(subset(nucleosome_content, select=c(Organelle,Nucleosome_number,TotalNucleosome_number)), row.names=FALSE)

# display RNA content of MLOs
print.data.frame(subset(rna_content, select=c(Organelle,RNA_number,TotalRNA_number)), row.names=FALSE)

# display scaffold proteins in MLOs
print.data.frame(subset(nucleolus_gc_scaffolds, select=c(uniprotswissprot,external_gene_name,Number_per_cell,Enrichment_MLO,Number_per_MLO,Size_AlphaFold,Size_relaxed,Size_expanded)), row.names=FALSE)
print.data.frame(subset(nucleolus_dfc_scaffolds, select=c(uniprotswissprot,external_gene_name,Number_per_cell,Enrichment_MLO,Number_per_MLO,Size_AlphaFold,Size_relaxed,Size_expanded)), row.names=FALSE)
print.data.frame(subset(txn_scaffolds, select=c(uniprotswissprot,external_gene_name,Number_per_cell,Enrichment_MLO,Number_per_MLO,Size_AlphaFold,Size_relaxed,Size_expanded)), row.names=FALSE)
print.data.frame(subset(heterochromatin_scaffolds, select=c(uniprotswissprot,external_gene_name,Number_per_cell,Enrichment_MLO,Number_per_MLO,Size_AlphaFold,Size_relaxed,Size_expanded)), row.names=FALSE)
print.data.frame(subset(polycomb_scaffolds, select=c(uniprotswissprot,external_gene_name,Number_per_cell,Enrichment_MLO,Number_per_MLO,Size_AlphaFold,Size_relaxed,Size_expanded)), row.names=FALSE)

# display intermolecular distances
print(nucleolus_gc_distances_multi$`Distance, all proteins and RNAs`)
print(nucleolus_dfc_distances_multi$`Distance, all proteins and RNAs`)
print(txn_distances_multi$`Distance, all proteins and RNAs`)
print(heterochromatin_distances_multi$`Distance, all proteins and RNAs`)
print(polycomb_distances_multi$`Distance, all proteins and RNAs`)

# display scores
disp_scores <- rbind.data.frame(nucleolus_gc_scores,nucleolus_dfc_scores,txn_scores,heterochromatin_scores,polycomb_scores)
rownames(disp_scores) <- c("Nucleolus, GC","Nucleolus, DFC","Txn condensate, small","Heterochromatin focus","Polycomb body")
colnames(disp_scores) <- c("S_Dist,p","S_Dist,pr","S_Dist,prn")
print.data.frame(disp_scores)
