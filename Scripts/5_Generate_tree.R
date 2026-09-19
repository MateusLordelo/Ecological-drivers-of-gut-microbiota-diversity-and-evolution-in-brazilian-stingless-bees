library(phyloseq)
library(DECIPHER)
library(phangorn)
library(Biostrings)
library(ape)

# 1. Load the phyloseq object
setwd("C:/Users/mateu/OneDrive/Area_de_Trabalho/2025-12_Gustavo_UFV")
ps <- readRDS("./ASV_tables/For_use/my_phyloseq_object.rds")

# 2. Extract DNA sequences directly from the ps object
dna_sequences <- refseq(ps)

# 3. Align sequences
alignment <- DECIPHER::AlignSeqs(dna_sequences, anchor = NA)

# 4. Construct Phylogenetic Tree
phang.align <- phyDat(as.matrix(alignment), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

# 5. Optimize with Maximum Likelihood (GTR model)
fit <- pml(treeNJ, data = phang.align)
fitGTR <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, 
                    rearrangement = "stochastic", control = pml.control(trace = 0))
final_tree <- fitGTR$tree
# 5.5 Root the tree at the midpoint
final_tree <- phangorn::midpoint(final_tree)

# 6. Merge tree back into the phyloseq object
ps_final <- merge_phyloseq(ps, final_tree)

# 7. Export files
outdir <- "./ASV_tables/For_use" # Ensure this directory exists
if(!dir.exists(outdir)) dir.create(outdir)

# Save the Newick tree file
write.tree(final_tree, file.path(outdir, "asv_tree.nwk"))

# Save the native R tree object
saveRDS(final_tree, file.path(outdir, "asv_tree.rds"))

# Save the updated phyloseq object (now containing OTU, Taxa, Samples, Seqs, and Tree)
saveRDS(ps_final, file.path(outdir, "ps_final.rds"))

cat("All files exported successfully to:", outdir, "\n")
