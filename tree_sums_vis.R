
.packages <- c("plyr", "dplyr", "phyloseq", "treelapse")
sapply(.packages, require, character.only = TRUE)

pregnancy_path <- "http://statweb.stanford.edu/~susan/papers/Pregnancy/PregnancyClosed15.Rdata"
tmp <- tempfile()
download.file(pregnancy_path, tmp)
load(tmp)

PS <- PS %>%
  subset_samples(SubjectID %in% c("10101") &
                 BodySite == "Vaginal_Swab") %>%
  filter_taxa(function(x) sum(x > 1) > 0.01 * length(x), TRUE)
el <- taxa_edgelist(tax_table(PS))

times <- sort(runif(100))
tree_data <- tree_sum(gp_data, as.matrix(el), times = times, bandwidth = 0.1)

mtree_data <- melt(tree_data)
mtree_data$time <- times[mtree_data$Var1]
mtree_data$Var1 <- NULL
mtree_data$Var2 <- NULL

colnames(mtree_data) <- c("value", "unit", "time")
timebox_tree(mtree_data, taxa, "Bacteria", 700, 500)
