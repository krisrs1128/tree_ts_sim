
.packages <- c("plyr", "dplyr", "phyloseq", "treelapse", "simData")
sapply(.packages, require, character.only = TRUE)

process_tree_data <- function(tree_data, times) {
  mtree_data <- melt(tree_data)
  mtree_data$time <- times[mtree_data$Var1]
  mtree_data$Var1 <- NULL
  mtree_data$Var2 <- NULL
  colnames(mtree_data) <- c("value", "unit", "time")
  mtree_data
}

# generate tree data
pregnancy_path <- "http://statweb.stanford.edu/~susan/papers/Pregnancy/PregnancyClosed15.Rdata"
tmp <- tempfile()
download.file(pregnancy_path, tmp)
load(tmp)

PS <- PS %>%
  subset_samples(SubjectID %in% c("10101") &
                 BodySite == "Vaginal_Swab") %>%
  filter_taxa(function(x) sum(x > 1) > 0.01 * length(x), TRUE)
el <- taxa_edgelist(tax_table(PS))

# generate HMM
P <- matrix(runif(4 * 4), 4, 4) + 10 * diag(4)
P <- diag(1 / rowSums(P)) %*% P
obs_densities <- lapply(1:4, function(mu) { function(n) { rnorm(n, mu, 1) }})
names(obs_densities) <- 1:4

tree_data <- tree_sum(
  function(...) { t(hmm_data(...)) },
  as.matrix(el),
  P = P,
  obs_densities = obs_densities
)

timebox_tree(
  process_tree_data(tree_data, 1:100),
  data.frame(el, stringsAsFactors = FALSE),
  "Bacteria",
  700,
  500
)

# generate GP
tree_data <- tree_sum(gp_data, as.matrix(el), times = times, bandwidth = 0.1)
timebox_tree(
  process_tree_data(tree_data, times),
  data.frame(el, stringsAsFactors = FALSE),
  "Bacteria",
  700,
  500
)

# Generate state space model
pdf <- pdf_factory(1)


tree_data <- tree_sum(
  node_fun = function(...) {
    matrix(ssm_sample(...), ncol = 1)
  },
  el = as.matrix(el),
  pdf = pdf
)

timebox_tree(
  process_tree_data(tree_data, 1:100),
  data.frame(el, stringsAsFactors = FALSE),
  "Bacteria",
  700,
  500
)
