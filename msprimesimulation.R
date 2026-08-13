library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/Users/rohannahasselkus/.virtualenvs/msprime-env/bin/python")
msprime <- import("msprime")
tskit <- import("tskit")

N1 <- 10000
N2 <- 10000
m  <- 1e-4

demography <- msprime$Demography()
demography$add_population(name = "pop1", initial_size = N1)
demography$add_population(name = "pop2", initial_size = N2)
demography$set_symmetric_migration_rate(
  populations = list("pop1", "pop2"),
  rate = m
)

ts <- msprime$sim_ancestry(
  samples = list(pop1 = 2L, pop2 = 1L),
  demography = demography,
  ploidy = 1L,
  sequence_length = 1,
  recombination_rate = 0,
  record_migrations = TRUE,
  random_seed = 42L
)

tree <- ts$first()
cat(tree$draw_text(), "\n")

n_nodes <- ts$num_nodes
for (i in 0:(n_nodes - 1)) {
  node <- ts$node(i)
  cat(sprintf("Node %d | time = %.4f | population = %d | is_sample = %s\n",
              i, node$time, node$population, node$is_sample()))
}

internal_nodes <- setdiff(0:(n_nodes - 1), 0:(ts$num_samples - 1))
coal_times <- sort(sapply(internal_nodes, function(u) tree$time(u)))
cat("\nCoalescence times:", coal_times, "\n")

for (i in 0:(ts$num_samples - 1)) {
  cat(sprintf("Leaf %d -> population %d\n", i, ts$node(i)$population))
}

cat("\nNewick:", tree$as_newick(), "\n")

# --- Migration history ---
n_mig <- ts$num_migrations
cat("\nNumber of migration events:", n_mig, "\n")

if (n_mig > 0) {
  cat("\nMigration table:\n")
  cat(sprintf("%-6s %-6s %-8s %-8s %-10s %-10s\n",
              "Node", "Left", "Right", "Time", "From", "To"))
  for (i in 0:(n_mig - 1)) {
    mig <- ts$migration(i)
    cat(sprintf("%-6d %-6.2f %-8.2f %-8.4f %-10d %-10d\n",
                mig$node, mig$left, mig$right, mig$time,
                mig$source, mig$dest))
  }
} else {
  cat("\nNo migration events occurred in this simulation.\n")
}
