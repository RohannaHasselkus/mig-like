
# PACKAGES

library(expm)


# TREE CONSTRUCTION

create_tree <- function(edges, times, leaf_states, root) {
  list(
    edges = edges,
    times = times,
    leaf_states = leaf_states,
    root = root
  )
}


# MIGRATION MODEL

make_migration <- function(K, Q) {
  list(K = K, Q = Q)
}


# BACKWARD PASS (TREE LIKELIHOOD)

backward_pass <- function(tree, mig) {
  
  K <- mig$K
  edges <- tree$edges
  
  children <- split(edges[,2], edges[,1])
  nodes <- unique(as.vector(edges))
  internal <- setdiff(nodes, names(tree$leaf_states))
  
  L <- list()
  
  # leaves
  for (leaf in names(tree$leaf_states)) {
    v <- rep(0, K)
    v[tree$leaf_states[[leaf]]] <- 1
    L[[leaf]] <- as.numeric(v)
  }
  
  # internal nodes
  for (a in internal) {
    
    kids <- children[[a]]
    left <- L[[kids[1]]]
    right <- L[[kids[2]]]
    
    La <- numeric(K)
    
    for (k in 1:K) {
      La[k] <- (left[k] * right[k]) / (2 * 1000)
    }
    
    L[[a]] <- as.numeric(La)
  }
  
  L
}


# FORWARD PASS (TIME-DEPENDENT MIGRATION)

forward_pass <- function(tree, mig, time_grid) {
  
  Q <- mig$Q
  K <- nrow(Q)
  
  pi_t <- list()
  
  pi_current <- rep(1 / K, K)
  pi_t[[1]] <- as.numeric(pi_current)
  
  dt <- diff(time_grid)
  
  for (t in 2:length(time_grid)) {
    
    T <- expm(Q * dt[t - 1])
    pi_current <- as.numeric(pi_current %*% T)
    
    pi_t[[t]] <- as.numeric(pi_current)
  }
  
  pi_t
}


# COALESCENT RATE λ(t)

compute_lambda_t <- function(tree, pi_t, mig, time_grid) {
  
  edges <- tree$edges
  
  lambda_t <- numeric(length(time_grid))
  
  for (t_idx in seq_along(time_grid)) {
    
    lambda <- 0
    
    pi_vec <- as.numeric(pi_t[[t_idx]])
    
    for (i in 1:nrow(edges)) {
      
      for (k in 1:length(pi_vec)) {
        lambda <- lambda + pi_vec[k] * pi_vec[k]
      }
    }
    
    lambda_t[t_idx] <- lambda / (2 * 1000)
  }
  
  lambda_t
}


# INTEGRATION

integrate_lambda <- function(lambda_t, time_grid) {
  dt <- diff(time_grid)
  sum(lambda_t[-1] * dt)
}


# LIKELIHOOD

compute_likelihood <- function(L, lambda_int, root) {
  log(sum(as.numeric(L[[root]])) + 1e-12) - lambda_int
}


# FULL PIPELINE

run_glike_continuous <- function(tree, mig, time_grid) {
  
  L <- backward_pass(tree, mig)
  pi_t <- forward_pass(tree, mig, time_grid)
  lambda_t <- compute_lambda_t(tree, pi_t, mig, time_grid)
  lambda_int <- integrate_lambda(lambda_t, time_grid)
  
  logL <- compute_likelihood(L, lambda_int, tree$root)
  
  list(
    log_likelihood = logL,
    L = L,
    pi_t = pi_t,
    lambda_t = lambda_t,
    lambda_integral = lambda_int
  )
}


# DIAGNOSTICS

run_diagnostics <- function(tree, mig, time_grid) {
  
  mig2 <- mig
  mig2$Q <- mig$Q * 10
  
  pi1 <- forward_pass(tree, mig, time_grid)
  pi2 <- forward_pass(tree, mig2, time_grid)
  
  diff_pi <- 0
  
  for (t in seq_along(time_grid)) {
    diff_pi <- diff_pi + sum(abs(pi1[[t]] - pi2[[t]]))
  }
  
  cat("PI sensitivity:", diff_pi, "\n")
  
  L <- backward_pass(tree, mig)
  pi <- forward_pass(tree, mig, time_grid)
  
  lambda <- compute_lambda_t(tree, pi, mig, time_grid)
  
  tree_term <- log(sum(L[[tree$root]]) + 1e-12)
  lambda_term <- integrate_lambda(lambda, time_grid)
  
  cat("Tree term:", tree_term, "\n")
  cat("Lambda term:", lambda_term, "\n")
  cat("Ratio:", abs(tree_term / (lambda_term + 1e-12)), "\n")
}


# TOY EXAMPLE


K <- 2

Q <- matrix(c(
  -0.2, 0.2,
  0.2, -0.2
), 2, 2, byrow = TRUE)

mig <- make_migration(K, Q)

edges <- matrix(c(
  "root", "a",
  "root", "b"
), byrow = TRUE, ncol = 2)

times <- c(root = 2, a = 0, b = 0)

leaf_states <- list(a = 1, b = 2)

tree <- create_tree(edges, times, leaf_states, "root")

time_grid <- seq(0, 2, length.out = 50)

# RUN MODEL
res <- run_glike_continuous(tree, mig, time_grid)

res$log_likelihood

# RUN DIAGNOSTICS
run_diagnostics(tree, mig, time_grid)
