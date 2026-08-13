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

make_migration <- function(K, Q, N, pi) {
  list(
    K = K,
    Q = Q,
    N = N,
    pi = pi
  )
}


# PROPAGATE LIKELIHOOD UP A BRANCH

propagate_branch <- function(L_child, delta, Q) {
  
  T <- expm(Q * delta)
  
  as.numeric(T %*% L_child)
  
}


# BACKWARD PASS

backward_pass <- function(tree, mig) {
  
  K <- mig$K
  Q <- mig$Q
  N <- mig$N
  
  children <- split(tree$edges[,2], tree$edges[,1])
  
  L <- list()
  
  
  # recursive postorder traversal
  
  recurse <- function(node) {
    
    
    # leaf
    
    if (node %in% names(tree$leaf_states)) {
      
      v <- rep(0, K)
      
      v[tree$leaf_states[[node]]] <- 1
      
      L[[node]] <<- v
      
      return(v)
    }
    
    
    # internal node
    
    kids <- children[[node]]
    
    
    left_child <- recurse(kids[1])
    right_child <- recurse(kids[2])
    
    
    left <- propagate_branch(
      left_child,
      abs(tree$times[node] - tree$times[kids[1]]),
      Q
    )
    
    right <- propagate_branch(
      right_child,
      abs(tree$times[node] - tree$times[kids[2]]),
      Q
    )
    
    
    node_L <- numeric(K)
    
    
    for (k in 1:K) {
      
      node_L[k] <-
        (left[k] * right[k]) /
        (2 * N[k])
      
    }
    
    
    L[[node]] <<- node_L
    
    
    node_L
    
  }
  
  
  recurse(tree$root)
  
  L
}



# FORWARD PASS (CURRENT PLACEHOLDER)
# This remains the original migration propagation.
# It will need redesign later for lineage/pair states.

forward_pass <- function(tree, mig, time_grid) {
  
  Q <- mig$Q
  K <- nrow(Q)
  
  pi_current <- mig$pi
  
  pi_t <- list()
  
  pi_t[[1]] <- pi_current
  
  
  dt <- diff(time_grid)
  
  
  for (t in 2:length(time_grid)) {
    
    Tmat <- expm(Q * dt[t-1])
    
    pi_current <- as.numeric(pi_current %*% Tmat)
    
    pi_t[[t]] <- pi_current
    
  }
  
  
  pi_t
  
}



# CURRENT LAMBDA PLACEHOLDER

compute_lambda_t <- function(tree, pi_t, mig, time_grid) {
  
  
  lambda_t <- numeric(length(time_grid))
  
  
  for (i in seq_along(time_grid)) {
    
    lambda_t[i] <- 0
    
  }
  
  
  lambda_t
  
}



# INTEGRATION

integrate_lambda <- function(lambda_t, time_grid) {
  
  dt <- diff(time_grid)
  
  sum(lambda_t[-1] * dt)
  
}



# LIKELIHOOD

compute_likelihood <- function(L, lambda_int, root, pi) {
  
  L_peel <- sum(pi * L[[root]])
  
  log(L_peel + 1e-12) - lambda_int
  
}



# FULL PIPELINE

run_glike_continuous <- function(tree, mig, time_grid) {
  
  
  L <- backward_pass(tree, mig)
  
  pi_t <- forward_pass(tree, mig, time_grid)
  
  lambda_t <- compute_lambda_t(
    tree,
    pi_t,
    mig,
    time_grid
  )
  
  
  lambda_int <- integrate_lambda(
    lambda_t,
    time_grid
  )
  
  
  logL <- compute_likelihood(
    L,
    lambda_int,
    tree$root,
    mig$pi
  )
  
  
  list(
    log_likelihood = logL,
    L = L,
    pi_t = pi_t,
    lambda_t = lambda_t,
    lambda_integral = lambda_int
  )
  
}



# TOY EXAMPLE


K <- 2


Q <- matrix(
  c(
    -0.2, 0.2,
    0.2,-0.2
  ),
  2,
  2,
  byrow = TRUE
)


mig <- make_migration(
  K = K,
  Q = Q,
  N = c(1000,1000),
  pi = c(0.5,0.5)
)



edges <- matrix(
  c(
    "root","a",
    "root","b"
  ),
  ncol = 2,
  byrow = TRUE
)


times <- c(
  root = 2,
  a = 0,
  b = 0
)


leaf_states <- list(
  a = 1,
  b = 2
)


tree <- create_tree(
  edges,
  times,
  leaf_states,
  "root"
)


time_grid <- seq(
  0,
  2,
  length.out = 50
)



# RUN

res <- run_glike_continuous(
  tree,
  mig,
  time_grid
)


res$L[["root"]]

res$log_likelihood
