#' BITFAM main function. BITFAM will infer the transcription factor activities
#' from single cell RNA-seq data based on the ChIP-seq data
#'
#' @param data A matrix [M, N] of normalized single cell RNA-seq data.
#'        The rows are genes and the columns are cells
#'        If Seurat object, access by `data <- GetAssayData(cells)`
#'        Recommended to subset to top 2000-5000 variable genes.
#' @param network A binary network matrix [M, K] of prior knowledge.
#'        The rows are genes and the columns are TFs.
#'        Recommended to subset to TFs present in the data and
#'        filter TFs with few (e.g., < 10) target genes.
#' @param number of CPU cores. Default is 8
#' @param number of max iteration. Default is 8000
#' @param convergence tolerance on the relative norm of the objective.
#'        Default is 0.005
#' @return sampling results of TF inferred activities and TF-gene weights
#' @export
#' @import rstan
#' @import Seurat

BITFAM <- function(
    data,
    network,
    ncores = 8,
    iter = 8000,
    tol_rel_obj = 0.005) {
  # Check if the column names of data match the row names of network
  if (!all(rownames(data) == rownames(network))) {
    stop("Row names of data must match row names of network")
  }

  X <- t(as.matrix(data))
  N <- dim(X)[1]
  D <- dim(X)[2]
  K <- ncol(network)
  data_to_model <- list(N = N, D = D, K = K, X = X, Mask = network)


  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = ncores)

  set.seed(100)
  pca_beta_piror <- "
data {
int<lower=0> N; // Number of samples
int<lower=0> D; // The original dimension
int<lower=0> K; // The latent dimension
matrix[N, D] X; // The data matrix
matrix[D, K] Mask; // The binary mask of prior knowledge indicate the target of TFs
}

parameters {
matrix<lower=0, upper=1>[N, K] Z; // The latent matrix
matrix[D, K] W; // The weight matrix
real<lower=0> tau; // Noise term
vector<lower=0>[K] alpha; // ARD prior
}

transformed parameters{
matrix<lower=0>[D, K] t_alpha;
real<lower=0> t_tau;
for(wmd in 1:D){
for(wmk in 1:K){
t_alpha[wmd, wmk] = Mask[wmd, wmk] == 1 ? inv(sqrt(alpha[wmk])) : 0.01;
}
}
t_tau = inv(sqrt(tau));
}
model {
tau ~ gamma(1,1);
to_vector(Z) ~ beta(0.5, 0.5);
alpha ~ gamma(1e-3,1e-3);
for(d in 1:D){
for(k in 1:K){
W[d,k] ~ normal(0, t_alpha[d, k]);
}
}
to_vector(X) ~ normal(to_vector(Z*W'), t_tau);
} "

  m_beta_prior <- stan_model(model_code = pca_beta_piror)
  suppressWarnings(
    fit.vb <- vb(
      m_beta_prior,
      data = data_to_model,
      algorithm = "meanfield",
      iter = iter,
      output_samples = 300,
      tol_rel_obj = tol_rel_obj
    )
  )
  BITFAM_list <- list(
    Model = fit.vb,
    TF_used = colnames(network),
    Genes = rownames(data),
    Cell_names = colnames(data)
  )
  return(BITFAM_list)
}
