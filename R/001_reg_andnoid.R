library(tidyverse)
library(here)
library(rstan)

set.seed(3082)
set.seed(1450198)

big_run <- TRUE
if (big_run) {
    my_cores <- 4
    my_chains <- 4
    my_iter <- 2000
    my_warmup <- 1000
    my_thin <- 1
} else {
    my_cores <- 1
    my_chains <- 1
    my_iter <- 40
    my_warmup <- 20
    my_thin <- 1
}


here::i_am("R/001_reg_andnoid.R")

### Read in the post-stratification frame
psf <- readRDS(here::here("data",
                          "ward_joint.rds"))

### Check these are all factors
psf <- psf |>
    mutate_at(vars(c(geogcode, Age, Sex, Quals, Ethnicity,
                     Quals_by_ethnicity,
                     Quals_by_sex)),
              factor)

### Read in the individual data
ind <- readRDS(here::here("working",
                          "tidy_ind_data.rds"))


### Read in the aggregate-level predictors
aux <- readRDS(here::here("working", "tidy_aux_data.rds"))


### Merge on to both sets of data
ind <- left_join(ind, aux, by = "geogcode")
psf <- left_join(psf, aux, by = "geogcode")

### Make sure everything is a data frame, not a tibble
ind <- as.data.frame(ind)
psf <- as.data.frame(psf)
aux <- as.data.frame(aux)

### Explicitly factor ind$geogcode using psf$geogcode levels
psf$geogcode <- factor(psf$geogcode)
ind$geogcode <- factor(ind$geogcode,
                       levels = levels(psf$geogcode))
ind <- ind |>
    subset(!is.na(geogcode))

### And the same for quals_by_ethnicity
psf$Quals_by_ethnicity <- factor(psf$Quals_by_ethnicity)
ind$Quals_by_ethnicity <- factor(ind$Quals_by_ethnicity,
                       levels = levels(psf$Quals_by_ethnicity))

### Make sure at the variable levels are the same
indvars <- c("Sex", "Age", "Quals", "Ethnicity",
             "Quals_by_sex", "LAD23NM")

for (v in indvars) {
    if (!is.element(v, names(ind))) {
        stop(paste0(v, " not in ind data"))
    }
    if (!is.element(v, names(psf))) {
        stop(paste0(v, " not in post-strat data"))
    }
    ind[,v] <- factor(ind[,v],
                          levels = sort(unique(as.character(ind[,v]))))
    psf[,v] <- factor(psf[,v],
                      levels = sort(unique(as.character(psf[,v]))))
    ind_levels <- levels(ind[,v])
    ps_levels <- levels(psf[,v])
    if (!all(ind_levels == ps_levels)) {
        stop(paste0("Factor levels different for ", v))
    }
}


### Make sure the auxiliary data frame is organised in the same way as
### the const. var
aux.bak <- aux
aux <- aux |>
    mutate(geogcode = factor(as.character(geogcode),
                             levels = levels(psf$geogcode))) |>
    arrange(geogcode) |>
    dplyr::select(pct_18_to_19, pct_20_to_24,
                  cob_uk, cob_eu,
                  pct_60plus, pct_noquals, ## pct_l4plus,
                  pct_nopass,
                  pct_black, pct_asian, turn) |>
    mutate_all(coalesce, 0.0) |>
    mutate_all(scale)



stan_data <- list(N = nrow(ind),
                  P = nrow(psf),
                  D = ncol(aux),
                  no_id = ind$no_id,
                  reg = ind$reg,
                  X = aux,
                  J_age = max(as.numeric(factor(ind$Age))),
                  J_sex = max(as.numeric(factor(ind$Sex))),
                  J_quals = max(as.numeric(factor(ind$Quals))),
                  J_ethnicity = max(as.numeric(factor(ind$Ethnicity))),
                  J_quals_by_ethnicity = max(as.numeric(psf$Quals_by_ethnicity)),
                  J_quals_by_sex = max(as.numeric(factor(ind$Quals_by_sex))),
                  J_consty = max(as.numeric(factor(psf$geogcode))),
                  J_lad = max(as.numeric(factor(psf$LAD23NM))),
                  G = 7, ## Need to check this later
                  idx_age = as.numeric(factor(ind$Age)),
                  idx_sex = as.numeric(factor(ind$Sex)),
                  idx_quals = as.numeric(factor(ind$Quals)),
                  idx_ethnicity = as.numeric(factor(ind$Ethnicity)),
                  idx_quals_by_ethnicity = as.numeric(ind$Quals_by_ethnicity),
                  idx_quals_by_sex = as.numeric(factor(ind$Quals_by_sex)),
                  idx_consty = as.numeric(ind$geogcode),
                  idx_lad = as.numeric(factor(ind$LAD23NM)),
                  psidx_age = as.numeric(factor(psf$Age)),
                  psidx_sex = as.numeric(factor(psf$Sex)),
                  psidx_quals = as.numeric(factor(psf$Quals)),
                  psidx_ethnicity = as.numeric(factor(psf$Ethnicity)),
                  psidx_quals_by_ethnicity = as.numeric(factor(psf$Quals_by_ethnicity)),
                  psidx_quals_by_sex = as.numeric(factor(psf$Quals_by_sex)),
                  psidx_consty = as.numeric(factor(psf$geogcode)),
                  psidx_lad = as.numeric(factor(psf$LAD23NM)),
                  psfweight = psf$count)

### Remove the psf from memory (it's in the stan data object)
consty_levels <- levels(factor(psf$geogcode))
### rm(psf)



stan_code <- "

data {
  int<lower=0> N;
  int<lower=0> P; // ?
  int<lower=1> D; // number of area-level predictors
  int<lower=0,upper=1> no_id[N];
  int<lower=0,upper=1> reg[N];
  int<lower=1> J_age; // number of factor levels for each
  int<lower=1> J_sex;
  int<lower=1> J_quals;
  int<lower=1> J_ethnicity;
  int<lower=1> J_quals_by_ethnicity;
  int<lower=1> J_quals_by_sex;
  int<lower=1> J_consty;
  int<lower=1> J_lad;
  int<lower=1> G;   // Number of random effects
  matrix[J_consty, D] X;  // area-level covariates
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1, upper=J_sex> idx_sex[N];
  int<lower=1, upper=J_ethnicity> idx_ethnicity[N];
  int<lower=1, upper=J_quals> idx_quals[N];
  int<lower=1, upper=J_quals_by_ethnicity> idx_quals_by_ethnicity[N];
  int<lower=1, upper=J_quals_by_sex> idx_quals_by_sex[N];
  int<lower=1, upper=J_consty> idx_consty[N];
  int<lower=1, upper=J_lad> idx_lad[N];
  int<lower=1, upper=J_age> psidx_age[P];
  int<lower=1, upper=J_sex> psidx_sex[P];
  int<lower=1, upper=J_quals> psidx_quals[P];
  int<lower=1, upper=J_ethnicity> psidx_ethnicity[P];
  int<lower=1, upper=J_quals_by_ethnicity> psidx_quals_by_ethnicity[P];
  int<lower=1, upper=J_quals_by_sex> psidx_quals_by_sex[P];
  int<lower=1, upper=J_consty> psidx_consty[P];
  int<lower=1, upper=J_lad> psidx_lad[P];
  int<lower=1> psfweight[P];

}
transformed data {
  matrix[J_consty, D] Q_ast;
  matrix[D, D] R_ast;
  matrix[D, D] R_ast_inverse;
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(J_consty - 1);
  R_ast = qr_thin_R(X) / sqrt(J_consty - 1);
  R_ast_inverse = inverse(R_ast);
}
parameters {
  vector[D] beta_noid;
  real alpha_noid_raw;
  real alpha_reg_raw;
  vector[J_age] eta_noid_age;
  vector[J_sex] eta_noid_sex;
  vector[J_quals] eta_noid_quals;
  vector[J_ethnicity] eta_noid_ethnicity;
  vector[J_quals_by_ethnicity] eta_noid_quals_by_ethnicity;
  vector[J_quals_by_sex] eta_noid_quals_by_sex;
  vector[J_consty] eta_noid_consty;
  vector[J_lad] eta_noid_lad;
  real<lower=0> sigma_noid_age;
  real<lower=0> sigma_noid_sex;
  real<lower=0> sigma_noid_quals;
  real<lower=0> sigma_noid_ethnicity;
  real<lower=0> sigma_noid_quals_by_ethnicity;
  real<lower=0> sigma_noid_quals_by_sex;
  real<lower=0> sigma_noid_consty;
  real<lower=0> sigma_noid_lad;
  vector<lower=0>[G] sigma_noid_inter_eqn;

  vector[D] beta_reg;
  vector[J_age] eta_reg_age;
  vector[J_sex] eta_reg_sex;
  vector[J_quals] eta_reg_quals;
  vector[J_ethnicity] eta_reg_ethnicity;
  vector[J_quals_by_ethnicity] eta_reg_quals_by_ethnicity;
  vector[J_quals_by_sex] eta_reg_quals_by_sex;
  vector[J_consty] eta_reg_consty;
  vector[J_lad] eta_reg_lad;
  real<lower=0> sigma_reg_age;
  real<lower=0> sigma_reg_sex;
  real<lower=0> sigma_reg_quals;
  real<lower=0> sigma_reg_ethnicity;
  real<lower=0> sigma_reg_quals_by_ethnicity;
  real<lower=0> sigma_reg_quals_by_sex;
  real<lower=0> sigma_reg_consty;
  real<lower=0> sigma_reg_lad;
  vector<lower=0>[G] sigma_reg_inter_eqn;

  real beta_has_id;

}
transformed parameters {
  vector[J_age] alpha_noid_age;
  vector[J_sex] alpha_noid_sex;
  vector[J_quals] alpha_noid_quals;
  vector[J_ethnicity] alpha_noid_ethnicity;
  vector[J_quals_by_ethnicity] alpha_noid_quals_by_ethnicity;
  vector[J_quals_by_sex] alpha_noid_quals_by_sex;
  vector[J_consty] alpha_noid_consty;
  vector[J_lad] alpha_noid_lad;

  vector[J_age] alpha_reg_age;
  vector[J_sex] alpha_reg_sex;
  vector[J_quals] alpha_reg_quals;
  vector[J_ethnicity] alpha_reg_ethnicity;
  vector[J_quals_by_ethnicity] alpha_reg_quals_by_ethnicity;
  vector[J_quals_by_sex] alpha_reg_quals_by_sex;
  vector[J_consty] alpha_reg_consty;
  vector[J_lad] alpha_reg_lad;


  alpha_noid_age = sigma_noid_age * eta_noid_age;
  alpha_noid_sex = sigma_noid_inter_eqn[1] * sigma_noid_sex * eta_noid_sex;
  alpha_noid_quals = sigma_noid_inter_eqn[2] * sigma_noid_quals * eta_noid_quals;
  alpha_noid_ethnicity = sigma_noid_inter_eqn[3] * sigma_noid_ethnicity * eta_noid_ethnicity;
  alpha_noid_quals_by_ethnicity = sigma_noid_inter_eqn[4] * sigma_noid_quals_by_ethnicity * eta_noid_quals_by_ethnicity;
  alpha_noid_quals_by_sex = sigma_noid_inter_eqn[5] * sigma_noid_quals_by_sex * eta_noid_quals_by_sex;
  alpha_noid_consty = sigma_noid_inter_eqn[6] * sigma_noid_consty * eta_noid_consty;
  alpha_noid_lad = sigma_noid_inter_eqn[7] * sigma_noid_lad * eta_noid_lad;

  alpha_reg_age = sigma_reg_age * eta_reg_age;
  alpha_reg_sex = sigma_reg_inter_eqn[1] * sigma_reg_sex * eta_reg_sex;
  alpha_reg_quals = sigma_reg_inter_eqn[2] * sigma_reg_quals * eta_reg_quals;
  alpha_reg_ethnicity = sigma_reg_inter_eqn[3] * sigma_reg_ethnicity * eta_reg_ethnicity;
  alpha_reg_quals_by_ethnicity = sigma_reg_inter_eqn[4] * sigma_reg_quals_by_ethnicity * eta_reg_quals_by_ethnicity;
  alpha_reg_quals_by_sex = sigma_reg_inter_eqn[5] * sigma_reg_quals_by_sex * eta_reg_quals_by_sex;
  alpha_reg_consty = sigma_reg_inter_eqn[6] * sigma_reg_consty * eta_reg_consty;
  alpha_reg_lad = sigma_reg_inter_eqn[7] * sigma_reg_lad * eta_reg_lad;

}
model {
  to_vector(beta_noid) ~ std_normal();
  for (i in 2:J_age) {
    eta_noid_age[i] ~ normal(eta_noid_age[i - 1], 1);
  }
  sum(eta_noid_age) ~ normal(0, 0.01 * J_age); // constraint so we can write likelihood for rw(1).
  to_vector(eta_noid_sex) ~ std_normal();
  to_vector(eta_noid_quals) ~ std_normal();
  to_vector(eta_noid_ethnicity) ~ std_normal();
  to_vector(eta_noid_quals_by_ethnicity) ~ std_normal();
  to_vector(eta_noid_quals_by_sex) ~ std_normal();
  to_vector(eta_noid_consty) ~ std_normal();
  to_vector(eta_noid_lad) ~ std_normal();
  sigma_noid_age ~ std_normal();
  sigma_noid_sex ~ std_normal();
  sigma_noid_quals ~ std_normal();
  sigma_noid_ethnicity ~ std_normal();
  sigma_noid_quals_by_ethnicity ~ std_normal();
  sigma_noid_quals_by_sex ~ std_normal();
  sigma_noid_consty ~ std_normal();
  sigma_noid_lad ~ std_normal();
  alpha_noid_raw ~ std_normal();
  sigma_noid_inter_eqn ~ std_normal();

  to_vector(beta_reg) ~ std_normal();
  for (i in 2:J_age) {
    eta_reg_age[i] ~ normal(eta_reg_age[i - 1], 1);
  }
  sum(eta_reg_age) ~ normal(0, 0.01 * J_age);
  to_vector(eta_reg_sex) ~ std_normal();
  to_vector(eta_reg_quals) ~ std_normal();
  to_vector(eta_reg_ethnicity) ~ std_normal();
  to_vector(eta_reg_quals_by_ethnicity) ~ std_normal();
  to_vector(eta_reg_quals_by_sex) ~ std_normal();
  to_vector(eta_reg_consty) ~ std_normal();
  to_vector(eta_reg_lad) ~ std_normal();
  sigma_reg_age ~ std_normal();
  sigma_reg_sex ~ std_normal();
  sigma_reg_quals ~ std_normal();
  sigma_reg_ethnicity ~ std_normal();
  sigma_reg_quals_by_ethnicity ~ std_normal();
  sigma_reg_quals_by_sex ~ std_normal();
  sigma_reg_consty ~ std_normal();
  sigma_reg_lad ~ std_normal();
  sigma_reg_inter_eqn ~ std_normal();
  beta_has_id ~ std_normal();

  vector[N] mu_reg = rep_vector(0.0, N);
  vector[N] mu_noid = rep_vector(0.0, N);
  vector[J_consty] consty_part_reg;
  vector[J_consty] consty_part_noid;
  consty_part_noid = Q_ast * beta_noid;
  consty_part_reg = Q_ast * beta_reg;

  for (n in 1:N) {
    
    mu_noid[n] = consty_part_noid[idx_consty[n]] +
  + alpha_noid_raw + alpha_noid_age[idx_age[n]]
  + alpha_noid_sex[idx_sex[n]]
  + alpha_noid_quals[idx_quals[n]]
  + alpha_noid_ethnicity[idx_ethnicity[n]]
  + alpha_noid_quals_by_ethnicity[idx_quals_by_ethnicity[n]]
  + alpha_noid_quals_by_sex[idx_quals_by_sex[n]]
  + alpha_noid_lad[idx_lad[n]]
  + alpha_noid_consty[idx_consty[n]];

    mu_reg[n] = consty_part_reg[idx_consty[n]]
  + alpha_reg_raw + beta_has_id * (1 - no_id[n])
  + alpha_reg_age[idx_age[n]]
  + alpha_reg_sex[idx_sex[n]]
  + alpha_reg_quals[idx_quals[n]]
  + alpha_reg_ethnicity[idx_ethnicity[n]]
  + alpha_reg_quals_by_ethnicity[idx_quals_by_ethnicity[n]]
  + alpha_reg_quals_by_sex[idx_quals_by_sex[n]]
  + alpha_reg_lad[idx_lad[n]]
  + alpha_reg_consty[idx_consty[n]];

  }

  no_id ~ bernoulli_logit(mu_noid);
  reg ~ bernoulli_logit(mu_reg);
}
generated quantities {
  int consty_reg_with_id[J_consty];
  int consty_reg_no_id[J_consty];
  int consty_noid[J_consty];


  for (n in 1:J_consty) {
    consty_noid[n] = 0;
    consty_reg_with_id[n] = 0;
    consty_reg_no_id[n] = 0;
  }

  ## begin a big local block not saved
  { 
  int psf_reg_with_id[P];
  int psf_reg_no_id[P];
  int psf_noid[P];

  vector[P] pmu_noid = rep_vector(0.0, P);
  vector[P] pmu_reg_with_id = rep_vector(0.0, P);
  vector[P] pmu_reg_no_id = rep_vector(0.0, P);
  vector[J_consty] consty_part_reg;
  vector[J_consty] consty_part_noid;
  consty_part_noid = Q_ast * beta_noid;
  consty_part_reg = Q_ast * beta_reg;

  for (p in 1:P) {
    pmu_noid[p] = consty_part_noid[psidx_consty[p]] +
  + alpha_noid_raw + alpha_noid_age[psidx_age[p]]
  + alpha_noid_sex[psidx_sex[p]]
  + alpha_noid_quals[psidx_quals[p]]
  + alpha_noid_ethnicity[psidx_ethnicity[p]]
  + alpha_noid_quals_by_ethnicity[psidx_quals_by_ethnicity[p]]
  + alpha_noid_quals_by_sex[psidx_quals_by_sex[p]]
  + alpha_noid_lad[psidx_lad[p]]
  + alpha_noid_consty[psidx_consty[p]];

    pmu_reg_no_id[p] = consty_part_reg[psidx_consty[p]]
  + alpha_reg_raw
  + alpha_reg_age[psidx_age[p]]
  + alpha_reg_sex[psidx_sex[p]]
  + alpha_reg_quals[psidx_quals[p]]
  + alpha_reg_ethnicity[psidx_ethnicity[p]]
  + alpha_reg_quals_by_ethnicity[psidx_quals_by_ethnicity[p]]
  + alpha_reg_quals_by_sex[psidx_quals_by_sex[p]]
  + alpha_reg_lad[psidx_lad[p]]
  + alpha_reg_consty[psidx_consty[p]];

    pmu_reg_with_id[p] = pmu_reg_no_id[p] + beta_has_id;

  psf_noid[p] = binomial_rng(psfweight[p], inv_logit(pmu_noid[p]));
  psf_reg_with_id[p] = binomial_rng(psfweight[p] - psf_noid[p],
    inv_logit(pmu_reg_with_id[p]));
  psf_reg_no_id[p] = binomial_rng(psf_noid[p],
    inv_logit(pmu_reg_no_id[p]));
  
  }

  for (p in 1:P) {
      consty_noid[psidx_consty[p]] = consty_noid[psidx_consty[p]] +
                               psf_noid[p];
      consty_reg_with_id[psidx_consty[p]] = consty_reg_with_id[psidx_consty[p]] +
                               psf_reg_with_id[p];
      consty_reg_no_id[psidx_consty[p]] = consty_reg_no_id[psidx_consty[p]] +
                               psf_reg_no_id[p];

  }
  }
}

"

samples <- stan(model_code = stan_code,
                data = stan_data,
                iter = my_iter,
                warmup = my_warmup,
                thin = my_thin,
                core = my_cores,
                chain = my_chains,
                init = 0,
                control = list(adapt_delta = 0.95,
                               max_treedepth = 12),
                pars = c("alpha_noid_raw",
                         "alpha_reg_raw",
                         "beta_noid",
                         "beta_has_id",
                         "consty_reg_with_id", "consty_reg_no_id",
                         "consty_noid"
                         ))


vi_levels <- levels(factor(ind$vi))

outfile <- paste0("mrp_samples_bivreg_",
                  Sys.Date(),
                  ".RData")

save(samples, vi_levels, consty_levels, aux, ind,
     file = here::here("working",
                       outfile))

## ### Post-process seat counts
## counts_vote_and_noid <- rstan::extract(samples,
##                                  pars = "consty_vote_and_noid",
##                                  permuted = TRUE)[[1]]

## counts_vote <- rstan::extract(samples,
##                                  pars = "consty_vote",
##                                  permuted = TRUE)[[1]]

## ### What does this imply for turnout?
## pop_totals <- psf |>
##     group_by(geogcode) |>
##     summarize(tot = sum(count),
##               .groups = "drop")
## pop_totals <- pop_totals$tot[match(consty_levels, pop_totals$geogcode)]

## pop_totals <- matrix(pop_totals,
##                      nrow = nrow(counts_vote),
##                      ncol = length(pop_totals),
##                      byrow = TRUE)

## foo <- colMeans(counts_vote / pop_totals)
