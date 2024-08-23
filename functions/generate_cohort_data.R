## autocorrelation function
ar1_cor <- function(n, rho) {
    exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                        (1:n - 1))
    rho^exponent
}

## generate random components
gen_data_comps <- function(n_obs, n_cohorts, expos_mean, expos_SD, expos_rho, prop_cohorts, sigsq_y_true, expos_names_prefix = "expos") {

    ## generate exposure
    cor_mat <- ar1_cor(n_cohorts, expos_rho)
    rownames(cor_mat) <- colnames(cor_mat) <- paste0(expos_names_prefix, 1:n_cohorts)
    expos <- MASS::mvrnorm(n_obs, mu = rep(expos_mean, n_cohorts), Sigma = expos_SD^2*cor_mat)

    ## generate cohort indicator
    ns_cohorts <- round(prop_cohorts*n_obs)
    n_prior_expos <- sample(rep(1:n_cohorts, ns_cohorts))

    ## dataset components
    dset <- tibble(
        study_id = 1:n_obs,
        eps = rnorm(n_obs, 0, sqrt(sigsq_y_true)),
        n_prior_expos = n_prior_expos,
        as_tibble(expos) ## add these to be use in imputation model
    )
    dset
}

## generate a dataset from the random components
gen_data_from_comps <- function(df_comps, DL_coef, expos_names_prefix = "expos") {
    ## Generate outcome data from components
    expos_df <- df_comps %>% select(starts_with(expos_names_prefix))
    expos_names <- colnames(expos_df)
    expos_df_true <- df_comps %>% arrange(study_id)
    expos_mat <- data.matrix(expos_df)
    mu <- as.vector(expos_mat %*% DL_coef)

    ## set up observed exposure data by setting to missing if n_prior_expos is smaller than the lag
    expos_df_obs <- expos_df_true %>%
        pivot_longer(all_of(expos_names), names_to = "lag", values_to = "expos") %>%
        mutate(expos = ifelse(as.numeric(gsub(expos_names_prefix, "", lag)) > n_prior_expos, NA, expos)) %>%
        pivot_wider(names_from = "lag", values_from = "expos") %>%
        arrange(study_id)

    df <- expos_df_obs %>%
        mutate(mu = mu,
               y = mu + eps)
    attr(df, "expos_names") <- expos_names
    df
}
