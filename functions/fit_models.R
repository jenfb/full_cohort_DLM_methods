## Fit unconstrained DLM using linear exposure-response function at each lag
fit_unconstr_dlm <- function(expos_var_names, data, form0) {
    form <- paste0(form0, " + ", paste(expos_var_names, collapse = " + "))
    fit0 <- lm(form, data = data)
    vc <- vcov(fit0)
    est_lags <- coef(fit0)[expos_var_names]
    vc_lags <- vc[expos_var_names, expos_var_names]
    cc <- rep(1, length(expos_var_names))
    cumm <- as.vector(cc %*% est_lags)
    cumm_se <- sqrt(as.vector((cc %*% vc_lags) %*% (cc)))
    cumm_est <- tibble(
        lag = NA,
        Estimate = cumm,
        SE = cumm_se,
        t = Estimate/SE,
        p = 2*(1 - pnorm(abs(t)))
    )
    summ <- summary(fit0)$coef %>%
        as_tibble(rownames = "variable") %>%
        filter(variable %in% expos_var_names) %>%
        transmute(
            lag = sapply(variable, function(x) which(expos_var_names == x)),
            Estimate,
            SE = `Std. Error`,
            t = `t value`,
            p = `Pr(>|t|)`
        ) %>%
        bind_rows(
            cumm_est
        )
}

## Fit constrained DLM with linear exposure-response function at each lag, using the `dlnm` R package
fit_DLM_constrained <- function(df_spl, max_lag, data, expos_var_names, form0) {
    expos_mat <- as.matrix(data[, expos_var_names])
    cb <- dlnm::crossbasis(expos_mat, lag = c(1,max_lag), argvar = list("lin"), arglag = list(fun = "ns", df = df_spl))
    model <- lm(as.formula(paste0(form0, " + cb")), data = data)
    pred <- dlnm::crosspred(cb, model, at = c(0, 1))
    fit <- tibble(
        lag = 1:max_lag %>% c(NA),
        Estimate = pred$matfit["1", ] %>% c(pred$allfit["1"]),
        SE = pred$matse["1", ] %>% c(pred$allse["1"])
    ) %>%
        mutate(t = Estimate/SE,
               p = 2*(1 - pnorm(abs(t))))

    fit
}

## For observed full-cohort data model:
  ## Set up model matrix with interaction terms
  ## data: should have a column 'n_prior_expos' indicating how many years of prior exposure history the person has
set_up_lagged_dset <- function(data, expos_names_prefix = "expos", n_cohorts = 12) {
    expos_var_names <- paste0(expos_names_prefix, 1:n_cohorts)

    ## create cohort variables
    anal2 <- data %>%
        mutate(cohort = factor(n_prior_expos, levels = n_cohorts:1)) %>%
        mutate(
            prior1 = n_prior_expos >= 1,
            prior2 = n_prior_expos >= 2,
            prior3 = n_prior_expos >= 3,
            prior4 = n_prior_expos >= 4,
            prior5 = n_prior_expos >= 5,
            prior6 = n_prior_expos >= 6,
            prior7 = n_prior_expos >= 7,
            prior8 = n_prior_expos >= 8,
            prior9 = n_prior_expos >= 9,
            prior10 = n_prior_expos >= 10,
            prior11 = n_prior_expos >= 11,
            prior12 = n_prior_expos >= 12
        )

    ## create interaction variables
    sep <- ""
    for (i in seq_along(expos_var_names)) {
        nm0 <- expos_var_names[i]
        lag <- gsub(expos_names_prefix, "", nm0) %>% as.numeric()
        vsel <- paste0("prior", lag)
        nm1 <- paste0(nm0, sep, vsel)
        anal2[[nm1]] <- with(anal2, as.numeric(ifelse(get(vsel), get(nm0), 0)))
    }
    expos_var_names_anal <- paste0(expos_var_names, sep, "prior", seq_along(expos_var_names))
    attr(anal2, "expos_var_names_anal") <- expos_var_names_anal
    anal2
}

## Fit unconstrained DLM using linear exposure-response function at each lag and combine across multiply-imputed datasets
fit_unconstr_dlm_MI <- function(imputed_data, expos_var_names, form0) {
    form <- paste0(form0, " + ", paste(expos_var_names, collapse = " + "))
    fitimp <- with(imputed_data, lm(as.formula(form)))
    pooled <- pool(fitimp)
    s <- summary(pooled) %>% as_tibble()
    fit_lag <- s %>%
        filter(term %in% expos_var_names) %>%
        transmute(
            lag = match(term, expos_var_names),
            Estimate = estimate,
            SE = std.error,
            t = statistic,
            p = p.value,
            df
        )
    ## overall estimate
    f <- expression({
        fit0 <- lm(as.formula(form))
        vc <- vcov(fit0)
        est_lags <- coef(fit0)[expos_var_names]
        vc_lags <- vc[expos_var_names, expos_var_names]
        cc <- rep(1, length(expos_var_names))
        cumm <- as.vector(cc %*% est_lags)
        cumm_se <- sqrt(as.vector((cc %*% vc_lags) %*% (cc)))
        tibble(estimate = cumm, std.error = cumm_se)
    })
    fitimp <- with(imputed_data, f)
    fitimp_output <- bind_rows(fitimp$analyses)
    res <- try(miceafter::pool_scalar_RR(est = fitimp_output$estimate, se = fitimp_output$std.error, dfcom = unique(pooled$pooled$dfcom)),  silent = TRUE)
    if (!inherits(res, "try-error")) {
        df <- if(is.null(res$v_adj)) {Inf} else {res$v_adj}
        fitcum <- res %$% tibble(lag = NA,
                                 Estimate = pool_est,
                                 SE = pool_se,
                                 t = Estimate/SE,
                                 df = df,
                                 p = 2*(1 - pt(abs(t), df))
        )
    } else {
        fitcum <- tibble(lag = NA,
                         Estimate = NA,
                         SE = NA,
                         t = NA,
                         df = NA,
                         p = NA)
    }

    ## prep output
    fit <- bind_rows(
        fit_lag,
        fitcum
    )
    fit
}

## Fit constrained DLM with linear exposure-response function at each lag, using the `dlnm` R package and combine across multiply-imputed datasets
fit_DLM_constrained_MI <- function(df_spl, max_lag, imputed_data, expos_var_names, form0) {
    ## Obtain pooled DLM coefficients
    fitimp <- with(imputed_data, {
        expos_mat <- sapply(expos_var_names, function(x) get(x))
        cb <- crossbasis(expos_mat, lag = c(1,max_lag), argvar = list("lin"), arglag = list(fun = "ns", df = df_spl))
        lm(as.formula(paste0(form0, " + cb")))
    })
    pooled <- pool(fitimp)
    ests <- lapply(fitimp$analyses, function(x) coef(x))
    vcs <- lapply(fitimp$analyses, function(x) vcov(x))
    pooled_dlm_cfs <- mitools::MIcombine(ests, vcs)

    ## Desired contrasts
    cb_pred <- crossbasis(rbind(rep(0, max_lag), rep(1, max_lag)), lag = c(1,max_lag), argvar = list("lin"), arglag = list(fun = "ns", df = df_spl))
    sel <- grep("cbv1", names(pooled_dlm_cfs$coefficients))
    pred <- crosspred(cb_pred, coef = pooled_dlm_cfs$coefficients[sel], vcov = pooled_dlm_cfs$variance[sel,sel], at = c(0, 1))

    ## prep output
    fit <- tibble(lag = c(1:max_lag, NA),
                  Estimate = c(pred$matfit["1", ], pred$allfit["1"]),
                  SE = c(pred$matse["1", ], pred$allse["1"]),
                  t = Estimate/SE,
                  p = 2*(1 - pnorm(abs(t))), ## use large sample
    )
    fit
}

## Fit constrained DLM with nonlinear exposure-response function at each lag, using the `dlnm` R package and combine across multiply-imputed datasets
## In this code, the exposure-response function is modeled using natural cubic splines with degrees of freedom given by `df_splX`
    ## To estimate lagged and cumulative effects, this code compares all exposures when their value is equal to 1 versus when their value is equal to 0
    ## to use different values instead of 0 and 1, either the function can be modified, or the exposures can all be transformed (centered and scaled) so that 0 corresponds to a different value (e.g., 25th percentile) and 1 corresponds to a second value (e.g., 75th percentile)
fit_DLM_constrained_MI_splX <- function(df_spl, max_lag, imputed_data, expos_var_names, form0, df_splX = 3, cen = 0) {

    comp <- complete(imputed_data)
    expos_comp <- comp[, expos_var_names]
    bs <- ns(data.matrix(expos_comp), df = df_splX)
    bs_info <- list(Boundary.knots = attr(bs, "Boundary.knots"),
                    knots = attr(bs, "knots"))

    ## Obtain pooled DLM coefficients
    fitimp <- with(imputed_data, {
        expos_mat <- sapply(expos_var_names, function(x) get(x))
        #cb <- crossbasis(expos_mat, lag = c(1,max_lag), argvar = list("lin"), arglag = list(fun = "ns", df = df_spl))
        cb <- crossbasis(expos_mat, lag = c(1,max_lag), argvar = list(fun = "ns", knots = bs_info$knots, Boundary.knots = bs_info$Boundary.knots), arglag = list(fun = "ns", df = df_spl))
        lm(as.formula(paste0(form0, " + cb")))
    })
    pooled <- pool(fitimp)
    ests <- lapply(fitimp$analyses, function(x) coef(x))
    vcs <- lapply(fitimp$analyses, function(x) vcov(x))
    pooled_dlm_cfs <- MIcombine(ests, vcs)

    ## Desired contrasts
    cb_pred <- crossbasis(rbind(rep(0, max_lag), rep(1, max_lag)), lag = c(1,max_lag), argvar = list(fun = "ns", knots = bs_info$knots, Boundary.knots = bs_info$Boundary.knots), arglag = list(fun = "ns", df = df_spl))
    sel <- grep("cbv", names(pooled_dlm_cfs$coefficients))
    pred <- crosspred(cb_pred, coef = pooled_dlm_cfs$coefficients[sel], vcov = pooled_dlm_cfs$variance[sel,sel], at = c(0, 1), cen = cen)

    ## prep output
    fit <- tibble(lag = c(1:max_lag, NA),
                  Estimate = c(pred$matfit["1", ], pred$allfit["1"]),
                  SE = c(pred$matse["1", ], pred$allse["1"]),
                  t = Estimate/SE,
                  p = 2*(1 - pnorm(abs(t))), ## use large sample
    )
    fit
}
