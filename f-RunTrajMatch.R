####################################################################################################
# Run trajectory matching for the coupling model 
####################################################################################################

FitOptim <- function(sim = sim_l$clim,
                     all_pars_nm = c("p", "R0", "alpha", "rho_mean", "rho_k", "e_Te", "e_RH"),
                     asum_e_RH = -0.2, asum_e_Te = -0.2, asum_coup = 0.01, covars = covars) {
  
  # Args:
  # sim: n simulations with df named CC_obs1 and CC_obs2 with columns per simulation (list)
  # all_pars_nm: names of the parameters to estimate (vector)
  # asum_e_RH: assumption and starting value on the RH effect (numeric)
  # asum_e_Te: assumption and starting value on the Te effect (numeric)
  # asum_coup: assumption and starting value on the coupling effect (0,1) (numeric)
  # covars: covars for the model (data frame)
  # Returns: Table with mle for the parameters (data frame)
  
  n_rep <- length(sim$CC_obs1)
  
  PompMod <- CreateMod(covars = covars)
  coef(PompMod)[names(parms)] <- parms
  
  # Assumption
  pomp::coef(PompMod, pars = c("e_Te", "e_RH")) <- c(asum_e_Te, asum_e_RH)
  pomp::coef(PompMod, pars = c("p")) <- asum_coup
  
  # Get initial values 
  base_pars_trans <- pomp::coef(PompMod, transform = TRUE)
  
  fits <- list()
  for(s in 1:n_rep) { print(s)
    # Assign data to pomp object
    PompMod@data[1, ] <- as.numeric(unlist(sim$CC_obs1[, s]))
    PompMod@data[2, ] <- as.numeric(unlist(sim$CC_obs2[, s]))
    # Create objective function
    LL_fun_cur <- traj_objfun(data = PompMod, est = all_pars_nm)
    parnames(LL_fun_cur) <- all_pars_nm
    # Run optimization
    fits[[s]] <- try(
      mle2(minuslogl = LL_fun_cur, 
           start = base_pars_trans[all_pars_nm], 
           optimizer = "optim", 
           method = "Nelder-Mead", 
           control = list(maxit = 1e9)))
  } 
  
  base_pars_df <- data.frame(par = all_pars_nm, 
                             value = unname(base_pars_trans[all_pars_nm]))
  
  # Remove estimations that failed,, if any
  fits_class <- sapply(fits, class)
  if(any(fits_class == "try-error")) fits <- fits[-which(fits_class == "try-error")]
  
  # Extracts information on model fit
  fits_info <- data.frame(sim_no = seq_along(fits), 
                          counts = map_dbl(.x = fits, .f = ~ .x@details$counts[1]),
                          ll = -map_dbl(.x = fits, .f = ~ .x@details$value), 
                          conv = map_int(.x = fits, .f = ~ .x@details$convergence))
  
  # Remove fits that didn't converge, if any
  if(any(fits_info$conv < 0)) {
    fits <- fits[-which(fits_info$conv < 0)]
  }
  
  # Extract MLE and SE of estimates
  pars_mle_se <- purrr::map(.x = fits, 
                            .f = ~ .x %>% summary() %>% coef() %>% as.data.frame() %>% rownames_to_column())
  
  pars_mle_se <- pars_mle_se %>% 
    bind_rows(.id = "sim_no") %>% 
    rename("par" = "rowname", 
           "mle" = "Estimate", 
           "se" = "Std. Error") %>% 
    select(sim_no, par, mle, se) %>% 
    left_join(y = base_pars_df %>% rename("true" = "value")) %>% 
    mutate(sim_no = as.integer(sim_no)) %>% 
    left_join(y = fits_info)
  
  # Replace NaN SE with NAs
  if(any(is.nan(pars_mle_se$se))) {
    pars_mle_se$se[is.nan(pars_mle_se$se)] <- NA
  }
  
  return(pars_mle_se)
}

#######################################################################################################
# End
#######################################################################################################
