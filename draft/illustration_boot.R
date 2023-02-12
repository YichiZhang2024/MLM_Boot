

# Setup ------------------------------------------------------------------------

library(lme4)
library(bootmlm)
library(lmeresampler)
library(boot)
library(MuMIn)
library(tidyverse)

# load data
temp <- tempfile()
download.file("https://www.stats.ox.ac.uk/~snijders/mlbook2_r_dat.zip", temp)
sj_dat <- read.table(unz(temp, "mlbook2_r.dat"), header = TRUE)

# unconditional model
mod_0 <- lmer(langPOST ~ (1 | schoolnr), data = sj_dat)
# final model
mod_f <- lmer(
  langPOST ~ IQ_verb*ses + sex + sch_iqv*sch_ses + (IQ_verb | schoolnr), 
  data = sj_dat
)


# Define functions -------------------------------------------------------------

# Function for obtaining fixed effects and R^2
fixef_rsq <- function(mod) {
  # fixed effect estimates
  fixef <- mod@beta
  # R^2 estimate
  rsq <- MuMIn::r.squaredGLMM(mod)[1]
  c(fixef = fixef, rsq = rsq)
}

# Function for obtaining intraclass correlation
icc <- function(mod) {
  # ICC estimate
  vc <- as.data.frame(VarCorr(mod))
  tau00sq <- vc[1, "vcov"]
  sigmasq <- vc[2, "vcov"]
  icc_est <- tau00sq / (tau00sq + sigmasq)
  # ICC variance with delta method
  dG_tau00sq <- sigmasq / (tau00sq + sigmasq)^2
  dG_sigmasq <- - tau00sq / (tau00sq + sigmasq)^2
  grad <- c(dG_tau00sq, dG_sigmasq)
  vcov_ranef <- vcov_vc(mod, sd_cor = FALSE)
  icc_var <- t(grad) %*% vcov_ranef %*% grad
  c(icc_est = icc_est, icc_var = icc_var)
}

# Function for obtaining random effects
raneff <- function(mod) {
  raneff <- as.data.frame(VarCorr(mod))$vcov
  c(raneff = raneff)
}

# Parametric bootstrap ---------------------------------------------------------

# Fixed effects and Rsq
boo_par_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, type = "parametric")
saveRDS(boo_par_fr, "draft/illustration_boo_res/boo_par_fr.rds")

# ICC
boo_par_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "parametric")
saveRDS(boo_par_icc, "draft/illustration_boo_res/boo_par_icc.rds")

# Random effects
boo_par_ran <- bootstrap_mer(mod_f, raneff, nsim = 1999L, type = "parametric")
saveRDS(boo_par_ran, "draft/illustration_boo_res/boo_par_ran.rds")

# Residual bootstrap -----------------------------------------------------------

boo_res_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, type = "residual")
boo_cgr_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, 
                            type = "residual_cgr")
saveRDS(boo_res_fr, "draft/illustration_boo_res/boo_res_fr.rds")
saveRDS(boo_cgr_fr, "draft/illustration_boo_res/boo_cgr_fr.rds")

boo_res_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "residual")
boo_cgr_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "residual_cgr")
saveRDS(boo_res_icc, "draft/illustration_boo_res/boo_res_icc.rds")
saveRDS(boo_cgr_icc, "draft/illustration_boo_res/boo_cgr_icc.rds")

boo_res_ran <- bootstrap_mer(mod_f, raneff, nsim = 1999L, type = "residual")
boo_cgr_ran <- bootstrap_mer(mod_f, raneff, nsim = 1999L, type = "residual_cgr")
saveRDS(boo_res_ran, "draft/illustration_boo_res/boo_res_ran.rds")
saveRDS(boo_cgr_ran, "draft/illustration_boo_res/boo_cgr_ran.rds")

# Influence values
inf_val_fr <- lapply(1:9, function(i) {
  empinf_mer(mod_f, fixef_rsq, index = i)
})
saveRDS(inf_val_fr, "draft/illustration_boo_res/inf_val_fr.rds")

inf_val_icc <- empinf_mer(mod_0, icc, index = 1)
saveRDS(inf_val_icc, "draft/illustration_boo_res/inf_val_icc.rds")

inf_val_ran <- lapply(1:4, function(i) {
  empinf_mer(mod_f, raneff, index = i)
})
saveRDS(inf_val_ran, "draft/illustration_boo_res/inf_val_ran.rds")

# Case bootstrap ---------------------------------------------------------------

boo_cas_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, type = "case")
boo_cas1_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, type = "case", 
                             lv1_resample = TRUE)
saveRDS(boo_cas_fr, "draft/illustration_boo_res/boo_cas_fr.rds")
saveRDS(boo_cas1_fr, "draft/illustration_boo_res/boo_cas1_fr.rds")

boo_cas_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "case")
boo_cas1_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "case", 
                              lv1_resample = TRUE)
saveRDS(boo_cas_icc, "draft/illustration_boo_res/boo_cas_icc.rds")
saveRDS(boo_cas1_icc, "draft/illustration_boo_res/boo_cas1_icc.rds")

boo_cas_ran <- bootstrap_mer(mod_f, raneff, nsim = 1999L, type = "case")
boo_cas1_ran <- bootstrap_mer(mod_f, raneff, nsim = 3999L, type = "case", 
                              lv1_resample = TRUE)
saveRDS(boo_cas_ran, "draft/illustration_boo_res/boo_cas_ran.rds")
saveRDS(boo_cas1_ran, "draft/illustration_boo_res/boo_cas1_ran.rds")

# with lmeresampler
boo_cas1_lr_icc <- case_bootstrap(mod_0, .f = icc, 
                                  B = 1999L, resample = c(TRUE, TRUE))
saveRDS(boo_cas1_lr_icc, "draft/illustration_boo_res/boo_cas1_lr_icc.rds")

# Wild bootstrap ---------------------------------------------------------------

boo_wil_fr <- wild_bootstrap(mod_f, .f = fixef_rsq, B = 1999L)
saveRDS(boo_wil_fr, "draft/illustration_boo_res/boo_wil_fr.rds")

boo_wil_icc <- wild_bootstrap(mod_0, .f = icc, B = 1999L)
saveRDS(boo_wil_icc, "draft/illustration_boo_res/boo_wil_icc.rds")

boo_wil_ran <- wild_bootstrap(mod_f, .f = raneff, B = 1999L)
saveRDS(boo_wil_ran, "draft/illustration_boo_res/boo_wil_ran.rds")


# Load boot results ------------------------------------------------------------

boo_par_fr <- readRDS("draft/illustration_boo_res/boo_par_fr.rds")
boo_par_icc <- readRDS("draft/illustration_boo_res/boo_par_icc.rds")
boo_par_ran <- readRDS("draft/illustration_boo_res/boo_par_ran.rds")

boo_res_fr <- readRDS("draft/illustration_boo_res/boo_res_fr.rds")
boo_cgr_fr <- readRDS("draft/illustration_boo_res/boo_cgr_fr.rds")
boo_res_icc <- readRDS("draft/illustration_boo_res/boo_res_icc.rds")
boo_cgr_icc <- readRDS("draft/illustration_boo_res/boo_cgr_icc.rds")
boo_res_ran <- readRDS("draft/illustration_boo_res/boo_res_ran.rds")
boo_cgr_ran <- readRDS("draft/illustration_boo_res/boo_cgr_ran.rds")

inf_val_fr <- readRDS("draft/illustration_boo_res/inf_val_fr.rds")
inf_val_icc <- readRDS("draft/illustration_boo_res/inf_val_icc.rds")
inf_val_ran <- readRDS("draft/illustration_boo_res/inf_val_ran.rds")

boo_cas_fr <- readRDS("draft/illustration_boo_res/boo_cas_fr.rds")
boo_cas1_fr <- readRDS("draft/illustration_boo_res/boo_cas1_fr.rds")
boo_cas_icc <- readRDS("draft/illustration_boo_res/boo_cas_icc.rds")
boo_cas1_icc <- readRDS("draft/illustration_boo_res/boo_cas1_icc.rds")
boo_cas_ran <- readRDS("draft/illustration_boo_res/boo_cas_ran.rds")
boo_cas1_ran <- readRDS("draft/illustration_boo_res/boo_cas1_ran.rds")

boo_wil_fr <- readRDS("draft/illustration_boo_res/boo_wil_fr.rds")
boo_wil_icc <- readRDS("draft/illustration_boo_res/boo_wil_icc.rds")
boo_wil_ran <- readRDS("draft/illustration_boo_res/boo_wil_ran.rds")


# Confidence intervals ---------------------------------------------------------

# Parametric bootstrap
boo_par_fr_ci <- lapply(1:9, function(i) {
  boot.ci(boo_par_fr, type = c("norm", "basic", "perc"), index = i)
})
boo_par_icc_ci <- boot.ci(boo_par_icc, 
                          type = c("norm", "basic", "perc", "stud"))
boo_par_ran_ci <- lapply(1:4, function(i) {
  boot.ci(boo_par_ran, type = c("norm", "basic", "perc"), index = i)
})

# Residual bootstrap
boo_res_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_res_fr, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_fr[[i]])
})
boo_res_icc_ci <- boot::boot.ci(boo_res_icc, 
                                type = c("norm", "basic", "perc", "bca", "stud"), 
                                L = inf_val_icc)
boo_res_ran_ci <- lapply(1:4, function(i) {
  boot::boot.ci(boo_res_ran, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_ran[[i]])
})

# Cases bootstrap (resampling only clusters)
boo_cas_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_cas_fr, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_fr[[i]])
})
boo_cas_icc_ci <- boot::boot.ci(boo_cas_icc, 
                                type = c("norm", "basic", "perc", "bca", "stud"), 
                                L = inf_val_icc)
boo_cas_ran_ci <- lapply(1:4, function(i) {
  boot::boot.ci(boo_cas_ran, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_ran[[i]])
})

# Cases bootstrap (resampling both clusters and individuals)
boo_cas1_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_cas1_fr, type = c("norm", "basic", "perc", "bca"),
                index = i, L = inf_val_fr[[i]])
})
boo_cas1_icc_ci <- boot::boot.ci(boo_cas1_icc,
                                 type = c("norm", "basic", "perc", "bca", "stud"),
                                 L = inf_val_icc)
boo_cas1_ran_ci134 <- lapply(c(1, 3, 4), function(i) {
  boot::boot.ci(boo_cas1_ran, type = c("norm", "basic", "perc", "bca"),
                index = i, L = inf_val_ran[[i]])
})
boo_cas1_ran_ci2 <- lapply(c(2), function(i) {
  boot::boot.ci(boo_cas1_ran, type = c("norm", "basic", "perc"),
                index = i, L = inf_val_ran[[i]])
})
boo_cas1_ran_ci <- list(boo_cas1_ran_ci134[[1]], boo_cas1_ran_ci2[[1]], 
                        boo_cas1_ran_ci134[[2]], boo_cas1_ran_ci134[[3]])

# Wild bootstrap
boo_wil_fr_ci <- confint(boo_wil_fr, method = c("norm", "basic", "perc", "bca"))
boo_wil_icc_ci <- confint(boo_wil_icc, method = c("norm", "basic", "perc", "bca"))
boo_wil_ran_ci <- confint(boo_wil_ran, method = c("norm", "basic", "perc", "bca"))

# Logit transformation
boo_par_rsq_ci <- boot.ci(boo_par_fr, h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc"))
boo_res_rsq_ci <- boot.ci(boo_res_fr, L = inf_val_fr[[9]], h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc", "bca"))
boo_cas_rsq_ci <- boot.ci(boo_cas_fr, L = inf_val_fr[[9]], h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc", "bca"))
boo_cas1_rsq_ci <- boot.ci(boo_cas1_fr, L = inf_val_fr[[9]], h = qlogis, 
                           hdot = function(x) 1 / (x - x^2), 
                           hinv = plogis, 
                           index = 9, type = c("norm", "basic", "perc", "bca"))

# Baseline CIs
base_ci <- confint(profile(mod_f, prof.scale = "varcov"))

# Baseline standard errors
sum_f <- summary(mod_f)
fixef_se <- sum_f$coefficients[, "Std. Error"]
ranef_se <- sqrt(diag(vcov_vc(mod_f, sd_cor = FALSE)))[c(1, 3, 2, 4)]
base_se <- c(fixef_se, ranef_se)

# Make tables ------------------------------------------------------------------

# get bias-corrected estimate
get_est <- function(out) {
  if (any(class(out) == "boot")) {
    boot_mean <- apply(out$t, 2, mean)
    bce <- out$t0 * 2 - boot_mean
  } else if (any(class(out) == "lmeresamp")) {
    bce <- unname(2 * out$stats$observed - out$stats$rep.mean)
  }
  sprintf(fmt = "%.2f", round(bce, 2))
}
# get bootstrap SE
get_se <- function(out) {
  if (any(class(out) == "boot")) {
    bse <- apply(out$t, 2, sd)
  } else if (any(class(out) == "lmeresamp")) {
    bse <- unname(out$stats$se)
  }
  sprintf(fmt = "%.2f", round(bse, 2))
}
comma <- function(x) {
  if (is.null(x)) NULL
  else gsub(" ", "", format(sprintf(fmt = "%.2f", round(x, digits = 2)), 
                            big.mark = ","))
}
# print confidence interval
print_ci <- function(boo_ci, type) {
  paste0("(", paste(comma(tail(boo_ci[[type]][1, ], 2L)), collapse = ", "), ")")
}
# get confidence interval
get_ci <- function(boo_ci, type) {
  if (any(class(boo_ci) == "data.frame")) {
    # output from confint, lmersampler 
    if (type == "normal") type_ci = "norm"
    if (type == "percent") type_ci = "perc"
    if (type == "basic") type_ci = "basic"
    if (type == "bca") return("()")
    ci_df <- subset(boo_ci, type == type_ci, select = c(lower, upper)) |>
      round(2L)
    apply(t(ci_df), 2, function(x) {
      paste0("(", paste(sprintf(fmt = "%.2f", x), collapse = ", "), ")")
    })
  } else if (class(boo_ci) == "bootci") {
    # output from boot.ci (list of parameters)
    print_ci(boo_ci, type)
  } else {
    # outcome from boot.ci (single parameter)
    unlist(lapply(boo_ci, function(x) { print_ci(x, type) }))
  }
}

base_ci_str <- apply(t(base_ci |> round(2)), 2, 
                     function(x) {
                       paste0("(", paste(sprintf(x, fmt = "%.2f"), collapse = ", "), ")")
                     })

orig_est <- c(boo_par_fr$t0, NA, boo_par_icc$t0[1], boo_par_ran$t0) |> 
  round(2) |> sprintf(fmt = "%.2f")
orig_ci <- c(base_ci_str[5:12], NA, NA, NA, base_ci_str[c(1, 3, 2, 4)])
orig_se <- c(base_se[1:8], NA, NA, NA, base_se[9:12]) |> 
  round(2) |> sprintf(fmt = "%.2f")
boo_est <- t(data.frame(
  par = c(get_est(boo_par_fr), NA, get_est(boo_par_icc)[1], 
          get_est(boo_par_ran)), 
  res = c(get_est(boo_res_fr), NA, get_est(boo_res_icc)[1], 
          get_est(boo_res_ran)), 
  wild = c(get_est(boo_wil_fr), NA, get_est(boo_wil_icc)[1], 
           get_est(boo_wil_ran)), 
  cas = c(get_est(boo_cas_fr), NA, get_est(boo_cas_icc)[1], 
          get_est(boo_cas_ran)), 
  cas1 = c(get_est(boo_cas1_fr), NA, get_est(boo_cas1_icc)[1], 
           get_est(boo_cas1_ran))
))
boo_se <- t(data.frame(
  par = c(get_se(boo_par_fr), NA, get_se(boo_par_icc)[1], 
          get_se(boo_par_ran)), 
  res = c(get_se(boo_res_fr), NA, get_se(boo_res_icc)[1], 
          get_se(boo_res_ran)), 
  wild = c(get_se(boo_wil_fr), NA, get_se(boo_wil_icc)[1], 
           get_se(boo_wil_ran)), 
  cas = c(get_se(boo_cas_fr), NA, get_se(boo_cas_icc)[1], 
          get_se(boo_cas_ran)), 
  cas1 = c(get_se(boo_cas1_fr), NA, get_se(boo_cas1_icc)[1], 
           get_se(boo_cas1_ran))
))
norm_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "normal"), get_ci(boo_par_rsq_ci, "normal"), 
          get_ci(boo_par_icc_ci, "normal"), get_ci(boo_par_ran_ci, "normal")), 
  res = c(get_ci(boo_res_fr_ci, "normal"), get_ci(boo_res_rsq_ci, "normal"), 
          get_ci(boo_res_icc_ci, "normal"), get_ci(boo_res_ran_ci, "normal")), 
  wild = c(get_ci(boo_wil_fr_ci, "normal"), "()", 
           get_ci(boo_wil_icc_ci, "normal")[1], 
           get_ci(boo_wil_ran_ci, "normal")), 
  cas = c(get_ci(boo_cas_fr_ci, "normal"), get_ci(boo_cas_rsq_ci, "normal"), 
          get_ci(boo_cas_icc_ci, "normal"), get_ci(boo_cas_ran_ci, "normal")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "normal"), get_ci(boo_cas1_rsq_ci, "normal"), 
           get_ci(boo_cas1_icc_ci, "normal"), 
           get_ci(boo_cas1_ran_ci, "normal")
  )
))
basic_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "basic"), get_ci(boo_par_rsq_ci, "basic"), 
          get_ci(boo_par_icc_ci, "basic"), get_ci(boo_par_ran_ci, "basic")), 
  res = c(get_ci(boo_res_fr_ci, "basic"), get_ci(boo_res_rsq_ci, "basic"), 
          get_ci(boo_res_icc_ci, "basic"), get_ci(boo_res_ran_ci, "basic")), 
  wild = c(get_ci(boo_wil_fr_ci, "basic"), "()", 
           get_ci(boo_wil_icc_ci, "basic")[1], 
           get_ci(boo_wil_ran_ci, "basic")), 
  cas = c(get_ci(boo_cas_fr_ci, "basic"), get_ci(boo_cas_rsq_ci, "basic"), 
          get_ci(boo_cas_icc_ci, "basic"), get_ci(boo_cas_ran_ci, "basic")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "basic"), get_ci(boo_cas1_rsq_ci, "basic"), 
           get_ci(boo_cas1_icc_ci, "basic"), 
           get_ci(boo_cas1_ran_ci, "basic")
  )
))
perc_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "percent"), get_ci(boo_par_rsq_ci, "percent"), 
          get_ci(boo_par_icc_ci, "percent"), get_ci(boo_par_ran_ci, "percent")), 
  res = c(get_ci(boo_res_fr_ci, "percent"), get_ci(boo_res_rsq_ci, "percent"), 
          get_ci(boo_res_icc_ci, "percent"), get_ci(boo_res_ran_ci, "percent")),
  wild = c(get_ci(boo_wil_fr_ci, "percent"), "()", 
           get_ci(boo_wil_icc_ci, "percent")[1], 
           get_ci(boo_wil_ran_ci, "percent")), 
  cas = c(get_ci(boo_cas_fr_ci, "percent"), get_ci(boo_cas_rsq_ci, "percent"), 
          get_ci(boo_cas_icc_ci, "percent"), get_ci(boo_cas_ran_ci, "percent")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "percent"), get_ci(boo_cas1_rsq_ci, "percent"), 
           get_ci(boo_cas1_icc_ci, "percent"), 
           get_ci(boo_cas1_ran_ci, "percent")
  )
))
stud_ci <- t(data.frame(
  par = c(rep("()", each = 10), get_ci(boo_par_icc_ci, "student"), 
          rep("()", each = 4)), 
  res = c(rep("()", each = 10), get_ci(boo_res_icc_ci, "student"), 
          rep("()", each = 4)), 
  wild = "()", 
  cas = c(rep("()", each = 10), get_ci(boo_cas_icc_ci, "student"), 
          rep("()", each = 4)), 
  cas1 = c(rep("()", each = 10), get_ci(boo_cas1_icc_ci, "student"), 
           rep("()", each = 4))
))
bca_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "bca"), get_ci(boo_par_rsq_ci, "bca"), 
          get_ci(boo_par_icc_ci, "bca"), get_ci(boo_par_ran_ci, "bca")), 
  res = c(get_ci(boo_res_fr_ci, "bca"), get_ci(boo_res_rsq_ci, "bca"), 
          get_ci(boo_res_icc_ci, "bca"), get_ci(boo_res_ran_ci, "bca")), 
  wild = "()", 
  cas = c(get_ci(boo_cas_fr_ci, "bca"), get_ci(boo_cas_rsq_ci, "bca"), 
          get_ci(boo_cas_icc_ci, "bca"), get_ci(boo_cas_ran_ci, "bca")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "bca"), get_ci(boo_cas1_rsq_ci, "bca"), 
           get_ci(boo_cas1_icc_ci, "bca"), 
           get_ci(boo_cas1_ran_ci[1:3], "bca"), "()"
  )
))

boo_tab <- rbind(
  orig_est, boo_est, orig_se, boo_se, orig_ci, 
  norm_ci, basic_ci, perc_ci, stud_ci, bca_ci
)
colnames(boo_tab) <- c(names(fixef(mod_f)), "rsq", "rsq_trans", "icc", 
                       "tau00sq", "tau11sq", "tau01", "sigmasq")
rownames(boo_tab) <- NULL
boo_names <- c("Parametric", "Residual", "Wild", "Cases (level-2)", 
               "Cases (both levels)")
ci_types <- c("Normal", "Basic", "Percentile", "Studentized", "BCa")
compare_boo <- data.frame(
  stat = c("Est.", rep("Bias-Corrected Est.", 5), "Asymptotic SE", 
           rep("Bootstrap SE", 5), "Likelihood-Based CI", 
           rep(ci_types, each = 5)), 
  boo_type = c(" ", boo_names, " ", boo_names, " ", rep(boo_names, 5)), 
  boo_tab
) %>%
  mutate_all(~ ifelse(. %in% c("()", "NA") | is.na(.), "--", .))

saveRDS(compare_boo, "draft/illustration_boo_res/compare_boo.rds")
