

# Setup ------------------------------------------------------------------------

library(lme4)
library(bootmlm)
library(lmeresampler)
library(boot)
library(MuMIn)
library(tidyverse)

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
  vcov_ranef <- vcov_vc(mod_0, sd_cor = FALSE)
  icc_var <- t(grad) %*% vcov_ranef %*% grad
  c(icc_est = icc_est, icc_var = icc_var)
}

# Parametric bootstrap ---------------------------------------------------------

# Fixed effects and Rsq
boo_par_fr <- bootstrap_mer(mod_f, fixef_rsq, nsim = 1999L, type = "parametric")
saveRDS(boo_par_fr, "draft/illustration_boo_res/boo_par_fr.rds")

# ICC
boo_par_icc <- bootstrap_mer(mod_0, icc, nsim = 1999L, type = "parametric")
saveRDS(boo_par_icc, "draft/illustration_boo_res/boo_par_icc.rds")

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

# Influence values
inf_val_fr <- lapply(1:9, function(i) {
  empinf_mer(mod_f, fixef_rsq, index = i)
})
saveRDS(inf_val_fr, "draft/illustration_boo_res/inf_val_fr.rds")

inf_val_icc <- empinf_mer(mod_0, icc, index = 1)
saveRDS(inf_val_icc, "draft/illustration_boo_res/inf_val_icc.rds")

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

# with lmeresampler
boo_cas1_lr_icc <- case_bootstrap(mod_0, .f = icc, 
                                  B = 1999L, resample = c(TRUE, TRUE))
saveRDS(boo_cas1_lr_icc, "draft/illustration_boo_res/boo_cas1_lr_icc.rds")

# Wild bootstrap ---------------------------------------------------------------

boo_wil_fr <- wild_bootstrap(mod_f, .f = fixef_rsq, B = 1999L)
saveRDS(boo_wil_fr, "draft/illustration_boo_res/boo_wil_fr.rds")

boo_wil_icc <- wild_bootstrap(mod_0, .f = icc, B = 1999L)
saveRDS(boo_wil_icc, "draft/illustration_boo_res/boo_wil_icc.rds")


# Load boot results ------------------------------------------------------------

boo_par_fr <- readRDS("illustration_boo_res/boo_par_fr.rds")
boo_par_icc <- readRDS("illustration_boo_res/boo_par_icc.rds")

boo_res_fr <- readRDS("illustration_boo_res/boo_res_fr.rds")
boo_cgr_fr <- readRDS("illustration_boo_res/boo_cgr_fr.rds")
boo_res_icc <- readRDS("illustration_boo_res/boo_res_icc.rds")
boo_cgr_icc <- readRDS("illustration_boo_res/boo_cgr_icc.rds")

inf_val_fr <- readRDS("illustration_boo_res/inf_val_fr.rds")
inf_val_icc <- readRDS("illustration_boo_res/inf_val_icc.rds")

boo_cas_fr <- readRDS("illustration_boo_res/boo_cas_fr.rds")
boo_cas1_fr <- readRDS("illustration_boo_res/boo_cas1_fr.rds")
boo_cas_icc <- readRDS("illustration_boo_res/boo_cas_icc.rds")
boo_cas1_icc <- readRDS("illustration_boo_res/boo_cas1_icc.rds")

boo_wil_fr <- readRDS("illustration_boo_res/boo_wil_fr.rds")
boo_wil_icc <- readRDS("illustration_boo_res/boo_wil_icc.rds")


# Confidence intervals ---------------------------------------------------------

# Parametric bootstrap
boo_par_fr_ci <- lapply(1:9, function(i) {
  boot.ci(boo_par_fr, type = c("norm", "basic", "perc"), index = i)
})
boo_par_icc_ci <- boot.ci(boo_par_icc, type = c("norm", "basic", "perc", "stud"))

# Residual bootstrap
boo_res_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_res_fr, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_fr[[i]])
})
boo_res_icc_ci <- boot::boot.ci(boo_res_icc, 
                                type = c("norm", "basic", "perc", "bca", "stud"), 
                                L = inf_val_icc)

# Cases bootstrap (resampling only clusters)
boo_cas_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_cas_fr, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_fr[[i]])
})
boo_cas_icc_ci <- boot::boot.ci(boo_cas_icc, 
                                type = c("norm", "basic", "perc", "bca", "stud"), 
                                L = inf_val_icc)

# Cases bootstrap (resampling both clusters and individuals)
boo_cas1_fr_ci <- lapply(1:9, function(i) {
  boot::boot.ci(boo_cas1_fr, type = c("norm", "basic", "perc", "bca"), 
                index = i, L = inf_val_fr)
})
boo_cas1_icc_ci <- boot::boot.ci(boo_cas1_icc, 
                                 type = c("norm", "basic", "perc", "bca", "stud"), 
                                 L = inf_val_icc)

# Wild bootstrap
boo_wil_fr_ci <- confint(boo_wil_fr, method = c("norm", "basic", "perc", "bca"))
boo_wil_icc_ci <- confint(boo_wil_icc, method = c("norm", "basic", "perc", "bca"))

# Logit transformation
boo_par_rsq_ci <- boot.ci(boo_par_fr, h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc"))
boo_res_rsq_ci <- boot.ci(boo_res_fr, L = inf_val_fr, h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc", "bca"))
boo_cas_rsq_ci <- boot.ci(boo_cas_fr, L = inf_val_fr, h = qlogis, 
                          hdot = function(x) 1 / (x - x^2), 
                          hinv = plogis, 
                          index = 9, type = c("norm", "basic", "perc", "bca"))
boo_cas1_rsq_ci <- boot.ci(boo_cas1_fr, L = inf_val_fr, h = qlogis, 
                           hdot = function(x) 1 / (x - x^2), 
                           hinv = plogis, 
                           index = 9, type = c("norm", "basic", "perc", "bca"))


# Make tables ------------------------------------------------------------------

# get bias-corrected estimate
get_est <- function(out) {
  if (any(class(out) == "boot")) {
    boot_mean <- apply(out$t, 2, mean)
    out$t0 * 2 - boot_mean
  } else if (any(class(out) == "lmeresamp")) {
    unname(2 * out$stats$observed - out$stats$rep.mean)
  }
}
# get bootstrap SE
get_se <- function(out) {
  if (any(class(out) == "boot")) {
    apply(out$t, 2, sd)
  } else if (any(class(out) == "lmeresamp")) {
    unname(out$stats$se)
  }
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
      paste0("(", paste(x, collapse = ", "), ")")
    })
  } else if (class(boo_ci) == "bootci") {
    # output from boot.ci (list of parameters)
    print_ci(boo_ci, type)
  } else {
    # outcome from boot.ci (single parameter)
    unlist(lapply(boo_ci, function(x) { print_ci(x, type) }))
  }
}

orig_est <- c(boo_par_fr$t0, NA, boo_par_icc$t0[1]) |> round(2)
boo_est <- t(data.frame(
  par = c(get_est(boo_par_fr), NA, get_est(boo_par_icc)[1]), 
  wild = c(get_est(boo_wil_fr), NA, get_est(boo_wil_icc)[1]), 
  res = c(get_est(boo_res_fr), NA, get_est(boo_res_icc)[1]), 
  cas = c(get_est(boo_cas_fr), NA, get_est(boo_cas_icc)[1]), 
  cas1 = c(get_est(boo_cas1_fr), NA, get_est(boo_cas1_icc)[1])
)) |> round(2)
boo_se <- t(data.frame(
  par = c(get_se(boo_par_fr), NA, get_se(boo_par_icc)[1]), 
  wild = c(get_se(boo_wil_fr), NA, get_se(boo_wil_icc)[1]), 
  res = c(get_se(boo_res_fr), NA, get_se(boo_res_icc)[1]), 
  cas = c(get_se(boo_cas_fr), NA, get_se(boo_cas_icc)[1]), 
  cas1 = c(get_se(boo_cas1_fr), NA, get_se(boo_cas1_icc)[1])
)) |> round(2)
norm_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "normal"), get_ci(boo_par_rsq_ci, "normal"), 
          get_ci(boo_par_icc_ci, "normal")), 
  wild = c(get_ci(boo_wil_fr_ci, "normal"), "()", 
           get_ci(boo_wil_icc_ci, "normal")[1]), 
  res = c(get_ci(boo_res_fr_ci, "normal"), get_ci(boo_res_rsq_ci, "normal"), 
          get_ci(boo_res_icc_ci, "normal")), 
  cas = c(get_ci(boo_cas_fr_ci, "normal"), get_ci(boo_cas_rsq_ci, "normal"), 
          get_ci(boo_cas_icc_ci, "normal")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "normal"), get_ci(boo_cas1_rsq_ci, "normal"), 
           get_ci(boo_cas1_icc_ci, "normal"))
))
basic_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "basic"), get_ci(boo_par_rsq_ci, "basic"), 
          get_ci(boo_par_icc_ci, "basic")), 
  wild = c(get_ci(boo_wil_fr_ci, "basic"), "()", 
           get_ci(boo_wil_icc_ci, "basic")[1]), 
  res = c(get_ci(boo_res_fr_ci, "basic"), get_ci(boo_res_rsq_ci, "basic"), 
          get_ci(boo_res_icc_ci, "basic")), 
  cas = c(get_ci(boo_cas_fr_ci, "basic"), get_ci(boo_cas_rsq_ci, "basic"), 
          get_ci(boo_cas_icc_ci, "basic")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "basic"), get_ci(boo_cas1_rsq_ci, "basic"), 
           get_ci(boo_cas1_icc_ci, "basic"))
))
perc_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "percent"), get_ci(boo_par_rsq_ci, "percent"), 
          get_ci(boo_par_icc_ci, "percent")), 
  wild = c(get_ci(boo_wil_fr_ci, "percent"), "()", 
           get_ci(boo_wil_icc_ci, "percent")[1]), 
  res = c(get_ci(boo_res_fr_ci, "percent"), get_ci(boo_res_rsq_ci, "percent"), 
          get_ci(boo_res_icc_ci, "percent")),
  cas = c(get_ci(boo_cas_fr_ci, "percent"), get_ci(boo_cas_rsq_ci, "percent"), 
          get_ci(boo_cas_icc_ci, "percent")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "percent"), get_ci(boo_cas1_rsq_ci, "percent"), 
           get_ci(boo_cas1_icc_ci, "percent"))
))
stud_ci <- t(data.frame(
  par = c(rep("()", each = 10), get_ci(boo_par_icc_ci, "student")), 
  wild = "()", 
  res = c(rep("()", each = 10), get_ci(boo_res_icc_ci, "student")), 
  cas = c(rep("()", each = 10), get_ci(boo_cas_icc_ci, "student")), 
  cas1 = c(rep("()", each = 10), get_ci(boo_cas1_icc_ci, "student"))
))
bca_ci <- t(data.frame(
  par = c(get_ci(boo_par_fr_ci, "bca"), get_ci(boo_par_rsq_ci, "bca"), 
          get_ci(boo_par_icc_ci, "bca")), 
  wild = "()", 
  res = c(get_ci(boo_res_fr_ci, "bca"), get_ci(boo_res_rsq_ci, "bca"), 
          get_ci(boo_res_icc_ci, "bca")), 
  cas = c(get_ci(boo_cas_fr_ci, "bca"), get_ci(boo_cas_rsq_ci, "bca"), 
          get_ci(boo_cas_icc_ci, "bca")), 
  cas1 = c(get_ci(boo_cas1_fr_ci, "bca"), get_ci(boo_cas1_rsq_ci, "bca"), 
           get_ci(boo_cas1_icc_ci, "bca"))
))

boo_tab <- rbind(
  orig_est, boo_est, boo_se, norm_ci, basic_ci, perc_ci, stud_ci, bca_ci
)
colnames(boo_tab) <- c(names(fixef(mod_f)), "rsq", "rsq_trans", "icc")
rownames(boo_tab) <- NULL
compare_boo <- data.frame(
  stat = c("Original Est.", 
           rep(c("Bias-Corrected Est.", "Bootstrap SE", "Normal", 
                 "Basic", "Percentile", "Studentized", "BCA"), each = 5)), 
  boo_type = c(" ", rep(c("Parametric", "Wild", "Residual", "Cases (level-2)", 
                          "Cases (both levels)"), 7)), 
  boo_tab
) %>%
  mutate_all(~ ifelse(. == "()" | is.na(.), "--", .))

saveRDS(compare_boo, "draft/illustration_boo_res/compare_boo.rds")
