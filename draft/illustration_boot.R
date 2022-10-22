

# Setup ------------------------------------------------------------------------

library(lme4)
library(bootmlm)
library(lmeresampler)
library(MuMIn)

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
  fixef <- mod@beta
  rsq <- MuMIn::r.squaredGLMM(mod)[1]
  c(fixef = fixef, rsq = rsq)
}
# Function for obtaining intraclass correlation
icc <- function(mod) {
  1 / (1 + mod@theta^(-2))
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

# Wild bootstrap ---------------------------------------------------------------

boo_wil_fr <- wild_bootstrap(mod_f, .f = fixef_rsq, B = 1999L)
saveRDS(boo_wil_fr, "draft/illustration_boo_res/boo_wil_fr.rds")

boo_wil_icc <- wild_bootstrap(mod_0, .f = icc, B = 1999L)
saveRDS(boo_wil_icc, "draft/illustration_boo_res/boo_wil_icc.rds")


