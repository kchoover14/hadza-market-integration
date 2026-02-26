######################## HADZA MARKET INTEGRATION -- REVISED ANALYSIS
# Berbesque & Hoover (2025)
# Journal of Development Studies
# DOI: https://doi.org/10.1080/00220388.2025.2555188
#
# This script consolidates data preparation, exploratory visualization,
# descriptive statistics, linear mixed-effects model, sensitivity analysis,
# and Bayesian confirmation for the portfolio page reanalysis.
# No analytical changes from the published paper.
#
# Data sourcing:
#   data2-hadza-depersonalized.csv -- depersonalized LEH field data
#   data3-hadza-analysis.csv -- analytical dataset derived from data2
#
# Output figures (saved as PNG):
#   revised-fig1.png -- two-panel time series (all individuals + by sex)
#   revised-fig2.png -- LMM predicted hypoplasia rates by birth year and sex
#   Supplemental figures saved as sfig1-sfig8 (not embedded in portfolio page)


######################## LIBRARIES

library(dplyr)        # data manipulation
library(ggplot2)      # data visualization
library(cowplot)      # arrange multi-panel figures
library(psych)        # descriptive statistics
library(car)          # model diagnostics, Levene test
library(lme4)         # linear mixed-effects models
library(lmerTest)     # significance testing for LMM fixed effects
library(effectsize)   # effect size estimation
library(effects)      # visualize model interactions
library(MuMIn)        # marginal and conditional R-squared
library(lattice)      # random effects dotplot
library(brms)         # Bayesian regression models
library(bayesplot)    # posterior predictive checks and MCMC diagnostics


######################## DATA PREPARATION

hadza = read.csv("data2-hadza-depersonalized.csv", stringsAsFactors = TRUE)

# create string factor for biological sex
hadza = rename(hadza, sexn = sex)
hadza$sex = with(hadza, ifelse(sexn > 1, "Female", "Male"))

# derive year of birth from age when data were collected
# if year of birth is already present, retain it; otherwise use derived value
hadza$yob2 = hadza$year - hadza$age
hadza = hadza |> mutate(yob = coalesce(yob, yob2))

# compute hypoplasia rate: total LEH count / total teeth observed
hadza$totalTeeth = rowSums(!is.na(hadza[, c("h_rxi1", "h_rxi2", "h_rxc",
                                            "h_lxi1", "h_lxi2", "h_lxc",
                                            "h_rni1", "h_lni1", "h_rni2",
                                            "h_lni2", "h_rnc", "h_lnc")]))
hadza$totalHypoplasia = rowSums(hadza[, c("h_rxi1", "h_rxi2", "h_rxc",
                                          "h_lxi1", "h_lxi2", "h_lxc",
                                          "h_rni1", "h_lni1", "h_rni2",
                                          "h_lni2", "h_rnc", "h_lnc")],
                                na.rm = TRUE)
hadza$hypoplasiaRate = with(hadza, totalHypoplasia / totalTeeth)
hadza = subset(hadza, !is.na(hypoplasiaRate))

# reduce to analytical variables
hadza = select(hadza, id, age, yob, sex, sexn,
               totalHypoplasia, totalTeeth, hypoplasiaRate)

# save analytical dataset
write.csv(hadza, "data3-hadza-analysis.csv", row.names = FALSE, quote = FALSE)


######################## DESCRIPTIVE STATISTICS

hadzaDesc = select(hadza, id, sex, age, yob, totalHypoplasia,
                   totalTeeth, hypoplasiaRate)

desc_sample = describe(hadzaDesc, na.rm = TRUE)
write.csv(desc_sample, "results-descriptives-sample.csv",
          quote = FALSE, row.names = TRUE, na = "")

desc_sex = describeBy(hadzaDesc, group = "sex", na.rm = TRUE, mat = TRUE)
write.csv(desc_sex, "results-descriptives-sex.csv",
          quote = FALSE, row.names = TRUE, na = "")

# tests for sex differences in hypoplasia count and rate
t.test(hadza$totalHypoplasia ~ hadza$sex)
kruskal.test(hadza$totalHypoplasia ~ hadza$sex)
t.test(hadza$hypoplasiaRate ~ hadza$sex)
kruskal.test(hadza$hypoplasiaRate ~ hadza$sex)

# variance
car::leveneTest(hypoplasiaRate ~ sex, data = hadza)


######################## SUPPLEMENTAL FIGURE -- DATA DISTRIBUTION
# sfig1: saved to repo but not embedded in portfolio page

means = hadza |>
  group_by(sex) |>
  summarize(Mean = mean(hypoplasiaRate, na.rm = TRUE))

sa = ggplot(hadza, aes(x = hypoplasiaRate)) +
  geom_histogram(aes(y = after_stat(density)),
                 colour = "black", fill = "white", bins = 5) +
  geom_vline(aes(xintercept = mean(hypoplasiaRate)),
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_density(alpha = 0.2) +
  annotate("text", x = mean(hadza$hypoplasiaRate), y = 0.55,
           label = round(mean(hadza$hypoplasiaRate), 2),
           color = "black", size = 4, hjust = -0.1) +
  theme_classic()

sb = ggplot(hadza, aes(x = hypoplasiaRate, color = sex)) +
  geom_histogram(aes(y = after_stat(density)),
                 fill = "white", position = "dodge", bins = 5) +
  geom_vline(data = means, aes(xintercept = Mean, color = sex),
             linetype = "dashed", linewidth = 1, show.legend = FALSE) +
  geom_density(alpha = 0.2, show.legend = FALSE) +
  geom_text(data = means,
            aes(x = Mean, y = 0.6, label = round(Mean, 2), color = sex),
            position = position_jitter(width = 0.15, height = 0.15),
            show.legend = FALSE, size = 4) +
  scale_color_viridis_d(end = 0.6) +
  theme_classic()

sfig1 = plot_grid(sa, sb, nrow = 2, labels = c("A", "B"),
                  rel_heights = c(1, 1))
ggsave("sfig1-distribution.png", plot = sfig1,
       width = 6, height = 8, dpi = 300)

rm(sa, sb, sfig1)


######################## FIGURE 1 -- TIME SERIES

a = ggplot(hadza, aes(x = yob, y = hypoplasiaRate)) +
  ggtitle("A") +
  stat_smooth(method = "lm", show.legend = FALSE) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1925, 2000, 5)) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14)) +
  xlab("") + ylab("Hypoplasia Rate")

b = ggplot(hadza, aes(x = yob, y = hypoplasiaRate,
                      fill = sex, colour = sex, shape = sex)) +
  ggtitle("B") +
  stat_smooth(method = "lm", show.legend = FALSE) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = seq(1925, 2000, 5)) +
  scale_color_viridis_d(begin = 0.15, end = 0.65) +
  scale_fill_viridis_d(begin = 0.15, end = 0.65) +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black", size = 14)) +
  xlab("Year of Birth") + ylab("Hypoplasia Rate")

fig1 = plot_grid(a, b, nrow = 2, labels = NULL, rel_heights = c(1, 1.25))
ggsave("revised-fig1.png", plot = fig1,
       width = 8, height = 10, dpi = 300)

rm(a, b, fig1, means)


######################## LINEAR MIXED-EFFECTS MODEL

# compute tooth offset to standardize for varying teeth observed
hadza$tooth_offset = log(hadza$totalTeeth / 12)

# center year of birth on 1975 (biologically meaningful turning point)
hadza$yob_centered_1975 = hadza$yob - 1975

# LMM: hypoplasia rate ~ birth year (centered 1975) * sex + offset + (1|id)
lmer_1975 = lme4::lmer(hypoplasiaRate ~ yob_centered_1975 * sex +
                         offset(tooth_offset) + (1 | id), data = hadza)

lmer_1975_test = lmerTest::lmer(hypoplasiaRate ~ yob_centered_1975 * sex +
                                  offset(tooth_offset) + (1 | id), data = hadza)
summary(lmer_1975_test)

# effect sizes for fixed effects
effectsize(lmer_1975_test)

# marginal (fixed) and conditional (fixed + random) R-squared
r.squaredGLMM(lmer_1975)

# post-1975 mean comparison by sex
hadza |>
  mutate(birth_era = ifelse(yob < 1975, "Before 1975", "After 1975")) |>
  group_by(birth_era, sex) |>
  summarise(mean_hypoplasia_rate = mean(hypoplasiaRate, na.rm = TRUE),
            .groups = "drop")


######################## FIGURE 2 -- LMM PREDICTED EFFECTS

eff = effect("yob_centered_1975:sex", lmer_1975)

png("revised-fig2.png", width = 8, height = 4, units = "in", res = 300)
plot(eff,
     xlab = "Year of Birth (centered at 1975)",
     ylab = "Predicted Hypoplasia Rate",
     main = "")
dev.off()

rm(eff)


######################## SUPPLEMENTAL FIGURES -- MODEL DIAGNOSTICS
# sfig2-sfig8: saved to repo but not embedded in portfolio page

# residual plot
png("sfig2-residuals.png", width = 6, height = 5, units = "in", res = 300)
plot(resid(lmer_1975) ~ fitted(lmer_1975),
     main = "Residuals vs Fitted",
     xlab = "Fitted values", ylab = "Residuals")
abline(h = 0, lty = 2)
dev.off()

# QQ plot of residuals with confidence interval
png("sfig3-qq.png", width = 6, height = 5, units = "in", res = 300)
car::qqPlot(resid(lmer_1975), main = "QQ Plot of Residuals")
dev.off()

# random effects distribution
re_df = as.data.frame(ranef(lmer_1975)$id)
re_df$id = rownames(re_df)
merged_df = merge(re_df, hadza, by = "id")

sfig4 = ggplot(re_df, aes(x = `(Intercept)`)) +
  geom_histogram(bins = 20, fill = "lightblue", color = "black") +
  labs(title = "Distribution of Random Effects (Intercept)",
       x = "Random Effect Value", y = "Frequency") +
  theme_classic()
ggsave("sfig4-re-histogram.png", plot = sfig4,
       width = 6, height = 5, dpi = 300)

sfig5 = ggplot(re_df, aes(x = `(Intercept)`)) +
  geom_density(fill = "lightblue", color = "black") +
  labs(title = "Density Plot of Random Effects (Intercept)",
       x = "Random Effect Value", y = "Density") +
  theme_classic()
ggsave("sfig5-re-density.png", plot = sfig5,
       width = 6, height = 5, dpi = 300)

# random effects boxplot with outlier label
outlier_value = boxplot.stats(re_df$`(Intercept)`)$out
outlier_row = which(re_df$`(Intercept)` == outlier_value)
outlier_id = rownames(re_df)[outlier_row]

sfig6 = ggplot(re_df, aes(y = `(Intercept)`)) +
  geom_boxplot(fill = "lightblue", color = "black",
               outlier.colour = "lightcoral", outlier.shape = 8,
               outlier.size = 3) +
  labs(title = "Boxplot of Random Effects (Intercept)",
       y = "Random Effect Value") +
  theme_classic() +
  geom_text(data = re_df[outlier_row, ],
            aes(x = 1, y = `(Intercept)`, label = outlier_id),
            nudge_y = 0.02, nudge_x = -0.975)
ggsave("sfig6-re-boxplot.png", plot = sfig6,
       width = 6, height = 5, dpi = 300)

# random effects by sex
sfig7 = ggplot(merged_df, aes(x = `(Intercept)`, color = sex)) +
  geom_density() +
  scale_color_viridis_d(end = 0.6) +
  labs(title = "Random Effects by Sex",
       x = "Random Effect Value", y = "Density") +
  theme_classic()
ggsave("sfig7-re-by-sex.png", plot = sfig7,
       width = 6, height = 5, dpi = 300)

rm(re_df, merged_df, sfig4, sfig5, sfig6, sfig7,
   outlier_value, outlier_row, outlier_id)


######################## SENSITIVITY ANALYSIS -- OUTLIER REMOVAL

# model summary with confidence intervals
summary_model = summary(lmer_1975_test)
coef_model = summary_model$coefficients
ci_lower = coef_model[, "Estimate"] - 1.96 * coef_model[, "Std. Error"]
ci_upper = coef_model[, "Estimate"] + 1.96 * coef_model[, "Std. Error"]
ci_full = cbind(coef_model, CI_lower = ci_lower, CI_upper = ci_upper)

# remove individual 450 (row 34)
hadza_no_450 = hadza[hadza$id != "450", ]
lmer_no_450 = lme4::lmer(hypoplasiaRate ~ yob * sex +
                           offset(tooth_offset) + (1 | id), data = hadza_no_450)
coef_no_450 = summary(lmer_no_450)$coefficients
ci_no_450 = cbind(coef_no_450,
                  CI_lower = coef_no_450[, "Estimate"] - 1.96 * coef_no_450[, "Std. Error"],
                  CI_upper = coef_no_450[, "Estimate"] + 1.96 * coef_no_450[, "Std. Error"])

# remove individual 16 (row 40)
hadza_no_16 = hadza[hadza$id != "16", ]
lmer_no_16 = lme4::lmer(hypoplasiaRate ~ yob * sex +
                          offset(tooth_offset) + (1 | id), data = hadza_no_16)
coef_no_16 = summary(lmer_no_16)$coefficients
ci_no_16 = cbind(coef_no_16,
                 CI_lower = coef_no_16[, "Estimate"] - 1.96 * coef_no_16[, "Std. Error"],
                 CI_upper = coef_no_16[, "Estimate"] + 1.96 * coef_no_16[, "Std. Error"])

# remove individual 861 (row 16)
hadza_no_861 = hadza[hadza$id != "861", ]
lmer_no_861 = lme4::lmer(hypoplasiaRate ~ yob * sex +
                           offset(tooth_offset) + (1 | id), data = hadza_no_861)
coef_no_861 = summary(lmer_no_861)$coefficients
ci_no_861 = cbind(coef_no_861,
                  CI_lower = coef_no_861[, "Estimate"] - 1.96 * coef_no_861[, "Std. Error"],
                  CI_upper = coef_no_861[, "Estimate"] + 1.96 * coef_no_861[, "Std. Error"])

# print sensitivity results
ci_full; ci_no_450; ci_no_16; ci_no_861

rm(hadza_no_450, hadza_no_16, hadza_no_861,
   lmer_no_450, lmer_no_16, lmer_no_861,
   coef_no_450, coef_no_16, coef_no_861,
   ci_no_450, ci_no_16, ci_no_861,
   ci_full, ci_lower, ci_upper, coef_model, summary_model)


######################## BRMS -- BAYESIAN CONFIRMATION

# compare LMM with and without random intercept
lmer_1975_fixed = lm(hypoplasiaRate ~ yob_centered_1975 * sex +
                       offset(tooth_offset), data = hadza)
anova(lmer_1975, lmer_1975_fixed)
# likelihood ratio test: chi-squared = 0.46, df = 1, p = 0.5
# no significant difference -- brms model uses fixed effects only

# derive weakly informative priors from LMM without random effect
lmer_summary = summary(lmer_1975_fixed)
lmer_coef = lmer_summary$coefficients
lmer_sd = lmer_coef[, "Std. Error"]

intercept_mean = lmer_coef["(Intercept)", "Estimate"]
intercept_sd = lmer_sd["(Intercept)"]
yob_centered_mean = lmer_coef["yob_centered_1975", "Estimate"]
yob_centered_sd = lmer_sd["yob_centered_1975"]
sexMale_mean = lmer_coef["sexMale", "Estimate"]
sexMale_sd = lmer_sd["sexMale"]
yob_sexMale_mean = lmer_coef["yob_centered_1975:sexMale", "Estimate"]
yob_sexMale_sd = lmer_sd["yob_centered_1975:sexMale"]

prior_list = c(
  prior(normal(intercept_mean, intercept_sd), class = "Intercept"),
  prior(normal(yob_centered_mean, yob_centered_sd), class = "b", coef = "yob_centered_1975"),
  prior(normal(sexMale_mean, sexMale_sd), class = "b", coef = "sexMale"),
  prior(normal(yob_sexMale_mean, yob_sexMale_sd), class = "b", coef = "yob_centered_1975:sexMale"),
  prior(normal(0, 0.2), class = "sigma", lb = 0)
)

stanvars_list = stanvar(intercept_mean, name = "intercept_mean") +
  stanvar(intercept_sd, name = "intercept_sd") +
  stanvar(yob_centered_mean, name = "yob_centered_mean") +
  stanvar(yob_centered_sd, name = "yob_centered_sd") +
  stanvar(sexMale_mean, name = "sexMale_mean") +
  stanvar(sexMale_sd, name = "sexMale_sd") +
  stanvar(yob_sexMale_mean, name = "yob_sexMale_mean") +
  stanvar(yob_sexMale_sd, name = "yob_sexMale_sd")

brms_model = brm(
  hypoplasiaRate ~ yob_centered_1975 * sex + offset(tooth_offset),
  data = hadza,
  prior = prior_list,
  stanvars = stanvars_list,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  iter = 10000,
  warmup = 5000
)

summary(brms_model)
posterior_summary(brms_model, variables = "sigma")
sd(hadza$hypoplasiaRate) # for context


######################## SUPPLEMENTAL FIGURES -- BRMS DIAGNOSTICS
# sfig8-sfig12: saved to repo but not embedded in portfolio page

y = hadza$hypoplasiaRate
y_rep = posterior_predict(brms_model, draws = 500)

png("sfig8-ppc-density.png", width = 7, height = 5, units = "in", res = 300)
ppc_dens_overlay(y, y_rep[1:100, ])
dev.off()

png("sfig9-ppc-hist.png", width = 7, height = 5, units = "in", res = 300)
ppc_hist(y, y_rep[1:10, ])
dev.off()

png("sfig10-ppc-scatter.png", width = 7, height = 5, units = "in", res = 300)
ppc_scatter_avg(y, y_rep)
dev.off()

png("sfig11-trace.png", width = 10, height = 6, units = "in", res = 300)
mcmc_trace(brms_model,
           pars = c("b_Intercept", "b_yob_centered_1975",
                    "b_sexMale", "b_yob_centered_1975:sexMale", "sigma"))
dev.off()

png("sfig12-acf.png", width = 10, height = 6, units = "in", res = 300)
mcmc_acf(brms_model,
         pars = c("b_Intercept", "b_yob_centered_1975",
                  "b_sexMale", "b_yob_centered_1975:sexMale", "sigma"))
dev.off()
