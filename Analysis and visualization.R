library(tidyverse)
library(deSolve)
library(phaseR)
library(ggplot2)
library(foreach)
library(doSNOW)
library(doParallel)
library(viridis)
library(ggpmisc)
library(car)
library(patchwork)
library(mgcv)
library(gratia)
library(gt)
library(broom)


data_HOI_1 <- read.csv('full_results_norm_and_lognorm.csv')
data_HOI_2 <- read.csv('new_combinations_norm_and_lognorm.csv')
data_HOI_full <- bind_rows(data_HOI_1, data_HOI_2)

coexistence_comparison <- data_HOI_full %>%
  mutate(model = str_replace_all(model, " - ", "_"),
         propn = n_end / n_start) %>%
  group_by(model, A_mean, A_sd, B_mean, B_sd, B_diag) %>%
  summarize(mean_propn = mean(propn, na.rm = TRUE),
            # sd_propn = sd(propn, na.rm = TRUE),
            stability_rate = mean(stability_result, na.rm = TRUE),
            .groups = 'drop')

# stability analysis under t-test
cc_stability <- coexistence_comparison %>%
  pivot_wider(names_from = model,
              values_from = c(mean_propn, stability_rate),
              names_sep = '_')

t_test_lognorm_stability <- t.test(cc_stability$stability_rate_HOI_lognorm, 
                                   cc_stability$stability_rate_Pairwise, 
                                   paired = TRUE,
                                   alternative = "two.sided")
print(t_test_lognorm_stability)

t_test_norm_stability <- t.test(cc_stability$stability_rate_HOI_norm, 
                                cc_stability$stability_rate_Pairwise, 
                                paired = TRUE,
                                alternative = "two.sided")
print(t_test_norm_stability)

coexistence_comparison$model <- as.factor(coexistence_comparison$model)
coexistence_comparison$model <- relevel(coexistence_comparison$model, ref = 'Pairwise')

stability_full_model <- lm(stability_rate ~ (A_mean + A_sd + B_mean + B_sd + B_diag) * model,
                           data = coexistence_comparison)
summary(stability_full_model)

anova_stability_full_model <- Anova(stability_full_model, type="III")
anova_stability_full_model

stability_model_mean <- lm(stability_rate ~ (A_mean + B_mean) * model,
                           data = coexistence_comparison)
summary(stability_model_mean)


stability_model_mean_coeffs <- tidy(stability_model_mean, conf.int = TRUE)


data_stability_model_mean_coeffs <- stability_model_mean_coeffs %>%
  filter(term != '(Intercept)') %>%
  mutate(term = fct_recode(term,
                           'Pairwise interaction strength (A_mean)' = 'A_mean',
                           'HOI Interaction strength (B_mean)' = 'B_mean',
                           'Competitive HOI model (HOI_lognorm)' = 'modelHOI_lognorm',
                           'Competitive & facilitative HOI model (HOI_norm)' = 'modelHOI_norm',
                           'A_mean x HOI_lognorm' = 'A_mean:modelHOI_lognorm',
                           'A_mean x HOI_norm' = 'A_mean:modelHOI_norm',
                           'B_mean x HOI_lognorm' = 'B_mean:modelHOI_lognorm',
                           'B_mean x HOI_norm' = 'B_mean:modelHOI_norm')) %>%
  mutate(term = fct_reorder(term, estimate))

plot_stability_model_mean_coeffs <- ggplot(data_stability_model_mean_coeffs,
                                           aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = "grey70") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.3, linewidth = 0.8) +
  geom_point(size = 4, color = "darkblue") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  labs(x = 'Coefficient (Effect on stability rate)',
       y = 'Parameter / Interaction Term') +
  theme_minimal(base_size = 24)

print(plot_stability_model_mean_coeffs)

# creating dataset for survival analysis
HOI_dataset <- coexistence_comparison %>%
  pivot_wider(names_from = model,
              values_from = c(mean_propn, stability_rate),
              names_sep = '_')%>%
  filter(stability_rate_Pairwise >= 0.9,
         stability_rate_HOI_norm >= 0.9,
         stability_rate_HOI_lognorm >= 0.9) %>%
  mutate(HOI_norm_effect = mean_propn_HOI_norm - mean_propn_Pairwise,
         HOI_norm_relative = mean_propn_HOI_norm / mean_propn_Pairwise,
         
         HOI_lognorm_effect = mean_propn_HOI_lognorm - mean_propn_Pairwise,
         HOI_lognorm_relative = mean_propn_HOI_lognorm / mean_propn_Pairwise,
         
         norm_promotes = HOI_norm_effect > 0,
         lognorm_promotes = HOI_lognorm_effect > 0)

# survival analysis under t-test
t_test_lognorm <- t.test(HOI_dataset$mean_propn_HOI_lognorm,
                         HOI_dataset$mean_propn_Pairwise,
                         paired = TRUE, 
                         alternative = "two.sided")
print(t_test_lognorm)

t_test_norm <- t.test(HOI_dataset$mean_propn_HOI_norm, 
                      HOI_dataset$mean_propn_Pairwise, 
                      paired = TRUE,
                      alternative = "two.sided")
print(t_test_norm)

t_test_HOI <- t.test(HOI_dataset$mean_propn_HOI_lognorm, 
                     HOI_dataset$mean_propn_HOI_norm, 
                     paired = TRUE,
                     alternative = "two.sided")
print(t_test_HOI)

#linear model analysis for HOI model
HOI_lm <- HOI_dataset %>%
  mutate(log_ratio_mean = log(B_mean / A_mean),
         log_sd_mean = log(B_sd / A_sd))


model_lm_lognorm <- lm(HOI_lognorm_effect ~ log_ratio_mean + log_sd_mean + B_diag, data = HOI_lm)
Anova(model_lm_lognorm, type = "III")
summary(model_lm_lognorm )

model_lm_norm <- lm(HOI_norm_effect ~ log_ratio_mean + log_sd_mean + B_diag, data = HOI_lm)
Anova(model_lm_norm, type = "III")
summary(model_lm_norm )

#logistic model for promoting rate analysis
logitmodel_HOI_lognorm <- glm(lognorm_promotes ~ A_mean + A_sd + B_mean + B_sd + B_diag,
                              data = HOI_dataset, 
                              family = binomial(link = "logit"))
summary(logitmodel_HOI_lognorm)

Anova(logitmodel_HOI_lognorm, type = "III")

logitmodel_HOI_norm <- glm(norm_promotes ~ A_mean + A_sd + B_mean + B_sd + B_diag,
                           data = HOI_dataset, 
                           family = binomial(link = "logit"))
summary(logitmodel_HOI_norm)

Anova(logitmodel_HOI_norm, type = "III")


# visualization for survival analysis
common_theme <- theme_classic(base_size = 20) + 
  theme(strip.text = element_text(size = 14, face = 'bold'),
        plot.margin = margin(5, 5, 5, 5))

## under B_diag = 0 situation
HOI_dataset_0 <- HOI_dataset %>% 
  filter(B_diag == 0)

lm_HOI_0 <- HOI_dataset_0 %>%
  mutate(log_ratio_mean = log(B_mean / A_mean),
         sd_ratio = B_sd / A_sd) %>%
  mutate(sd_ratio_group = cut(sd_ratio,
                              breaks = quantile(sd_ratio, probs = c(0, 0.33, 0.66, 1.0), na.rm = TRUE),
                              labels = c('Low ratio', 'Medium ratio', 'High ratio'),
                              include.lowest = TRUE)) %>%
  pivot_longer(cols = c(mean_propn_Pairwise, mean_propn_HOI_norm, mean_propn_HOI_lognorm),
               names_to = 'model',
               values_to = 'propn') %>%
  mutate(model = factor(model, 
                        levels = c('mean_propn_Pairwise', 'mean_propn_HOI_norm', 'mean_propn_HOI_lognorm'),
                        labels = c('Pairwise', 'HOI_norm', 'HOI_lognorm')))

plot_lm_HOI_0 <- ggplot(lm_HOI_0, aes(x = log_ratio_mean, y = propn, color = model)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = T, linewidth = 1.2) +
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~ sd_ratio_group, nrow = 1) +
  labs(x = 'Relative strength',
       y = 'Living proportion',
       tag = bquote(B[iii] == 0),
       color = 'Model') +
  common_theme +
  scale_color_manual(values = c('Pairwise' = "red",
                                'HOI_lognorm' = 'green',
                                'HOI_norm' = 'black'))

print(plot_lm_HOI_0)

## under B_diag = 1.5 situation
HOI_dataset_1.5 <- HOI_dataset %>% 
  filter(B_diag == 1.5)

lm_HOI_1.5 <- HOI_dataset_1.5 %>%
  mutate(log_ratio_mean = log(B_mean / A_mean),
         sd_ratio = B_sd / A_sd) %>%
  mutate(sd_ratio_group = cut(sd_ratio,
                              breaks = quantile(sd_ratio, probs = c(0, 0.33, 0.66, 1.0), na.rm = TRUE),
                              labels = c('Low ratio', 'Medium ratio', 'High ratio'),
                              include.lowest = TRUE)) %>%
  pivot_longer(cols = c(mean_propn_Pairwise, mean_propn_HOI_norm, mean_propn_HOI_lognorm),
               names_to = 'model',
               values_to = 'propn') %>%
  mutate(model = factor(model, 
                        levels = c('mean_propn_Pairwise', 'mean_propn_HOI_norm', 'mean_propn_HOI_lognorm'),
                        labels = c('Pairwise', 'HOI_norm', 'HOI_lognorm')))

plot_lm_HOI_1.5 <- ggplot(lm_HOI_1.5, aes(x = log_ratio_mean, y = propn, color = model)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = T, linewidth = 1.2) +
  geom_point(size = 0.2, alpha = 0.3) +
  facet_wrap(~ sd_ratio_group, nrow = 1) +
  labs(x = 'Relative strength',
       tag = bquote(B[iii] == 1.5),
       color = 'Model') +
  common_theme +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank()) +
  scale_color_manual(values = c('Pairwise' = "red",
                                'HOI_lognorm' = 'green',
                                'HOI_norm' = 'black'))

print(plot_lm_HOI_1.5)

# creating final plot of gam regression on HOI models
final_plot_gam_seq <- (plot_lm_HOI_0 + plot_lm_HOI_1.5) 

final_plot_gam <- final_plot_gam_seq +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')

print(final_plot_gam)

# creating heatmap for promoting situations
HOI_dataset_0_new   <- HOI_dataset_0   %>% mutate(B_diag_level = 'B_iii = 0')
HOI_dataset_1.5_new <- HOI_dataset_1.5 %>% mutate(B_diag_level = 'B_iii = 1.5')
HOI_dataset_new <- bind_rows(HOI_dataset_0_new, HOI_dataset_1.5_new)

x_breaks <- seq(min(HOI_dataset_new$A_mean), max(HOI_dataset_new$A_mean), length.out = 11)
y_breaks <- seq(min(HOI_dataset_new$B_mean), max(HOI_dataset_new$B_mean), length.out = 11)

HOI_new_grid_summary  <- HOI_dataset_new %>%
  mutate(A_mean_grid = cut(A_mean, breaks = x_breaks, include.lowest = TRUE, labels = FALSE),
         B_mean_grid = cut(B_mean, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)) %>%
  group_by(A_mean_grid, B_mean_grid, B_diag_level) %>%
  summarise(promotion_rate = mean(lognorm_promotes, na.rm = TRUE),
            .groups = 'drop') %>%
  filter(!is.na(A_mean_grid) & !is.na(B_mean_grid))

HOI_new_grid <- expand.grid(center_A_mean = (x_breaks[-length(x_breaks)] + x_breaks[-1]) / 2,
                            center_B_mean = (y_breaks[-length(y_breaks)] + y_breaks[-1]) / 2) %>%
  mutate(A_mean_grid = as.integer(cut(center_A_mean, breaks = x_breaks, include.lowest = TRUE, labels = FALSE)),
         B_mean_grid = as.integer(cut(center_B_mean, breaks = y_breaks, include.lowest = TRUE, labels = FALSE)))

HOI_new_grid_data <- HOI_new_grid %>%
  crossing(B_diag_level = c('B_iii = 1.5', 'B_iii = 0')) %>%
  left_join(HOI_new_grid_summary, by = c('A_mean_grid', 'B_mean_grid', 'B_diag_level'))

facet_labels <- as_labeller(c(`B_iii = 0` = "B[iii]==0", `B_iii = 1.5` = "B[iii]==1.5"), label_parsed)

plot_mean_HOI_promote_lognorm <- ggplot(HOI_new_grid_data,
                                        aes(x = center_A_mean, y = center_B_mean, fill = promotion_rate)) + 
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(name = "Promotion\nRate",
                       colors = c("red", "white", "blue"),
                       values = scales::rescale(c(0, 0.1, 1)),
                       labels = scales::percent,
                       limits = c(0, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = 'Strength of Pairwise interaction', y = 'Strength of HOI interaction') +
  common_theme +
  coord_equal() +
  facet_wrap(~ B_diag_level, labeller = facet_labels) +
  theme(strip.text = element_text(size = 22))

print(plot_mean_HOI_promote_lognorm)