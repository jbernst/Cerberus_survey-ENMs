#############################################################
#                                                           #
# Cerberus schneiderii linear models and significance tests #
#                                                           #
#############################################################
# Linear models and significance tests for Cerberus schneiderii in the Philippines
# Code by: Justin M. Bernstein
# Date: 2/8/2024


# Load libraries
library(tidyverse)
library(tidyquant)
library(ggdist)
library(ggthemes)

# set working directory
setwd("D:/Documents/Publications/Cerberus-Autecology/")

# import data. Note, 0, 0.5, and 1 represent the following:
# 0   = >5 km from mangrove
# 0.5 = 1-5 km from mangrove
# 1   = <1 km from mangrove 

cerb <- read.csv("cerb_coordinates_suitability_min-head.csv")
View(cerb)

### Set up a linear model for habitat suitability and the distance of occurrence records
# to the coast

lm.coast <- lm(suitability~km2coast, data = cerb)
summary(lm.coast)
# Call:
#   lm(formula = suitability ~ km2coast, data = cerb)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.59071 -0.08183  0.03548  0.14733  0.49014 
# 
# Coefficients:
#              Estimate Std.  Error    t value   Pr(>|t|)    
# (Intercept)    0.769337   0.016014  48.042    <2e-16 ***
#   km2coast    -0.008992   0.001973  -4.557   8.66e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1913 on 217 degrees of freedom
# Multiple R-squared:  0.08734,	Adjusted R-squared:  0.08313 
# F-statistic: 20.77 on 1 and 217 DF,  p-value: 8.663e-06

rsq.coast <- summary(lm.coast)$r.squared

#create scatterplot with fitted regression line
ggplot(lm.coast, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb$km2coast), y = max(cerb$suitability), 
           label = paste("R-squared =", round(rsq.coast, 3)), hjust = 1, vjust = 1)

### Set up a linear model for habitat suitability and the distance of occurrence records
# to the mangroves

lm.mang <- lm(suitability~km2man, data = cerb)
summary(lm.mang)

# Call:
#   lm(formula = suitability ~ km2man, data = cerb)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62583 -0.08889  0.04616  0.18734  0.29890 
# 
# Coefficients:
#               Estimate   Std. Error  t value    Pr(>|t|)    
# (Intercept)   0.753641    0.017927   42.040    <2e-16 ***
#   km2man      -0.004333   0.001891   -2.292     0.0229 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1979 on 217 degrees of freedom
# Multiple R-squared:  0.02364,	Adjusted R-squared:  0.01914 
# F-statistic: 5.253 on 1 and 217 DF,  p-value: 0.02287

rsq.mang <- summary(lm.mang)$r.squared

#create scatterplot with fitted regression line
ggplot(lm.mang, aes(x = km2man, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb$km2man), y = max(cerb$suitability), 
           label = paste("R-squared =", round(rsq.mang, 3)), hjust = 1, vjust = 1)

# Raincloud plot to show suitability estimates in relation to mangrove categories.

cerb %>% 
  filter(mangroves %in% c(0, 0.5, 1)) %>% 
  ggplot(aes(x = factor(mangroves), y = suitability, fill = factor(mangroves))) +
  labs(fill = "Mangrove Presence") +
  scale_fill_manual(values=c("#ccbc6e", "#56B4E9", "#2da146")) +
  # add half-violin from {ggdist} package
  stat_halfeye(
    # adjust bandwidth
    adjust = 0.65,
    # move to the right
    justification = -0.2,
    # remove the slub interval
    .width = 0,
    point_colour = NA,
  ) +
  geom_boxplot(
    width = 0.12,
    # removing outliers
    outlier.color = NA,
    alpha = 0.5
  ) +
  stat_dots(
    # ploting on left side
    side = "left",
    # adjusting position
    justification = 1.1,
    # adjust grouping (binning) of observations
    binwidth = 0.01
  ) +
  # Themes and Labels
  theme_tq() +
  labs(
    title = "Suitability to Mangrove Plots",
    x = "Mangrove Category",
    y = "Habitat Suitability", 
    
  ) +
  coord_flip()

### Test for significant between suitability values in mangrove categories of 0, 0.5, and 1

# Make separate mangrove datasets
mangrove_0 <- subset(cerb, mangroves == 0)
mangrove_05 <- subset(cerb, mangroves == 0.5)
mangrove_1 <- subset(cerb, mangroves == 1)

# Make one for combining categories 0.5 and 1 using rbind()

mangroves_05_1 <- rbind(mangrove_05, mangrove_1)

# Perform the t-tests
t_test_0x05  <- t.test(mangrove_0$suitability, mangrove_05$suitability)
t_test_0x1   <- t.test(mangrove_0$suitability, mangrove_1$suitability)
t_text_05x1  <- t.test(mangrove_05$suitability, mangrove_1$suitability)
t_test_0_051 <- t.test(mangrove_0$suitability, mangroves_05_1$suitability)

mangrove.color <- c("#ccbc6e", "#56B4E9", "#2da146", "black")

boxplot(
  mangrove_0$suitability,
  mangrove_05$suitability,
  mangrove_1$suitability,
  mangroves_05_1$suitability,
  col = c("#ccbc6e", "#56B4E9", "#2da146", "#008080"),
  names = c("Inland", "Inland-Coastal", "Coastal", "Inland-Coastal+Coastal"),
  cex.axis = 0.9
)

# View the results

# inland vs. inland-coastal
print(t_test_0x05) # t = -2.6144, df = 105.05, p-value = 0.01025

# inland vs. coastal
print(t_test_0x1) # t = -6.6093, df = 167.91, p-value = 4.908e-10

# inland-coastal vs. coastal
print(t_text_05x1) # t = -3.1486, df = 88.965, p-value = 0.002234

# inland vs. (inland-coastal + coastal)
print(t_test_0_051) # t = -4.879, df = 150.5, p-value = 2.689e-06

cerb.labeled <- cerb %>%
  mutate(habitat.type = case_when(
    mangroves == 0.0 ~ "coastal",
    mangroves == 0.5 ~ "coastal-inland",
    mangroves == 1.0 ~ "inland"
  ))

aov.one.way <- aov(suitability ~ habitat.type, data = cerb.labeled)

summary(aov.one.way)
#                Df   Sum Sq  Mean     Sq F value   Pr(>F)    
# habitat.type   2    1.478   0.7390   22.09        1.87e-09 ***
# Residuals      216  7.226   0.0335                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

TukeyHSD(aov.one.way)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = suitability ~ habitat.type, data = cerb.labeled)
# 
# $habitat.type
#                        diff       lwr         upr       p adj
# coastal-inland-coastal 0.08699741 0.008763378 0.1652314 0.0251146
# inland-coastal         0.18500734 0.119318862 0.2506958 0.0000000
# inland-coastal-inland  0.09800993 0.018667890 0.1773520 0.0109214

########################################################################################
##
## Below we will remove the samples that had extreme distances from the coast (> 30 km2)
##
########################################################################################


cerb.no.out <- read.csv("cerb_coordinates_suitability_min-head.csv")
cerb.no.out <- cerb.no.out %>%
  filter(km2coast < 30)

tail(cerb.no.out)


lm.coast.no.out <- lm(suitability~km2coast, data = cerb.no.out)
summary(lm.coast.no.out)
# Call:
#   lm(formula = suitability ~ km2coast, data = cerb.no.out)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.55799 -0.08702  0.01866  0.12676  0.31460 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.792578   0.016119  49.169  < 2e-16 ***
#   km2coast    -0.015267   0.002319  -6.584 3.45e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1832 on 215 degrees of freedom
# Multiple R-squared:  0.1678,	Adjusted R-squared:  0.1639 
# F-statistic: 43.35 on 1 and 215 DF,  p-value: 3.454e-10

rsq.coast.no.out <- summary(lm.coast.no.out)$r.squared

#create scatterplot with fitted regression line
pdf("km2coast-x-suitability_lm_no-outlier.pdf")
ggplot(lm.coast.no.out, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb.no.out$km2coast), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(rsq.coast.no.out, 3)), hjust = 1, vjust = 1)
dev.off()

### Set up a linear model for habitat suitability and the distance of occurrence records
# to mangroves

lm.mang.no.out <- lm(suitability~km2man, data = cerb.no.out)
summary(lm.mang.no.out)

# Call:
#   lm(formula = suitability ~ km2man, data = cerb.no.out)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.62076 -0.08283  0.03303  0.16972  0.25591 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.767886   0.018620  41.241  < 2e-16 ***
#   km2man      -0.007184   0.002185  -3.288  0.00118 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1959 on 215 degrees of freedom
# Multiple R-squared:  0.04788,	Adjusted R-squared:  0.04345 
# F-statistic: 10.81 on 1 and 215 DF,  p-value: 0.001178

rsq.mang.no.out <- summary(lm.mang.no.out)$r.squared

#create scatterplot with fitted regression line
# pdf("km2man-x-suitability_lm_no-outlier.pdf")
ggplot(lm.mang.no.out, aes(x = km2man, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb.no.out$km2man), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(rsq.mang.no.out, 3)), hjust = 1, vjust = 1)
# dev.off()


cerb.labeled.no.out <- cerb.no.out %>%
  mutate(habitat.type = case_when(
    mangroves == 0.0 ~ "coastal",
    mangroves == 0.5 ~ "coastal-inland",
    mangroves == 1.0 ~ "inland"
  ))

aov.one.way.no.out <- aov(suitability ~ habitat.type, data = cerb.labeled.no.out)

summary(aov.one.way.no.out)
#                Df   Sum Sq  Mean Sq  F value   Pr(>F)    
# habitat.type   2    1.543   0.7715   23.18     7.76e-10 ***
# Residuals      214  7.124   0.0333                     
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

TukeyHSD(aov.one.way.no.out)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = suitability ~ habitat.type, data = cerb.labeled.no.out)
# 
# $habitat.type
#                        diff       lwr        upr       p adj
# coastal-inland-coastal 0.09205783 0.01370806 0.1704076 0.0165781
# inland-coastal         0.19006776 0.12417719 0.2559583 0.0000000
# inland-coastal-inland  0.09800993 0.01885417 0.1771657 0.0107022


# Make separate mangrove datasets
mangrove_0 <- subset(cerb.no.out, mangroves == 0)
mangrove_05 <- subset(cerb.no.out, mangroves == 0.5)
mangrove_1 <- subset(cerb.no.out, mangroves == 1)

# Make one for combining categories 0.5 and 1 using rbind()

mangroves_05_1 <- rbind(mangrove_05, mangrove_1)

# Perform the t-tests
t_test_0x05  <- t.test(mangrove_0$suitability, mangrove_05$suitability)
t_test_0x1   <- t.test(mangrove_0$suitability, mangrove_1$suitability)
t_text_05x1  <- t.test(mangrove_05$suitability, mangrove_1$suitability)
t_test_0_051 <- t.test(mangrove_0$suitability, mangroves_05_1$suitability)

mangrove.color <- c("#ccbc6e", "#56B4E9", "#2da146", "black")

# Plot boxplots
# pdf("Habitat-Boxplots_no-outliers.pdf")
boxplot(
  mangrove_0$suitability,
  mangrove_05$suitability,
  mangrove_1$suitability,
  mangroves_05_1$suitability,
  col = c("#ccbc6e", "#56B4E9", "#2da146", "#008080"),
  names = c("Inland", "Inland-Coastal", "Coastal", "Inland-Coastal+Coastal"),
  cex.axis = 0.9
)
# dev.off()





########################################################################################
## 
## TRYING OTHER MODELS
##
########################################################################################


# no interaction term (i.e., distance from mangroves does not depend on distance from coast)
lm.coast.mang.no.out <- lm(suitability ~ km2man + km2coast, data = cerb.no.out)
summary(lm.coast.mang.no.out)

# Call:
# lm(formula = suitability ~ km2man + km2coast, data = cerb.no.out)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.53701 -0.10799  0.00424  0.12736  0.29739 
# 
# Coefficients:
#                Estimate Std. Error     t value    Pr(>|t|)    
# (Intercept)    0.777574     0.017322   44.890   < 2e-16 ***
#   km2man       0.006905     0.003086   2.237      0.0263 *  
#   km2coast    -0.021185     0.003504  -6.047      6.51e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1815 on 214 degrees of freedom
# Multiple R-squared:  0.1868,	Adjusted R-squared:  0.1792 
# F-statistic: 24.58 on 2 and 214 DF,  p-value: 2.456e-10


# with an interaction term (i.e., distance from mangroves depends on distance from coast)
lm.coast.mang.no.out.int <- lm(suitability ~ km2man * km2coast, data = cerb.no.out)
summary(lm.coast.mang.no.out.int)
# Call:
#   lm(formula = suitability ~ km2man * km2coast, data = cerb.no.out)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.52052 -0.09738  0.01797  0.13296  0.29108 
# 
# Coefficients:
#                    Estimate Std. Error      t value   Pr(>|t|)    
# (Intercept)        0.7970802     0.0200845  39.686   < 2e-16 ***
#   km2man           0.0054887     0.0031583   1.738    0.0837 .  
# km2coast          -0.0324793     0.0069255  -4.690    4.88e-06 ***
#   km2man:km2coast  0.0007624     0.0004041   1.887    0.0606 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 0.1804 on 213 degrees of freedom
# Multiple R-squared:  0.2002,	Adjusted R-squared:  0.1889 
# F-statistic: 17.77 on 3 and 213 DF,  p-value: 2.482e-10


# plot the linear model without the interaction term
rsq.man.coast.no.out <- summary(lm.coast.mang.no.out)$r.squared

pdf("km2coast+km2mang-x-suitability_lm_no-outlier.pdf")
ggplot(lm.coast.mang.no.out, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb.no.out$km2coast), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(rsq.man.coast.no.out, 3)), hjust = 1, vjust = 1)
dev.off()

# plot the linear model with the interaction term
rsq.man.coast.no.out.int <- summary(lm.coast.mang.no.out.int)$r.squared

pdf("km2coast-interaction-km2mang-x-suitability_lm_no-outlier.pdf")
ggplot(lm.coast.mang.no.out.int, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm") + 
  annotate("text", x = max(cerb.no.out$km2coast), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(rsq.man.coast.no.out.int, 3)), hjust = 1, vjust = 1)
dev.off()

##
## Can try other models too, but while the R squared may increase, they do not make much sense biologically
##
# Fit Generalized Additive Models (GAM) 
gam.coast.mang.no.out <- gam(suitability ~ s(km2coast) + km2man, data = cerb.no.out)


summary(gam.coast.mang.no.out)
rsq.man.coast.no.out.gam <- summary(gam.coast.mang.no.out)$r.sq

# Plot the GAM model
ggplot(cerb.no.out, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "gam", formula = y ~ s(x), color = "red") +
  annotate("text", x = max(cerb.no.out$km2coast), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(summary(gam.coast.mang.no.out)$r.sq, 3)), 
           hjust = 1, vjust = 1) +
  labs(title = "Generalized Additive Model",
       x = "km2coast",
       y = "Suitability")



# Fit an exponential model
lm.coast.mang.no.out.exp <- lm(log(suitability) ~ km2coast + km2man, data = cerb.no.out)

summary(lm.coast.mang.no.out.exp)
rsq.man.coast.no.out.exp <- summary(lm.coast.mang.no.out.exp)$r.squared


# Plot the quadratic model
ggplot(cerb.no.out, aes(x = km2coast, y = suitability)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ exp(x), color = "purple") +
  annotate("text", x = max(cerb.no.out$km2coast), y = max(cerb.no.out$suitability), 
           label = paste("R-squared =", round(summary(lm.coast.mang.no.out.exp)$r.squared, 3)), 
           hjust = 1, vjust = 1) +
  labs(title = "Exponential Regression",
       x = "km2coast",
       y = "Suitability")

