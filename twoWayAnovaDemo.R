

#  demonstration of using 2-factor ANOVA
# Firstly with only additive effects
# Secondly after including an interaction.


groupsA <- c("me1", "me2", "me3")
treatment <- c("control", "KO")


# simulate some data with no interaction
set.seed(489723)    # there is a small chance the null data will lead to a significant result, if so, we can change the seed to reproducably produce "non-significant" data

statsData1 <- data.frame(group = rep(groupsA, each=100),
                        treatment="control",
                        metric = rnorm(300, mean=.07, sd=.001))

# add a small amount to set me2
statsData1[statsData1$group == "me2", "metric"] <- statsData1[statsData1$group == "me2", "metric"]  + rnorm(100, mean=.01, sd=.001)
# add a small amount to set me3
statsData1[statsData1$group == "me3", "metric"] <- statsData1[statsData1$group == "me3", "metric"]  + rnorm(100, mean=.025, sd=.001)
# check the data:-
boxplot(metric ~ group, data=statsData1)
# now copy and shift the means by a fixed amount (additive model, no interaction)
statsData2 <- statsData1
statsData2$treatment <- "KO"
statsData2$metric <- statsData2$metric + rnorm(300, mean=.02, sd=.001)

# combine both data-sets into one data.frame.
statsDataAll <- rbind(statsData1, statsData2)
boxplot(metric ~ treatment + group , data=statsDataAll , col=c("palegreen", "lightblue"))

twoWayAdd <- lm(metric ~ treatment * group , data=statsDataAll)
summary(twoWayAdd)
#  Should have the following ouput.       "<-  my annotation"
# # Call:
# lm(formula = metric ~ treatment * group, data = statsDataAll)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0049193 -0.0008217  0.0000274  0.0008668  0.0036702 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           0.0701374  0.0001363 514.656   <2e-16 ***       <- is the mean of all data significantly different from zero.  Boring but true.
#   treatmentKO           0.0200769  0.0001927 104.172   <2e-16 ***     <- is the KO treatment signficantly different from the control.  YES
#   groupme2              0.0101376  0.0001927  52.600   <2e-16 ***     <- is group me2 signficantly different from me1 - YES
#   groupme3              0.0248428  0.0001927 128.900   <2e-16 ***     <- is group me3 signficantly different from me1  - YES
#   treatmentKO:groupme2 -0.0001879  0.0002726  -0.689    0.491         <- is the effect of KO dependent on me2   (an interaction )  - NO 
# treatmentKO:groupme3 -0.0001694  0.0002726  -0.622    0.534           <- is the effect of KO dependent on me3   (an interaction )  - NO 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.001363 on 594 degrees of freedom
# Multiple R-squared:  0.991,	Adjusted R-squared:  0.9909               <- what proportion of variance is explained by the model. Should be close to 1.0 here as we determined the main effects with very low random scatter.
# F-statistic: 1.311e+04 on 5 and 594 DF,  p-value: < 2.2e-16

#  Section 2 - data with an interaction
#  NOW edit the data to include an interaction.

# now copy and shift the means by a different amount for each group ( interaction)
statsData3 <- statsData2
statsData3$treatment <- "KO"
statsData3[statsData3$group == "me2", "metric"] <- statsData3[statsData3$group == "me2", "metric"]  + rnorm(100, mean=.01, sd=.001)
# we already added an amount to all the KO samples in statsData2.  
# Don't do that again
#statsData3$metric <- statsData3$metric + rnorm(300, mean=.02, sd=.001)
# Instead, add a different amount to each group WITHIN the KO group
statsData3[statsData3$group == "me2", "metric"] <- statsData3[statsData3$group == "me2", "metric"]  + rnorm(100, mean=.01, sd=.001)
# add a small amount to set me3
statsData3[statsData3$group == "me3", "metric"] <- statsData3[statsData3$group == "me3", "metric"]  + rnorm(100, mean=.025, sd=.001)

# combine both data-sets into one data.frame.
statsDataInt <- rbind(statsData1, statsData3)
boxplot(metric ~ treatment + group , data=statsDataInt , col=c("palegreen", "lightblue"))
# they should now all show differences, but the differences should change for each group.

twoWayInt <- lm(metric ~ treatment * group , data=statsDataInt)
summary(twoWayInt)
# this time, the interaction effects "treatmentKO:groupme2" and "treatmentKO:groupme3" should also be significant.










# Reading the two CSV files
wild_type <- read.csv("/Users/maryyuhanna/desktop/simulation_results_ks_100.csv")
ks_sim <- read.csv("/Users/maryyuhanna/desktop/simulation_results_ksKMT2D_100.csv")

# Combine the two dataframes if they have the same structure
combined_results_df <- rbind(wild_type, ks_sim)




# Check the structure of the combined dataframe
str(combined_results_df)

# View the first few rows
head(combined_results_df)




# Add a 'Simulation' column to identify the simulation type
wild_type$Simulation <- "WildType"
ks_sim$Simulation <- "KMT2D_Knockout"

# Combine the dataframes
combined_results_df <- rbind(wild_type, ks_sim)



# Check the structure and first few rows of the combined dataframe
str(combined_results_df)
head(combined_results_df)





# Assuming that your combined_results_df contains the relevant data for the analysis

# Two-way ANOVA for Specificity
specificity_anova <- lm(Specificity ~ Layer * Simulation, data = combined_results_df)
summary(specificity_anova)

# Two-way ANOVA for Sensitivity
sensitivity_anova <- lm(Sensitivity ~ Layer * Simulation, data = combined_results_df)
summary(sensitivity_anova)

# Two-way ANOVA for Accuracy
accuracy_anova <- lm(Accuracy ~ Layer * Simulation, data = combined_results_df)
summary(accuracy_anova)

# Two-way ANOVA for Percentage
percentage_anova <- lm(Percentage ~ Layer * Simulation, data = combined_results_df)
summary(percentage_anova)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("MotifDb")



citation("GenomicRanges")

