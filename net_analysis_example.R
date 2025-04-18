
##############################################################################################
##############################################################################################
##############################
############################## Exploring public risk perception of multiple hazards 
############################## through network analysis
##############################
############################## April 2025
##############################
############################## Blanca Paniello-Castillo  
##############################
##############################################################################################
##############################################################################################

# This is an example code for the period 0 (August 2020) of Italy for the perception of likelihood. 
# The same code could be ran changing the variables for other periods and perceptions. 
# The survey data is available open access on Zenodo: https://zenodo.org/records/5653322

#####################
# Install Packages
#####################

install.packages("qgraph")
install.packages("IsingFit")
install.packages("bootnet")
install.packages("usethis")
install.packages("network")

################################################################################
# Clean up Environment
################################################################################

rm(list=ls(all=TRUE))  # Clear entire R workspace

################################################################################
# Load Data
################################################################################

library(readxl)

# change for your path
# Load dataset from Excel file with explicit column type specification
data <- read_excel("you_path/data.xlsx")  

View(data)

################################################################################
# Load Required Libraries
################################################################################

library("usethis")
library("mgm")          # For Mixed Graphical Models
library("devtools")     # To install packages from GitHub
library("ggplot2")      # Plotting
library("bootnet")      # Main package for network estimation and stability
library("IsingFit")     # For Ising model fitting (binary networks)
library("qgraph")       # For network visualization
library("dplyr")        # Data manipulation
library("networktools") # For stability measures

################################################################################
# Create Subsets 
################################################################################

# Rename column if needed (fix encoding issue)
names(data)[names(data) == "ï..country"] <- "country"

# Select variable names for different types of hazards
select_like <- c("like_ep", "like_fl", "like_dr", "like_wf", "like_eq", 
                 "like_ta", "like_dv", "like_ec", "like_cc")

################################################################################
# Examine Missing Values
################################################################################

# Summaries of variables to inspect missingness
summary(data[select_like])

################################################################################
########################## ITALY - PERIOD 0 ANALYSIS ###########################
################################################################################

# Subset Italy Period 0 data and select 'likelihood' variables
data_like <- data[select_like]

# Estimate the network using EBICglasso and Spearman correlation
pred_it0_like_net <- estimateNetwork(data_like, 
                                     default = "EBICglasso", 
                                     corMethod = "spearman", 
                                     tuning = 0.5, 
                                     weighted = TRUE,
                                     threshold = TRUE)

# Define labels for nodes in the network
Names_like_it0 <- c("Epidemics", "Floods", "Droughts", "Wildfires", "Earthquakes", 
                    "Terrorist Attack", "Domestic Violence", "Economic Crisis", "Climate Change")
Labels_it0 <- c("Ep", "Fl", "Dr", "Wf", "Eq", "TA", "DV", "EC", "CC")

# Remove NAs for predictability analysis
pred_it0_like <- data_like %>% 
  dplyr::select(all_of(select_like)) %>%
  na.omit()

# Specify variable types and levels
pred_it0_like_types <- rep("g", 9)      # 'g' = Gaussian
pred_it0_like_levels <- rep("1", 9)     # Continuous data

# Fit MGM model (used for predictability)
set.seed(1234)
predict_like_it0 <- mgm(pred_it0_like,
                        type = pred_it0_like_types, 
                        level = pred_it0_like_levels,
                        labels = Labels_it0,
                        lambdaSel = "EBIC", 
                        lambdaGam = 0.5)

# Compute Predictability (R² values)
pred_it0_like_net_final <- predict(predict_like_it0, pred_it0_like, 
                                   errorCon = c("RMSE", "R2"),
                                   errorCat = c("CC", "nCC"), 
                                   tvMethod = "weighted")

# R² values for each node
pieit0_like <- pred_it0_like_net_final$errors$R2
mean(pieit0_like)  # Mean predictability

# Plot the network
plot1_it0_like <- plot(pred_it0_like_net, layout = "spring", labels = Labels_it0, theme = "colorblind", 
                       nodeNames = Names_like_it0, label.cex = 1.3, legend.cex = .38, 
                       repulsion = 1, 
                       pie = pieit0_like, pieBorder = 0.25, pieColor = "#88CCEE",
                       title = "Likelihood - Italy Period 0")

################################################################################
# Network Stability and Centrality Measures
################################################################################

# Small-worldness index
set.seed(1234)
smallworldness(plot1_it0_like)

# Plot centrality (influence and strength)
centralityPlot(plot1_it0_like, include = c("ExpectedInfluence", "Strength"), weighted = TRUE,
               orderBy = "ExpectedInfluence")

# Compute expected influence
inf_it0_like <- expectedInf(plot1_it0_like)

# Stepwise influence results
inf_it0_like$step1
inf_it0_like$step2

# Plot influence with z-score scaling
plot(inf_it0_like, zscore=TRUE)
plot(inf_it0_like, order="value", zscore=TRUE)

################################################################################
# Bootstrapping for Stability (No Replication)
################################################################################

# Non-parametric bootstrapping (2500 samples)
boot1_norep <- bootnet(pred_it0_like_net, nBoots = 2500)

# Plot edge-weight CIs
plot(boot1_norep, labels=TRUE, order="sample")

# Edge difference test
plot(boot1_norep, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")

# Centrality stability plots
plot(boot1_norep, "strength")

################################################################################
# Bootstrapping with Case Dropping
################################################################################

# Case-dropping bootstrap for robustness of centrality metrics
boot2_norep <- bootnet(pred_it0_like_net, 
                       nBoots = 25, 
                       type = "case", 
                       statistics = c("betweenness", "strength", "closeness"), 
                       nonPositiveDefinite = "continue")

# Plot results
plot(boot2_norep)
corStability(boot2_norep)

# High-powered case-dropping bootstrapping (2500 samples, 9 cores)
bootnet_case_droppingDV <- bootnet(pred_it0_like_net, 
                                   nBoots = 2500,
                                   type = "case",
                                   nCores = 9,
                                   statistics = c("strength", "expectedInfluence"))

plot(bootnet_case_droppingDV, "all")
corStability(bootnet_case_droppingDV)

################################################################################
# End
