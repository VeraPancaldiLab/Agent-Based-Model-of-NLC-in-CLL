library(tidyverse)
library(clusterSim)
library(lattice)
library(cowplot)
library(ggplot2)
library(writexl)
library(numbers)
library(corrplot)
library(reshape2)
library(ggpubr)

################# OUTPUT FILE MANIPULATION ##############
# Load the outputs
parameters_agg <-
  read.table("outputs_ABM_eqFit_1.txt", header = T, sep = ",")

# Remove the duplicates
parameters_agg_rmdups <- parameters_agg %>% unique()

################# FITNESS EVOLUTION ##############

# Plot fitness evolution (max and min value) over the generations

source("~/Documents/ABM/ABM_Netlogo/Analyse_outputs/R_Scripts/fitness_evolution.R")

svg("viability_fitness_evolution.svg", width = 10, height = 6)
ggplot(Viability_evolution, aes(x = Generation, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  xlab("Generation") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

svg("concentration_fitness_evolution.svg", width = 10, height = 6)
ggplot(Concentration_evolution, aes(x = Generation, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  xlab("Generation") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 45))
dev.off()

################# PARAMETER CORRELATION ##############

parameters_agg$evolution.generation <- NULL
parameters_agg$evolution.samples <- NULL

my_data <- as.data.frame(parameters_agg)

svg("Corr_plot_parameters.svg", width = 7, height = 6)
corrplot(
  cor(my_data[, 1:(ncol(my_data)-2)]),
  method = "color",
  type = "full", #upper",
  addCoef.col = "black",
  order = "hclust",
  hclust.method = "complete",
  tl.srt = 45,
  tl.cex = 0.8,
  tl.col = 'black',
  number.cex = 0.65,
  mar = c(0, 0, 0, 0),
  addrect = 3,
  rect.col = 'black',
  rect.lwd = 3
)
dev.off()

################# FITNESS CORRELATION ##############

corr_coef <- cor(my_data) #%>% as.data.frame()

corr_coef_ord <- corr_coef[order(abs(corr_coef[, ncol(corr_coef) - 1])),]


corrplot_obj <- melt(corr_coef_ord[, (ncol(corr_coef_ord)-1):ncol(corr_coef_ord)])
colnames(corrplot_obj) <- c("Parameter", "Fitness", "value")

svg("Corr_plot_fitnesses.svg", width = 7, height = 7)
ggplot(corrplot_obj, aes(Parameter, round(value, digits = 2), fill = Fitness)) +
  facet_wrap( ~ Fitness, nrow = 1) + #place the factors in separate facets
  geom_bar(stat = "identity") + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal
  xlab("Parameter") +
  ylab("Correlation Coefficient") + #improve y-axis label
  theme_bw(base_size = 12) +
  geom_text(aes(label = round(value, digits = 2)), color="black", size = 3) #show the values in the bars
dev.off()

################# PARETO FRONT DISTRIBUTION ##############

pareto_front <-
  read.table("pareto_eqFit1_100.txt", header = T, sep = ",")

## VIOLIN PLOTS WITH GGPLOT2

# Create ggplot object and use factor() to group the variables

ggplot_obj <- melt(pareto_front[, 3:length(pareto_front)])
svg("pareto_front_distribution.svg", width = 8.3, height = 8)
ggplot(ggplot_obj, aes(factor(variable), value, fill = variable)) +
  geom_violin() +
  geom_point() + 
  geom_boxplot(alpha = 0.3,
               color = "black",
               width = .1) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  facet_wrap(~ variable, scales = "free", nrow = 5)
dev.off()

# Load the second output
parameters_agg_2 <-
  read.table("outputs_ABM_10_14.txt", header = T, sep = ",")

source("corrplot_difference.R")
