library(tidyverse)
library(ggplot2)
library(reshape2)


pareto_front_1 <-
  read.table("pareto10_15.txt", header = T, sep = ",")
pareto_front_2 <-
  read.table("pareto10_16.txt", header = T, sep = ",")
pareto_front_3 <-
  read.table("pareto10_17.txt", header = T, sep = ",")
pareto_front_4 <-
  read.table("pareto10_18.txt", header = T, sep = ",")

## VIOLIN PLOTS WITH GGPLOT2

pareto_front_1$delta_fitness_via <- NULL
pareto_front_1$delta_fitness_conc <- NULL

pareto_front_2$delta_fitness_via <- NULL
pareto_front_2$delta_fitness_conc <- NULL

pareto_front_3$delta_fitness_via <- NULL
pareto_front_3$delta_fitness_conc <- NULL

pareto_front_4$delta_fitness_via <- NULL
pareto_front_4$delta_fitness_conc <- NULL

pareto_front_1$Output <- "Output_1"
pareto_front_2$Output <- "Output_2"
pareto_front_3$Output <- "Output_3"
pareto_front_4$Output <- "Output_4"


several_paretos <- rbind(pareto_front_1, pareto_front_2, pareto_front_3, pareto_front_4)
# Create ggplot object and use factor() to group the variables

ggplot_obj <- melt(several_paretos)


ggplot(ggplot_obj, aes(factor(variable), value, fill = Output)) +
  geom_violin() +
  geom_point(position = position_dodge(width = 0.9)) + 
  geom_jitter(position = position_dodge(width = 0.9)) +
  geom_boxplot(position = position_dodge(width = 0.9),
               alpha = 0.3,
               color = "black",
               width = .1) +
  labs(x = "", y = "", title = "Parameter distribution in the Pareto fronts 10_15-18") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    #    legend.title = element_blank(),
    #    legend.position = "none"
  ) +
  facet_wrap(~ variable, scales = "free", nrow = 5)
