######################
# Plot the Viability fitness evolution
via_max_vals <-
  parameters_agg %>% group_by(evolution.generation) %>% top_n(1, fitnessVia)
via_max <-
  as.data.frame(cbind(
    via_max_vals$evolution.generation,
    via_max_vals$fitnessVia
  ))
colnames(via_max) <- c("Generation", "Viability_fitness_max_value")

via_min_vals <-
  parameters_agg %>% group_by(evolution.generation) %>% top_n(-1, fitnessVia)
via_min <-
  as.data.frame(cbind(
    via_min_vals$evolution.generation,
    via_min_vals$fitnessVia
  ))
colnames(via_min) <- c("Generation", "Viability_fitness_min_value")

via_ranges <- cbind(via_max, via_min$Viability_fitness_min_value)
colnames(via_ranges) <-
  c("Generation", "Viability fitness max value", "Viability fitness min value")
Viability_evolution <- via_ranges %>% melt(id = "Generation")

ggplot(Viability_evolution, aes(x = Generation, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  xlab("Generation") +
  ylab("Value") +
  ggtitle("Viability fitness evolution") +
  theme(axis.text.x = element_text(angle = 45))


# Plot the Concentration fitness evolution
conc_max_vals <-
  parameters_agg %>% group_by(evolution.generation) %>% top_n(1, fitnessConc)

conc_max <-
  as.data.frame(cbind(
    conc_max_vals$evolution.generation,
    conc_max_vals$fitnessConc
  ))
colnames(conc_max) <- c("Generation", "conc_max_value")
conc_max <- conc_max %>% distinct(Generation, .keep_all = TRUE)

conc_min_vals <-
  parameters_agg %>% group_by(evolution.generation) %>% top_n(-1, fitnessConc)
conc_min <-
  as.data.frame(cbind(
    conc_min_vals$evolution.generation,
    conc_min_vals$fitnessConc
  ))
colnames(conc_min) <- c("Generation", "conc_min_value")

conc_ranges <- cbind(conc_max, conc_min$conc_min_value)
colnames(conc_ranges) <-
  c("Generation", "Concentration fitness max value", "Concentration fitness min value")
Concentration_evolution <- conc_ranges %>% melt(id = "Generation")

ggplot(Concentration_evolution, aes(x = Generation, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  xlab("Generation") +
  ylab("Value") +
  ggtitle("Concentration fitness evolution") +
  theme(axis.text.x = element_text(angle = 45))
