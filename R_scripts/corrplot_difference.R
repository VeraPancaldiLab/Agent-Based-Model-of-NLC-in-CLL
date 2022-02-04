
# Rename parameters (should be same as the output 1)
colnames(parameters_agg_2)[which(colnames(parameters_agg_2) == "adherenceTime")] <-
  "DifferentiationTime"
colnames(parameters_agg_2)[which(colnames(parameters_agg_2) == "adherenceInitStd")] <-
  "DifferentiationInitStd"

parameters_agg_2$evolution.generation <- NULL
parameters_agg_2$evolution.samples <- NULL

parameters_agg_rmdubs_2 <- parameters_agg_2 %>% distinct()

my_data_2 <- as.data.frame(parameters_agg_rmdubs_2)
# Plot the matrix of the differences between 2 correlation matrixes 
corrplot(
  cor(my_data) - cor(my_data_2),
  method = "color",
  type = "upper",
  addCoef.col = "black",
  tl.srt = 45,
  tl.cex = 0.9,
  tl.col = 'black',
  #  main = "Correlation matrix output_ABM_10_5 (30 days)",
  #  mar = c(0, 0, 2, 0),
  # order = 'hclust',
  addrect = 3,
  #  diag = FALSE,
  rect.col = 'black',
  rect.lwd = 3
)