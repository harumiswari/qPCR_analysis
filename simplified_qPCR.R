#library(tidyqpcr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(rstatix)
library(multcomp)
library(ggsignif)


rm(list = ls())
cat("\014") 

qpcr_data <- read.csv("C:/Users/harum/OneDrive/Documents/LAB/mhpQPCR/qpcrcsv/ckmmbqpcr.csv",  stringsAsFactors = FALSE)
head(qpcr_data)

tidy_data <- separate(qpcr_data, Well, into = c("row", "column"), 
                      sep = 1, convert = TRUE)
#tidy_data <- filter(tidy_data, Content != c("CoSCR-4", "CMV-3"))

summarised_data<- tidy_data %>%
  mutate(Cq = as.numeric(Cq)) %>%
  group_by(Content, Target) %>%
  summarise(mean_Cq = mean(Cq))

summarised_data <- summarised_data %>%
  separate(Content, into = c("Content", "replicate"), sep = "-")

test_data <- summarised_data %>%
  filter(Target == "F9")

ref_data <- summarised_data %>%
  filter(Target == "GAPDH") %>%
  rename("ref_Cq" = "mean_Cq")

combined_data <- left_join(test_data, ref_data, by = c("Content", "replicate"))

combined_data <- mutate(combined_data, delta_Cq = mean_Cq - ref_Cq)

treatment_summary <- combined_data %>%
  group_by(Content) %>%
  summarise(mean_delta_Cq = mean(delta_Cq))

mean_control <- filter(treatment_summary, Content == "AAVPrimSCR") %>% pull(mean_delta_Cq)

combined_data <- combined_data %>% 
  mutate(delta_delta_Cq = delta_Cq - mean_control)

combined_data <- combined_data %>%
  mutate(rel_conc = 2^-delta_delta_Cq)

combined_data <- combined_data %>%
  group_by(Content) %>%
  mutate(mean_rel_conc = mean(rel_conc), sd_conc=sd(rel_conc))

#$combined_data$lower <- combined_data$mean_rel_conc - combined_data$sd_conc
#combined_data$upper <- combined_data$mean_rel_conc + combined_data$sd_conc
missing_values <- na.omit(combined_data)

# Calculate SEM for each Content group
sem_data <- combined_data %>%
  group_by(Content) %>%
  summarize(
    n = n(),
    mean_rel_conc = mean(rel_conc),
    sd_conc=sd(rel_conc),
    se = sd_conc/ sqrt(n)
  )


sem_data$lower <- sem_data$mean_rel_conc - sem_data$se
sem_data$upper <- sem_data$mean_rel_conc + sem_data$se

# Ensure data is ungrouped
combined_data <- ungroup(combined_data)

# Perform t-tests
t_test_results <- combined_data %>%
  t_test(rel_conc ~ Content, ref.group = "AAVPrimSCR") %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance() %>%
  mutate(y.position = max(combined_data$rel_conc, na.rm = TRUE) * 1.1)

# Plot the data with error bars and t-test results
ggplot(combined_data, aes(x = Content, y = mean_rel_conc, fill = Content, group = Content)) +
  geom_col(position = "dodge") +
  geom_point(aes(y = rel_conc), position = position_dodge(width = 0.9)) +
  geom_errorbar(data = sem_data, aes(ymin = mean_rel_conc, ymax = upper), position = position_dodge(width = 0.9)) +
  scale_y_log10("Relative concentration") +
  scale_x_discrete("") +
  ggtitle("F9 expression") +
  theme_classic(base_size = 12)
  #stat_pvalue_manual(t_test_results, label = "p.adj.signif", tip.length = 0.01)


ggplot(combined_data, aes(x = Target.x, y = mean_rel_conc, fill = Content, group = Content)) +
  geom_col(position = "dodge") +
  geom_point(aes(y = rel_conc), position = position_dodge(width = 0.9)) +
  scale_y_continuous("Relative concentration", labels = scales::percent) +
  scale_fill_manual("", values = c("CoSCR" = "brown", "CKPL" = "lightpink", "MB" = "green", "CMV" = "blue", "CK8" = "gray", "R26" = "orange")) +
  scale_x_discrete("") +
  ggtitle("F9 expression") +
  coord_cartesian(ylim = c(1, 15)) +  # Use coord_cartesian to set limits without removing data
  theme_classic(base_size = 20)

combined_data$Target.x <- factor(combined_data$Target.x)
# Perform ANOVA
anova_result <- aov(delta_delta_Cq ~ Content, data = combined_data)
tukey_result <- TukeyHSD(anova_result)
summary(anova_result)

df <- as.data.frame(combined_data[combined_data$Content == 'CoSCR',])
dunn <- rstatix::dunn_test(df, mean_rel_conc~Content, detailed = FALSE)
plot1 <- ggbarplot(combined_data, x = "Content", y = "rel_conc", fill = "Target.x", add = "mean_sd", position = position_dodge())+
  stat_compare_means(method="anova", label.x = 1.5, label.y = 6)+
  facet_grid(~Target.x)+
  stat_compare_means(label = "p.signif", method="t.test", ref.group = "Control")+
  stat_pvalue_manual(dunn, label = "p.adj.signif", y.position = 20, step.increase = 0.2, size = 6)+
  scale_y_continuous("Relative concentration") +
  theme(legend.position = "right")
plot1