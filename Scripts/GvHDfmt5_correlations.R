### Description: Correlations                             ###
### Author: Anna Li, date: 15062026                       ###
### Email: anna.li@adelaide.edu.au                        ###

## Load Libraries
library(ggplot2)          # Required for graphing
library(ggpubr)           # Required for stat_compare_means when graphing shannon
library(rmcorr)           #for repeated measures correlations
library(dplyr)
library(tidyr)
library(lme4)

## Set up folders
project_dir <- here::here()
results_dir <- file.path(project_dir, "Results")

#use Recipient_data.df from script GvHDfmt5_stats-plots.R for relevant samples 
#(no cecal samples, no donor samples, no mice died from procedure = correct timepoints)

#calculate change in each factor to our Recipient_data.df dataframe (lag model?)
delta_df <- Recipient_data.df %>%
  arrange(Mouse_ID, Timepoint_ordered) %>%
  group_by(Mouse_ID) %>%
  mutate(
    delta_shannon = Shannon - lag(Shannon),
    delta_weight = Corrected_weight - lag(Corrected_weight),
    delta_dist = dist_to_donor_a - lag(dist_to_donor_a),
    delta_GvHD = GvHD_score - lag(GvHD_score)
  ) %>%
  filter(!is.na(delta_dist))
#remember that delta_df has taken out S0 as a consequence

#Taking a look at the data overall
## mean trajectories of variables (does not factor group)
Recipient_data.df %>%
  group_by(Timepoint_ordered) %>%
  summarise(
    mean_weight = mean(Corrected_weight, na.rm = TRUE),
    mean_GvHD = mean(GvHD_score, na.rm = TRUE),
    mean_shannon = mean(Shannon, na.rm = TRUE),
    mean_distance = mean(dist_to_donor_a, na.rm = TRUE)
  )

# run simple correlation stat test to decide what to plot
#Q: Are donor similarity and GvHD score correlated/associated overall
cor.test(
  Recipient_data.df$dist_to_donor_a,
  Recipient_data.df$GvHD_score,
  method = "spearman"
)

#run repeated measures correlation test
#Q:for a given mouse, does distance to donor and GvHD score associate with each other?
#testing in-between mice

rm<-rmcorr(
  participant = Mouse_ID,
  measure1 = dist_to_donor_a,
  measure2 = GvHD_score,
  dataset = Recipient_data.df
)
rm
plot(rm)
#but these results suggest as dist to donor decreases, GvHD score increases (accumulation / time may need to be factored)

# now run fixed mixed model - does association remain after adjusting for time?
#Q: at a given timepoint, are mice that are closer to the donor in distance have lower GvHD scores?

lmer_GvHDdist<- lmer(
  GvHD_score ~ dist_to_donor_a + Timepoint_ordered +
    (1|Mouse_ID),
  data = Recipient_data.df
)

summary(lmer_GvHDdist)
#A: Not significant...so accumulation does need to be factored

#Interrogate two ways:
# 1) within and 2) between mice effects
# 2) Calculates each mouse's average distance to donor as an overall value 
#This accounts for distances changing over time individually as more of a single accumulative readout

Recipient_data.df <- Recipient_data.df %>%
  group_by(Mouse_ID) %>%
  mutate(
    mean_dist = mean(dist_to_donor_a, na.rm = TRUE),
    within_dist = dist_to_donor_a - mean_dist
  ) %>%
  ungroup()

#Q: When a mouse becomes more donor-like than its own average, does GvHD decrease?
#Q: Does mean distance of a mouse from the donor associate with GvHD score?

library(lmerTest)
m <- lmer(
  GvHD_score ~ within_dist +
    mean_dist +
    Timepoint_ordered +
    (1|Mouse_ID),
  data = Recipient_data.df
)

summary(m)

#CONCLUSIONS
#within_dist - estimate = 0.0227, p = 0.333 
#--> when a mouse becomes more donor like than its own average it is not associated with a lower GvHD score
# mean_dist - estimate = 0.0695, p = 0.077
#-> almost significant
#--> mice that are generally further away from the donor microbiome tend to have higher GvHD scores
#---> overall between-mouse differences matter more than within-mouse fluctuations in distance

#new framework of questioning: does donor similarity at earlier timepoint predict future GvHD severity?
#delta_future_GvHD is the change in GvHD score from one timepoint to the next in order
future_df <- Recipient_data.df %>%
  arrange(Mouse_ID, Timepoint_ordered) %>%
  group_by(Mouse_ID) %>%
  mutate(
    next_GvHD = lead(GvHD_score),
    delta_future_GvHD = next_GvHD - GvHD_score
  ) %>%
  filter(!is.na(delta_future_GvHD))

#subset of relevant times if necessary====
future_subset <- future_df %>%
  filter(Timepoint_ordered %in% c("D0", "D15", "F7", "F14", "F21"))
#====

#where delta_future_GvHD = next GvHD score -(subtract)- current GvHD score
#positive = disease worsens
#negative = disease improves

#Q: does microbiome distance at the current timepoint predict how much the GvHD score changes by the next timepoint?
library(lmerTest)

delta_model <- lmer(
  delta_future_GvHD ~ dist_to_donor_a +
    (1|Mouse_ID),
  data = future_df
)

summary(delta_model)

#range of distances:
summary(Recipient_data.df$dist_to_donor_a)

#CONCLUSION
#distance_to_donor_a - estimate: 0.06421, p=0.000251***
#->great distance gives positive value 
#--> mice that are further from the donor microbiome tend to experience larger increases in GvHD score by the next timepoint
#---> mice that are more donor like tend to show less subsequent worsening of GvHD

#----PLOTTING----

#future_df
future_plot<-ggplot(future_df, aes(x= dist_to_donor_a,
                      y= delta_future_GvHD))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")+
  geom_point(alpha = 0.7)+
  geom_smooth(method = "lm", se=TRUE)+
  theme_classic(base_size = 14)+
  labs(x = "Aitchison distance to LPJ microbiome",
       y = expression(Delta*"GvHD score (next - current)"))
future_plot

future_plot_group<-ggplot(future_df, aes(x= dist_to_donor_a,
                                   y= delta_future_GvHD,
                                   color = Group))+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")+
  geom_point(alpha = 0.7)+
  geom_smooth(method = "lm", se=TRUE)+
  scale_color_manual(values = group_cols) +
  theme_classic(base_size = 14)+
  labs(x = "Aitchison distance to LPJ microbiome",
       y = expression(Delta*"GvHD score (next - current)"))
future_plot_group

#----save plots----
save_plot <- function(plot_obj, filename, width = 8, height = 6) {
  ggsave(
    filename = file.path(results_dir, filename),
    plot = plot_obj,
    width = width,
    height = height
  )
}

save_plot(future_plot, "future_plot.pdf",width = 4, height = 4)
save_plot(future_plot_group, "future_plot_coloured.pdf",width = 6, height = 4)
