# ------------------------------------------------------------------------------
# Visualise results
# ------------------------------------------------------------------------------

# load required packages -------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)

# read in results df -----------------------------------------------------------
df_wave1 <- readRDS("inst/extdata/results/wave1_results.rds")
df_wave2 <- readRDS("inst/extdata/results/wave2_results.rds")
df_wave3 <- readRDS("inst/extdata/results/wave3_results.rds")

# create a single data frame
df_all <- bind_rows(df_wave1, df_wave2, df_wave3, .id = "wave") %>%
  mutate(target_variable = factor(target_variable,
                                  levels = c("inc infection", "inc hosp", "occ hosp", 
                                             "inc icu", "occ icu", "inc death")),
         sample = factor(sample),
         wave = factor(wave))

# Line plot --------------------------------------------------------------------
# plot each sample individually
# first sum over age groups for each target_variable
df_sum_ag <- df_all %>%
  group_by(wave, scenario_id, target_variable, date, sample) %>%
  summarise(sum = sum(value)) %>%
  ungroup()

# plot
p_lines <- ggplot(data = df_sum_ag, aes(x = date, y = sum, group = sample, 
                                     color = scenario_id)) +
  geom_line() +
  facet_gird(target_variable~wave)

p_lines

# Average with confidence bounds -----------------------------------------------
df_summary <- df_all %>%
  group_by(wave,scenario_id, target_variable, date) %>%
  summarise(mean  = median(sum),
            q025 = quantile(sum, probs = 0.025),
            q25  = quantile(sum, probs = 0.25),
            q75  = quantile(sum, probs = 0.75),
            q975 = quantile(sum, probs = 0.975)
            ) %>%
  select(date, scenario_id, target_variable, mean:q975) 

# plot 
p_ribbon <- ggplot(data = df_summary %>%
                   filter(target_variable %in% c("occ hosp"),
                          scenario_id != "A-Wave1"
                  #        #date < as.Date("2022-10-01")
                         ), # "inc hosp", "inc icu", "inc death"
                 aes(x = date, y = mean, color = scenario_id, fill = scenario_id)) +
  geom_line() +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.1, color = NA) +
  #geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
  xlab("Date") + 
  ylab("Mean value") +
  scale_x_date(date_breaks = "1 month", date_labels = "%d %b %Y") +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.box="vertical",
        axis.title=element_text(size=14,face="bold")) +
  guides(fill=guide_legend("Scenario ID"), colour = guide_legend("Scenario ID")) +
  facet_grid(.~target_variable, scales = "free") #+
  # annotate("rect", xmin = as.Date("2022-09-22"), xmax = as.Date("2022-12-15"), ymin = 0, ymax = 200000, 
  #          alpha = .5)
  #geom_vline(xintercept = as.Date("2022-09-15"), linetype = "dashed", color = "grey70")
p_ribbon


# save outputs -----------------------------------------------------------------
ggsave(filename = "/rivm/s/ainsliek/results/scenario_hub/round2/case_plot_round2.jpg", 
       plot = p_ribbon,
       units = "in", height = 8, width = 13, dpi = 300)

