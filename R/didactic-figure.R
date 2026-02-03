## making committee meeting didactic figure 
library(tidyverse)
theme_set(theme_bw())

## source function for generating noise 
source("R/functions/generate_noise_plotting.R")

## generate two time series 
beta_0 = generate_noise(beta = 0,
                        p = 0,
                        n_ts = 1, 
                        L1 = 500,
                        L2 = 365*5,
                        sigma = 0.5)
beta_1 = generate_noise(beta = 1,
                       p = 0,
                       n_ts = 1, 
                       L1 = 500,
                       L2 = 365*5,
                       sigma = 0.5)

## save data 
# saveRDS(beta_0, "outputs/figures/didactic/beta_0.rds")
# saveRDS(beta_1, "outputs/figures/didactic/beta_1.rds")

beta_0 = readRDS( "outputs/figures/didactic/beta_0.rds")
beta_1 = readRDS( "outputs/figures/didactic/beta_1.rds")

beta_0_df = data.frame(time = 1:1500, signal =  unlist(beta_0[[1]])[501:2000], beta = 0)
beta_1_df = data.frame( time = 1:1500, signal = unlist(beta_1[[1]])[501:2000], beta = 1)

ts_plot = beta_0_df %>%
  rbind(beta_1_df) %>%
  mutate(beta = as.character(beta)) %>%
  filter(time <= 365*5) %>%
  ggplot(aes(x = time, y = signal, colour = beta)) +
  geom_line(linewidth = 0.25) +
  theme(panel.grid = element_blank(), 
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  scale_colour_manual(values = c("grey75", "black")) +
  theme(legend.position = "none") +
  labs(x = "Time (days)", y = "Environment")


## plot spectra
spectra_0 = as.data.frame(beta_0[[3]])
spectra_1 = as.data.frame(beta_1[[3]])

spectra_0 %>%
  ggplot(aes(x = freq, y = power)) + 
  geom_line() +
  scale_y_log10() + 
  scale_x_log10() + 
  geom_smooth(method = "lm") 
  
spectrum_plot = spectra_1 %>%
  ggplot(aes(x = freq, y = power)) + 
  geom_line(aes(x = freq, y = power), data = spectra_0, colour = "grey75", linewidth = 0.5) +
  geom_smooth(method = "lm", se = F, colour = "grey75",
              aes(x = freq, y = power), data = spectra_0,  linewidth = 0.6) +
  geom_line(linewidth = 0.5,  alpha = 0.7) +
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_smooth(method = "lm", se = F, colour = "black", linewidth = 0.5) +
  labs(x = "Frequency", y = "Power") +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted")  +
  theme(panel.grid = element_blank(), 
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) 

ggsave(ts_plot, width = 4, height = 2, path = "outputs/figures/didactic", filename = "ts_plot.png")
ggsave(spectrum_plot, width = 4, height = 3, path = "outputs/figures/didactic", filename = "spectrum_plot.png")


## make plot for synchrony
beta_1_sync = generate_noise(beta = 1.5,
                        p = 0,
                        n_ts = 2, 
                        L1 = 500,
                        L2 = 365*5,
                        sigma = 0.5)
## save data 
# saveRDS(beta_1_sync, "outputs/figures/didactic/beta_1_sync.rds")

beta_1_sync = readRDS( "outputs/figures/didactic/beta_1_sync.rds")

ts = beta_1_sync[[1]]
first = unlist(ts[[1]])
second = unlist(ts[[2]])


beta_1_sync_df1 = data.frame(time = 1:1500, signal =  first[501:2000], beta = 0)
beta_1_sync_df2 = data.frame( time = 1:1500, signal = second[501:2000], beta = 1)


beta_1_sync_df1 %>%
  mutate(beta = as.character(beta)) %>%
  filter(time <= 365) %>%
  ggplot(aes(x = time, y = signal)) +
  geom_line(linewidth = 0.25, colour = "black") +
  geom_line(data = beta_1_sync_df2 %>% filter(time <= 365), aes(x = time, y = signal), colour = "red") +
  theme(panel.grid = element_blank(), 
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none") +
  labs(x = "Time (days)", y = "Environment")



beta_1_async = generate_noise(beta = 1.5,
                              p = 0.9,
                              n_ts = 2, 
                              L1 = 500,
                              L2 = 365*5,
                              sigma = 0.5)


## save data 
# saveRDS(beta_1_async, "outputs/figures/didactic/beta_1_async.rds")

beta_1_async = readRDS( "outputs/figures/didactic/beta_1_async.rds")

ts = beta_1_async[[1]]
first = unlist(ts[[1]])
second = unlist(ts[[2]])

beta_1_async_df1 = data.frame(time = 1:1500, signal =  first[501:2000], beta = 0)
beta_1_async_df2 = data.frame(time = 1:1500, signal = second[501:2000], beta = 1)


synchrony_plot <- beta_1_async_df1 %>%
  mutate(beta = as.character(beta)) %>%
  filter(time <= 365) %>%
  ggplot(aes(x = time, y = signal)) +
  geom_line(linewidth = 0.25, colour = "black") +
  geom_line(data = beta_1_async_df2 %>% filter(time <= 365), aes(x = time, y = signal+0.5), colour = "#1875a5ff") +
  theme(panel.grid = element_blank(), 
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none") +
  labs(x = "Time (days)", y = "Environment")


beta_1_sync_df1 <- filter(beta_1_sync_df1, time<= 365) 
beta_1_sync_df1_new = beta_1_sync_df1

beta_1_sync_df1_new$signal<- c(beta_1_sync_df1$signal[201:365], beta_1_sync_df1$signal[1:200]) 

asynchrony_plot = beta_1_async_df1 %>%
  mutate(beta = as.character(beta)) %>%
  filter(time <= 365) %>%
  ggplot(aes(x = time, y = signal)) +
  geom_line(linewidth = 0.25, colour = "black") +
  geom_line(data = beta_1_sync_df1_new %>% filter(time <= 365), aes(x = time, y = signal +1), colour = "#4fb2e5ff") +
  theme(panel.grid = element_blank(), 
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")) +
  theme(legend.position = "none") +
  labs(x = "Time (days)", y = "Environment")


ggsave(asynchrony_plot, width = 3, height = 2, path = "outputs/figures/didactic", filename = "asynchrony_plot.png")
ggsave(synchrony_plot, width = 3, height = 2, path = "outputs/figures/didactic", filename = "synchrony_plot.png")

