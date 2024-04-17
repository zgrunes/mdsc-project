
library(tidyverse)
library(ggplot2)
library(ggthemes)

carbonDF <- as.data.frame(read.table("data/co2_monthly.txt", skip = 158, header = TRUE)) %>% 
  rename(co2 = "value") %>% 
  select(year, month, co2) %>% 
  filter(co2 > 0)

ggplot(carbonDF, aes(x = year, y = co2)) +
  geom_smooth(color = "black",
              se = FALSE) +
  scale_x_continuous(breaks = seq(1975, 2025, 10),
                     limits = c(1973, 2027)) +
  theme_few()

