
library(tidyverse)
library(ggplot2)
library(ggthemes)

methaneDF <- as.data.frame(read.table("data/ch4_monthly.txt", skip = 146, header = TRUE)) %>% 
  rename(ch4 = "value") %>% 
  select(year, month, ch4) %>% 
  filter(ch4 > 0)

ggplot(methaneDF, aes(x = year, y = ch4)) +
  geom_smooth(color = "black",
              se = FALSE) +
  scale_x_continuous(breaks = seq(1985, 2025, 10),
                     limits = c(1983, 2027)) +
  theme_few()

