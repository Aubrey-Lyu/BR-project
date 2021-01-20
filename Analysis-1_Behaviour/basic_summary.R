library(tidyverse)
data_trial = read_csv('/media/dian/D/data/Binocular_Rivalry/tables/table_trial.txt')
data_trial <- data_trial %>% filter(!is.nan(ALT_DELAY) & !is.nan(ALT_UNTIL))
data_green <- data_trial %>% filter(ALT_PERCEPT==1)
data_red <- data_trial %>% filter(ALT_PERCEPT==2)
data_mixed <- data_trial %>% filter(ALT_PERCEPT==3)

summary(data_mixed[,5]>=0.3 & data_mixed[,6]>=0.5)