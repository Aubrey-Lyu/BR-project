library(tidyverse)

setwd('/media/dian/D/data/Binocular_Rivalry/info')
a <- read_delim('Meta_GlobalQuality_ArtRepair.txt', delim  = ' ')

u <- read_delim('Meta_GlobalQuality_unrepaired.txt', delim = ' ')
# strip white space
a <- apply(a,2,function(x)gsub('\\s+', '',x))
a <- as.data.frame(apply(a,2,function(x)gsub('\\s+', '',x)))
colnames(a) <- trimws(colnames(a))

u <- apply(u,2,function(x)gsub('\\s+', '',x))
u <- as.data.frame(apply(u,2,function(x)gsub('\\s+', '',x)))
colnames(u) <- trimws(colnames(u))

a$type = 'Repaired'
u$type = 'Unrepaired'

d <- rbind(a,u)

u_t <- u %>% filter(Region=='total')
a_t <- a %>% filter(Region=='total')

# 2 key stats to determine whether to use artrepair
RSTD_total <- (as.numeric(u_t$Std) - as.numeric(a_t$Std))/as.numeric(u_t$Std)
RSTD <- (as.numeric(u$Std) - as.numeric(a$Std))/as.numeric(u$Std)

key_stat1_total <- u_t %>% select(Sub, image, Region) %>% mutate(RSTD_total = RSTD_total)
key_stat1 <- u %>% select(Sub, image, Region) %>% mutate(RSTD = RSTD)

u_Res_t <- u_t %>% filter(image=='ResMS')
a_Res_t <- a_t %>% filter(image=='ResMS')
RResMS_total <- (as.numeric(u_Res_t$Mean) - as.numeric(a_Res_t$Mean))/as.numeric(u_Res_t$Mean)
key_stat2_total <- u_Res_t %>% select(Sub, image, Region) %>% mutate(RResMS_total = RResMS_total)

# use art repair
a1_t <- key_stat1_total[key_stat1_total$RSTD_total>0,]
a2_t <- key_stat2_total[key_stat2_total$RResMS_total>0,]
# omit art repair
o1_t <- key_stat1_total[key_stat1_total$RSTD_total<0,]
o2_t <- key_stat2_total[key_stat2_total$RResMS_total<0,]

u_Res <- u %>% filter(image=='ResMS')
a_Res <- a %>% filter(image=='ResMS')
RResMS <- (as.numeric(u_Res$Mean) - as.numeric(a_Res$Mean))/as.numeric(u_Res$Mean)
key_stat2 <- u_Res %>% select(Sub, image, Region) %>% mutate(RResMS = RResMS)

a1 <- key_stat1[key_stat1$RSTD>0,]
a2 <- key_stat2[key_stat2$RResMS>0,]
o1 <- key_stat1[key_stat1$RSTD<0,]
o2 <- key_stat2[key_stat2$RResMS<0,]




