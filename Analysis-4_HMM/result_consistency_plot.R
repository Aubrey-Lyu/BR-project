# note: not completed; may not be necessary
library(R.matlab)
library(gpplot2)
library(hrbrthemes)

setwd('/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses')
matVars <- readMat('examine_state_correlation.mat')
r_mean_diff <- data.frame(correlations = unlist(matVars['r.mean.diff']))

# plot
ggplot(r_mean_diff,aes(x=correlations)) + 
  geom_histogram( aes(y = ..density..), color="#69b3a2", fill="WHITE")+
  #geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6)+
  geom_density(fill="#69b3a2", color="#69b3a2", alpha=0.5)+
  geom_vline(aes(xintercept=0, color='red'),
             linetype="dashed")+
  annotate(geom="label", x=-0.015, y=30, label="reference line",
             color="black") +
  annotate(geom="text", x=-0.1, y=10, label="One-sample T-test: \n t = -137.07, df = 22, p-value < 0.00",
           color="black")+
  theme_ipsum()+
  theme(legend.position="none",
        axis.text = element_text(size=12),
       ) +
  labs(x='Correlation difference to the right match of states') +
  xlim(c(-0.3,0))