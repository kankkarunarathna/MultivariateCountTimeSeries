####################################################### 
#This is the scrip contain necessary codes for producing
# models and figures in the manuscript.


####################################################### 
#Install mvgam package from GitHub repository 
#remotes::install_github('nicholasjclark/mvgam', force=TRUE)



####Load necessary libraries####
# This tutorial relies on the following packages:
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidybayes)       # Tidy manipulation / plots of posterior draws
library(patchwork)



# Load the pre-prepared Portal rodent abundance data
portal_ts <- read.csv('Data/portal_data.csv')


# Inspect the data structure
dplyr::glimpse(portal_ts)


# Adding 'series' column that acts as a factor 
portal_ts %>%
dplyr::mutate(series = as.factor(species)) -> portal_ts
dplyr::glimpse(portal_ts)


# see the series level
#levels(portal_ts$series)


# Figure 1
# Figure 1_(a)
png(file="Figures/Figure_1(a).png", res=400, units='in', width=6, height=5)
plot_mvgam_series(data = portal_ts, y = 'captures', series = 'all')                   
dev.off()


# Figure 1_(b) : only for the first plot
png(file="Figures/Figure_1(b).png", res=400, units='in', width=6, height=5)
plot_mvgam_series(data = portal_ts, y = 'captures', series = 1)
dev.off



# Figure S1
png(file='Figures/Figure_S1_1.png', res=400, units='in', width=6, height=3)
par(mar = c(1, 1, 1, 1))
# For species DM
portal_ts %>% 
  dplyr::filter(species == 'DM') %>%
  ggplot(aes(x = mintemp, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'DM',
       y = "log(captures)", 
       x = 'Minimum temperature') +
  portal_ts %>% 
  dplyr::filter(species == 'DM') %>%
  ggplot(aes(x = ndvi_ma12, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')
dev.off()

# For species DO
png(file='Figures/Figure_S1_2.png', res=400, units='in', width=6, height=3)
par(mar = c(1, 1, 1, 1))
portal_ts %>% 
  dplyr::filter(species == 'DO') %>%
  ggplot(aes(x = mintemp, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'DO',
       y = "log(captures)", 
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'DO') %>%
  ggplot(aes(x = ndvi_ma12, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')
dev.off()

# For species PB
png(file='Figures/Figure_S1_3.png', res=400, units='in', width=6, height=3)
portal_ts %>% 
  dplyr::filter(species == 'PB') %>%
  ggplot(aes(x = mintemp, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'PB',
       y = "log(captures)",
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'PB') %>%
  ggplot(aes(x = ndvi_ma12, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')
dev.off()

# For species PP
png(file='Figures/Figure_S1_4.png', res=400, units='in', width=6, height=3)
portal_ts %>% 
  dplyr::filter(species == 'PP') %>%
  ggplot(aes(x = mintemp, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'PP',
       y = "log(Captures)", 
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'PP') %>%
  ggplot(aes(x = ndvi_ma12, y = log(captures))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')
dev.off()



# split the data into training and testing folds to evaluate predictions
data_train <- portal_ts %>%
  dplyr::filter(time <= 68)
data_test <- portal_ts %>%
  dplyr::filter(time > 68)





####################################################### 
#### Fitting Models ####

# Model 1
mod1 <- mvgam(captures ~ 
               
             # Hierarchical intercepts capture variation in average Captures
             s(series, bs = 're') +
               
             # A shared smooth of minimum temperature
             s(mintemp, k = 8) +
              
             # allowing each species' response to mintemp to vary from the shared smooth
             s(mintemp, series, bs = 'sz', k = 8) - 1,
             
             # Define training data
             data = data_train,
             
             # Define test data
             newdata = data_test,
             
             # Poisson family for observations
             family = poisson(),
             
             # cmdstanr as backend
             backend = 'cmdstanr')

# saving model output
save(mod1, file = 'Model outputs/mod1.rda')

# loading model output
load('Model outputs/mod1.rda')

# Obtaining model summary
summary(mod1)

# Obtaining model's Stan code
code(mod1)



# Figure S2.2
png(file = 'Figures/Figure_S2.png', res = 400,
    units = 'in', width = 6, height =4 )
par(mar = c(1, 1, 1, 1))
mcmc_plot(mod1, type = 'rhat_hist')
dev.off()


# Figure S2.3
png(file = 'Figures/Figure_S3.png', res = 400,
    units = 'in', width = 8, height = 5)
par(mar = c(1, 1, 1, 1))
mcmc_plot(mod1, type = 'trace')
dev.off()


#Figure S2.4
png(file = 'Figures/Figure_S4.png', res = 400,
    units = 'in', width = 6, height = 4)
par(mar = c(1, 1, 1, 1))
pairs(mod1, variable = c('mean(series)', 'sd(series)'))
dev.off()



# Figure S2.5
# Plot the hierarchical intercepts
png(file = 'Figures/Figure_S5.png', res = 400,
    units = 'in', width =5, height =3)
par(mar = c(1, 1, 1, 1))
plot(mod1, type = 're')
dev.off()



# Figure S2.6
# Plot the hierarchical smooth components with the S3 'plot' function
png(file = 'Figures/Figure_S6.png', res = 400,
    units = 'in', width = 7, height = 4)
par(mar = c(1, 1, 1, 1))
plot(mod1, type = 'smooths')
dev.off()


# Figure 2
# Figure 2(a)
png(file = 'Figures/Figure_2(a).png', res = 400, units = 'in', width = 6, height = 4)
par(mar = c(1, 1, 1, 1))
plot_predictions(mod1, 
                 condition = c('mintemp', 'series', 'series'),
                 points = 0.5,
                 rug = TRUE) +
  theme(legend.position = 'none') +
  labs(y = 'Captures', x = 'Minimum temperature')
dev.off()




# Figure S2.7
png(file = 'Figures/Figure_S7.png', res = 400,
    units = 'in', width = 6, height = 4)
par(mar = c(1, 1, 1, 1))
plot_predictions(mod1, 
                 condition = c('mintemp', 'series', 'series'),
                 type = 'link') +
  theme(legend.position = 'none') +
  labs(y = 'log(captures)', x = 'Minimum temperature')
dev.off()




# Figure S2.8
png(file = 'Figures/Figure_S8.png', res = 400,
    units = 'in', width = 5, height = 4)
par(mar = c(1, 1, 1, 1))
plot_slopes(mod1, variable = 'mintemp',
            condition = c('series', 'series'),
            type = 'link') +
  theme(legend.position = 'none') +
  labs(y = 'log(Captures)', x = 'Series')
labs(y = 'log(captures)', x = 'Minimum temperature')
dev.off()




# Figure 2(b)
# Figure in evaluation
png(file = 'Figures/Figure_2(b).png', res = 400,
    units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
post_contrasts <- comparisons(mod1,
                              
                              # Set the contrast to compute
                              variables = 
                                list(mintemp = c(-1.75, 1.75)),
                              
                              # Compute contrast for each level of
                              # the 'series' variable
                              newdata = 
                              datagrid(series = unique)) %>%

    
# Take posterior draws
posteriordraws()


# posterior draw  
post_contrasts %>%
  ggplot(aes(x = draw)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  
  # Use the stat_halfeye function from {tidybayes} for a nice visual
  stat_halfeye(fill = "#C79999", alpha = 0.75) +
  facet_wrap(~ series, scales = 'free') +
  labs(x = "Change in Captures if mintemp is high vs low",
       y = "Density")+
  theme(legend.position = 'none')
dev.off()



# Residuals plots
par(mar = c(1, 1, 1, 1))
plot(mod1, type = 'residuals', series = 1)
plot(mod1, type = 'residuals', series = 2)
plot(mod1, type = 'residuals', series = 3)
plot(mod1, type = 'residuals', series = 4)



# obtaining forecasts
forecast(mod1, newdata = data_test)




##################################################
##### Model 2 #####

mod2 <- mvgam(captures ~ 
               
               # Hierarchical intercepts capture variation in average Captures
               s(series, bs = 're') +
                
               # shared interactions of mintemp and NDVI through hierarchical tensor products
               te(mintemp, ndvi_ma12, k = c(3, 5)) + 
               
               # within species interactions of mintemp and NDVI through hierarchical tensor products 
               te(mintemp, ndvi_ma12, by = series, k = c(3, 5), m = 1) - 1,
              
             # Define training data
             data = data_train,
             
             # Define test data
             newdata = data_test,
             
             # Poisson family for observations
             family = poisson(),
             
             # cmdstanr as backend
             backend = 'cmdstanr')


# saving model output
save(mod2, file = 'Model outputs/mod2.rda')

# loading model output
load('Model outputs/mod2.rda')

# Obtaining model summary
summary(mod2, include_betas = FALSE) # with suppressed some coefficients
summary(mod2, include_betas = TRUE)  # with all coefficients


# Sampling diagnostics
mcmc_plot(mod2, type = 'rhat_hist')

# The hierarchical smooth components are now a series of bivariate
# nonlinear interaction smooths  : 

# Figure S2.9:
png(file = 'Figures/Figure_S9.png', res = 400,
    units = 'in', width = 7, height = 6)
par(mar = c(1, 1, 1, 1))
plot(mod2, type = 'smooths')
dev.off()


# Figure S2.10:
png(file = 'Figures/Figure_S10.png', res = 400, units = 'in', width = 8, height = 4)
par(mar = c(1, 1, 1, 1))
gratia::draw(mod2$mgcv_model)
dev.off()




# Figure 3
png(file = 'Figures/Figure_3.png', res = 400, units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
plot_predictions(mod2, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'), type = 'link') +
                labs(y = 'log(Captures)', x = 'NDVI moving average')
dev.off()




# Figure S11
png(file = 'Figures/Figure_S11.png', res = 400, units = 'in', width = 7, height = 4)
par(mar = c(1, 1, 1, 1))
plot_predictions(mod2, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'),
                 points = 0.5, conf_level = 0.5) +
                labs(y = 'Captures', x = 'NDVI moving average')

dev.off()



# Figure 
# par(mar = c(1, 1, 1, 1))
# layout(matrix(1:4, nrow = 2, byrow = TRUE))
# ppc(mod1, type = 'hist', series = 1)
# title(main = 'Model 1')
# ppc(mod1, type = 'pit', series = 1)
# 
# ppc(mod2, type = 'hist', series = 1)
# title(main = 'Model 2')
# ppc(mod2, type = 'pit', series = 1)
# layout(1)



# Posterior predictions still leave a lot to be desired
par(mar = c(1, 1, 1, 1))
plot(mod2, type = 'forecast', series = 1)
plot(mod2, type = 'forecast', series = 2)
plot(mod2, type = 'forecast', series = 3)
plot(mod2, type = 'forecast', series = 4)

#
par(mar = c(1, 1, 1, 1))
plot(mod2, type = 'residuals', series = 1)
plot(mod2, type = 'residuals', series = 2)
plot(mod2, type = 'residuals', series = 3)
plot(mod2, type = 'residuals', series = 4)








##################################################
##### Model 3 #####

mod3 <- mvgam(captures ~ -1,
                
              # GAM components now moved to the latent process model
              trend_formula = ~ 
                s(trend, bs = 're') +
                
                # shared interactions of mintemp and NDVI through hierarchical tensor products
                te(mintemp, ndvi_ma12, k = c(3, 5)) +
                
                # contribution of interactions of mintemp and NDVI on trend through hierarchical tensor products
                te(mintemp, ndvi_ma12, by = trend, k = c(3, 5), m = 1),
              
              # define trend by AR1
              trend_model = 'AR1',
              
              # Set priors
              priors = c(prior(exponential(1), class = sigma),
                         prior(normal(0.5, 0.25), class = ar1,
                             lb = -1, ub = 1)),
              
              # Define training data
              data = data_train,
              
              # Define test data
              newdata = data_test,
              
              # Poisson family for observations
              family = poisson(),
              
              # cmdstanr as backend
              backend = 'cmdstanr'
              )

# saving model output
save(mod3, file = 'Model outputs/mod3.rda')


# loading model output
load('Model outputs/mod3.rda')

# obtain Stan code
code(mod3)

# obtain summary without beta etimates
summary(mod3, include_betas = FALSE)


# obtain summary with beta etimates
summary(mod3, include_betas = TRUE)




# Figure S2.12
png(file = 'Figures/Figure_S12.png', res = 400,
    units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
plot(mod3, type = 'smooths', trend_effects = TRUE)
dev.off()




# Figure 4
# Figure 4(a)  : 
par(mar = c(1, 1, 1, 1))
mcmc_plot(mod3, variable = 'ar1', regex = TRUE, type = 'areas')


# Following was used in article
m3a<-  MCMCvis::MCMCchains(mod3$model_output, 'ar1')
m3a<-as.data.frame(m3a)

m3a1<-ggplot(m3a, aes(x=m3a[,1])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3a)[1]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3a2<-ggplot(m3a, aes(x=m3a[,2])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3a)[2]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3a3<-ggplot(m3a, aes(x=m3a[,3])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3a)[3]) +
  theme(axis.ticks=element_blank(),axis.text.y =element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3a4<-ggplot(m3a, aes(x=m3a[,4])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3a)[4], ) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))


png(file = 'Figures/Figure_4(a).png', res = 400, units = 'in', width = 8, height = 7)
m3a1+m3a2+m3a3+m3a4+plot_layout(ncol =1 )
dev.off()




# Figure 4(b) : without formatting 
par(mar = c(1, 1, 1, 1))
mcmc_plot(mod3, variable = c('sigma[1]','sigma[2]','sigma[3]','sigma[4]') , regex = FALSE, type = 'areas')


# Following was used for formatting.
m3s<-  MCMCvis::MCMCchains(mod3$model_output, 'sigma')
m3s<-as.data.frame(m3s)

m3p1<-ggplot(m3s, aes(x=m3s[,1])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3s)[1]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3p2<-ggplot(m3s, aes(x=m3s[,2])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3s)[2]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3p3<-ggplot(m3s, aes(x=m3s[,3])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3s)[3]) +
  theme(axis.ticks=element_blank(),axis.text.y =element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

m3p4<-ggplot(m3s, aes(x=m3s[,4])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(colnames(m3s)[4]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))


png(file = 'Figures/Figure_4(b).png', res = 400,
    units = 'in', width = 8, height = 7)
m3p1+m3p2+m3p3+m3p4+plot_layout(ncol =1 )
dev.off()










##################################################
##### Model 4 #####

mod4 <- mvgam(captures ~ -1,
              
              # GAM components now moved to the latent process model
              trend_formula = ~
                
                # models than we have series)
                s(trend, bs = 're') +

                # covariates
                te(mintemp, ndvi_ma12, k = c(3, 5)) +
                
                
                te(mintemp, ndvi_ma12, by = trend, k = c(3, 5), m = 1),
              
              # define trend by Vector Autoregressive model with correlated errors
              trend_model = 'VAR1cor',
              
              # Prior for process error variance components 
              # LKJ prior on the error correlations)
              priors = prior(exponential(1), class = sigma),
              
              # Define training data
              data = data_train,
              
              # Define test data
              newdata = data_test,
              
              # Poisson family for observations
              family = poisson(),
              
              # cmdstanr as backend
              backend = 'cmdstanr')


# saving model output
save(mod4, file = 'Model outputs/mod4.rda')

# loading model output
load('Model outputs/mod4.rda')


# the full process error covariance matrix ('Sigma')
summary(mod4, include_betas = TRUE)



# Figure 5 
par(mar = c(1, 1, 1, 1))
A_pars <- matrix(NA, nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    A_pars[i, j] <- paste0('A[', i, ',', j, ']')
  }
}
mcmc_plot(mod4, 
          variable = as.vector(t(A_pars)), 
          type = 'hist')



# For formatting Figure 5, used the following.

A_pars <- matrix(NA, nrow = 4, ncol = 4)
nms=c('DM', 'DO', 'PB', 'PP')
for(i in 1:4){
  for(j in 1:4){
    A_pars[i, j] <- paste0(nms[j], ',', nms[i])
  }
}


s1<-  MCMCvis::MCMCchains(mod4$model_output, 'A')
colnames(s1)<-c(as.vector(A_pars))
cm<-colnames(s1)
ss<-as.data.frame(s1)


p1<-ggplot(ss, aes(x=ss[,1])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[1]) +
theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
    
p2<-ggplot(ss, aes(x=ss[,2])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[2]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p3<-ggplot(ss, aes(x=ss[,3])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[3]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p4<-ggplot(ss, aes(x=ss[,4])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[4]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p5<-ggplot(ss, aes(x=ss[,5])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[5]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p6<-ggplot(ss, aes(x=ss[,6])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[6]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p7<-ggplot(ss, aes(x=ss[,7])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[7]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

p8<-ggplot(ss, aes(x=ss[,8])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[8]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p9<-ggplot(ss, aes(x=ss[,9])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[9]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p10<-ggplot(ss, aes(x=ss[,10])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[10]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p11<-ggplot(ss, aes(x=ss[,11])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[11]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p12<-ggplot(ss, aes(x=ss[,12])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[12]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p13<-ggplot(ss, aes(x=ss[,13])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[13]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

p14<-ggplot(ss, aes(x=ss[,14])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[14]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p15<-ggplot(ss, aes(x=ss[,15])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[15]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
  
p16<-ggplot(ss, aes(x=ss[,16])) + 
geom_histogram(fill="#A25050",colour="#630000")+
ggtitle(cm[16]) +
theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

layout(1)
png(file = 'Figures/Figure_5.png', res = 400,
    units = 'in', width = 8, height = 7)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+plot_layout(ncol = 4)
dev.off()




# Figure s2.13

# without formatting
Sigma_pars <- matrix(NA, nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    Sigma_pars[i, j] <- paste0('Sigma[', i, ',', j, ']')
  }
}
mcmc_plot(mod4, 
          variable = as.vector(t(Sigma_pars)), 
          type = 'hist')


# with formatting
s2<-  MCMCvis::MCMCchains(mod4$model_output, 'Sigma')
colnames(s2)<-c(as.vector(A_pars))
cm<-colnames(s2)
ss2<-as.data.frame(s2)


pp1<-ggplot(ss2, aes(x=ss2[,1])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[1]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp2<-ggplot(ss2, aes(x=ss2[,2])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[2]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp3<-ggplot(ss2, aes(x=ss2[,3])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[3]) +
  theme(axis.ticks=element_blank(),axis.text.y =element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp4<-ggplot(ss2, aes(x=ss2[,4])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[4]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp5<-ggplot(ss2, aes(x=ss2[,5])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[5]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp6<-ggplot(ss2, aes(x=ss2[,6])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[6]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp7<-ggplot(ss2, aes(x=ss2[,7])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[7]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))
pp8<-ggplot(ss2, aes(x=ss2[,8])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[8]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp9<-ggplot(ss2, aes(x=ss2[,9])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[9]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp10<-ggplot(ss2, aes(x=ss2[,10])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[10]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp11<-ggplot(ss2, aes(x=ss2[,11])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[11]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp12<-ggplot(ss2, aes(x=ss2[,12])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[12]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp13<-ggplot(ss2, aes(x=ss2[,13])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[13]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp14<-ggplot(ss2, aes(x=ss2[,14])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[14]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp15<-ggplot(ss2, aes(x=ss2[,15])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[15]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))

pp16<-ggplot(ss2, aes(x=ss2[,16])) + 
  geom_histogram(fill="#A25050",colour="#630000")+
  ggtitle(cm[16]) +
  theme(axis.ticks=element_blank(),axis.text.y=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0.5))


png(file = 'Figures/Figure_S13.png', res = 400,
    units = 'in', width = 8, height = 7)
pp1+pp2+pp3+pp4+pp5+pp6+pp7+pp8+pp9+pp10+pp11+pp12+pp13+pp14+pp15+pp16+plot_layout(ncol = 4)
dev.off()






# Figure 6
# Figure 6(a)
png(file = 'Figures/Figure_6(a).png', res = 400,units = 'in', width = 8, height = 7)
par(mar = c(5, 4, 3, 2))
layout(matrix(1:4, nrow = 2, byrow = FALSE))
ppc(mod1, type = 'hist', series = 1)
title(main = 'Model 1')
ppc(mod1, type = 'pit', series = 1)
ppc(mod4, type = 'hist', series = 1)
title(main = 'Model 4')
ppc(mod4, type = 'pit', series = 1)
layout(1)
dev.off()

#Figure 6(b)
png(file = 'Figures/Figure_6(b).png', res = 400, units = 'in', width = 8, height = 7)
par(mar = c(5, 5, 3, 1))
layout(matrix(1:2, nrow = 1, byrow = TRUE))
ppc(mod1, type = 'density', series = 1)
title(main = 'Model 1')
ppc(mod4, type = 'density', series = 1)
title(main = 'Model 4')
layout(1)
dev.off()

#Figure 6(c)
png(file = 'Figures/Figure_6(c).png', res = 400,units = 'in', width = 8, height = 7)
par(mar = c(4, 4.5, 2, 2))
plot(mod1, type = 'forecast', series = 1,main='Model1')
layout(1)
dev.off()

#Figure 6(d)
png(file = 'Figures/Figure_6(d).png', res = 400, units = 'in', width = 8, height = 7)
par(mar = c(4, 4.5, 2, 2))
plot(mod4, type = 'forecast', series = 1,main='Model4')
layout(1)
dev.off()

#Figure 6(e)
png(file = 'Figures/Figure_6(e).png', res = 400, units = 'in', width = 8, height = 7)
diff_scores <- score(forecast(mod4), score = 'energy')$all_series$score -
  score(forecast(mod1), score = 'energy')$all_series$score
par(mar = c(4.5, 4.5, 1, 1))
plot(diff_scores, pch = 16, col = 'darkred', 
     ylim = c(-1*max(abs(diff_scores), na.rm = TRUE),
              max(abs(diff_scores), na.rm = TRUE)),
     bty = 'l',
     xlab = 'Forecast horizon',
     ylab = expression(Energy[model_4]~-~Energy[model_1]))
title("Difference in Energy Values")
box(bty = 'l', lwd = 2)
abline(h = 0, lty = 'dashed', lwd = 2)
dev.off()





# Figure S11
png(file = 'Figures/Figure_1.png', res = 400,
    units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
plot(mod4, type = 'forecast', series = 1)
plot(mod4, type = 'forecast', series = 2)
plot(mod4, type = 'forecast', series = 3)
plot(mod4, type = 'forecast', series = 4)


layout(matrix(1:4, nrow = 4)) 
par(mar = c(3.75, 4.5, 1, 1))
for (i in 1:4) { 
  plot(mod4, type = 'forecast', series = i) 
}
layout(1)






##################################################
##### Model 5 #####

mod5 <- mvgam(captures ~ 
                 
                 # Hierarchical intercepts capture variation in average Captures
                 s(series, bs = 're') +
                 
                 # interactions of mintemp and NDVI through Hierarchical tensor products 
                 te(mintemp, ndvi_ma12, k = c(3, 5)) +
                
                # interactions of mintemp and NDVI withing each species - through Hierarchical tensor products 
                 te(mintemp, ndvi_ma12, by = series, k = c(3, 5), m = 1) - 1,
               
               # define trend witg Gaussian process
               trend_model = 'GP',
               use_lv = TRUE,
               n_lv = 2,
               
              # Define training data
               data = data_train,
              
              # Define test data
               newdata = data_test,
              
              # Poisson family for observations
               family = poisson(),
              
              # cmdstanr as backend
               backend = 'cmdstanr')


# saving model output
save(mod5, file = 'Model outputs/mod5.rda')

# loading model output
load('Model outputs/mod5.rda')

# The usual summaries
summary(mod5, include_betas = FALSE)



# Figure S2.14
png(file = 'Figures/Figure_S14.png', res = 400, units = 'in', width = 8, height = 7)
plot_mvgam_factors(mod5)
dev.off()




# Figure S2.15
png(file = 'Figures/Figure_S15.png', res = 400, units = 'in', width = 8, height = 7)
correlations <- lv_correlations(object = mod5)
mean_correlations <- correlations$mean_correlations
mean_correlations[upper.tri(mean_correlations)] <- NA
mean_correlations <- data.frame(mean_correlations)
mean_correlations %>%
  dplyr::add_rownames("series1") %>%
  tidyr::pivot_longer(-c(series1), 
                      names_to = "series2", 
                      values_to = "Correlation") -> mean_correlations
ggplot(mean_correlations,
       aes(x = series1, y = series2)) + 
  geom_tile(aes(fill = Correlation)) +
  scale_fill_gradient2(low = "darkred", 
                       mid = "white", 
                       high = "darkblue",
                       midpoint = 0,
                       breaks = seq(-1,1,length.out = 5),
                       limits = c(-1, 1),
                       name = 'Trend\ncorrelation') + 
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()





# Figure S2.16
png(file = 'Figures/Figure_S16.png', res = 400, units = 'in', width = 8, height = 7)
# dynamic trend components
plot(mod5, type = 'trend', series = 1)
dev.off()



# Figure S2.17
# Figure S2.17(a)
png(file = 'Figures/Figure_S17(a) .png', res = 400, units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
plot(mod1, type = 'residuals', series = 3)
title("Model 1 ")
dev.off()


# Figure S2.17(b)
png(file = 'Figures/Figure_S17(b).png', res = 400,
    units = 'in', width = 8, height = 7)
par(mar = c(1, 1, 1, 1))
plot(mod4, type = 'residuals', series = 3)
title("Model 4")
dev.off()



# Table S3.3

fc1 <- forecast(mod1, newdata = data_test)
fc2 <- forecast(mod2, newdata = data_test)
fc3 <- forecast(mod3, newdata = data_test)
fc4 <- forecast(mod4, newdata = data_test)
fc5 <- forecast(mod5, newdata = data_test)


energy_fc1 <- score(fc1, score = 'energy')
energy_fc2 <- score(fc2, score = 'energy')
energy_fc3 <- score(fc3, score = 'energy')
energy_fc4 <- score(fc4, score = 'energy')
energy_fc5 <- score(fc5, score = 'energy')
energy_fc1$all_series
energy_fc2$all_series
energy_fc3$all_series
energy_fc4$all_series
energy_fc5$all_series





