library(SLGP)
library(tidyr)
library(dplyr)

library(ggplot2); library(ggpubr); library(viridis); library(forcats)

have_quantreg <- requireNamespace("quantreg", quietly = TRUE)
have_goftest  <- requireNamespace("goftest",  quietly = TRUE)

utils::packageVersion("SLGP")
sessionInfo()

## >>> Chunk of code 1 (paper, Section 4.1) >>>
library(SLGP); data("quakes"); df <- quakes
range_response <- c(40, 680); range_x <- c(165, 190)
modelPrior <- slgp(depth~long, # Formula to specify predictors VS response
                   data=df, method="Prior",
                   basisFunctionsUsed = "RFF", interpolateBasisFun="WNN",
                   hyperparams = list(lengthscale=c(0.15, 0.15), sigma2=1),
                   sigmaEstimationMethod = "heuristic", # rewrites sigma2
                   predictorsLower= range_x[1], predictorsUpper= range_x[2],
                   responseRange= c(40, 680), seed=1, opts=list(ndraws = 3),
                   opts_BasisFun = list(nFreq=200, MatParam=5/2))
## <<< end chunk of code 1 <

## >>> Chunk of code 2 (paper, Section 4.1) >>>
print(modelPrior)
## <<< end chunk of code 2 <

## >>> To produce Figure 1 (paper, Section 1) >>>
df <- quakes %>%
  mutate(long_bin = cut(long, breaks = seq(165, 190, by = 2.5), include.lowest = FALSE)) %>%
  group_by(long_bin) %>%
  mutate(long_bin = paste0(long_bin, "\nn=", n()))%>%
  ungroup()%>%
  mutate(long_bin = factor(long_bin, 
                           levels = sort(unique(long_bin), decreasing = FALSE))) %>%
  data.frame()

# Scatterplot: long vs depth
scatter_plot <- ggplot(df, aes(x = long, y = depth)) +
  geom_point(alpha = 0.5, color = "grey", pch=20) +
  labs(y = "Hypocentral depth (km)",
       x = "Longitude (°)",
       title = "Scatterplot of the data") +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)

# Histogram: Distribution of depth by 'long' bin
hist_plot <- ggplot(df, aes(x = depth)) +
  geom_histogram(mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 40),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(sides = "b", color = "navy", alpha = 0.5)+
  facet_wrap(~ long_bin, scales = "free_y", nrow=2) +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density", 
       title = "Histogram of 'depth' by 'longitude' group") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.015))
ggarrange(scatter_plot, hist_plot, ncol = 2, nrow = 1,
          widths = c(0.3, 0.7))
# ggsave("./FiguresQuake/scatter.pdf", width=10, height=3.5, scale=1.0)
## <<< End of code to produce Figure 1 <

## >>> Chunk of code 3 (paper, Section 4.1) >>>
plot(modelPrior, draw = c(1:3), panels = TRUE)
## <<< end chunk of code 3 <

## >>> Chunk of code 4 (paper, Section 4.1) >>>
dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 1), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")
predPrior <- predict(modelPrior, newdata = dfGrid, type="density")
## <<< end chunk of code 4 <


## >>> To produce Figure 6 (paper, Section 4.1) >>>
colnames(predPrior) <- c("long", "depth", paste0("Draw from the prior n°", seq(3)))

predPrior <- predPrior %>%  pivot_longer(-c("long", "depth"))
scale_factor <- 300

ggplot()  +
  labs(title = "Prior SLGP samples of conditional depth densities, visualised across slices",
       x = "Hypocentral depth (km)",
       y = "Longitude (°)")+
  theme_bw()+
  geom_ribbon(data=predPrior,
              mapping=aes(x=depth, ymax=scale_factor*value+long, 
                          ymin=long, group=-long, fill=long),
              col="grey", alpha=0.9)+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(title = "Indexing variable: longitude",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  coord_flip()+
  facet_grid(.~name)
# ggsave(paste0("./FiguresQuake/ribbonsPrior",  ".pdf"), width=10, height=4)
## <<< End of code to produce Figure 6 <


## >>> To produce Figure 5 (paper, Section 3.1) >>>
u <- seq(-1, 1, 0.001)
du <- diff(u[1:2])
df <- data.frame()
n_draws <- 3
for(n_poly in c(10, 100, 1000)){
  # Define the Legendre polynomials
  basis_fun <- matrix(1, nrow=n_poly, ncol=length(u))
  basis_fun[2, ] <- u
  ## recurrence (n+1) P_{n+1} = (2n+1) u P_n - n P_{n-1}.
  for(i in seq(3, n_poly)){
    n <- i - 2
    basis_fun[i, ] <- ((2*n+1)*u*basis_fun[i-1, ] - n*basis_fun[i-2, ])/(n+1) 
  }
  # Orthonormalize !
  # basis_fun <- basis_fun * sqrt((2 * (seq_len(n_poly) - 1) + 1) / 2)
  set.seed(1)
  epsilon <- matrix(rnorm(n_draws * n_poly), nrow = n_draws)
  Z <- epsilon %*% basis_fun
  SLGP <- t(apply(Z, 1, function(z) {
    e <- exp(z - max(z))
    e / (sum(e) * du)}))
  
  df <- rbind(df, data.frame(x = rep(u, times = n_draws),
                             value = as.vector(t(SLGP)),
                             draw = rep(1:n_draws, each = length(u)),
                             p = n_poly))
}
df$panel <- factor(paste0("Number of basis functions:", df$p),
                   levels = paste0("Number of basis functions:", c(10, 100, 1000)))
df$draw  <- factor(paste0("Draw from the basis function n°", df$draw))
ggplot(df, aes(x = x, y = value, colour = draw)) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ panel, nrow = 1) +
  coord_cartesian(ylim = c(0, 5)) +
  labs(x = "Response", y = "Probability density", colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave(paste0("./FiguresQuake/badBasisFun",  ".pdf"), width=10, height=4)
## <<< End of code to produce Figure 5 <


## >>> Chunk of code 5 (paper, Section 4.1) >>>
common_args <- list( formula = depth ~ long, data = df, method = "Prior",
                     basisFunctionsUsed = "RFF", interpolateBasisFun = "WNN",
                     hyperparams = list(lengthscale = c(0.15, 0.15), sigma2 = 1),
                     sigmaEstimationMethod = "heuristic", opts=list(ndraws=1),
                     predictorsLower = range_x[1], predictorsUpper = range_x[2],
                     responseRange = range_response, seed = 1)
# One model per smoothness: Matérn 1/2 (exponential), 3/2 and Inf (Gaussian)
priors <- lapply(c(1/2, 3/2, Inf), function(nu)
  do.call(slgp, c(common_args,  list(opts_BasisFun = list(nFreq = 200, MatParam = nu)))))
## <<< end chunk of code 5 <

## >>> To produce Figure 7 (paper, Section 4.1) >>>
# With Exponential = Matérn 1/2 kernel
modelPrior2 <- priors[[1]]
# With Matérn 3/2 kernel
modelPrior3 <-  priors[[2]]
# With Gaussian = Matérn Inf kernel
modelPrior4 <-  priors[[3]]

dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 1), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")

predPrior <- rbind(predict(modelPrior2, newdata = dfGrid),
                   predict(modelPrior3, newdata = dfGrid),
                   predict(modelPrior4, newdata = dfGrid))
predPrior$Kernel <- c(sapply(seq(3), function(i){rep(i, nrow(dfGrid))}))
colnames(predPrior) <- c("long", "depth", "value", "Kernel")

predPrior$Kernel <- factor(c("Exponential", "Matérn 3/2", "Gaussian")[predPrior$Kernel], 
                           levels=c("Exponential", "Matérn 3/2", "Gaussian"))
scale_factor <- 300
ggplot()  +
  labs(x = "Hypocentral depth (km)",
       y = "Longitude (°)",
       title = paste("Prior SLGP samples, visualised across slices.",
                     "SLGPs have varying degrees of smoothness.")) +
  theme_bw()+
  geom_ribbon(data=predPrior,
              mapping=aes(x=depth, ymax=scale_factor*value+long, 
                          ymin=long, group=-long, fill=long),
              col="grey", alpha=0.9)+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(title = "Indexing variable: Longitude",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  coord_flip()+
  facet_grid(.~Kernel)
#ggsave(paste0("./FiguresQuake/ribbonsPriorOthers",  ".pdf"), width=10, height=4)
## <<< End of code to produce Figure 7 <


## >>> Chunk of code 6 (paper, Section 4.1) >>>
modelMAP <- update(modelPrior, newdata = df, method = "MAP")
## <<< end chunk of code 6 <

## >>> Chunk of code 7 (paper, Section 4.1) >>>
summary(modelMAP)                   
## <<< end chunk of code 7 <

plot(modelMAP, draw = "mean", panels = TRUE)

## >>> To produce Figure 8 (paper, Section 4.1) >>>
dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 1), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")
pred <- predict(modelMAP, newdata= dfGrid)

scale_factor <- 300
ggplot()  +
  labs(x = "Hypocentral depth (km)",
       y = "Longitude (°)")+
  theme_bw()+
  geom_point(data=df,
             mapping=aes(x = depth, y = long), alpha = 0.5, 
             pch=20, color = "navy")+
  geom_ribbon(data=pred,
              mapping=aes(x=depth, ymax=scale_factor*pdf_1+long, 
                          ymin=long, group=-long, fill=long),
              col="grey", alpha=0.8)+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(title = "Indexing variable: Longitude",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  coord_flip()
# ggsave(paste0("./FiguresQuake/ribbonsMAP",  ".pdf"), width=10, height=4)
## <<< End of code to produce Figure 8 <

## >>> To produce Figure 9 (paper, Section 4.1) >>>
selected_values <- c(167, 180, 185)
gap <- 0.5
df_filtered <- df %>%
  mutate(interval=findInterval(long, c(0, 
                                       selected_values[1]-gap, 
                                       selected_values[1]+gap, 
                                       selected_values[2]-gap, 
                                       selected_values[2]+gap, 
                                       selected_values[3]-gap, 
                                       selected_values[3]+gap)))%>%
  filter(interval %in% c(2, 4, 6))%>%
  group_by(interval)%>%
  mutate(category = paste0("long close to ", c("", selected_values[1],
                                               "", selected_values[2],
                                               "", selected_values[3])[interval], 
                           "\nn=", n()))
names <- sort(unique(df_filtered$category))
dfGrid <- data.frame(expand.grid(selected_values, 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")
predMAP <- predict(modelMAP, newdata = dfGrid)
colnames(predMAP) <- c("long", "depth", "MAP estimator")
predMAP <- predMAP %>%
  pivot_longer(-c("long", "depth"))
predMAP$category <-ifelse(predMAP$long==selected_values[1], names[1],
                          ifelse(predMAP$long==selected_values[2], names[2], names[3]))

ggplot(mapping=aes(x = depth)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 20),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=predMAP, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=predMAP, mapping=aes(y=value, group=name, col=name), lwd=1.1)+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density",
       title = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP-MAP estimators at bins centers") +
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.02))
# ggsave(paste0("./FiguresQuake/histMAP",  ".pdf"), width=10, height=3.5)
## <<< End of code to produce Figure 9 <

## >>> To produce Figure A-3 (appendix, Section C-2) >>>
set.seed(1)
newOrder <- sample(seq(nrow(df)))

run_long <- FALSE ### Remove if ready to run the benchmark [! Long runtime]

### Create the synthetic datasets
if(run_long){
  
  dinvgamma <- function(x, alpha=4.5, beta=0.35) {
    ifelse(x<=0, 0, (beta^alpha / gamma(alpha)) * x^(-alpha - 1) * exp(-beta / x))
  }
  plot(dinvgamma, from=0, to=1)
  if(file.exists("OptimisingLenQuakes.RData")){
    load(file="OptimisingLenQuakes.RData")
    
  }else{
    lengthscale_grid <- seq(0.025, 0.5, 0.025)
    df_res <- data.frame(expand.grid(lengthscale_grid, lengthscale_grid))
    colnames(df_res) <- c("l_long", "l_depth")
    df_res$logPostData <- NaN
    df_res$logPostData2 <- NaN
  }
  
  
  df_res$logPrior <- log(dinvgamma(df_res$l_long))+log(dinvgamma(df_res$l_depth))
  for(i in which(round(df_res$l_long*100)%%5==0&
                 round(df_res$l_depth*100)%%5==0)){
    l1 <- df_res$l_long[i]
    l2 <- df_res$l_depth[i]
    if(is.na(df_res$logPostData[i])){
      
      modelMAPtemp <- update(modelMAP, 
                             newdata = df, 
                             method="MAP", 
                             hyperparams=list(sigma2=1, 
                                              lengthscale=c(l1, l2)),
                             sigmaEstimationMethod="heuristic")
      df_res$logPostData[i] <- modelMAPtemp@logPost
      if(sum(!is.na(df_res$logPost))%%10==0){
        save(df_res, file="OptimisingLenQuakes.RData")
      }
    }
    if(is.na(df_res$logPostData2[i])){
      modelMAPtemp <- update(modelMAP, 
                             newdata = df[newOrder[1:750], 
                                          c("long", "depth")], 
                             method="MAP", 
                             hyperparams=list(sigma2=1, 
                                              lengthscale=c(l1, l2)),
                             sigmaEstimationMethod="heuristic")
      temp <- predict(modelMAPtemp, newdata = df[newOrder[-c(1:750)], 
                                                 c("long", "depth") ])
      df_res$logPostData2[i] <- sum(log(temp$pdf_1))
      if(sum(!is.na(df_res$logPostData2))%%10==0){
        save(df_res, file="OptimisingLenQuakes.RData")
      }
    }
  }
  save(df_res, file="OptimisingLenQuakes.RData")
}
load(file="OptimisingLenQuakes.RData")

min_done <- df_res[which.min(-df_res$logPostData2-df_res$logPrior), ]

df_res %>%
  dplyr::filter(!is.na(logPostData2))%>%
  ggplot(mapping=aes(x=l_long, y=l_depth))+
  geom_tile(mapping=aes(fill=-logPostData2-logPrior))+
  geom_contour(mapping=aes(z=-logPostData2-logPrior), col="black", bins=20)+
  geom_point(data=min_done, col="red")+
  theme_bw()+
  scale_fill_viridis(direction=-1, 
                     guide = guide_colourbar(barwidth = 15))+
  labs(fill="log Posterior")+
  theme(legend.position = "bottom",
        legend.direction="horizontal")+
  ggtitle("Trained on 75% of dataset, evaluated on remaining 25%")+
  xlab("Lengthscale for 'long'")+
  ylab("Lengthscale for 'depth'")

# ggsave(paste0("./FiguresQuake/optimLandscape",  ".pdf"), width=5, height=4)
## <<< End of code to produce Figure A-3 <


## >>> To produce Figure A-4 (appendix, Section C-2) >>>
modelMAP2 <- update(modelPrior, 
                    newdata = df, 
                    method="MAP", 
                    hyperparams=list(sigma2=1, 
                                     lengthscale=c(0.01, 0.01)),
                    sigmaEstimationMethod="heuristic")
modelMAP3 <-  update(modelPrior,
                     newdata = df, 
                     method="MAP", 
                     hyperparams=list(sigma2=1, 
                                      lengthscale=c(0.5, 0.5)),
                     sigmaEstimationMethod="heuristic")
modelMAP4 <-  update(modelPrior,
                     newdata = df, 
                     method="MAP", 
                     hyperparams=list(sigma2=1, 
                                      lengthscale=c(min_done$l_depth, min_done$l_long)),
                     sigmaEstimationMethod="heuristic")
gc()
dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 1), 
                                 seq(range_response[1], range_response[2],, 101)))

colnames(dfGrid) <- c("long", "depth")
pred1 <- predict(modelMAP, newdata = dfGrid)
pred1$Len <- "Heuristic lengthscale"
pred2 <- predict(modelMAP2, newdata = dfGrid)
pred2$Len <- "Short lengthscale"
pred3 <- predict(modelMAP3, newdata = dfGrid)
pred3$Len <- "Long lengthscale"
pred4 <- predict(modelMAP4, newdata = dfGrid)
pred4$Len <- "Optimised lengthscale"
pred <- rbind(pred1, pred2, pred3, pred4)
scale_factor <- 100
ggplot()  +
  labs(x = "Hypocentral depth (km)",
       y = "Longitude (°)") +
  theme_bw()+
  geom_point(data=df,
             mapping=aes(x = depth, y = long), alpha = 0.5, color = "navy", pch=20)+
  geom_ribbon(data=pred,
              mapping=aes(x=depth, ymax=scale_factor*pdf_1+long, 
                          ymin=long, group=-long, fill=long),
              col="grey", alpha=0.8)+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(title = "Indexing variable: Longitude",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  facet_grid(Len~.)+
  coord_flip()

# ggsave(paste0("./FiguresQuake/lengthscaleRibbons",  ".pdf"),  width=10, height=10)
## <<< End of code to produce Figure A-4 <


## >>> To produce Figures A-5 (appendix, Section C-2) >>>
dfGrid <- data.frame(expand.grid(selected_values, 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("long", "depth")

pred1 <- predict(modelMAP, newdata = dfGrid)
pred1$Len <- "Heuristic lengthscale"
pred2 <- predict(modelMAP2, newdata = dfGrid)
pred2$Len <- "Short lengthscale"
pred3 <- predict(modelMAP3, newdata = dfGrid)
pred3$Len <- "Long lengthscale"
pred4 <- predict(modelMAP4, newdata = dfGrid)
pred4$Len <- "Optimised lengthscale"
pred <- rbind(pred1, pred2, pred3, pred4)

colnames(pred) <- c("long", "depth", "MAP estimator", "Len")
pred <- pred%>%
  pivot_longer(-c("long", "depth", "Len"))
pred$category <-ifelse(pred$long==selected_values[1], names[1],
                       ifelse(pred$long==selected_values[2], names[2], names[3]))

ggplot(mapping=aes(x = depth)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 20),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=pred, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=pred, mapping=aes(y=value, group=name, col=name), lwd=1.1)+
  facet_grid(Len ~ category, scales = "free_y") +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density",
       title = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP-MAP estimators at bins centers") +
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.02))

# ggsave(paste0("./FiguresQuake/lengthscaleHist",  ".pdf"), width=10, height=10)
rm(modelMAP2, modelMAP3, modelMAP4, pred1, pred2, pred3, pred4); gc()
## <<< End of code to produce Figure A-5 <

## >>> Chunk of code 8, produces Figure 10 (paper, Section 4.1) >>>
modelLaplace <- update(modelMAP, newdata = df, method="Laplace")
summary(modelLaplace)
plot(modelLaplace, draw = c("mean", 1:10), panels = TRUE)
## <<< end chunk of code 8, produces Figure 10  <

options(mc.cores = max(1L, parallel::detectCores() - 1L))
rstan::rstan_options(auto_write = TRUE)  
## >>> Chunk of code 9 (paper, Section 4.1) >>>
### ! This can take 5-10 minutes on a laptop !
modelMCMC <- update(modelMAP, newdata = df, method="MCMC",
                    opts = list(stan_chains=2, stan_iter=1000), seed=1)
## <<< end chunk of code 9 <

summary(modelMCMC)

## >>> To produce Figure 2a (paper, Section 1) >>>
dfGrid <- data.frame(expand.grid(seq(3), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("ID", "depth")
dfGrid$long <- selected_values[dfGrid$ID]
pred <- predict(modelMCMC, newdata = dfGrid)
pred$meanpdf <- rowMeans(pred[, -c(1:3)])
pred$category <- names[pred$ID]

set.seed(1)
selected_cols <- sample(seq(1000), size=10, replace=FALSE)
df_plot <- pred %>%
  dplyr::select(c("long", "depth", "category", 
                  paste0("pdf_", selected_cols)))%>%
  pivot_longer(-c("long", "depth", "category"))


ggplot(mapping=aes(x = depth)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(40, 680, 20),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=df_plot, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=pred, mapping=aes(y=meanpdf, group=category), color = "red")+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Hypocentral depth (km)",
       y = "Probability density",
       title = "Binned 'depth' histograms by 'long' (width = 1) VS SLGP-MCMC estimators at bins centers") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.02))
# ggsave(paste0("./FiguresQuake/histMCMC",  ".pdf"), width=10, height=3.5)
## <<< End of code to produce Figure 2a <



## >>> Chunk of code 10 (paper, Section 4.1) >>>
dfX <- data.frame(long=seq(range_x[1], range_x[2], 0.1))
predMean <- predict(modelMCMC, type= "moments", newdata = dfX,
                    power=c(1), centered=FALSE)
# Centered moments: variance, Kurtosis and Skewness
predVar <- predict(modelMCMC, type= "moments", newdata = dfX,
                   power=c(2, 3, 4), centered=TRUE)
## <<< end chunk of code 10 <

## >>> To produce Figure 2b (paper, Section 1) >>>
pred <- rbind(predMean, predVar)
pred <- pred %>%
  pivot_longer(-c("long", "power")) %>%
  pivot_wider(values_from = value, names_from = power) %>%
  mutate(`3` = `3` / `2`^1.5,   # skewness
         `4` = `4` / `2`^2,     # kurtosis
         `2` = sqrt(`2`)) %>%   # s.d., after standardising
  pivot_longer(-c("long", "name"), names_to = "power") %>%
  data.frame()

pred$power <- factor(c("Expected value", 
                       "Standard deviation",
                       "Skewness","Kurtosis")[as.numeric(pred$power)],
                     levels=c("Expected value", "Standard deviation",
                              "Skewness", "Kurtosis"))
df_plot <- pred %>%
  group_by(long, power)%>%
  summarise(q10 = quantile(value, probs=c(0.1)),
            q50 = quantile(value, probs=c(0.5)),
            q90 = quantile(value, probs=c(0.9)),
            mean = mean(value), .groups="keep")%>%
  ungroup() # summarise uncertainty


ggplot(df_plot, mapping=aes(x = long, group=power)) +
  geom_ribbon(mapping = aes(ymin=q10, ymax=q90),
              alpha = 0.25, lty=2, col="black", fill="cornflowerblue")+
  geom_line(mapping=aes(y=q50))+
  facet_wrap(.~power,
             scales = "free", nrow=1)+
  labs(x =  "Longitude (°)", 
       y = "Moment value") +
  theme_bw()+
  coord_cartesian(xlim=range_x)
ggsave(paste0("./FiguresQuake/MomentsMCMC",  ".pdf"), width=8, height=3)
## <<< End of code to produce Figure 2b <

## >>> Chunk of code 11 (paper, Section 4.1) >>>
pred <- predict(modelMCMC, type = "quantiles", newdata = dfX,
                probs = c(5, 25, 50, 75, 95)/100)
## <<< end chunk of code 11 <

## >>> To produce Figure 2c (paper, Section 1) >>>
df_plot <- pred %>%
  pivot_longer(-c("long", "probs")) %>%
  group_by(long, probs) %>%
  summarise(q10 = quantile(value, 0.1), q50 = quantile(value, 0.5),
            q90 = quantile(value, 0.9), .groups = "drop") %>%
  mutate(probs = factor(paste0("Quantile: ", 100 * probs, "%"),
                        levels = paste0("Quantile: ", c(5, 25, 50, 75, 95), "%")))

# Comparison 1: classical quantile regression (cubic polynomial)
library(quantreg)
long_center <- mean(df$long)
long_scale <- sd(df$long)
df$long2 <- (df$long-long_center)/long_scale
dfX$long2 <- (dfX$long-long_center)/long_scale

fit_q <- quantreg::rq(depth ~ poly(long2, degree=3, raw=TRUE), 
                      data = df, tau = c(5, 25, 50, 75, 95)/100, method="br")
pred_q <- predict(fit_q, newdata = dfX)
pred$qQuantileReg <- c(t(pred_q))

# Comparison 2: empirical quantiles on bins, kept only when n >= 10
bin_w <- 1  
min_n  <- 10

emp <- do.call(rbind, lapply(dfX$long, function(x0) {
  in_win <- df$depth[df$long >= x0 & df$long <= x0 + bin_w]
  ok <- length(in_win) >= min_n
  data.frame(long  = x0,
             n     = length(in_win),
             probs =  c(5, 25, 50, 75, 95)/100,
             q     = if (ok) as.numeric(quantile(in_win, probs =  c(5, 25, 50, 75, 95)/100))
             else rep(NA_real_, 5))}))
emp$probs <- factor(paste0("Quantile: ", 100*emp$probs, "%"),
                    levels = paste0("Quantile: ", c(5, 25, 50, 75, 95), "%"))

plot_slgp <- ggplot(df_plot, aes(x = long)) +
  geom_ribbon(aes(ymin = q10, ymax = q90, fill = probs, group = probs),
              alpha = 0.2, colour = NA) +
  geom_line(aes(y = q50, col = probs, group = probs, lty = "Model (SLGP or QR)"),
            lwd = 0.75) +
  geom_line(data = emp, aes(y = q, col = probs, group = probs, lty = "Empirical"),
            lwd = 0.9) +
  labs(x = "Longitude [°E]", y = "Hypocentre depth [km]",
       col = "Quantile levels", fill = "Quantile levels",
       lty = "Quantile estimation method") +
  theme_bw() +
  coord_cartesian(xlim = range_x, ylim = range_response) +
  scale_linetype_manual(values = c("Empirical" = 2,
                                   "Model (SLGP or QR)" = 1))

pred$probs <- factor(paste0("Quantile: ", 100*pred$probs, "%"),
                     levels = paste0("Quantile: ", c(5, 25, 50, 75, 95), "%"))
plot_qr <- ggplot() +
  geom_line(data = pred, mapping = aes(x = long, y = qQuantileReg, 
                                       lty = "Model (SLGP or QR)",
                                       col = probs, group = probs), lwd = 0.75) +
  geom_line(data = emp, 
            mapping = aes(x = long, y = q, col = probs, group = probs,
                          lty = "Empirical"), lwd = 0.9) +
  labs(x = "Longitude [°E]", y = "Hypocentre depth [km]",
       col = "Quantile levels", lty = "Quantile estimation method") +
  theme_bw() +
  coord_cartesian(xlim = range_x, ylim = range_response) +
  scale_linetype_manual(values = c("Empirical" = 2,
                                   "Model (SLGP or QR)" = 1))

plot_hist <- ggplot(df, mapping = aes(x = long)) +
  geom_histogram(col = "navy", fill = "grey", alpha = 0.3,
                 breaks = seq(range_x[1], range_x[2], 1)) +
  theme_minimal() +
  labs(x = NULL, y = "Sample\ncount")

library(ggpubr)
p_with_marginals <- ggarrange(plot_hist + 
                                labs(title = "SLGP vs empirical quantiles"),
                              plot_hist+ 
                                labs(title = "Quantile regression from quantreg vs empirical quantiles"),
                              plot_slgp,
                              plot_qr,
                              ncol = 2,
                              nrow = 2,
                              heights = c(1.5, 4),
                              common.legend = TRUE,
                              legend = "bottom",
                              align = "hv")

print(p_with_marginals)
# ggsave(p_with_marginals, file = "./FiguresQuake/QuantilesMAP2.pdf", width = 8, height = 5, scale=1.5)
## <<< End of code to produce Figure 2c <

## >>> To produce Figure 11 (paper, Section 4.2) >>>
range_x   <- c(165, 190)
mag_levels <- seq(4.0, 6.4, by = 0.1)  # observed discrete support


scatter_mag <- ggplot(df, aes(x = long, y = mag)) +
  geom_point(size = 1, alpha = 0.4, colour = "navy") +
  coord_cartesian(xlim = range_x, ylim=c(0, 6.4)) +
  labs(x = "Longitude [°E]", y = "Magnitude",
       title = "Observed magnitudes") +
  theme_bw()

df$long_bin <- cut(df$long, breaks = seq(165, 190, by = 2.5))
hist_mag <- ggplot(df, aes(x = mag)) +
  geom_histogram(aes(y = after_stat(density)),
                 breaks = seq(0.0, 6.5, by = 0.1),
                 fill = "darkgrey", colour = "grey50", linewidth = 0.2) +
  facet_wrap(~ long_bin, nrow = 2) +
  labs(x = "Magnitude", y = "Probability density",
       title = "Magnitude distribution by longitude band") +
  theme_bw()

ggarrange(scatter_mag, hist_mag, ncol = 2, widths = c(0.35, 0.65))
# ggsave("./FiguresQuake/scatterDiscrete.pdf", width=10, height=3.5)
## <<< End of code to produce Figure 11 <


## >>> Chunk of code 12, produces figure 12 (paper, Section 4.2) >>>
trend_fun <- function(df){return(ifelse(df$mag < 4, -10, 0))}
modelMAGdisc <- slgp(mag~long,
                     data=df, method="MAP", discrete=TRUE, nIntegral = 66,
                     basisFunctionsUsed = "RFF", interpolateBasisFun="WNN",
                     hyperparams = list(lengthscale=c(0.15, 0.15), sigma2=1),
                     sigmaEstimationMethod = "heuristic", # rewrites sigma2
                     predictorsLower= range_x[1], predictorsUpper= range_x[2],
                     responseRange= c(0, 6.5), trend=trend_fun, seed=1,
                     opts_BasisFun = list(nFreq=250, MatParam=5/2))
plot(modelMAGdisc, draw = "mean", panels = TRUE)
# dev.copy(pdf, "./FiguresQuake/defaultMAPDiscrete.pdf", width = 10, height = 3.5)
# dev.off()
## <<< end chunk of code 12, produces figure 12  <

## >>> To produce Figure 13 (paper, Section 4.2) >>>
range_mag_full <- c(0, 6.5)   # FULL modelled range, including the unobserved part
K_full <- 66                  

selected_values <- c(167, 181, 185)
gap <- 1.0
df_filtered <- df %>%
  mutate(interval=findInterval(long, c(0, 
                                       selected_values[1]-gap, 
                                       selected_values[1]+gap, 
                                       selected_values[2]-gap, 
                                       selected_values[2]+gap, 
                                       selected_values[3]-gap, 
                                       selected_values[3]+gap)))%>%
  filter(interval %in% c(2, 4, 6))%>%
  group_by(interval)%>%
  mutate(category = paste0("long close to ", c("", selected_values[1],
                                               "", selected_values[2],
                                               "", selected_values[3])[interval], 
                           "\nn=", n()))
names <- sort(unique(df_filtered$category))
dfGrid <- data.frame(expand.grid(selected_values, 
                                 seq(range_mag_full[1], range_mag_full[2],, K_full)))
colnames(dfGrid) <- c("long", "mag")
predMAP <- predict(modelMAGdisc, newdata = dfGrid, discrete=TRUE, nIntegral = K_full)
colnames(predMAP) <- c("long", "mag", "MAP estimator")
predMAP <- predMAP%>%
  pivot_longer(-c("long", "mag"))
predMAP$category <-ifelse(predMAP$long==selected_values[1], names[1],
                          ifelse(predMAP$long==selected_values[2], names[2], names[3]))

ggplot(mapping=aes(x = mag)) +
  geom_histogram(df_filtered,
                 mapping=aes(y = after_stat(density)/10, group=category),
                 closed = "left", binwidth = 0.1, boundary = 0,
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_step(data=predMAP, mapping=aes(y=value, group=name, col=name), lwd=1.1)+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Magnitude (-)",
       y = "Probability density",
       title = "Binned 'mag' histograms by 'long' (width = 2) VS SLGP-MAP estimators at bins centers")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        legend.title = element_blank())+
  coord_cartesian(xlim=range_mag_full,
                  ylim=c(0, 0.15))

# ggsave(paste0("./FiguresQuake/histMAPdisc",  ".pdf"), width=10, height=3.5)
## <<< End of code to produce Figure 13 <


## >>> Chunk of code 13 (paper, Section 4.3) >>>
set.seed(1); id_train <- sample(seq_len(nrow(df)), size = 750)
df_train <- df[id_train, ]; df_test  <- df[-id_train, ]
modelMAP2D <- slgp(depth ~ long + lat,
                   data = df_train, method = "MAP",
                   basisFunctionsUsed = "RFF", interpolateBasisFun = "WNN",
                   sigmaEstimationMethod = "heuristic", seed = 1,
                   predictorsLower = c(165,-40), predictorsUpper = c(190, -10),
                   responseRange = range_response,
                   opts_BasisFun = list(nFreq = 250, MatParam = 5/2))
summary(modelMAP2D)
## <<< end chunk of code 13 <

## >>> Chunk of code 14 (paper, Section 4.3) >>>
pit_pred <- predict(modelMAP2D, type = "cdf", newdata = df_test)
ks.test(pit_pred$cdf_1, "punif", 0, 1) # D = 0.0625, p = 0.282
library(goftest); cvm.test(pit_pred$cdf_1) # omega2 = 0.3467, p = 0.100
## <<< end chunk of code 14 <


## >>> To produce Figure 14 (paper, Section 4.2) >>>
pit_hist <- ggplot(pit_pred, aes(x = cdf_1)) +
  geom_histogram(breaks = seq(0, 1, by = 0.1),
                 fill = "darkgrey", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = nrow(df_test) / 10,
             linetype = 2, linewidth = 0.8, colour = "navy") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(x = "PIT", y = "Count", title = "PIT histogram") +
  theme_bw()

n_pit <- length(pit_pred$cdf_1)
qq_df <- data.frame(theoretical = (seq_len(n_pit) - 0.5) / n_pit,   # uniform plotting positions
                    empirical   = sort(pit_pred$cdf_1))

pit_qq <- ggplot(qq_df, aes(x = theoretical, y = empirical)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "navy") +
  geom_point(size = 1, alpha = 0.6, colour = "grey20") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Theoretical uniform quantile", y = "Empirical PIT quantile",
       title = "Uniform Q-Q plot") +
  theme_bw()

ggarrange(pit_hist, pit_qq, ncol = 2, widths = c(0.5, 0.5))
# ggsave("./FiguresQuake/pit2D.pdf", width = 9, height = 3.0)
## <<< End of code to produce Figure 14 <


## >>> To produce Figure 15 (paper, Section 4.2) >>>
range_xlong <- c(165, 190) 
range_xlat <- c(-40, -10) 

map_plot <- ggplot(mapping = aes(x = long, y = lat, fill = depth)) +
  geom_point(data = df_train, mapping = aes(col = "Train", pch = "Train"),
             size = 1.5, alpha = 0.8) +
  geom_point(data = df_test,  mapping = aes(col = "Test",  pch = "Test"),
             size = 2, alpha = 0.8) +
  scale_fill_viridis_c(name = "Hypocentre\ndepth [km]",
                       limits = range_response, option = "viridis") +
  scale_shape_manual(values = c("Train" = 21, "Test" = 22)) +
  scale_color_manual(values = c("Train" = "grey", "Test" = "red")) +
  coord_cartesian(xlim = range_xlong, ylim = range_xlat) +
  labs(x = "Longitude [°E]", y = "Latitude [°]",
       title = "Observed hypocentre depths",
       col = "Data split", shape = "Data split") +
  theme_bw()

long_min <- c(166, 180, 178, 182, 180, 178)
long_max <- c(170, 184, 182, 186, 184, 182)
lat_min  <- c(-15, -20, -25, -25, -30, -35)
lat_max  <- c(-10, -15, -20, -20, -25, -30)

bins <- data.frame(bin_id = seq_along(long_min),
                   long_min, long_max, lat_min, lat_max) %>%
  mutate(bin_name = paste0("Bin (", letters[bin_id], ")"))

df_train$set <- "Train"; df_test$set <- "Test"
df_split <- bind_rows(df_train, df_test) %>%
  mutate(set = factor(set, levels = c("Train", "Test")))

df_selected <- do.call(rbind, lapply(seq_len(nrow(bins)), function(i) {
  tmp <- df_split %>%
    filter(long >  bins$long_min[i], long <= bins$long_max[i],
           lat  >  bins$lat_min[i],  lat  <= bins$lat_max[i])
  tmp$bin_id <- bins$bin_id[i]; tmp$bin_name <- bins$bin_name[i]
  tmp}))

bin_counts <- df_selected %>%
  count(bin_id, bin_name, set) %>%
  tidyr::pivot_wider(names_from = set, values_from = n, values_fill = 0)

bins <- bins %>%
  left_join(bin_counts, by = c("bin_id", "bin_name")) %>%
  mutate(facet_label = paste0(bin_name, "\n",
                              "long in (", long_min, ", ", long_max, "]\n",
                              "lat in (", lat_min, ", ", lat_max, "]\n",
                              "nTrain = ", Train, ", nTest = ", Test))

df_selected <- df_selected %>%
  left_join(bins %>% select(bin_id, facet_label), by = "bin_id") %>%
  mutate(facet_label = factor(facet_label, levels = bins$facet_label))

bin_centers <- bins %>%
  mutate(long = (long_min + long_max) / 2, lat = (lat_min + lat_max) / 2)

depth_grid <- seq(range_response[1], range_response[2], length.out = 101)
dfGrid <- merge(bin_centers %>% select(bin_id, bin_name, facet_label, long, lat),
                data.frame(depth = depth_grid))

predMAP <- predict(modelMAP2D, newdata = dfGrid %>% select(long, lat, depth))
predMAP <- cbind(dfGrid, predMAP %>% select(starts_with("pdf_")))


hist_plotMAP <- ggplot(mapping = aes(x = depth)) +
  geom_histogram(data = df_selected, aes(y = after_stat(density), 
                                         fill = set, colour = set),
                 breaks = seq(range_response[1], range_response[2], by = 20),
                 position = "identity", alpha = 0.35, linewidth = 0.3) + 
  geom_rug(data = df_selected, sides = "b", colour = "navy", alpha = 0.5) +
  geom_line(data = predMAP,
            aes(y = pdf_1, group = "MAP estimator", col = "MAP estimator"),
            linewidth = 1.1) +
  scale_fill_manual(values = c("Train" = "grey50", "Test" = "navy")) +
  scale_colour_manual(values = c("Train" = "grey30", "Test" = "navy", 
                                 "MAP estimator" = "tomato"))+
  guides(colour = guide_legend(
    override.aes = list(fill = c("grey50", "navy", NA),
                        colour = c("grey30", "navy", "tomato"),
                        linetype = c(0, 0, 1))),
    fill = "none")+
  facet_wrap(~ facet_label, scales = "free_y", nrow = 2) +
  labs(x = "Hypocentre depth [km]", y = "Probability density",
       title = "Empirical depth distributions vs SLGP estimates at bin centres") +
  coord_cartesian(xlim = range_response) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank())


map_plot2 <- map_plot +
  geom_rect(data = bins,
            aes(xmin = long_min, xmax = long_max,
                ymin = lat_min, ymax = lat_max),
            inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.8) +
  geom_label(data = bins,
             aes(x = (long_min + long_max) / 2, y = (lat_min + lat_max) / 2,
                 label = bin_name),
             inherit.aes = FALSE, fontface = "bold", size = 3,
             fill = "white", alpha = 0.8)

ggarrange(map_plot2, hist_plotMAP, ncol = 2, nrow = 1, widths = c(0.45, 0.55))
# ggsave("./FiguresQuake/scatter2D.pdf", width = 10, height = 4, scale=1.1)
## <<< End of code to produce Figure 15 <

## >>> Chunk of code 16 (paper, Section 4.3) >>>
grid2D_cdf <- expand.grid(long = seq(165, 190), lat = seq(-40, -10), depth=300)
pred_cdf <- predict(modelMAP2D, type = "cdf", newdata = grid2D_cdf)
grid2D_cdf$exceed <- 1 - pred_cdf$cdf_1
## <<< end chunk of code 16 <

## >>> To produce Figure 3 (paper, Section 1) >>>
threshold <- 300
df_train2 <- df_train %>%
  mutate(depth_class = ifelse(depth > threshold,
                              "Observed depth > 300 km",
                              "Observed depth <= 300 km"))

ggplot() +
  geom_raster(data = grid2D_cdf,
              aes(x = long, y = lat, fill = exceed),
              interpolate = TRUE) +
  geom_point(data = df_train2,
             aes(x = long, y = lat, shape = depth_class, colour = depth_class),
             inherit.aes = FALSE,
             size = 1.7,
             stroke = 0.6) +
  scale_shape_manual(name = "Observed events",
                     values = c("Observed depth > 300 km" = 0,
                                "Observed depth <= 300 km" = 1)) +
  scale_colour_manual(name = "Observed events",
                      values = c("Observed depth > 300 km" = "grey40",
                                 "Observed depth <= 300 km" = "grey80")) +
  scale_fill_viridis_c(name = paste0("P(depth > ", threshold, " km)"),
                       limits = c(0, 1),
                       option = "plasma") +
  coord_cartesian(xlim = range_xlong,
                  ylim = range_xlat) +
  labs(x = "Longitude [°E]",
       y = "Latitude [°]",
       title = paste0("Probability of an event deeper than ", threshold, " km")) +
  theme_bw()  +
  theme(legend.position = "right",
        legend.key = element_rect(fill = "white", colour = NA))

# ggsave("./FiguresQuake/exceed2D.pdf", width = 6, height = 4)
## <<< End of code to produce Figure 3 <


## >>> Chunk of code 15 (paper, Section 4.3) >>>
sims <- simulate(modelMAP2D, nsim = 1, newdata = df[, c("long", "lat")])
## <<< end chunk of code 15 <

## >>> To produce Figure 16 (paper, Section 4.2) >>>
df_obs <- df %>%
  mutate(event_type = "Events to train model")

df_sim <- sims %>%
  mutate(event_type = "Simulated new events at the same locations")

df_plot <- bind_rows(df_obs, df_sim) %>%
  mutate(event_type = factor(event_type,
                             levels = c("Events to train model",
                                        "Simulated new events at the same locations")))

ggplot(df_plot, aes(x = long, y = lat, fill = depth)) +
  geom_point(size = 1.5, alpha = 0.8, pch = 21, col = "grey", alpha=0.8) +
  scale_fill_viridis_c(name = "Hypocentre\ndepth [km]",
                       limits = range_response,
                       option = "viridis") +
  coord_cartesian(xlim = range_xlong, ylim = range_xlat) +
  facet_wrap(~ event_type, nrow = 1) +
  ggtitle("Observed and simulated events from the fitted SLGP") +
  labs(x = "Longitude [°E]",
       y = "Latitude [°]") +
  theme_bw()
# ggsave("FiguresQuake/simulate.pdf", width = 8, height = 4)
## <<< End of code to produce Figure 16 <


## >>> To produce the benchmark Sections 5.1 - 5.3 >>>
### Setup
modelRef <- modelMAP

n_list <- c(10, 25, 50, 
            100, 250, 500, 
            1000, 2500, 5000, 
            10000, 25000, 50000,
            100000)

config1 <- expand.grid(1, c(25, 50, 100, 200, 300),
                       seq(25), 0, 0.15, "Unif")
config2 <- expand.grid(2, 200, 
                       1, seq(25), c(0.01, 0.05, 0.15, 0.25, 0.5), "Unif")
config3 <- expand.grid(3, 200, 
                       1, seq(25), 0.15, c("Unif", "Clustered", 
                                           "Clustered2",  "Clustered3",
                                           "Hole-y"))
colnames(config1) <-
  colnames(config2) <-
  colnames(config3) <- c("Scenario", "nFreq", "SeedFreq", 
                         "SeedSamples", "lengthscale", "Sampling")
config1$SeedSamples <- config1$SeedFreq
config <- rbind(config1, config2, config3)
config <- config[, c("Scenario", "Sampling", "SeedSamples",  
                     "nFreq", "SeedFreq", "lengthscale")]
run_long <- FALSE ### Remove if ready to run the benchmark [! Long runtime]

### Create the synthetic datasets
if(run_long){
  n <- max(n_list)
  for(strat in unique(config$Sampling)){
    for(SeedSamples in sort(unique(config$SeedSamples))){
      #cat("Strategy ", strat, " Seed ", SeedSamples, "\n")
      title <- paste0("./resQuake/samp_", strat, "_num_", SeedSamples, ".Rdata")
      samp <- data.frame()
      if(strat == "Hole-y"){
        title2 <- paste0("./resQuake/samp_", "Unif", 
                         "_num_", SeedSamples, ".Rdata")
        load(title2)
        samp <- samp[abs(samp$long-172.5)>=2.5 & abs(samp$long-185)>=1, ]
        samp <- samp[1:50000, ]      
        save(samp, file=title)
      }
      if(!(file.exists(title))){
        for(slice in seq(25)){  
          cat("Strategy ", strat, " Seed ", 
              SeedSamples, "slice", slice, "/25\n")
          current_seed <- SeedSamples*1000+slice
          set.seed(current_seed)
          if(strat == "Unif"){
            newX <- data.frame(long=runif(n/25, range_x[1], range_x[2]))
            nsamp=1
          }
          if(strat == "Clustered"){
            x <- sample(seq(range_x[1], range_x[2],, 21), 
                        size=n/25, replace=T)
            t <- table(x)
            newX <- data.frame(long=as.numeric(names(t)))
            nsamp <- as.vector(t)
          }
          if(strat == "Clustered2"){
            x <- sample(seq(range_x[1], range_x[2],, 11), 
                        size=n/25, replace=T)
            t <- table(x)
            newX <- data.frame(long=as.numeric(names(t)))
            nsamp <- as.vector(t)
          }
          if(strat == "Clustered3"){
            x <- sample(seq(range_x[1], range_x[2],, 5), 
                        size=n/25, replace=T)
            t <- table(x)
            newX <- data.frame(long=as.numeric(names(t)))
            nsamp <- as.vector(t)
          }
          temp <- simulate(modelMAP, 
                           newdata = newX, 
                           seed=current_seed*2,
                           nsim=nsamp, 
                           interpolateBasisFun = "WNN")
          samp <- rbind(samp, temp[sample(seq(nrow(temp))), ])
        }
        save(samp, file=title)
        gc()
      }
    }
  }
}

### Run the evaluations
if(run_long){
  dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2],, 101), 
                                   seq(range_response[1], range_response[2],, 101)))
  colnames(dfGrid) <- c("long", "depth")
  dt <- diff(sort(unique(dfGrid$depth))[1:2])
  predRef <- predict(modelMAP, newdata = dfGrid)
  predRefCDF <- predict(modelMAP, newdata = dfGrid, type = "cdf")
  predRefCDF <- round(predRefCDF, 15)
  df_compute <- predRef
  colnames(df_compute)[3] <- "TruePDF"
  df_compute$TrueCDF <- predRefCDF$cdf_1
  df_compute$TruePDFsafe <- df_compute$TruePDF + 1e-15
  rm(predRefCDF, predRef)
  
  pb <- txtProgressBar(min = 0, max = nrow(config), style = 3)
  
  
  for(i in seq(nrow(config))){
    Scenario <- config$Scenario[i]
    nFreq <- config$nFreq[i]
    SeedFreq <- config$SeedFreq[i]
    SeedSamples <- config$SeedSamples[i]
    lengthscale <- config$lengthscale[i]
    Sampling <- config$Sampling[i]
    
    title <- paste0("./resQuake/samp_", Sampling, "_num_",
                    SeedSamples, ".Rdata")
    if(file.exists(title)){
      
      load(file=title)
      
      set.seed(SeedFreq)
      title2 <- paste0("./resQuake/SLGP_len_", 100*lengthscale,
                       "_nFreq_", nFreq,
                       "_rep_", SeedFreq, ".Rdata")
      if(!file.exists(title2)){
        modelPrior <- slgp(depth~long, 
                           data=samp[1:n_list[1], ],
                           method="none", 
                           basisFunctionsUsed = "RFF",
                           interpolateBasisFun="WNN", 
                           hyperparams = list(lengthscale=c(lengthscale,
                                                            lengthscale), 
                                              sigma2=1), 
                           sigmaEstimationMethod = "heuristic", 
                           predictorsLower= c(range_x[1]),
                           predictorsUpper= c(range_x[2]),
                           responseRange= range_response,
                           opts_BasisFun = list(nFreq=nFreq,
                                                MatParam=5/2),
                           seed=SeedFreq)
        save(modelPrior, file=title2)
        cat("\nCreated and saved model ", title2, "\n")
      }else{
        load(title2)
      }
      title3 <- paste0("./resQuake/strat_", Sampling,
                       "_SLGP_nFreq_", nFreq,
                       "_seedFreq_", SeedFreq, 
                       "_len_", lengthscale*100, "%",
                       "_rep_", SeedSamples, ".Rdata")
      
      if(!file.exists(title3)){
        d_list <- data.frame()
        if(Sampling=="Hole-y"){
          n_list2 <- n_list[1:(length(n_list)-1)]
        }else{
          n_list2 <- n_list
        }
        for(j in seq_along(n_list2)){
          start_time <- Sys.time()
          cat(".")
          n <- n_list2[j]
          modelCurrent <- update(modelPrior, newdata = samp[1:n, ], method="MAP",
                                 sigmaEstimationMethod = "heuristic")
          predTemp <- predict(modelCurrent, newdata = dfGrid)
          predTempCDF <- predict(modelCurrent, type = "cdf", newdata = dfGrid)
          predTempCDF <- round(predTempCDF, 15)
          df_compute$estPDF <- predTemp$pdf_1
          df_compute$estPDFsafe <- df_compute$estPDF + 1e-15
          df_compute$estCDF <- predTempCDF$cdf_1
          
          # ggplot(df_compute, mapping=aes(x=long, y=depth, fill=TruePDF-estPDF))+
          #   geom_raster()+theme_bw()+scale_fill_viridis()+#limits = c(-0.12, 0.12))+        
          #   labs(fill="pdf")
          
          distances <- df_compute %>%
            group_by(long)%>%
            summarise(H=sqrt(0.5* dt*
                               sum((sqrt(estPDF) - sqrt(TruePDF))^2)),
                      TV = 0.5 * sum(abs(estPDF - TruePDF)) * dt,
                      KL_TE = sum(TruePDFsafe * 
                                    log(TruePDFsafe / estPDFsafe)) * dt,
                      KL_ET = sum(estPDFsafe * 
                                    log(estPDFsafe / TruePDFsafe)) * dt,
                      W1 = sum(abs(TrueCDF - estCDF)) * dt,
                      KS = max(abs(TrueCDF - estCDF)),
                      CvM =  sum((TrueCDF - estCDF)^2 * TruePDF) * dt,
                      CvMsym =  sum((TrueCDF - estCDF)^2) * dt)%>%
            data.frame()
          distances$n <- n
          end_time <- Sys.time()
          diff_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
          distances$time <- diff_minutes 
          
          d_list <- rbind(d_list, distances)
        }
        
        save(d_list, file=title3)
        cat(title3, "(",round(sum(d_list$time)/101, 2),  "m)\n")
        setTxtProgressBar(pb, i)  # Update progress
        gc()
      }
    }
    
  }
  close(pb)
}

### Assemble results
df_res <- data.frame()

for(i in seq(nrow(config))){
  Scenario <- config$Scenario[i]
  nFreq <- config$nFreq[i]
  SeedFreq <- config$SeedFreq[i]
  SeedSamples <- config$SeedSamples[i]
  lengthscale <- config$lengthscale[i]
  Sampling <- config$Sampling[i]
  title3 <- paste0("./resQuake/strat_", Sampling,
                   "_SLGP_nFreq_", nFreq,
                   "_seedFreq_", SeedFreq, 
                   "_len_", lengthscale*100, "%",
                   "_rep_", SeedSamples, ".Rdata")
  
  if(file.exists(title3)){
    load(file=title3)
    
    temp <- data.frame(Scenario=Scenario,
                       nFreq=nFreq,
                       SeedFreq=SeedFreq,
                       SeedSamples=SeedSamples,
                       lengthscale=lengthscale,
                       Sampling=Sampling)
    if(!is.null(d_list$nFreq)){ d_list$nFreq <- NULL}
    temp <- cbind(temp, d_list)
    df_res <- rbind(df_res, temp)
  }
}
df_res <- df_res %>%
  group_by(Scenario, nFreq, SeedFreq, SeedSamples, lengthscale, Sampling)%>%
  mutate(time=ifelse(var(time)<=1e-15, NA, time)) %>%
  ungroup()%>%
  data.frame()
## <<< End of code to produce the benchmark Sections 5.1 - 5.3  <


## >>> To produce Figure 17 (paper, Section 5.1) >>>
library(forcats)
df_res %>%
  filter(Scenario == 1 ) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "lengthscale", 
                   "Sampling", "time"))%>%
  pivot_longer(-c("nFreq", "SeedSamples", "long", "n")) %>%
  filter(as.character(name) %in% c("KL_ET", "H", "TV", "CvM")) %>%
  # rename facets
  mutate(name = fct_recode(name,
                           "Hellinger" = "H",
                           "Total Variation" = "TV",
                           "KL True -> Est" = "KL_TE",
                           "Kullback-Leibler (pred, ref)" = "KL_ET",
                           "Wasserstein-1" = "W1",
                           "Kolmogorov–Smirnov" = "KS",
                           "Cramér–von Mises" = "CvM",
                           "CvM Symmetric" = "CvMsym"),
  )%>%
  group_by(nFreq, SeedSamples, name, n)%>%
  summarise(value=sum(dx*value), .groups="keep")%>%
  ungroup()%>%
  group_by(nFreq, name, n) %>%
  summarise(meanV=mean(value),
            q25=quantile(value, probs=0.25),
            q50=quantile(value, probs=0.5),
            q75=quantile(value, probs=0.75), .groups="keep")%>%
  ggplot(aes(x=n, group=nFreq, 
             col=as.factor(nFreq), fill=as.factor(nFreq)))+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, lty=2)+
  geom_line(mapping= aes(y=meanV))+
  theme_bw()+
  facet_wrap(name~., scales="free", ncol=4)+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  scale_x_log10()+
  xlab("Size of the training sample")+
  ylab("Integrated dissimilarities")+
  labs(fill = "Number of frequencies",
       col= "Number of frequencies")+
  theme(legend.direction = "horizontal", legend.position = "bottom")

# ggsave(paste0("./FiguresQuake/benchmarkS1",  ".pdf"), width=10, height=3., scale=0.8)
## <<< End of code to produce Figure 17 <

## >>> To produce Figure 18 (paper, Section 5.2) >>>
df_res %>%
  filter(Scenario == 2 ) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "Sampling", "time"))%>%
  pivot_longer(-c("lengthscale", "SeedSamples", "long", "n"))%>%
  filter(as.character(name) %in% c("KL_ET", "H", "TV", "CvM")) %>%
  # rename facets
  mutate(name = fct_recode(name,
                           "Hellinger" = "H",
                           "Total Variation" = "TV",
                           "KL(True, Est)" = "KL_TE",
                           "Kullback-Leibler (pred, ref)" = "KL_ET",
                           "Wasserstein-1" = "W1",
                           "Kolmogorov–Smirnov" = "KS",
                           "Cramér–von Mises" = "CvM",
                           "CvM Symmetric" = "CvMsym"),
  )%>%
  group_by(lengthscale, SeedSamples, name, n)%>%
  summarise(value=sum(dx*value), .groups="keep")%>%
  ungroup()%>%
  group_by(lengthscale, name, n) %>%
  summarise(meanV=mean(value),
            q25=quantile(value, probs=0.25),
            q75=quantile(value, probs=0.75), .groups="keep")%>%
  ggplot(aes(x=n, group=lengthscale, 
             col=as.factor(lengthscale), fill=as.factor(lengthscale)))+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, lty=2)+
  geom_line(mapping= aes(y=meanV))+
  theme_bw()+
  facet_wrap(name~., scales="free", ncol=4)+
  scale_x_log10()+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  xlab("Size of the training sample")+
  ylab("Integrated dissimilarities")+
  labs(fill = "lengthscales",
       col= "lengthscales")+
  theme(legend.direction = "horizontal", legend.position = "bottom")
# ggsave(paste0("./FiguresQuake/benchmarkS2",  ".pdf"), width=10, height=3., scale=0.8)
## <<< End of code to produce Figure 18 <

## >>> To produce Figure 19 (paper, Section 5.3) >>>
df_res %>%
  filter(Scenario == 3 ) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "lengthscale", "time"))%>%
  pivot_longer(-c("Sampling", "SeedSamples", "long", "n"))%>%
  filter(as.character(name) %in% c("KL_ET", "H", "TV", "CvM")) %>%
  mutate(name = fct_recode(name,
                           "Hellinger" = "H",
                           "Total Variation" = "TV",
                           "KL(True, Est)" = "KL_ET",
                           "Cramér–von Mises" = "CvM"),
         Sampling = fct_recode(Sampling,
                               "Uniform" = "Unif",
                               "Uniform with holes" = "Hole-y",
                               "Unif. on regular grid with 21pts" = "Clustered",
                               "Unif. on regular grid with 11pts" = "Clustered2",
                               "Unif. on regular grid with 05pts" = "Clustered3"),
  )%>%
  group_by(Sampling, SeedSamples, name, n)%>%
  summarise(value=sum(dx*value), .groups="keep")%>%
  ungroup()%>%
  group_by(Sampling, name, n) %>%
  summarise(meanV=mean(value),
            q25=quantile(value, probs=0.25),
            q75=quantile(value, probs=0.75), .groups="keep")%>%
  ggplot(aes(x=n, group=Sampling, 
             col=as.factor(Sampling), fill=as.factor(Sampling)))+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, lty=2)+
  geom_line(mapping= aes(y=meanV))+
  theme_bw()+
  facet_wrap(name~., ncol=4, scales="free")+
  scale_x_log10()+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  xlab("Size of the training sample")+
  ylab("Integrated dissimilarities")+
  theme(legend.direction = "horizontal", legend.position = "bottom")+
  labs(fill = "Sampling",
       col = "Sampling") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         fill = guide_legend(nrow = 2, byrow = TRUE)) 
# ggsave(paste0("./FiguresQuake/benchmarkS3",  ".pdf"), width=10, height=3.5, scale=0.8)
## <<< End of code to produce Figure 19 <

## >>> To produce Figure A-6 (appendix, Section D) >>>
temp <- df_res %>%
  filter(Scenario == 3 ) %>%
  filter(n %in% c(100, 500, 5000, 50000)) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "lengthscale", "time"))%>%
  pivot_longer(-c("SeedSamples", "long", "n", "Sampling"))%>%
  filter(name == "KL_TE" ) %>%
  group_by(name, n, long, Sampling) %>%
  summarise(meanV=mean(value),
            q75=quantile(value, probs=0.75),
            q50=quantile(value, probs=0.5),
            q25=quantile(value, probs=0.25), .groups="keep")%>%
  ungroup()%>%
  # rename facets
  mutate(name = "Kullback-Leibler (pred, ref)",
         Sampling = fct_recode(Sampling,
                               "Uniform" = "Unif",
                               "Uniform with holes" = "Hole-y",
                               "Unif. on regular grid with 21pts" = "Clustered",
                               "Unif. on regular grid with 11pts" = "Clustered2",
                               "Unif. on regular grid with 05pts" = "Clustered3"),
  )%>%
  data.frame()
temp2 <- temp %>%
  group_by(name) %>%
  summarise(ymin = 0,
            ymax = max(q75),
            Sampling = "Uniform with holes", .groups="keep")%>%
  data.frame()

ggplot(temp, aes(x=long))+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, 
              lty=2, fill="cornflowerblue", col="navy")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 170, xmax = 175, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 184, xmax = 186, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 05pts",
                               long = seq(range_x[1], range_x[2],, 5)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 11pts",
                               long = seq(range_x[1], range_x[2],, 11)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 21pts",
                               long = seq(range_x[1], range_x[2],, 21)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_line(mapping= aes(y=meanV), col="navy")+
  theme_bw()+
  facet_grid(Sampling~paste0("Sample size: ", n))+
  ylab("Kullback-Leibler divergence")+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme(legend.direction = "horizontal", legend.position = "bottom",
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

# ggsave(paste0("./FiguresQuake/benchmarkS3localKL",  ".pdf"), width=10, height=10, scale=0.975)
## <<< End of code to produce Figure A-6 <

## >>> To produce Figure A-7 (appendix, Section D) >>>
temp <- df_res %>%
  filter(Scenario == 3 ) %>%
  filter(n %in% c(100, 500, 5000, 50000)) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "lengthscale", "time"))%>%
  pivot_longer(-c("SeedSamples", "long", "n", "Sampling"))%>%
  filter(name == "H" ) %>%
  group_by(name, n, long, Sampling) %>%
  summarise(meanV=mean(value),
            q75=quantile(value, probs=0.75),
            q50=quantile(value, probs=0.5),
            q25=quantile(value, probs=0.25), .groups="keep")%>%
  ungroup()%>%
  # rename facets
  mutate(name = "Hellinger",
         Sampling = fct_recode(Sampling,
                               "Uniform" = "Unif",
                               "Uniform with holes" = "Hole-y",
                               "Unif. on regular grid with 21pts" = "Clustered",
                               "Unif. on regular grid with 11pts" = "Clustered2",
                               "Unif. on regular grid with 05pts" = "Clustered3"),
  )%>%
  data.frame()
temp2 <- temp %>%
  group_by(name) %>%
  summarise(ymin = 0,
            ymax = max(q75),
            Sampling = "Uniform with holes", .groups="keep")%>%
  data.frame()

ggplot(temp, aes(x=long))+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, 
              lty=2, fill="cornflowerblue", col="navy")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 170, xmax = 175, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 184, xmax = 186, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 05pts",
                               long = seq(range_x[1], range_x[2],, 5)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 11pts",
                               long = seq(range_x[1], range_x[2],, 11)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 21pts",
                               long = seq(range_x[1], range_x[2],, 21)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_line(mapping= aes(y=meanV), col="navy")+
  theme_bw()+
  facet_grid(Sampling~paste0("Sample size: ", n))+
  ylab("Hellinger distance")+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme(legend.direction = "horizontal", legend.position = "bottom",
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
# ggsave(paste0("./FiguresQuake/benchmarkS3localH",  ".pdf"), width=10, height=10, scale=0.975)
## <<< End of code to produce Figure A-7 <


## >>> To produce Figure A-8 (appendix, Section D) >>>
temp <- df_res %>%
  filter(Scenario == 3 ) %>%
  filter(n %in% c(100, 500, 5000, 50000)) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "lengthscale", "time"))%>%
  pivot_longer(-c("SeedSamples", "long", "n", "Sampling"))%>%
  filter(name == "CvM" ) %>%
  group_by(name, n, long, Sampling) %>%
  summarise(meanV=mean(value),
            q75=quantile(value, probs=0.75),
            q50=quantile(value, probs=0.5),
            q25=quantile(value, probs=0.25), .groups="keep")%>%
  ungroup()%>%
  # rename facets
  mutate(name = "Cramér–von Mises",
         Sampling = fct_recode(Sampling,
                               "Uniform" = "Unif",
                               "Uniform with holes" = "Hole-y",
                               "Unif. on regular grid with 21pts" = "Clustered",
                               "Unif. on regular grid with 11pts" = "Clustered2",
                               "Unif. on regular grid with 05pts" = "Clustered3"),
  )%>%
  data.frame()
temp2 <- temp %>%
  group_by(name) %>%
  summarise(ymin = 0,
            ymax = max(q75),
            Sampling = "Uniform with holes", .groups="keep")%>%
  data.frame()

ggplot(temp, aes(x=long))+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, 
              lty=2, fill="cornflowerblue", col="navy")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 170, xmax = 175, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 184, xmax = 186, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 05pts",
                               long = seq(range_x[1], range_x[2],, 5)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 11pts",
                               long = seq(range_x[1], range_x[2],, 11)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 21pts",
                               long = seq(range_x[1], range_x[2],, 21)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_line(mapping= aes(y=meanV), col="navy")+
  theme_bw()+
  facet_grid(Sampling~paste0("Sample size: ", n))+
  ylab("Cramér–von Mises value")+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme(legend.direction = "horizontal", legend.position = "bottom",
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

# ggsave(paste0("./FiguresQuake/benchmarkS3localCvM",  ".pdf"), width=10, height=10, scale=0.975)
## <<< End of code to produce Figure A-8 <


## >>> To produce Figure A-9 (appendix, Section D) >>>
temp <- df_res %>%
  filter(Scenario == 3 ) %>%
  filter(n %in% c(100, 500, 5000, 50000)) %>%
  dplyr::select(-c("Scenario", "SeedFreq", "nFreq", "lengthscale", "time"))%>%
  pivot_longer(-c("SeedSamples", "long", "n", "Sampling"))%>%
  filter(name == "TV" ) %>%
  group_by(name, n, long, Sampling) %>%
  summarise(meanV=mean(value),
            q75=quantile(value, probs=0.75),
            q50=quantile(value, probs=0.5),
            q25=quantile(value, probs=0.25), .groups="keep")%>%
  ungroup()%>%
  # rename facets
  mutate(name = "Total Variation",
         Sampling = fct_recode(Sampling,
                               "Uniform" = "Unif",
                               "Uniform with holes" = "Hole-y",
                               "Unif. on regular grid with 21pts" = "Clustered",
                               "Unif. on regular grid with 11pts" = "Clustered2",
                               "Unif. on regular grid with 05pts" = "Clustered3"),
  )%>%
  data.frame()
temp2 <- temp %>%
  group_by(name) %>%
  summarise(ymin = 0,
            ymax = max(q75),
            Sampling = "Uniform with holes", .groups="keep")%>%
  data.frame()

ggplot(temp, aes(x=long))+
  geom_ribbon(mapping=aes(ymin = q25, ymax=q75), alpha=0.1, 
              lty=2, fill="cornflowerblue", col="navy")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 170, xmax = 175, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_rect(data = temp2,
            mapping = aes(xmin = 184, xmax = 186, ymin = ymin, ymax = ymax),
            alpha = 0.2, fill = "grey", inherit.aes = FALSE,
            lty=2, col = "darkgrey")+
  geom_hline(yintercept=0, lty=2, col="black")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 05pts",
                               long = seq(range_x[1], range_x[2],, 5)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 11pts",
                               long = seq(range_x[1], range_x[2],, 11)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_vline(data = data.frame(Sampling="Unif. on regular grid with 21pts",
                               long = seq(range_x[1], range_x[2],, 21)),
             mapping = aes(xintercept = long),
             lty=2, col = "darkgrey")+
  geom_line(mapping= aes(y=meanV), col="navy")+
  theme_bw()+
  facet_grid(Sampling~paste0("Sample size: ", n))+
  ylab("Total Variation distance")+
  scale_y_continuous(transform = scales::pseudo_log_trans())+
  theme(legend.direction = "horizontal", legend.position = "bottom",
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

#ggsave(paste0("./FiguresQuake/benchmarkS3localTV",  ".pdf"), width=10, height=10, scale=0.975)
## <<< End of code to produce Figure A-9 <


## >>> To produce the benchmark Section 5.4 >>>
### Setup
timing_n     <- c(100, 1000, 10000)
timing_sch   <- c("nothing", "NN", "WNN")
timing_meth   <- c("MAP", "Laplace", "MCMC")
timing_reps  <- seq(10)
timings_nfreq <- c(25, 50, 100, 250)

config1 <- expand.grid(n = timing_n,
                       scheme      = timing_sch,
                       method      = "MAP",
                       rep         = timing_reps,
                       nfreq       = timings_nfreq,
                       stringsAsFactors = FALSE)
config2 <- expand.grid(n           = timing_n,
                       scheme      = "WNN",
                       method      = timing_meth,
                       rep         = timing_reps,
                       nfreq       = timings_nfreq,
                       stringsAsFactors = FALSE)
config_time <- unique(rbind(config1, config2))
nrow(config_time)

timing_file <- "./resQuake/timings.Rdata"
checkpoint_every <- 10

if (!file.exists(timing_file)) {
  t_list <- data.frame()
}else{
  load(timing_file)
}
if(nrow(t_list)< nrow(config_time)){
  pb <- txtProgressBar(min = 0, max = nrow(config_time), style = 3)
  since_save <- 0
  for (i in seq(nrow(t_list)+1, nrow(config_time))) {
    cfg <- config_time[i, ]
    load(file = paste0("./resQuake/samp_Unif_num_", cfg$rep, ".Rdata"))
    set.seed(cfg$rep)
    fit <- slgp(depth ~ long,
                data = samp[1:cfg$n, ],
                method =  cfg$method,
                basisFunctionsUsed  = "RFF",
                interpolateBasisFun = cfg$scheme,
                hyperparams = list(lengthscale = rep(0.15, 2),
                                   sigma2 = 1),
                sigmaEstimationMethod = "heuristic",
                predictorsLower = range_x[1], predictorsUpper = range_x[2],
                responseRange   = range_response,
                opts_BasisFun = list(nFreq = cfg$nfreq, MatParam = 5 / 2),
                seed =  cfg$rep)
    
    tm <- timing(fit)
    tm$ndistinct   <- nrow(unique(samp[1:cfg$n, "long", drop = FALSE]))
    tm$scheme      <- cfg$scheme
    tm$rep         <- cfg$rep
    t_list <- rbind(t_list, tm)
    
    cat(sprintf(
      "\n[%3d/%3d] %-7s %-7s n=%6d p=%4d  setup %6.2fs  est %7.2fs",
      i, nrow(config_time), cfg$method, cfg$scheme, tm$nobs, tm$p,
      tm$setup, tm$estimation))
    setTxtProgressBar(pb, i)
    rm(fit); gc()
    since_save <- since_save + 1L
    if (since_save >= checkpoint_every) {
      save(t_list, file = timing_file)
      since_save <- 0L
      cat(sprintf("\n    ... checkpoint saved (%d rows)\n", nrow(t_list)))
    }
  }
  close(pb)
  save(t_list, file = timing_file)
} 
## <<< End of code to produce the benchmark Section 5.4  <

## >>> To produce Figure 20 (paper, Section 5.4) >>>
rank_levels <- paste("p =", sort(unique(t_list$p)))

df_tot <- t_list %>%
  filter(method == "MAP") %>%
  mutate(scheme = factor(scheme, levels = c("nothing", "NN", "WNN")),
         p_lab  = factor(paste("p =", p), levels = rank_levels))

ggplot(df_tot, aes(x = nobs, y = total, colour = scheme, fill = scheme)) +
  stat_summary(fun = median, geom = "line", linewidth = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.5) +
  stat_summary(fun.min = function(z) quantile(z, 0.1),
               fun.max = function(z) quantile(z, 0.9),
               geom = "ribbon", alpha = 0.2, colour = NA) +
  facet_wrap(~ p_lab, nrow = 1) +
  scale_x_log10(breaks = c(100, 1000, 10000),
                labels = c("100", "1000", "10000")) +
  scale_y_log10() +
  labs(x = "Training sample size", y = "Total fitting time [s]",
       colour = "Integral approximation scheme",
       fill   = "Integral approximation scheme") +
  theme_bw() + theme(legend.position = "bottom")
# ggsave("./FiguresQuake/timingTotalScheme.pdf", width = 10, height = 3.5)
## <<< End of code to produce Figure 20 <

## >>> To produce Figure A-10 (appendix, Section E) >>>
n_levels <- paste("n =", sort(unique(t_list$nobs)))

df_ratio <- t_list %>%
  filter(scheme == "WNN") %>%
  mutate(ratio  = estimation / total,
         method = factor(method, levels = c("MAP", "Laplace", "MCMC")),
         n_lab  = factor(paste("n =", nobs), levels = n_levels))

ggplot(df_ratio, aes(x = nobs, y = total, colour = method, fill = method)) +
  stat_summary(fun = median, geom = "line", linewidth = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.5) +
  stat_summary(fun.min = function(z) quantile(z, 0.1),
               fun.max = function(z) quantile(z, 0.9),
               geom = "ribbon", alpha = 0.2, colour = NA) +
  facet_wrap(~ paste("p =", p), nrow = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Training sample size", y = "Total fitting time [s]",
       colour = "Estimation method", fill = "Estimation method") +
  theme_bw() + theme(legend.position = "bottom")
# ggsave("./FiguresQuake/timingMethod.pdf", width = 10, height = 3.5)
## <<< End of code to produce Figure A-10 <


## >>> To produce Figure A-11 (appendix, Section E) >>>
n_levels <- paste("n =", sort(unique(t_list$nobs)))

df_ratio <- t_list %>%
  filter(scheme == "WNN") %>%
  mutate(ratio  = estimation / total,
         method = factor(method, levels = c("MAP", "Laplace", "MCMC")),
         n_lab  = factor(paste("n =", nobs), levels = n_levels))

ggplot(df_ratio, aes(x = p, y = ratio, colour = method, fill = method)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  stat_summary(fun = median, geom = "line", linewidth = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.5) +
  stat_summary(fun.min = function(z) quantile(z, 0.1),
               fun.max = function(z) quantile(z, 0.9),
               geom = "ribbon", alpha = 0.2, colour = NA) +
  facet_wrap(~ n_lab, nrow = 1) +
  scale_x_log10(breaks = sort(unique(t_list$p))) +
  scale_y_log10() +
  labs(x = "Rank p", y = "Estimation time / Total time",
       colour = "Estimation method", fill = "Estimation method") +
  theme_bw() + theme(legend.position = "bottom")
# ggsave("./FiguresQuake/timingRatioMethod.pdf", width = 10, height = 3.5)
## <<< End of code to produce Figure A-11 <


