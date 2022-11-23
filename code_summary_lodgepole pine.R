---
title: "code summary for manuscript Predicting fundamental climate niches of forest trees based on species occurrence data"
author: "Yueru Zhao"
date: '2022-11-22'
note: 'to avoid redundancy, this file only contains code for lodgepole pine, because the code for Douglas-fir will be exactly the same but with different dataset'
---
  
setwd("C:/Users/zyrkim/Desktop/PhD/data")
library(rJava)
library(plyr) 
library(dplyr)
library(Matrix)
library(sp)
library(raster)
library(lme4)
library(nlme)
library(bit64)
library(data.table)
library(gridExtra)
library(Matrix)
library(rsq)
library(r2glmm)
library(MuMIn)
library(raster)
library(CEMT)
library(gbm)
library(extraTrees)
library(climatenaAPI)
library(randomForest)
library(xgboost)
library(maps)
library(mapdata)
library(dismo)
library(ggplot2)
library(caret)
library(biomod2)
library(verification)
library(enmSdm)
library(mgcv)
library(sf);library(ROCR);library(cvAUC);library(blockCV)

######selecting clim variables for each algorithm
dat_original <- read.csv("p_a_data/pa_wo_dup_yk1.csv")#native observations
dat_p <- subset(dat_original, pa==1)
dat_a<- subset(dat_original, pa==0)
#random forest
mf <- mcmfRF2(dat_p, dat_a, nr=1, varList=varList_rf, yCol='PINUCON', reg=F, nTree=100, nForest=10)
sort(mf$importance)
#MaxEnt
for(i in 1:10){ 
  dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)),])
  x <- dat_sample[,varList_max];x
  p <- dat_sample[,c("PINUCON")] #make p/a data a vector, without log lat.
  Maxent <- maxent(x, p=p, silent=TRUE);
  imp_max <- var.importance(Maxent)#var importance rank
  sorT(imp_max,byCols=c(2))
}
#GAM
dat2 <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)),])
g1 <- bam(PINUCON~s(mat, bs="tp")+s(map, bs="tp")+s(mwmt, bs="tp"), family=binomial, data=dat2, nthreads=6)
summary(g1)
######

######set up variable list
varList <- c('mat','mwmt','mcmt','td', 'map','msp', 'ahm', 'dd5', 'nffd', 'bffp', 'effp', 'ffp', 'pas', 'emt', 'ext', 'eref', 'cmd', 'tave_wt', 'tave_sp', 'tave_sm', 'tave_at', 'tmax_wt', 'tmax_sp', 'tmax_sm', 'tmax_at', 'tmin_wt', 'tmin_sp', 'tmin_sm', 'tmin_at', 'ppt_wt', 'ppt_sp', 'ppt_sm', 'ppt_at', 'dd_0_wt', 'dd_0_sp', 'dd_0_sm', 'dd_0_at', 'dd5_wt', 'dd5_sp', 'dd5_sm', 'dd5_at', 'dd18_wt', 'dd18_sp', 'dd18_sm', 'dd18_at', 'nffd_wt', 'nffd_sp', 'nffd_sm', 'nffd_at', 'pas_wt', 'pas_sp', 'pas_sm', 'pas_at', 'eref_wt', 'eref_sp', 'eref_sm', 'eref_at', 'cmd_wt', 'cmd_sp', 'cmd_sm', 'cmd_at', 'dd_18', 'dd18', 'dd_0')#full list - all annual and seasonal clim var
varList_rf <- c("ppt_sm","ppt_at","ppt_sp","map", "msp","ppt_wt","td","eref_at", "tmax_wt","emt","dd_0_sp","tmin_wt","dd_0_at","ahm","pas_wt", "dd_0","tmin_at","dd_0_wt","eref_sm","tmin_sp","pas_at", "bffp","cmd_sm", "tave_sp", "eref", "pas","dd_18","effp","tmin_sm","ext","mwmt","nffd","tmax_at", "cmd_sp",  "mcmt",    "ffp","tave_at", "tave_wt", "pas_sp", "dd5_sm", "tmax_sm", "nffd_at", "tmax_sp", "cmd")#for random forest - lodgepole pine
varList_max <- c('pas_at','eref_at','dd5_wt','ahm','dd18_sm','cmd_at','td','eref_wt','ppt_wt','tmax_at','pas_sp','dd_0_sp','tmax_sm','ppt_at','nffd_sm','dd18_sp','tmin_sm', 'nffd_sp', 'ppt_sm', 'pas_wt', 'dd18_at')#for maxent - lodgepole pine
stk_w <- rasterStack('World data/Normal_1961_1990SY/',varList, rType='tif', vConvert=TRUE);stk_w#Raster stack for world current period- data generated using ClimateNA
stk_w_45 <- rasterStack('World data/15GCM-Ensemble_rcp45_2055SY/',varList, rType='tif', vConvert=TRUE);stk_w_45;names(stk_w_45)<- tolower(names(stk_w_45))#world-future period
stk_wna <- rasterStack('WNA mapping/climatedata/1961_1990/',varList,rType='tif', vConvert=TRUE);names(stk_wna)<- tolower(names(stk_wna))#raster stack for western north america current period
stk_wna_45 <- rasterStack('WNA mapping/climatedata/8GCMs_ensemble_ssp245_2041-2070MSY/',varList,rType='tif', vConvert=TRUE);names(stk_wna_45)<- tolower(names(stk_wna_45))#raster stack for western north america future period
######

######generating dataset used for modelling
dat_original <- read.csv("p_a_data/pa_wo_dup_yk1.csv")#native observations
dat_p <- subset(dat_original, pa==1)
######genertaing pseudo-absence data for native range - same for both species
xyv <- rasterToPoints(stk_wna);hd(xyv) #convert raster to points
xyv2 <- xyv[sample(1:nrow(xyv),120000,replace=FALSE),];head(xyv2);dim(xyv2) #sampling
dat <- data.frame(xyv2);hd(dat) #make it as a dataframe
dat$pa <- 0;hd(dat) #add a new variable
# combine the presence and absence data, and remove overlapping absence data points---
colnames(dat_p)[3] <- "y";colnames(dat_p)[4] <- "x";dat_p$pa <- '1'
cmb <- rbind(dat_p[,names(dat)], dat); hd(cmb)
cmb$x2 <- round(cmb$x,1); #make a new variable x2 from rounded lon
cmb$y2 <- round(cmb$y,1);head(cmb) #make a new variable x2 from rounded lat
#drop the data points falling into the same grid at 0.1 degree (10 km)
cmb2 <- cmb[!(duplicated(cmb[,c('x2','y2')])|duplicated(cmb[,c('x2','y2')])), ];hd(cmb2)
xy_n <- subset(cmb2,pa==0);hd(xy_n) #re-generate the absence data
cmb3 <- rbind(dat_p[,names(dat)],xy_n[,names(dat)]);hd(cmb3); #combine the presence data with the new absence data
write.csv(cmb3,"p_a_data/Pl_all_abs.csv")
#######
dat_full <- read.csv("p_a_data/Pl_all_abs.csv")
names(dat_original[,2:91]) <- tolower(names(dat_original[,2:91]))
dat_full_ab <- subset(dat_full,pa==0);nrow(dat_full_ab)
dat2 <- bind(dat_full_ab, dat_original)
dat_p <- subset(dat2, pa==1)
dat_a <- subset(dat2, pa==0)

mp <- raster("WNA mapping/climatedata/1961_1990/mat.tif")
pt <- st_as_sf(dat2, coords = c("x", "y"), crs = crs(mp))
plot(mp)
points(dat2$x, dat2$y, pch = ".")

#####looking for best pr/ab ratio - correspond to figure2 in manuscript
SDMs <- c("rf","gam","maxent")
SS <- array(NA,c(10,4,3,5), list(c(seq(from = 0.1, to = 1, by = 0.1)), c("sensitivity", "specificity", "S+S", "AUC"), c(SDMs),c(seq(from = 1, to = 5, by = 1))))

sb <- spatialBlock(speciesData = pt,
                   species = "lpp",
                   rasterLayer = mp,
                   theRange = 1000000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE,
                   xOffset = 0, # shift the blocks horizontally
                   yOffset = 0)
sb$plots + geom_sf(data = pt, alpha = 0.5)
dat2$fold <- sb$foldID;folds <- dat2$fold
for (k in 1:5){#repeat CV for each of the 5 fold
  train.data <- dat2[which(folds != k),]# training set indices
  test.data <- dat2[which(folds == k),]
  dat_p <- subset(train.data, pa == 1)
  dat_a <- subset(train.data, pa == 0)
  for (i in seq(from = 0.1, to = 1, by = 0.1)){
    dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)*i),])
    for (SDM in SDMs){
      if(SDM == "maxent"){
        x <- dat_sample[,varList_max] #train.data[,varList_max]
        p <- dat_sample[,c("pa")] #make p/a data a vector, without log lat.
        Maxent <- maxent(x, p=p, silent=TRUE)
        pr_m <- predict(Maxent, test.data, "prob")
        predictions <- pr_m
      } 
      else if(SDM == "rf"){
        for(f in 1:10){ 
          library("foreach"); library("doSNOW");library(dismo);library(randomForest)
          registerDoSNOW(makeCluster(7, type="SOCK"))
          x2 <- dat_sample[,varList_rf]
          #class
          y2 <- as.factor(as.matrix(dat_sample[,'pa']))
          rf2 <- foreach(ntree = rep(floor(nTree/7),7), .combine = randomForest::combine, .packages = "randomForest") %dopar% randomForest(x2, y2, ntree = ntree, importance=T)
          if(f==1){rfC <- rf2}else{rfC <- randomForest::combine(rfC,rf2)}
        }
        gc()
        mf <- rfC
        pr_rf <- predict(mf, test.data[,varList_rf], "prob")
        predictions <- pr_rf[,2]
      }
      else if(SDM == "gam"){    ####GAM
        gc()
        library(mgcv)
        g1 <- bam(pa~s(mat, bs="tp")+ s(mwmt, bs="tp")+ s(map, bs="tp"), family=binomial, data=dat_sample, nthreads=6, discrete = TRUE)
        pr_g <- predict(g1, test.data)
        pr_g[pr_g > 0]<- 0
        predictions <- exp(pr_g)
      }
      SS[i*10,4,c(SDM),k] <- AUC(predictions,test.data$pa)
      predictions[predictions<0.5] <- 0; predictions[predictions>=0.5] <- 1
      tab.out <- table(as.numeric(test.data$pa), as.numeric(predictions))
      a <- tryCatch(tab.out["1", "1"], error=function(e) 0)%>%as.numeric()
      b <- tryCatch(tab.out["0", "1"], error=function(e) 0)%>%as.numeric()
      c <- tryCatch(tab.out["1", "0"], error=function(e) 0)%>%as.numeric()
      d <- tryCatch(tab.out["0", "0"], error=function(e) 0)%>%as.numeric()
      n <- a + b + c + d
      SS[i*10,1,c(SDM),k] <- a /(a+c)
      SS[i*10,2,c(SDM),k] <- d/(b+d)
      SS[i*10,3,c(SDM),k] <- a /(a+c) + d/(b+d)
    }  
  }
}

###average each fold
for(SDM in SDMs){
  for (i in seq(1:4)){
    SS_cv[,i,c(SDM)] <- rowMeans(SS[,i,c(SDM),])
  }
}
save(SS_cv,file = "LPP_ratio_CV.RData")
#plot for Figure2 in manuscript - lodgepole pine
for(SDM in SDMs){
  dd <- SS_cv[,1:2,c(SDM)]%>%as.data.frame()
  dd$ratio <- seq(0.1, 1, 0.1)
  #dd$tss <- dd$sensitivity + dd$specificity - 1
  #plot(dd$ratio, dd$tss)
  dd_m <- reshape2::melt(dd, id.vars = "ratio", variable.name = 'metric')
  ggplot(dd_m, aes(ratio,value)) + geom_line(aes(colour = metric), size = 1.5) + ggtitle(paste0(SDM)) +  xlab("p/a ratio")+ ylab("Value")+ scale_color_manual(values=c("black", "red"))+
    theme(legend.position= "right", 
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.text.x = element_text(size=14, color = 'black'),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text(size=14, color = 'black'),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.key=element_blank(),
          legend.key.width=unit(1.5, "cm")) + #xlim(0.75, 1)+ #ylim(0,1)+
    scale_x_continuous(breaks=seq(from = 0.2, to = 1, by = 0.2), labels=c("1:0.2", "1:0.4", "1:0.6", "1:0.8", "1:1"))+
    #scale_x_continuous(breaks=seq(from = 1:0.1, to = 1:1, by = 0.2))+
    scale_y_continuous(limits = c(0.6, 1),)
  
}


######spatial block cross-validation to test each model's classification power, correspond to Figure3 
#array for storing 5 fold and 10 CV result
SDMs <- c("rf","gam","maxent","ensemble")
nTree <- 100
# arrray for each metrics
arr_hss <- array(NA,c(9, 5, 4,10),list(c(seq(from = 0.1, to = 0.9, by = 0.1)), c(seq(from = 1, to = 5, by = 1)), c(SDMs), c(seq(from = 1, to = 10, by = 1))))
arr_auc <- array(NA,c(1,5,4,10),list(c(1), c(seq(from = 1, to = 5, by = 1)), c(SDMs),c(seq(from = 1, to = 10, by = 1))));
arr_tss <- arr_pofd <- arr_pod <- arr_orss <- arr_hss

#array for storing 10 CV result
arr_hss_10 <- array(NA,c(9,1,4), list(c(seq(from = 0.1, to = 0.9, by = 0.1)), c(1), c(SDMs)))
arr_auc_10 <- array(NA,c(1,1,4), list(c(1), c(1), c(SDMs)))
arr_tss_10<- arr_pofd_10 <- arr_pod_10 <- arr_orss_10 <- arr_hss_10

#array for storing 5 fold and 10 CV result
for (i in 1:10){#repeat 10 times of the random CV fold assignment
  #i=1
  sb <- spatialBlock(speciesData = pt,
                     species = "lpp",
                     rasterLayer = mp,
                     theRange = 1000000, # size of the blocks
                     k = 5,
                     selection = "random",
                     iteration = 100, # find evenly dispersed folds
                     biomod2Format = TRUE,
                     xOffset = 0, # shift the blocks horizontally
                     yOffset = 0)
  sb$plots + geom_sf(data = pt, alpha = 0.5)
  dat2$fold <- sb$foldID;folds <- dat2$fold
  for (k in 1:5){#repeat CV for each of the 5 fold
    #k=5
    train.data <- dat2[which(folds != k),]# training set indices
    test.data <- dat2[which(folds == k),]
    dat_p <- subset(train.data, pa == 1)
    dat_a <- subset(train.data, pa == 0)
    for (SDM in SDMs){
      if(SDM == "maxent"){
        dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)*0.4),])
        x <- dat_sample[,varList_max] #train.data[,varList_max]
        p <- dat_sample[,c("pa")] #make p/a data a vector, without log lat.
        Maxent <- maxent(x, p=p, silent=TRUE)
        pr_m <- predict(Maxent, test.data, "prob")
        predictions <- pr_m
      } 
      else if(SDM == "rf"){
        dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)*0.5),])
        #mf <- mcmfRF2(dat_p, dat_a, nr=0.5, varList=varList_df, yCol='pa', reg=F, nTree=100, nForest=10)
        for(f in 1:10){ 
          library("foreach"); library("doSNOW");library(dismo);library(randomForest)
          registerDoSNOW(makeCluster(7, type="SOCK"))
          x2 <- dat_sample[,varList_rf]
          #class
          y2 <- as.factor(as.matrix(dat_sample[,'pa']))
          rf2 <- foreach(ntree = rep(floor(nTree/7),7), .combine = randomForest::combine, .packages = "randomForest") %dopar% randomForest(x2, y2, ntree = nTree, importance=T)
          if(f==1){rfC <- rf2}else{rfC <- randomForest::combine(rfC,rf2)}
        }
        mf <- rfC
        pr_rf <- predict(mf, test.data[,varList_rf], "prob")
        predictions <- pr_rf[,2]
      }
      else if(SDM == "ensemble"){
        pr_en <- rowMeans(cbind(pr_rf[,2],pr_m))
        predictions <- pr_en
      }
      else if(SDM == "gam"){    ####GAM
        gc()
        library(mgcv)
        dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)),])
        g1 <- bam(pa~s(mat, bs="tp")+ s(mwmt, bs="tp")+ s(map, bs="tp"), family=binomial, data=dat_sample, nthreads=6, discrete = TRUE)
        pr_g <- predict(g1, test.data)
        pr_g[pr_g > 0]<- 0
        predictions <- exp(pr_g)
      }
      arr_auc_lpp[1,k,c(SDM),i] <- AUC(predictions, test.data$pa)
      ############## contingency table metrics ------ error of integer flow might show up for Df data, in that case, use line193-207 to replace line174-187.
      for (j in c(1:9)){
        #SDM = 'gam'
        predictions_t <- predictions
        predictions_t[predictions_t<0.1*j] <- 0; predictions_t[predictions_t>=0.1*j] <- 1
        tbs <- table.stats(test.data$pa, predictions_t)
        #arr_acc_lpp[j,k,c(SDM),i] <- (tbs$tab[1,1]+tbs$tab[2,2])/sum(tbs$tab)
        arr_hss_lpp[j,k,c(SDM),i] <- tbs$HSS #Kappa
        #arr_pod_lpp[j,k,c(SDM),i] <- tbs$POD #sensitivity
        #arr_pofd_lpp[j,k,c(SDM),i] <- 1 - tbs$F #specificity
        arr_tss_lpp[j,k,c(SDM),i] <- tbs$POD + tbs$F - 1 #TSS
        arr_orss_lpp[j,k,c(SDM),i] <- tbs$orss #Yule's Q
        #arr_ts_lpp[j,k,c(SDM),i] <- table.stats(test.data$PINUCON, predictions_t)$TS
        #arr_sr_lpp[j,k,c(SDM),i] <- 1 - table.stats(test.data$PINUCON, predictions_t)$FAR
      }
    }  
  }
}

for(SDM in SDMs){
  arr_hss_10[,,c(SDM)] <- rowMeans(arr_hss_lpp[,,c(SDM),1])
  arr_orss_10[,,c(SDM)] <- rowMeans(arr_orss_lpp[,,c(SDM),1])
  arr_tss_10[,,c(SDM)] <- rowMeans(arr_tss_lpp[,,c(SDM),1])
}

metrics <- list(arr_hss_10,arr_tss_10, arr_orss_10)
metric_list = list()
for (l in 1:3){
  aaa <- data.frame(metrics[[l]][,,c("rf")], metrics[[l]][,,c("maxent")], metrics[[l]][,,c("gam")], metrics[[l]][,,c("ensemble")])
  names(aaa) = c("RF","MaxEnt","GAM","Ensemble")
  metric_list[[l]]=aaa
}
names(metric_list)<- c("HSS","TSS","ORSS")

###write all metric tables into csv
out_dir <- "p_a_data/acc_result/sbCV/new-pa-ratio/Pl_sbCV_three_model_"
for (i in 1:length(metric_list)){
  metric_list[[i]]$prob <- seq(from = 0.1, to = 0.9, by = 0.1)
  write.csv(metric_list[[i]], file = paste0(out_dir,names(metric_list[i]),".csv"), row.names = FALSE)
}

#ploting all metrics- Figure3 in manuscript
for (i in 1:(length(metric_list))){
  df <- melt(metric_list[[2]], id.vars = "prob", variable.name = 'model')
  p1 <- ggplot(df, aes(prob,value)) + geom_line(aes(colour = model), size = 1.5) +
    ggtitle(paste0(names(metric_list[2]))) + xlab("threshold")+ scale_color_manual(values=c("steelblue4", "aquamarine3", "azure3", "lightsalmon")) +
    theme(legend.position= "right", 
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          axis.text.x = element_text( size=14, color = 'black'),
          axis.title.x = element_text(size=14),
          axis.text.y = element_text( size=14, color = 'black'),
          axis.title.y = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.key=element_blank(),
          legend.key.width=unit(1.5, "cm")) + #xlim(0, 1)+ #ylim(0,1)+
    scale_x_continuous(breaks=seq(from = 0.1, to = 0.9, by = 0.1))
  print(p1)
}


######calculate global sensitivity and specificity for all models - Pl - correspond to Figure 4
######generate global pseudo-absence points
dat <- rasterToPoints(stk_w[[1]])%>%data.frame() #convert raster to points
colnames(a1)[1] <- "x";colnames(a1)[2] <- "y"
dat$pa <- 0
a1$pa <- 1
# combine the presence and absence data, and remove overlapping absence data points---
cmb <- rbind(a1, dat[,names(a1)]); hd(cmb)
cmb$x2 <- round(cmb$x,1); #make a new variable x2 from rounded lon
cmb$y2 <- round(cmb$y,1);head(cmb) #make a new variable x2 from rounded lat
#drop the data points falling into the same grid at 0.1 degree (10 km)
cmb2 <- cmb[!(duplicated(cmb[,c('x2','y2')])|duplicated(cmb[,c('x2','y2')])), ];hd(cmb2)
a_ab <- subset(cmb2,pa==0)[,1:2];#re-generate the absence data

######presence points
a0 <- read.csv("World data/all_wd_obs.csv")
#a0 <- subset(a,Longitude>-100)
a1 <- data.frame(a0$Longitude,a0$Latitude);colnames(a1) <- c("x", "y")
a1 <- na.omit(a1)

SDMs <- c("rf","gam","maxent","ensemble")
df_Pl <- data.frame(seq(from = 0.1, to = 1, by = 0.1), matrix(nrow = 10, ncol = 4));colnames(df_Pl) <- c("prob_Pl", "rf", "maxent", "gam", "ensemble")

dat_original <- read.csv("p_a_data/pa_wo_dup_yk1.csv")
names(dat_original[,2:91]) <- tolower(names(dat_original[,2:91]))
dat_full <- read.csv("p_a_data/Pl_all_abs.csv")#;subset(dat_lpp,pa==0)%>%nrow()#read.csv("p_a_data/Pl_all_abs.csv")##
dat_full_ab <- subset(dat_full,pa==0);nrow(dat_full_ab)
dat_lpp <- bind(dat_full_ab, dat_original)
dat_p <- subset(dat_lpp,pa==1)
dat_a <- subset(dat_lpp,pa==0)

for (SDM in SDMs){
  if(SDM == "maxent"){
    dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)),])
    x <- dat_sample[,varList_max]
    p <- dat_sample[,c("pa")] #make p/a data a vector, without log lat.
    Maxent <- maxent(x, p=p, silent=TRUE)
    predictions <- predict(stk_w,Maxent,type='prob',index=2)
  } 
  else if(SDM == "rf"){
    mf <- mcmfRF2(dat_p, dat_a, nr=1, varList=varList_rf, yCol='pa', reg=F, nTree=100, nForest=10)
    predictions <- predict(stk_w,mf,type='prob',index=2)
  }
  else if(SDM == "ensemble"){
    predictions <- calc(stack(pr_rf,pr_m), fun = mean)
  }
  else if(SDM == "gam"){
    gc()
    library(mgcv)
    dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)),])
    g1 <- bam(pa~s(mat, bs="tp")+ s(mwmt, bs="tp")+ s(map, bs="tp"), family=binomial, data=dat_sample, nthreads=6, discrete = TRUE)
    pr <- predict(stk_w, g1)
    pr[pr>0]<- 0
    pr_g <- exp(pr)
    predictions <- pr_g
    predictions <- raster("p_a_data/tif/adjusted_pa_ratio_pre/Pl_gam.tif")
  }
  for (i in seq(0,10,1)){
    model0 <- predictions
    model0[model0<=(i/10)]<-0;#plot(model0)
    pt <- raster::extract(model0, a_ab, method='simple', df=TRUE, cellnumbers=TRUE)
    df_Pl[i,SDM] <- nrow(pt[pt[,3] <= (i/10), ])/nrow(a_ab)#specificity
    #df_Pl[i,SDM] <- nrow(pt[pt[,3] <= (i/10), ])/nrow(a1)#sensitivity
  }
}  
plot(pr_g)

write.csv(df_Pl, 'p_a_data/threshold/df_Pl_overlap_rate_gam_updated_specificity.csv')#write sensitivity and specificity


#Pl
df_Pl_p <- read.csv('p_a_data/threshold/df_Pl_overlap_rate_gam_updated.csv')[,2:5]
df_Pl_a <- read.csv('p_a_data/threshold/df_Pl_overlap_rate_gam_updated_specificity.csv')[,2:5]

#df_Pl <- cbind(df_Pl_p, df_Pl_a)
Pl_a_m <- reshape2::melt(df_Pl_a, id.vars = "prob_Pl")
Pl_p_m <- reshape2::melt(df_Pl_p, id.vars = "prob_Pl")
Pl_m <- rbind(Pl_p_m,Pl_a_m)

###########plot figure4 in manuscript
ggplot(Pl_m, aes(x = prob_Pl, y= value)) +geom_line(aes(colour=variable,linetype = variable),size = 1) +
  xlab("threshold")+ ylab("Value")+
  scale_color_manual(values=c(rep(c("steelblue4", "aquamarine3", "lightsalmon"),2)))+
  scale_linetype_manual(values = c(rep("solid", 3), rep("longdash", 3))) +
  theme(legend.position= "right", 
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        axis.text.x = element_text( size=14, color = 'black'),
        axis.title.x = element_text(size=14),
        axis.text.y = element_text( size=14, color = 'black'),
        axis.title.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.key.width=unit(1.5, "cm"),
        aspect.ratio = 1/1
  ) + #xlim(0, 1)+ #ylim(0,1)+
  scale_x_continuous(breaks=seq(from = 0.1, to = 0.9, by = 0.1))


######final model (ensemble model) training - current and future prediction
##########lodgepole pine
dat_original <- read.csv("p_a_data/pa_wo_dup_yk1.csv")
names(dat_original[,2:91]) <- tolower(names(dat_original[,2:91]))
dat_full <- read.csv("p_a_data/Pl_all_abs.csv")#;subset(dat_lpp,pa==0)%>%nrow()#read.csv("p_a_data/Pl_all_abs.csv")##
dat_full_ab <- subset(dat_full,pa==0);nrow(dat_full_ab)
dat_lpp <- bind(dat_full_ab, dat_original)
dat_p <- subset(dat_lpp,pa==1)
dat_a <- subset(dat_lpp,pa==0)

#RF
#mf <- mcmfRF2(dat_p, dat_a, nr=1, varList=varList_rf, yCol='pa', reg=F, nTree=100, nForest=10)
for(f in 1:10){ 
  library("foreach"); library("doSNOW");library(dismo);library(randomForest)
  registerDoSNOW(makeCluster(7, type="SOCK"))
  dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)*0.5),])
  x2 <- dat_sample[,varList_rf]
  #class
  y2 <- as.factor(as.matrix(dat_sample[,'pa']))
  rf2 <- foreach(ntree = rep(floor(nTree/7),7), .combine = randomForest::combine, .packages = "randomForest") %dopar% randomForest(x2, y2, ntree = ntree, importance=T)
  if(f==1){rfC <- rf2}else{rfC <- randomForest::combine(rfC,rf2)}
}
gc()
mf <- rfC

#save(mf,file = "final_RF_og+full_Pl.RData")
#load("final_RF_og+full_Pl.RData")
rf_w <- predict(stk_w,mf,type='prob',index=2)
rf_w_45 <- predict(stk_w_45,mf,type='prob',index=2)
plot(rf_w);plot(rf_w_45)

######wna
rf_wna <- predict(stk_wna,mf,type='prob',index=2)
rf_wna_45 <- predict(stk_wna_45,mf,type='prob',index=2)
plot(rf_wna, xlim = c(-150,-105),ylim = c(40,65));plot(rf_wna_45)

#maxent    
dat_sample <- bind(dat_p,dat_a[sample(nrow(dat_a), nrow(dat_p)*0.3),])
x <- dat_sample[,varList_max] 
p <- dat_sample[,c("pa")]
Maxent <- maxent(x, p=p, silent=TRUE)
#save(Maxent,file = "final_MAX_og+full_Pl.RData")
#load()
max_w <- predict(Maxent, stk_w)
max_w_45 <- predict(Maxent, stk_w_45)
plot(max_w);plot(max_w_45)

ens_stk <- stack(max_w,rf_w);ens_stk_45 <- stack(max_w_45,rf_w_45)
en_w <- calc(ens_stk, fun = mean);en_w_45 <- calc(ens_stk_45, fun = mean)
plot(en_w);plot(en_w_45)

#writeRaster(en_w, "p_a_data/tif/full+og_ab_pred/Pl_ensemble_final.tif")
#writeRaster(en_w_45, "p_a_data/tif/full+og_ab_pred/Pl_ensemble_final_45.tif")
#randomForest()

