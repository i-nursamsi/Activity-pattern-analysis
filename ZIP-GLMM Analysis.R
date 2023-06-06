## RULES OF THUMB 
## Divide # of obs by 25 and that's the max # of variables in a model 
# Commpn palm civet has ~100 obs, so that means 5 total variables 
# Landscape and Effort = 2 variables, so that leaves 3 more covariates available 

setwd("C:/Users/ilyas/Dropbox/UQ summer 2020 species projects/Offline GLMMs and species projects raw datasets - summer species 2020/Offline_glmms_Sunda_pangolin_ilyas/GLMM_ZIP_No_Singapore")

# Load libraries and packages


library(dplyr); library(data.table); library(lme4); library(lmerTest); library(ggplot2); library(sjPlot)
library(plyr); library(MuMIn); library(GLMMadaptive);library(AER);library(DMwR2);
library(merTools); library(lattice); library(sjlabelled); library(ggeffects);library(here)

#====================   EXTRACT PRESENCE/ABSENCE DATA   ====================


# Extract occurrence records within ranges
All_sp_data = read.csv("Full_Data_by_Survey_n.0.NA_site2_20211025.csv")
my_species= 'Manis_javanica' #name your species
My_sp_data_with_NAs = All_sp_data[All_sp_data$Species == my_species,]
My_sp_data_NO_NAs=My_sp_data_with_NAs[!is.na(My_sp_data_with_NAs$records),];Data=My_sp_data_NO_NAs
dim(All_sp_data);dim(My_sp_data_with_NAs);dim(My_sp_data_NO_NAs);dim(Data)
Data$records = as.numeric(Data$records)

head(Data)
Data$RAI=Data$records/(Data$effort/100)
boxplot(Data$RAI~Data$country)
boxplot(Data$records~Data$country)
boxplot(Data$HumanFootprint20K~Data$country)
Data_No_SG=Data[Data$country!="Singapore",]

#xxx
# ggeffects crashes when logs and scales are INSIDE the regression, so transform them first:
Data_No_SG$HumanFootprint20K_scale=scale(Data_No_SG$HumanFootprint20K)
Data_No_SG$HumanFootprint20K_scale_log=scale(log(Data_No_SG$HumanFootprint20K+1))
Data_No_SG$AvgElevation20K_scale=scale(Data_No_SG$AvgElevation20K)
Data_No_SG$effort_log_scale=scale(log(Data_No_SG$effort+1))
Data_No_SG$OilPalm20K_scale=scale(Data_No_SG$OilPalm20K)
Data_No_SG$ForestCover20K_scale=scale(Data_No_SG$ForestCover20K) 
Data_No_SG$Forest_integrity20K_scale=scale(Data_No_SG$Forest_integrity20K)
Data_No_SG$Roughness20K_scale=scale(Data_No_SG$Roughness20K)
Data_No_SG$AvgSlope20K_scale=scale(Data_No_SG$AvgSlope20K)
Data_No_SG$Forest_intactness20K_scale=scale(Data_No_SG$Forest_intactness20K)
Data_No_SG$NightLights20K_scale=scale(Data_No_SG$NightLights20K)
Data_No_SG$HumanPop20K_log_scale=scale(log(Data_No_SG$HumanPop20K+1))
Data_No_SG$HumanPop20K_log=log(Data_No_SG$HumanPop20K+1)



# tests for best "null model" with effort and poisson or ZIP
# ZIP = zero-inflated poisson
reg_poisson0_No_SG=mixed_model(fixed = records ~ 1+ scale(effort), random = ~1 | Landscape, data = Data_No_SG, family = poisson(),zi_fixed = NULL, zi_random = NULL) #normal poisson
ZIP_No_SG=mixed_model(fixed = records ~ 1+ scale(effort), random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1, zi_random = NULL) #zero inflated poisson
reg_poisson_log_effort_No_SG=mixed_model(fixed = records ~ 1+ effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = poisson(),zi_fixed = NULL, zi_random = NULL) #normal poisson
ZIP_log_effort_No_SG=mixed_model(fixed = records ~ 1+ effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1, zi_random = NULL) #zero inflated poisson
ZIP_effort2_No_SG=mixed_model(fixed = records ~ 1+ scale(effort)+I(scale(effort)^2), random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1, zi_random = NULL) #zero inflated poisson, quadratic effort
ZIP_log_effort2_No_SG=mixed_model(fixed = records ~ 1+ effort_log_scale+I(effort_log_scale^2), random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1, zi_random = NULL) #zero inflated poisson, quadratic logged effort
#test effect of random effect in the ZIP part of the model
ZIP_log_effort_rand_zi_Landscape_No_SG=mixed_model(fixed = records ~ 1+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1, zi_random = ~1 | Landscape) #zero inflated poission w/ random effect on ZIP

# compile and assess AIC values
AICc_reg_poisson0l_No_SG=AICc(reg_poisson0_No_SG);AICc_ZIP=AICc(ZIP_No_SG);AICc_reg_poisson_log_effort=AICc(reg_poisson_log_effort_No_SG);AICc_ZIP_log_effort=AICc(ZIP_log_effort_No_SG);AICc_ZIP_effort2=AICc(ZIP_effort2_No_SG);AICc_ZIP_log_effort2=AICc(ZIP_log_effort2_No_SG);AICc_ZIP_log_effort_rand_zi_Landscape=AICc(ZIP_log_effort_rand_zi_Landscape_No_SG)
All_AIC_Zips_No_SG=c(AICc_reg_poisson0l_No_SG,AICc_ZIP,AICc_reg_poisson_log_effort,AICc_ZIP_log_effort,AICc_ZIP_effort2,AICc_ZIP_log_effort2,AICc_ZIP_log_effort_rand_zi_Landscape)
mods_No_SG=c('reg_poisson0l','ZIP','reg_poisson_log_effort','ZIP_log_effort','ZIP_effort2','ZIP_log_effort2','ZIP_log_effort_rand_zi_Landscape')
length(mods_No_SG);length(All_AIC_Zips_No_SG)
mod_sel_ZIP_No_SG=as.data.frame(cbind(mods_No_SG,All_AIC_Zips_No_SG))
colnames(mod_sel_ZIP_No_SG)=c("model","AICc")
mod_sel_ZIP_No_SG=mod_sel_ZIP_No_SG[order(mod_sel_ZIP_No_SG$AICc),]
mod_sel_ZIP_No_SG
summary(ZIP_log_effort_No_SG) #best null model, use this moving forward
write.csv(ZIP_log_effort_No_SG,"ZIP_model_selection_No_SG.csv")

################# Main covariate testing #######################
####
#### Single Variable Model Testing (all variables at 20Km radius) in order of AIC



Null=mixed_model(fixed = records ~ 1+ effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Elevation_No_SG=mixed_model(fixed = records ~ AvgElevation20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_No_SG=mixed_model(fixed = records ~ HumanFootprint20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_log_No_SG=mixed_model(fixed = records ~ HumanFootprint20K_scale_log + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
OilPalm_No_SG=mixed_model(fixed = records ~ OilPalm20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
TreeCover_No_SG=mixed_model(fixed = records ~ ForestCover20K_scale+ effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
ForestIntegrity_No_SG=mixed_model(fixed = records ~ Forest_integrity20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Roughness_No_SG=mixed_model(fixed = records ~ Roughness20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Slope_No_SG=mixed_model(fixed = records ~ AvgSlope20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanPop_log_No_SG=mixed_model(fixed = records ~ HumanPop20K_log_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1) 
ForestIntactness_No_SG=mixed_model(fixed = records ~ Forest_intactness20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
NightLights_No_SG=mixed_model(fixed = records ~ NightLights20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)

#### Nonlinear Univariate Model Testing

Elevation2_No_SG=mixed_model(fixed = records ~ AvgElevation20K_scale + I(AvgElevation20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint2_No_SG=mixed_model(fixed = records ~ HumanFootprint20K_scale + I(HumanFootprint20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_log2_No_SG=mixed_model(fixed = records ~ HumanFootprint20K_scale_log + I(HumanFootprint20K_scale_log^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
OilPalm2_No_SG=mixed_model(fixed = records ~ OilPalm20K_scale + I(OilPalm20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
TreeCover2_No_SG=mixed_model(fixed = records ~ ForestCover20K_scale+ I(scale(ForestCover20K)^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
ForestIntegrity2_No_SG=mixed_model(fixed = records ~ Forest_integrity20K_scale + I(Forest_integrity20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Roughness2_No_SG=mixed_model(fixed = records ~ Roughness20K_scale + I(Roughness20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Slope2_No_SG=mixed_model(fixed = records ~ AvgSlope20K_scale + effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanPop_log2_No_SG=mixed_model(fixed = records ~ HumanPop20K_log_scale + I(HumanPop20K_log_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1) 
ForestIntactness2_No_SG=mixed_model(fixed = records ~ Forest_intactness20K_scale + I(Forest_intactness20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
NightLights2_No_SG=mixed_model(fixed = records ~ NightLights20K + I(NightLights20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)


##### univariate model selection loop
# compile univariate models
mods=c('Null','Elevation','HumanFootprint','OilPalm','TreeCover','ForestIntegrity','Roughness','HumanPop_log','ForestIntactness','NightLights')
mods_list=list(Null,Elevation_No_SG,HumanFootprint_No_SG,OilPalm_No_SG,TreeCover_No_SG,ForestIntegrity_No_SG,Roughness_No_SG,HumanPop_log_No_SG,ForestIntactness_No_SG,NightLights_No_SG)
names(mods_list)=mods

# create dataframe to fill
colnamestoextract=c('Model','beta','P_value','AICc','LogLik')
empty=as.data.frame(matrix(NA,nrow=length(mods),ncol=5))
colnames(empty)=colnamestoextract
empty$Model=mods
head(empty)
# loop to pull model selection info

for(i in 1:length(mods)){
  temp=mods_list[empty$Model[i]][[1]] # this make sure it gets the right models even if the order changes
  empty$beta[i]=summary(temp)$coef_table[2,1]
  empty$P_value[i]=summary(temp)$coef_table[2,4]
  empty$AICc[i]=AICc(temp)
  empty$LogLik[i]=logLik(temp)[1]
  # df=df(temp) # would be great to pull df later
}
empty
# caculate some model weights
ModSel=empty
ModSel=ModSel[order(ModSel$AICc),]
ModSel$deltaAICc=ModSel$AICc-ModSel$AICc[1]
ModSel$relative_liklihood=exp(ModSel$deltaAICc*-0.5)
ModSel$weight=ModSel$relative_liklihood/sum(ModSel$relative_liklihood)
ModSel=ModSel[,colnames(ModSel)!='relative_liklihood']
ModSel[,c(2,4:dim(ModSel)[2])]=round(ModSel[,c(2,4:dim(ModSel)[2])],2)
ModSel[,3]=round(ModSel[,3],3);ModSel[,4:5]=round(ModSel[,4:5],1)

ModSel
# save it 
write.csv(ModSel,paste("UNIvariate model selection_No_SG",my_species,".csv"))

# now check to see if coefficients make sense

# BEST Univariate model is Forestintegrity
summary(ForestIntegrity_No_SG)
# Negative relationship fits our ecological knowledge that this is a edge/disturabance adapted species
summary(HumanPop_log2_No_SG)
# next best variable is HumanPop_log2, which is probably not worth including for two reasons:
# First, this capture similar information as HumanFootprint and is driven by Singapore
# Second, it's nonlinear meaning we'd need to add 2 more terms (HumanPop_log AND HumanPop_log2) and that's overfitting = too complicated


##################################  
########  OPTIONAL  ############
###   Mulivariate models #######
##################################  

# choose top 2 variables only, and only those variables you like
# can't put in highly correlated variables together, so not two human impact varibale
# only one of Human Footprint, Forest Integrity, Forest Intactness, Human population, etc 
# And only one of elevation, roughness
ForestIntactness_Elevation_No_SG = mixed_model(fixed = records ~ Forest_intactness20K_scale + AvgElevation20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
ForestIntegrity_Elevation_No_SG = mixed_model(fixed = records ~ Forest_integrity20K_scale + AvgElevation20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_Elevation_No_SG = mixed_model(fixed = records ~ HumanFootprint20K_scale + AvgElevation20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
ForestIntactness_Elev2_No_SG = mixed_model(fixed = records ~ Forest_intactness20K_scale + AvgElevation20K_scale+ I(AvgElevation20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
ForestIntegrity_Elev2_No_SG = mixed_model(fixed = records ~ Forest_integrity20K_scale + AvgElevation20K_scale+ I(AvgElevation20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_Elev2_No_SG = mixed_model(fixed = records ~ HumanFootprint20K_scale + AvgElevation20K_scale+ I(AvgElevation20K_scale^2)+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_Roughness_No_SG = mixed_model(fixed = records ~ HumanFootprint20K_scale + Roughness20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Forest_integrity20K_Roughness_No_SG = mixed_model(fixed = records ~ Forest_integrity20K_scale + Roughness20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
#HumanFootprint_Roughness = mixed_model(fixed = records ~ HumanFootprint20K_scale + Roughness20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data, family = zi.poisson(), zi_fixed = ~1)
Forest_integrity_Roughness_No_SG = mixed_model(fixed = records ~ Forest_integrity20K_scale + Roughness20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
#ForestIntactness_OilPalm = mixed_model(fixed = records ~ Forest_intactness20K_scale + OilPalm20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data, family = zi.poisson(), zi_fixed = ~1)
ForestIntegrity_OilPalm_No_SG = mixed_model(fixed = records ~ Forest_integrity20K_scale + OilPalm20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
HumanFootprint_OilPalm_No_SG = mixed_model(fixed = records ~ HumanFootprint20K_scale + OilPalm20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Elevation_OilPalm_No_SG = mixed_model(fixed = records ~ AvgElevation20K_scale + OilPalm20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)
Roughness_OilPalm_No_SG = mixed_model(fixed = records ~ Roughness20K_scale + OilPalm20K_scale+effort_log_scale, random = ~1 | Landscape, data = Data_No_SG, family = zi.poisson(), zi_fixed = ~1)

mods_MV=c('ForestIntegrity_Elevation','ForestIntegrity_Elevation','HumanFootprint_Elevation','ForestIntactness_Elev2','ForestIntegrity_Elev2','HumanFootprint_Elev2','Forest_integrity_Roughness','ForestIntegrity_OilPalm','HumanFootprint_OilPalm','Elevation_OilPalm','Roughness_OilPalm','Elevation2','HumanFootprint2','OilPalm2','TreeCover2','ForestIntegrity2','Roughness2','HumanPop_log2','ForestIntactness2','NightLights2')
mods_MV_list=list(ForestIntegrity_Elevation_No_SG,ForestIntegrity_Elevation_No_SG,HumanFootprint_Elevation_No_SG,ForestIntactness_Elev2_No_SG,ForestIntegrity_Elev2_No_SG,HumanFootprint_Elev2_No_SG,Forest_integrity_Roughness_No_SG,ForestIntegrity_OilPalm_No_SG,HumanFootprint_OilPalm_No_SG,Elevation_OilPalm_No_SG,Roughness_OilPalm_No_SG,Elevation2_No_SG,HumanFootprint2_No_SG,OilPalm2_No_SG,TreeCover2_No_SG,ForestIntegrity2_No_SG,Roughness2_No_SG,HumanPop_log2_No_SG,ForestIntactness2_No_SG,NightLights2_No_SG)
names(mods_MV_list)=mods_MV
#### re-run Model Selection for MV
##### Model selection loop

# create dataframe to fill
MVcolnamestoextract=c('Model','AICc','LogLik')
MVempty=as.data.frame(matrix(NA,nrow=length(mods_MV),ncol=3))
colnames(MVempty)=MVcolnamestoextract
MVempty$Model=mods_MV
head(MVempty)
# loop to pull model selection info
i=1
for(i in 1:length(mods_MV)){
  temp=mods_MV_list[MVempty$Model[i]][[1]] # this make sure it gets the right models even if the order changes
  #empty$beta[i]=summary(temp)$coef_table[2,1]
  #empty$P_value[i]=summary(temp)$coef_table[2,4]
  MVempty$AICc[i]=AICc(temp)
  MVempty$LogLik[i]=logLik(temp)[1]
  # df=df(temp) # would be great to pull df later
}

# caculate some model weights
ModSel_MV=MVempty
ModSel_full=rbind(ModSel[,c('Model','AICc','LogLik')],ModSel_MV)
ModSel_full=ModSel_full[order(ModSel_full$AICc),]
ModSel_full$deltaAICc=ModSel_full$AICc-ModSel_full$AICc[1]
ModSel_full$relative_liklihood=exp(ModSel_full$deltaAICc*-0.5)
ModSel_full$weight=ModSel_full$relative_liklihood/sum(ModSel_full$relative_liklihood)
ModSel_full=ModSel_full[,colnames(ModSel_full)!='relative_liklihood']
ModSel_full[,2:dim(ModSel_full)[2]]=round(ModSel_full[,2:dim(ModSel_full)[2]],2)
ModSel_full
# save it 
write.csv(ModSel_full,paste("non-linear & multivariate model selection_No_SG",my_species,".csv"))

# results for common palm civiet explained below
# removing any term increases AICc by >20 points, so we can be sure this is best model
# the elevation quadratic term is negative (hump-shaped relationship)
# this DOES make biological/ecological sense, since the species might mid-elevations
# Positive relationships HumanFootprint also fits our ecological knowledge that this is a edge/disturabance adapted species
# for the 2nd best model, the negative relationships Forest_Intactness also fits the edge/disturabance adapted species
# neither models flip coefficents when excluding Singapore



####################################################################
####  Plotting the effect sizes from UNIvariate models #########
####  holding everything else constant                    #########
####################################################################


########################## 
#DOES NOT WORK FOR ME, USE THE NEXT SECTION (NON-LINEAR UNIVARIATE PLOTTING) to plot linear univariate models
# Choose the model based on its AIC rank in your model selection:
ModSel
mod_rank=1
### extract best univariate model
plot_Uni_model_from_ModSel=function(mod_rank){
  mod_name=ModSel[mod_rank,1]
  mod=mods_list[empty$Model==mod_name][[1]] 
  covariate=rownames(summary(mod)$coef_table)[2]
  sum_mod=summary(mod)
  
  # get ZIP model preditions 
  mod_pred=ggpredict(mod, terms = covariate, condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale)))
  
  # only run this for quadratic
  mod=mods_MV_list$Elevation2
  covariate="AvgElevation20K_scale"
  mod_pred=ggpredict(mod, terms = "AvgElevation20K_scale [all]", condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale)))
  
  ############################
  #backtransform scaled x-axis
  OG_covariate=gsub("_scale","",covariate)
  Mean=mean(na.exclude(Data_No_SG[,OG_covariate])); Sd=sd(na.exclude(Data_No_SG[,OG_covariate]))
  mod_pred$covar_back_trans= (mod_pred$x*Sd) + Mean
  mod_pred_backtrans=as.data.frame(cbind(mod_pred[[1]],mod_pred[[2]],mod_pred[[3]],mod_pred[[4]],mod_pred[[5]],mod_pred[[6]],mod_pred[[7]]))
  names(mod_pred_backtrans)=colnames(mod_pred);
  #str(mod_pred_backtrans)
  
  #################
  ####### dashed line for non-significant
  if(summary(mod_pred_backtrans)$coef_table[2,"p-value"] > 0.05){
    line_type = "dashed"
    line_color = "gray40"}else{
      line_type = "solid"
      line_color = "black" }
  Data_No_SG$Cols = ifelse(Data_No_SG$records > 0, "blue", "darkred")
  Data_No_SG$dotshape = ifelse(Data_No_SG$records > 0, "19", "1")
  plot_name=
    ggplot(mod_pred_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
    geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
    geom_line(size=1) +#geom_line(linetype=line_type, size=1, color = line_color) +
    geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
    coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate]), max(Data_No_SG[,OG_covariate])))+
    scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
    labs(x = OG_covariate, y = "Number of captures + 1")+
    theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
          legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_y_log10() #+ scale_x_log10()
  print(sum_mod)
  return(plot_name)
  plot_name
}

dir.create(file.path(getwd(),paste("glmmZIPfigs",my_species)), recursive = TRUE)

print_it=function(x){
  plot_name
  ggsave(filename=paste(ModSel[x,1],my_species,".png"), path = paste("glmmZIPfigs",my_species), width = 10, height = 10, units = "cm",dpi=500,scale=1.1)
  ggsave(filename=paste(ModSel[x,1],my_species,"_wide.png"), path = paste("glmmZIPfigs",my_species),width = 12, height = 10, units = "cm",dpi=500)
}

#for the quadratic model
#ggsave(filename=paste("Elevation2",my_species,".png"), path = paste("glmmZIPfigs",my_species), width = 10, height = 10, units = "cm",dpi=500,scale=1.1)
#ggsave(filename=paste("Elevation2",my_species,"_wide.png"), path = paste("glmmZIPfigs",my_species), width = 12, height = 12, units = "cm",dpi=500,scale=1.1)

#look at one:
plot_Uni_model_from_ModSel(6)
#print all of them, and will automatticaly save into your Offline_glmm species folder
for(i in 1:dim(ModSel)[1]){
  plot_Uni_model_from_ModSel(i)
  print_it(i)
}
# they will stop printing once they hit an error at the Null model


############################################################
####  Plotting the effect sizes           ##################
####  from NON-LINEAR-UNIvariate models   ##################
####  holding everything else constant    ##################
############################################################
##Can also be used to plot LINEAR univariate models if edited accordingly
#choose nonlinear model based
ModSel_full
mod2=ForestIntegrity_No_SG

mod2_name= "ForestIntegrity"
covariate=rownames(summary(mod2)$coef_table)[2]
sum_mod2=summary(mod2)

# set Elevation to the minimum to prevent captures from blowing up
mod_pred2=ggpredict(mod2, terms = covariate, condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale)))
#plot(predicted~x,data=mod_pred,type='o',col='darkgreen',lwd=2,ylim=c(0,20),xlim=c(-1,2))

############################
#backtransform scaled x-axis
OG_covariate=gsub("_scale","",covariate)
Mean=mean(na.exclude(Data_No_SG[,OG_covariate])); Sd=sd(na.exclude(Data_No_SG[,OG_covariate]))
mod_pred2$covar_back_trans= (mod_pred2$x*Sd) + Mean
mod_pred2_backtrans=as.data.frame(cbind(mod_pred2[[1]],mod_pred2[[2]],mod_pred2[[3]],mod_pred2[[4]],mod_pred2[[5]],mod_pred2[[6]],mod_pred2[[7]]))
names(mod_pred2_backtrans)=colnames(mod_pred2);
#str(mod_pred_backtrans)

#################
####### dashed line for non-significant
if(summary(mod2)$coef_table[2,"p-value"] > 0.05){
  line_type = "dashed"
  line_color = "gray40"}else{
    line_type = "solid"
    line_color = "black" }
Data_No_SG$Cols = ifelse(Data_No_SG$records > 0, "blue", "darkred")
Data_No_SG$dotshape = ifelse(Data_No_SG$records > 0, "19", "1")
plot_mod2=
  ggplot(mod_pred2_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
  geom_line(linetype=line_type, size=1, color = line_color) +#geom_line(size=1) 
  geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
  coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate]), max(Data_No_SG[,OG_covariate])))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
  labs(x = "Forest integrity 20 km", y = "Number of captures + 1")+
  theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() #+ scale_x_log10()

plot_mod2


###########

ModSel_full
mod4=HumanFootprint_No_SG

mod4_name= "HumanFootprint"
covariate=rownames(summary(mod4)$coef_table)[2]
sum_mod4=summary(mod4)

# set Elevation to the minimum to prevent captures from blowing up
mod_pred4=ggpredict(mod4, terms = covariate, condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale)))
#plot(predicted~x,data=mod_pred,type='o',col='darkgreen',lwd=2,ylim=c(0,20),xlim=c(-1,2))

############################
#backtransform scaled x-axis
OG_covariate=gsub("_scale","",covariate)
Mean=mean(na.exclude(Data_No_SG[,OG_covariate])); Sd=sd(na.exclude(Data_No_SG[,OG_covariate]))
mod_pred4$covar_back_trans= (mod_pred4$x*Sd) + Mean
mod_pred4_backtrans=as.data.frame(cbind(mod_pred4[[1]],mod_pred4[[2]],mod_pred4[[3]],mod_pred4[[4]],mod_pred4[[5]],mod_pred4[[6]],mod_pred4[[7]]))
names(mod_pred4_backtrans)=colnames(mod_pred4);
#str(mod_pred_backtrans)

#################
####### dashed line for non-significant
if(summary(mod4)$coef_table[2,"p-value"] > 0.05){
  line_type = "dashed"
  line_color = "gray40"}else{
    line_type = "solid"
    line_color = "black" }
Data_No_SG$Cols = ifelse(Data_No_SG$records > 0, "blue", "darkred")
Data_No_SG$dotshape = ifelse(Data_No_SG$records > 0, "19", "1")
plot_mod4=
  ggplot(mod_pred4_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
  geom_line(linetype=line_type, size=1, color = line_color) +#geom_line(size=1) 
  geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
  coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate]), max(Data_No_SG[,OG_covariate])))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
  labs(x = "Human footprint 20 km", y = "Number of captures + 1")+
  theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() #+ scale_x_log10()

plot_mod4

########### Human footprint

ModSel_full
mod3=OilPalm_No_SG

mod3_name= "OilPalm"
covariate=rownames(summary(mod3)$coef_table)[2]
sum_mod3=summary(mod3)

# set Elevation to the minimum to prevent captures from blowing up
mod_pred3=ggpredict(mod3, terms = covariate, condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale)))
#plot(predicted~x,data=mod_pred,type='o',col='darkgreen',lwd=2,ylim=c(0,20),xlim=c(-1,2))

############################
#backtransform scaled x-axis
OG_covariate=gsub("_scale","",covariate)
Mean=mean(na.exclude(Data_No_SG[,OG_covariate])); Sd=sd(na.exclude(Data_No_SG[,OG_covariate]))
mod_pred3$covar_back_trans= (mod_pred3$x*Sd) + Mean
mod_pred3_backtrans=as.data.frame(cbind(mod_pred3[[1]],mod_pred3[[2]],mod_pred3[[3]],mod_pred3[[4]],mod_pred3[[5]],mod_pred3[[6]],mod_pred3[[7]]))
names(mod_pred3_backtrans)=colnames(mod_pred3);
#str(mod_pred_backtrans)

#################
####### dashed line for non-significant
if(summary(mod3)$coef_table[2,"p-value"] > 0.05){
  line_type = "dashed"
  line_color = "gray40"}else{
    line_type = "solid"
    line_color = "black" }
Data_No_SG$Cols = ifelse(Data_No_SG$records > 0, "blue", "darkred")
Data_No_SG$dotshape = ifelse(Data_No_SG$records > 0, "19", "1")
plot_mod3=
  ggplot(mod_pred3_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
  geom_line(linetype=line_type, size=1, color = line_color) +#geom_line(size=1) 
  geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
  coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate]), max(Data_No_SG[,OG_covariate])))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
  labs(x = "Oil palm 20 km", y = "Number of captures + 1")+
  theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() #+ scale_x_log10()

plot_mod3
############################################################
####  Plotting the effect sizes           ##################
####  from NON-LINEAR-MULTI-variate models   ##################
####  holding everything else constant    ##################
############################################################


# Choose model

### extract best univariate model
MVmod_name=ModSel_full[1,1]
MVmod=mods_MV_list[MVempty$Model==MVmod_name][[1]] 
covariate1=rownames(summary(MVmod)$coef_table)[2]
covariate2=rownames(summary(MVmod)$coef_table)[3]

# CONDITIONS: set Elevation to the minimum to prevent captures from blowing up
# probably set others to median? 
mod_pred_cov1=ggpredict(MVmod, terms = "HumanFootprint20K_scale [all]", condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale),AvgElevation20K_scale = min(Data_No_SG$AvgElevation20K_scale)))
mod_pred_cov2=ggpredict(MVmod, terms = "AvgElevation20K_scale [all]", condition = c(effort_log_scale = median(Data_No_SG$effort_log_scale),HumanFootprint20K_scale = median(Data_No_SG$HumanFootprint20K_scale)))
#plot(predicted~x,data=mod_pred,type='o',col='darkgreen',lwd=2,ylim=c(0,20),xlim=c(-1,2))

############################
#backtransform scaled x-axis
OG_covariate1=gsub("_scale","",covariate1); OG_covariate2=gsub("_scale","",covariate2)
Mean1=mean(na.exclude(Data_No_SG[,OG_covariate1])); Sd1=sd(na.exclude(Data_No_SG[,OG_covariate1]))
Mean2=mean(na.exclude(Data_No_SG[,OG_covariate2])); Sd2=sd(na.exclude(Data_No_SG[,OG_covariate2]))
mod_pred_cov1$covar_back_trans= (mod_pred_cov1$x*Sd1) + Mean1
mod_pred_cov2$covar_back_trans= (mod_pred_cov2$x*Sd2) + Mean2

mod_pred_cov1_backtrans=as.data.frame(cbind(mod_pred_cov1[[1]],mod_pred_cov1[[2]],mod_pred_cov1[[3]],mod_pred_cov1[[4]],mod_pred_cov1[[5]],mod_pred_cov1[[6]],mod_pred_cov1[[7]]))
names(mod_pred_cov1_backtrans)=colnames(mod_pred_cov1);
mod_pred_cov2_backtrans=as.data.frame(cbind(mod_pred_cov2[[1]],mod_pred_cov2[[2]],mod_pred_cov2[[3]],mod_pred_cov2[[4]],mod_pred_cov2[[5]],mod_pred_cov2[[6]],mod_pred_cov2[[7]]))
names(mod_pred_cov2_backtrans)=colnames(mod_pred_cov2);

#################
# aethetics 
line_type = "solid" # if non-significnat (unlikely) then manually change to  line_type = "dashed"
line_color = "black"  # if non-significnat (unlikely) then manually change to line_color = "gray40"
Data_No_SG$Cols = ifelse(Data_No_SG$records > 0, "blue", "darkred")
Data_No_SG$dotshape = ifelse(Data_No_SG$records > 0, "19", "1")

plot_MV1=
  ggplot(mod_pred_cov1_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
  geom_line(size=1) +#geom_line(linetype=line_type, size=1, color = line_color) +
  geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate1], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
  coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate1]), max(Data_No_SG[,OG_covariate1])))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
  labs(x = OG_covariate1, y = "Number of captures + 1")+
  theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() #+ scale_x_log10()

plot_MV1

plot_MV2=
  ggplot(mod_pred_cov2_backtrans, aes(y=predicted+1, x=covar_back_trans)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.15)+ # errors not working atm
  geom_line(size=1) +#geom_line(linetype=line_type, size=1, color = line_color) +
  geom_jitter(data = Data_No_SG, aes(x = Data_No_SG[,OG_covariate2], y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=1),cex=1.2,width = .05, height = .05)+
  coord_cartesian(ylim = c(.9, max(Data_No_SG$records)),xlim = c(min(Data_No_SG[,OG_covariate2]), max(Data_No_SG[,OG_covariate2])))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(1, 19))+
  labs(x = OG_covariate2, y = "Number of captures + 1")+
  theme(axis.title = element_text(size = 18,color='black'), axis.text = element_text(size = 14,color='black'),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() #+ scale_x_log10()

plot_MV2
# Elevation plot starts at 0 because contingent on HumFoot
# All the high human footprints are at lowest elevation (e.g. Singapore)
plot(Data_No_SG$Elevation20K_scale~Data_No_SG$HumanFootprint20K_scale)



#########################################
## Try logging the axes

xlab = "Human Footprint Index" #name your x-axis
ggplot(HFP2, aes(y=predicted+.1, x=HumanFootprint_back_trans+.1)) + 
  geom_ribbon(aes(ymin = conf.low+.1, ymax = conf.high+.1), alpha = 0.1)+ # errors not working atm
  geom_line(linetype=line_type, size=1, color = line_color) +
  geom_jitter(data = Data_No_SG, aes(x = HumanFootprint20K+.1, y = records+.1, pch=dotshape,colour=Cols, alpha =.5,size=2),width = .1, height = .1)+
  coord_cartesian(ylim = c(.1, 100),xlim = c(1, 80))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(16, 18))+
  labs(x = xlab, y = "Number of captures")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_y_log10() # scale_x_log10() 


############################
#backtransform scaled x-axis (just don't do this NightLights)
Mean=mean(na.exclude(Data_No_SG$AvgAvgElevation20K)); Sd=sd(na.exclude(Data_No_SG$AvgAvgElevation20K))
Elev$Elev_back_trans= (Elev$x*Sd) + Mean
plot(predicted~Elev_back_trans,data=Elev,type='l',col='darkgreen',lwd=2,ylim=c(0,20))
points(records~AvgAvgElevation20K,data=Data_No_SG, col='red')

dim(Elev)
Elev2=as.data.frame(cbind(Elev[[1]],Elev[[2]],Elev[[3]],Elev[[4]],Elev[[5]],Elev[[6]],Elev[[7]]))
names(Elev2)=colnames(Elev)
str(Elev2)
xlab="Elevation"

ggplot(Elev2, aes(y=predicted+1, x=Elev_back_trans+1)) + 
  geom_ribbon(aes(ymin = conf.low+1, ymax = conf.high+1), alpha = 0.1)+ # errors not working atm
  geom_line(linetype=line_type, size=1, color = line_color) +
  geom_jitter(data = Data_No_SG, aes(x = AvgAvgElevation20K, y = records+1, pch=dotshape,colour=Cols, alpha =.5,size=2),width = .1, height = .1)+
  #coord_cartesian(ylim = c(0, 20))+
  scale_color_identity()+ scale_shape_manual(name = "Captures:", values = c(16, 18))+
  labs(x = xlab, y = "Number of captures")+
  theme(axis.title = element_text(size = 14), axis.text = element_text(size = 14),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10() + scale_y_log10() 


#########################
### End here
######################


save.image(file = "GLMM_Offline_20220203.RData")


load("GLMM_Offline_20220203.RData")
