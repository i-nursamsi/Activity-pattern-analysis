#### Compare Activity Patterns of single-species in degraded vs non-degraded habitats. 
### Create a loop to test all relevant species
## ILYAS NURSAMSI



library(tidyverse)
library(activity)
library(overlap)
library(lubridate)
#library(chron)

setwd("C:/Users/Ilyas Nursamsi/Dropbox/UQ summer 2020 species projects/Activity pattern analyses - summer species 2020/Sunda Pangolin - Activity Pattern Analysis")

###### Import data and format times #####

# Import
caps = read_csv("ECL captures summer spp projects_20210322.csv")
caps$Species[caps$Species=='Manidae_spp']='Manis_javanica'
str(caps)
head(caps)
dim(caps) #62675 rows


### Scale and convert capture times to radians

# Extract the number of second of each capture
caps$Time <- hms(caps$Photo.Time)
caps$time.sec <- as.numeric(seconds(caps$Time))

#standardize times from 0-n
caps$time.std <- caps$time.sec - min(caps$time.sec)
range(caps$time.std)

#convert to radians by dividing the number of seconds in a day 
#in order for the data to scale properly from 0 to 1
#Then multiply by 2*pi
caps$Time.rad <- (caps$time.std/86400)*2*pi
range(caps$Time.rad)

# trim the fat, dont need dates either
names(caps)
caps =subset(caps, select = -c(camera_start.date,camera_end.date,
                               Photo.Date, Time, time.sec, time.std))
head(caps)
dim(caps) #62675, all good 



## Import covaraite data
covs = read.csv("ECL_metadata_cam_level_summer_spp_20210322.csv")
names(covs) #these are the same variables as used in the occupancy model analysis
# but these are at the camera level (i.e. not spatially re-sampled)
# Need a landscape covaraite as well (?)

str(covs) #looks good

# Make sure all camera names and survey names match between captures and covariates
setdiff(caps$camera_id, covs$camera_id) #Dan20 cams
setdiff(covs$camera_id, caps$camera_id) #no difference

setdiff(caps$survey_id, covs$survey_id) #Danum_2020
setdiff(covs$survey_id, caps$survey_id) #no difference

# Because we do not have Danum_2020 covariates calculated yet, remove them from the captures so the data matches
caps = caps[caps$survey_id != "Danum_2020",]
#you can re-run the setdiff() functions again to see that there is no difference. 

#merge caps and covs together so we can subset caps later on the basis of covs
#acts = merge(caps, covs, by = c("camera_id", "survey_id"))
#head(data.frame(acts)) #looks good!



### Inspect which degradation variables will be included in the loop
names(covs)

#HFP
hist(covs$human_footprint)
summary(covs$human_footprint) #median = 4, mean = 6.3, set cut off at 5? or use median
summary(acts$human_footprint) #everything is larger, b/c repeated values

##FLII
hist(covs$forest_integrity)
summary(covs$forest_integrity) #median = 8.7, mean = 7.42--> use median?
summary(acts$forest_integrity) #slightly different values

##edges
hist(covs$dist_to_edge)
summary(covs$dist_to_edge) #median = 2553, mean = 3292..--> use median?

## degraded land
hist(covs$degraded_forest_1km)
summary(covs$degraded_forest_1km) #median = 0, mean = 14 --> use median

## forest land
hist(covs$forest_cover_1km)
summary(covs$forest_cover_1km) #median = 100, mean = 84 --> use median





##### Activity comparison loop #####

## species projects, can update this vector w/ more species later
sp = c( "Manis_javanica")

res = list() #store results here

for(i in 1:length(sp)){ #repeat for each species 
  
  # specify species and subset acts for it
  s = sp[i]
  a = caps[caps$Species == s,]
  
  # subset covs
  c = covs[covs$camera_id %in% unique(a$camera_id),]
  
  ## Create degredation statuses
  c$HFP_status= "No_Humans"
  c$HFP_status[c$human_footprint > median(c$human_footprint)]= "Humans"
  
  c$FLII_status = "Low_Integrity"
  c$FLII_status[c$forest_integrity > median(c$forest_integrity)]= "High_Integrity"
  
  c$edge_status = "Edge"
  c$edge_status[c$dist_to_edge > median(c$dist_to_edge)]= "Interior"
  
  c$degraded_status = "Non_Degraded"
  c$degraded_status[c$degraded_forest_1km > median(c$degraded_forest_1km)]= "Degraded"
  
  c$forest_status = "Forest"
  c$forest_status[c$forest_cover_1km < median(c$forest_cover_1km)]= "Non_Forest"
  
  ## merge covs w/ caps
  d = merge(a,c, by = c("camera_id", "survey_id"))
  
  
  ### Determine if activity patterns are significantly different 
  ## First fit actmods
  if(nrow(d)> 100){ #if more than 100 captures, sample = "model"
    
    t1 = fitact(d$Time.rad[d$HFP_status == "No_Humans"],sample = "model", adj = 1, reps = 1000)
    t2 = fitact(d$Time.rad[d$HFP_status == "Humans"],sample = "model", adj = 1, reps = 1000)
    t3 = fitact(d$Time.rad[d$FLII_status == "Low_Integrity"],sample = "model", adj = 1, reps = 1000)
    t4 = fitact(d$Time.rad[d$FLII_status == "High_Integrity"],sample = "model", adj = 1, reps = 1000)
    t5 = fitact(d$Time.rad[d$edge_status == "Edge"],sample = "model", adj = 1, reps = 1000)
    t6 = fitact(d$Time.rad[d$edge_status == "Interior"],sample = "model", adj = 1, reps = 1000)
    t7 = fitact(d$Time.rad[d$degraded_status == "Non_Degraded"],sample = "model", adj = 1, reps = 1000)
    t8 = fitact(d$Time.rad[d$degraded_status == "Degraded"],sample = "model", adj = 1, reps = 1000)
    t9 = fitact(d$Time.rad[d$forest_status == "Forest"],sample = "model", adj = 1, reps = 1000)
    t10 = fitact(d$Time.rad[d$forest_status == "Non_Forest"],sample = "model", adj = 1, reps = 1000)
  } else{ #if less than 100 captures, sample = "data"
    
    t1 = fitact(d$Time.rad[d$HFP_status == "No_Humans"],sample = "data", adj = 1, reps = 1000)
    t2 = fitact(d$Time.rad[d$HFP_status == "Humans"],sample = "data", adj = 1, reps = 1000)
    t3 = fitact(d$Time.rad[d$FLII_status == "Low_Integrity"],sample = "data", adj = 1, reps = 1000)
    t4 = fitact(d$Time.rad[d$FLII_status == "High_Integrity"],sample = "data", adj = 1, reps = 1000)
    t5 = fitact(d$Time.rad[d$edge_status == "Edge"],sample = "data", adj = 1, reps = 1000)
    t6 = fitact(d$Time.rad[d$edge_status == "Interior"],sample = "data", adj = 1, reps = 1000)
    t7 = fitact(d$Time.rad[d$degraded_status == "Non_Degraded"],sample = "data", adj = 1, reps = 1000)
    t8 = fitact(d$Time.rad[d$degraded_status == "Degraded"],sample = "data", adj = 1, reps = 1000)
    t9 = fitact(d$Time.rad[d$forest_status == "Forest"],sample = "data", adj = 1, reps = 1000)
    t10 = fitact(d$Time.rad[d$forest_status == "Non_Forest"],sample = "data", adj = 1, reps = 1000)
    
  }
  ##Next, test for sig differences per status
  r1 = data.frame(compareAct(list(t1,t2)))
  r1$test = "HFP_status"
  
  r2 = data.frame(compareAct(list(t3,t4)))
  r2$test = "FLII_status"
  
  r3 = data.frame(compareAct(list(t5,t6)))
  r3$test = "edge_status"
  
  r4 = data.frame(compareAct(list(t7,t8)))
  r4$test = "degraded_status"
  
  r5 = data.frame(compareAct(list(t9,t10)))
  r5$test = "forest_status"
  
  r = rbind(r1,r2,r3,r4,r5)
  
  
  ### Add overlap for each comparison
  r$overlap = 0
  
  #HFP
  o = overlapEst(d$Time.rad[d$HFP_status == "No_Humans"],
                 d$Time.rad[d$HFP_status == "Humans"])
  
  if(length(d$Time.rad[d$HFP_status == "No_Humans"]) > 75 & 
     length(d$Time.rad[d$HFP_status == "Humans"]) > 75){
    
    r$overlap[r$test == "HFP_status"] = o[2] #Dhat 4 for big sample size
    
  } else{
    r$overlap[r$test == "HFP_status"] = o[1] #Dhat 1 for small sample size
  }
  
  #FLII
  o = overlapEst(d$Time.rad[d$FLII_status == "Low_Integrity"],
                 d$Time.rad[d$FLII_status == "High_Integrity"])
  
  if(length(d$Time.rad[d$FLII_status == "Low_Integrity"]) > 75 & 
     length(d$Time.rad[d$FLII_status == "High_Integrity"]) > 75){
    
    r$overlap[r$test == "FLII_status"] = o[2]
    
  } else{
    r$overlap[r$test == "FLII_status"] = o[1]
  }
  
  #edge
  o = overlapEst(d$Time.rad[d$edge_status == "Edge"],
                 d$Time.rad[d$edge_status == "Interior"])
  
  if(length(d$Time.rad[d$edge_status == "Edge"]) > 75 & 
     length(d$Time.rad[d$edge_status == "Interior"]) > 75){
    
    r$overlap[r$test == "edge_status"] = o[2]
    
  } else{
    r$overlap[r$test == "edge_status"] = o[1]
  }
  
  #degraded
  o = overlapEst(d$Time.rad[d$degraded_status == "Non_Degraded"],
                 d$Time.rad[d$degraded_status == "Degraded"])
  
  if(length(d$Time.rad[d$degraded_status == "Non_Degraded"]) > 75 & 
     length(d$Time.rad[d$degraded_status == "Degraded"]) > 75){
    
    r$overlap[r$test == "degraded_status"] = o[2]
    
  } else{
    r$overlap[r$test == "degraded_status"] = o[1]
  }
  
  #forest
  o = overlapEst(d$Time.rad[d$forest_status == "Forest"],
                 d$Time.rad[d$forest_status == "Non_Forest"])
  
  if(length(d$Time.rad[d$forest_status == "Forest"]) > 75 & 
     length(d$Time.rad[d$forest_status == "Non_Forest"]) > 75){
    
    r$overlap[r$test == "forest_status"] = o[2]
    
  } else{
    r$overlap[r$test == "forest_status"] = o[1]
  }
  
  
  #### Save results and data used
  temp = list("data" = d, "results"= r)
  
  res[[i]] = temp #storing both original data and results per species
  names(res)[i]= s #save name per species. 
  
  
  
} ## Takes about 3-4 hours to run completely 
rm(c,a,s,i,d,t1,t2,t3,t4,t5,t6,t7,
   t8,t9,t10,r1,r2,r3,r4,r5, o, temp, r)

## Inspect

head(res$Manis_javanica$data)



#### Save results per species ######


for(i in 1:length(res)){
  
  l = res[[i]]
  s = names(res)[i]
  
  d = l$data
  r = l$results
  
  ### save the data
  path = paste0("Loop_Output/Data/", s, "_Activity_Pattern_Analysis_Data.csv", sep = "" )
  
  write.csv(d, path, row.names = FALSE)
  
  
  ### Save the results
  path = paste0("Loop_Output/Results/", s, "_Activity_Analysis_Results.csv", sep = "")
  
  write.csv(r, path, row.names = FALSE)
  
  
}
rm(i, l,s,d,r, path)




###### Create figures for each comparison ######

plots = list() #store all plots here

for(i in 1:length(res)){
  
  l = res[[i]]
  s = names(res)[i]
  
  d = l$data
  
  p = list() #store all plots per species here
  
  ### HFP plot
  e1 = densityPlot(d$Time.rad[d$HFP_status == "No_Humans"])
  e1$HFP_status = "No_Humans"
  e2 = densityPlot(d$Time.rad[d$HFP_status == "Humans"])
  e2$HFP_status = "Humans"
  e = rbind(e1,e2)
  
  a=
    ggplot(e, aes(x = x, y = y, fill = HFP_status))+
    geom_line(size = 1)+
    geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
    geom_rug(data = d[d$HFP_status == "No_Humans",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "cornflowerblue", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    geom_rug(data = d[d$HFP_status == "Humans",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "coral", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("Midnight", "Sunrise", "Noon", "Sunset", "Midnight"), limits = c(0, 24))+
    coord_cartesian(ylim = c(0, max(e$y)))+ ## Need to limit y axis to positive values due to geom_rug's jitter
    annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    annotate("rect", xmin = 18, xmax = 24, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    labs(x = "\nTime of day", y = "Activity\n")+
    scale_fill_manual(name = "HFP Status:", values = c("cornflowerblue", "coral"))+ # Make sure to change the name to your coviaraite!
    theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
          legend.position = "top", legend.key = element_rect(fill = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  p[[1]] = a
  names(p)[1]= "HFP_status"
  
  
  
  ### FLII plot
  e1 = densityPlot(d$Time.rad[d$FLII_status == "Low_Integrity"])
  e1$FLII_status = "Low_Integrity"
  e2 = densityPlot(d$Time.rad[d$FLII_status == "High_Integrity"])
  e2$FLII_status = "High_Integrity"
  e = rbind(e1,e2)
  
  a=
    ggplot(e, aes(x = x, y = y, fill = FLII_status))+
    geom_line(size = 1)+
    geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
    geom_rug(data = d[d$FLII_status == "Low_Integrity",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "cornflowerblue", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    geom_rug(data = d[d$FLII_status == "High_Integrity",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "coral", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("Midnight", "Sunrise", "Noon", "Sunset", "Midnight"), limits = c(0, 24))+
    coord_cartesian(ylim = c(0, max(e$y)))+ ## Need to limit y axis to positive values due to geom_rug's jitter
    annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    annotate("rect", xmin = 18, xmax = 24, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    labs(x = "\nTime of day", y = "Activity\n")+
    scale_fill_manual(name = "FLII Status:", values = c("cornflowerblue", "coral"))+ # Make sure to change the name to your coviaraite!
    theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
          legend.position = "top", legend.key = element_rect(fill = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  p[[2]] = a
  names(p)[2]= "FLII_status"
  
  
  
  
  #### Edge plot
  e1 = densityPlot(d$Time.rad[d$edge_status == "Edge"])
  e1$edge_status = "Edge"
  e2 = densityPlot(d$Time.rad[d$edge_status == "Interior"])
  e2$edge_status = "Interior"
  e = rbind(e1,e2)
  
  a=
    ggplot(e, aes(x = x, y = y, fill = edge_status))+
    geom_line(size = 1)+
    geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
    geom_rug(data = d[d$edge_status == "Edge",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "cornflowerblue", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    geom_rug(data = d[d$edge_status == "Interior",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "coral", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
    scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("Midnight", "Sunrise", "Noon", "Sunset", "Midnight"), limits = c(0, 24))+
    coord_cartesian(ylim = c(0, max(e$y)))+ ## Need to limit y axis to positive values due to geom_rug's jitter
    annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    annotate("rect", xmin = 18, xmax = 24, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
    labs(x = "\nTime of day", y = "Activity\n")+
    scale_fill_manual(name = "Edge Status:", values = c("cornflowerblue", "coral"))+ # Make sure to change the name to your coviaraite!
    theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
          legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
          legend.position = "top", legend.key = element_rect(fill = NA),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  p[[3]] = a
  names(p)[3]= "edge_status"
  
  
  
  if(length(d$Time.rad[d$degraded_status == "Degraded"])> 1){ #bypass lack of dhole observations in degraded
    
    #### Degraded plots
    e1 = densityPlot(d$Time.rad[d$degraded_status == "Non_Degraded"])
    e1$degraded_status = "Non_Degraded"
    e2 = densityPlot(d$Time.rad[d$degraded_status == "Degraded"])
    e2$degraded_status = "Degraded"
    e = rbind(e1,e2)
    
    a=
      ggplot(e, aes(x = x, y = y, fill = degraded_status))+
      geom_line(size = 1)+
      geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
      geom_rug(data = d[d$degraded_status == "Non_Degraded",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "cornflowerblue", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
      geom_rug(data = d[d$degraded_status == "Degraded",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "coral", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("Midnight", "Sunrise", "Noon", "Sunset", "Midnight"), limits = c(0, 24))+
      coord_cartesian(ylim = c(0, max(e$y)))+ ## Need to limit y axis to positive values due to geom_rug's jitter
      annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
      annotate("rect", xmin = 18, xmax = 24, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
      labs(x = "\nTime of day", y = "Activity\n")+
      scale_fill_manual(name = "Degraded Status:", values = c("cornflowerblue", "coral"))+ # Make sure to change the name to your coviaraite!
      theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
            legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
            legend.position = "top", legend.key = element_rect(fill = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    p[[4]] = a
    names(p)[4]= "degraded_status"
  }
  
  
  
  
  if(length(d$Time.rad[d$forest_status == "Non_Forest"]) > 1){ #bypass lack of mongoose observations in Non-forest
    
    #### Forest plots
    e1 = densityPlot(d$Time.rad[d$forest_status == "Forest"])
    e1$forest_status = "Forest"
    e2 = densityPlot(d$Time.rad[d$forest_status == "Non_Forest"])
    e2$forest_status = "Non_Forest"
    e = rbind(e1,e2)
    
    a=
      ggplot(e, aes(x = x, y = y, fill = forest_status))+
      geom_line(size = 1)+
      geom_ribbon(aes(x = x, ymax = y), ymin = 0, alpha = 0.3) +
      geom_rug(data = d[d$forest_status == "Forest",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "cornflowerblue", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
      geom_rug(data = d[d$forest_status == "Non_Forest",], aes(x = hour(Photo.Time), y=0), position = "jitter", sides = "b", color = "coral", inherit.aes = FALSE)+ # this function uses data NOT in the ggplot function!
      scale_x_continuous(breaks = c(0, 6, 12, 18, 24), labels = c("Midnight", "Sunrise", "Noon", "Sunset", "Midnight"), limits = c(0, 24))+
      coord_cartesian(ylim = c(0, max(e$y)))+ ## Need to limit y axis to positive values due to geom_rug's jitter
      annotate("rect", xmin = 0, xmax = 6, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
      annotate("rect", xmin = 18, xmax = 24, ymin = 0, ymax= max(e$y), color = "grey", alpha = .3)+ 
      labs(x = "\nTime of day", y = "Activity\n")+
      scale_fill_manual(name = "Forest Status:", values = c("cornflowerblue", "coral"))+ # Make sure to change the name to your coviaraite!
      theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14),
            legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14),
            legend.position = "top", legend.key = element_rect(fill = NA),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    p[[5]] = a
    names(p)[5]= "forest_status"
    
    
  }
  
  
  ### Save all plots and name per species
  plots[[i]]= p
  names(plots)[i]= s
  
  
  
  
}
rm(l,s,d,e1,e2,e,i, p, a)

## inspect
plots$Argus_sp.$HFP_status
plots$Tragulus_sp.$HFP_status
plots$Manis_javanica$degraded_status
plots$Catopuma_temminckii$edge_status
## looks good! 



### Save plots-

for(i in 1:length(plots)){
  
  a = plots[[i]]
  s = names(plots)[i]
  
  for(l in 1:length(a)){
    
    v = names(a)[l]
    p = a[[l]]
    
    path = paste0("Loop_Output/Figures/", s, "_", v, ".png", sep = "")
    
    ggsave(path, p, height = 5, width = 8, units = "in")
    
  }
} 
rm(a,s,v,p,path,l,i)
## And thats it! Nothing else to do I beleive. 

