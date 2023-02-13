# Activity-pattern-shifts-and-GLMM-analysis

ACTIVITY PATTERN SHIFTS

This analysis is aimed to investigate shifts in species' nocturnality impacted by human activities.

I used time-stamped detections from our new camera trapping sessions to investigate variability in the species's diel activity within the study areas. 
I tested if disturbances (Human Footprint Index and distance from forest edges) affect pangolin behaviour by testing for significant differences in diel activity patterns. 
Specifically, I split our time-stamped detections based on the median value of the disturbance variable and ran a bootstrap procedure to simulate 1,000 distributions of activity pattern data to conduct a Wald test using the function compareAct() in the R package 'activity' (Rowcliffe et al. 2014). 
When significant differences in activity patterns were detected, the coefficient of overlapping was calculated from the R package 'overlap' (Ridout & Linkie 2009). 
I also used the Schmid and Schmidt (2006) Dhat estimator to compute the coefficient of overlap for each type of disturbance variable.

GLMM ANALYSIS

This analysis is aimed to investigate landscape-scale habitat associations of mammal species across SoutEast Asia. 

I assessed landscape-scale associations between the capture rates of published and new camera trapping sessions and landscape-level anthropological and environmental variables using generalized linear mixed models (GLMMs). 
I used a zero-inflated Poisson (ZIP) distribution as we treated detections as count data and included sampling effort per camera trapping session (measured in trap nights) as a continuous fixed effect and landscape as a random effect. 
The sampling unit was a camera trapping session (one set of camera deployments from a single landscape), and all variables were calculated for the entire session. 
I choose to run GLMMs on the raw count data as an improvement on using linear mixed models with on relative abundance index (RAI, usually the number of independent captures per 100 trap nights) following Ash et al. (2020) while acknowledging that either approach does not account for variation in detection probability and thus do not linearly reflect true abundance (Sollmann et al. 2013). 
Therefore, in this analysis, I am implicitly assuming that detection probability does not vary between camera trapping sessions and acknowledge that this approach may introduce unexplained variation in captures owing to slight differences in equipment and deployment methodology between sessions. 
These sources of measurement error may reduce the chances of detecting significant 'true' relationships.

I used AICc model selection to test for relationships with three biophysical descriptors (latitude, annual precipitation, average elevation) and seven indicators of habitat degradation and anthropogenic impact: forest patch size, percentage of forest cover, human population density, nightlight intensity, and the Forest Landscape Integrity Index ('forest integrity' hereafter) and the Human Footprint Index (HFI). 
These variables describe the area within a 20 km radius around the centroid of each landscape (1256 km2) to account for the large areas covered by some camera trapping grids and consider the spatial scale relevant to persistence between landscapes. 
I also tested multivariate additive models after filtering highly correlated variables (| r | > 0.6) (Burnhan & Anderson 2002). 











