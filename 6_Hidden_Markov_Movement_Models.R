  #'  ======================================================
  #'  Code from Bassing et al. 2024 "Predator-prey space-use  
  #'  and landscape features influence animal movement  
  #'  behaviors in a large-mammal community". Ecology.
  #'  
  #'  Hidden Markov Movement Models 
  #'  Washington Predator-Prey Project
  #'  ======================================================
  #'  Script to run hidden Markov movement models for deer, elk, cougars, & wolves
  #'  for summer & winter, July 2018 - March 2021 in northeastern Washington. 
  #'  Data were collected & generously provided by WPPP collaborators including
  #'  M. Devivo, B. Kertson, T.Ganz, T.Roussin, L.Satterfield, & others. Code 
  #'  adapted from momentuHMM GitHub, J.Merkle, L.Satterfield, and R.Emmet.
  #'  
  #'  NOTE: First section of script (lines 18 - 131) provided for transparency and 
  #'  reproducibility but it will not run with data provided on Dryad (animal 
  #'  relocation data are sensitive; contact Director of the Science Division 
  #'  with the Washington Dept. of Fish and Wildlife at (360) 902-2515 for raw data). 
  #'  Data available on Dryad repository associated with this publication can be
  #'  loaded so script will run starting at line 133.
  #'  DOI: 10.5061/dryad.kh1893292
  #'  ======================================================
  
  #'  Clear memory
  rm(list=ls())

  #'  Load libraries
  library(momentuHMM)
  library(rgdal)
  library(ggplot2)
  library(tidyverse)

  #'  Load crwOut & covaraite data
  load("./Outputs/Telemetry_crwOut/crwOut_ALL.RData")    
  load("./Outputs/Telemetry_covs/spp_telem_covs.RData")
  
  #'  Merge datasets and create momentuHMMData object
  #'  Data merged and scaled by study area separately b/c different species collared
  #'  in each study area- can't test effect of study area-specific species
  #'  across both study areas.
  #'  OKANOGAN data sets
  spp_dataPrep_OK <- function(crwOut, telem_covs){
    #'  Merge crawlOut data with extracted covariate data
    crwlMerge <- crawlMerge(crwOut, telem_covs, Time.name = "time")
    #'  Make categorical variables factors
    crwlMerge$crwPredict$StudyArea <- as.factor(crwlMerge$crwPredict$StudyArea)
    crwlMerge$crwPredict$Sex <- as.factor(crwlMerge$crwPredict$Sex)
    crwlMerge$crwPredict$Season <- as.factor(crwlMerge$crwPredict$Season)
    crwlMerge$crwPredict$SnowCover <- as.factor(crwlMerge$crwPredict$SnowCover)
    crwlMerge$crwPredict$daytime <- as.factor(crwlMerge$crwPredict$daytime)
    #'  Standardize continuous variables
    crwlMerge$crwPredict$Dist2Road <- scale(crwlMerge$crwPredict$Dist2Road)
    crwlMerge$crwPredict$PercOpen <- scale(crwlMerge$crwPredict$PercOpen)
    crwlMerge$crwPredict$TRI <- scale(crwlMerge$crwPredict$TRI)
    crwlMerge$crwPredict$MD_RSF <- scale(crwlMerge$crwPredict$MD_RSF)
    crwlMerge$crwPredict$COUG_RSF <- scale(crwlMerge$crwPredict$COUG_RSF)
    crwlMerge$crwPredict$WOLF_RSF <- scale(crwlMerge$crwPredict$WOLF_RSF)
    crwlMerge$crwPredict$hour <- as.integer(crwlMerge$crwPredict$hour)
    crwlMerge$crwPredict$hour_fix <- as.integer(crwlMerge$crwPredict$hour_fix)
    crwlMerge$crwPredict$hour3 <- as.integer(crwlMerge$crwPredict$hour3)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen", 
                                                    "SnowCover", "TRI", "MD_RSF", 
                                                    "COUG_RSF", "WOLF_RSF",  
                                                    "hour", "hour_fix",
                                                    "hour3", "daytime", "Sex", 
                                                    "StudyArea", "Season"))  
    return(Data)
  }
  #'  Run season & species-specific data from the Okanogan through prep function
  #'  Warnings are due to missing data for interpolated locations. prepData 
  #'  command automatically fills in values with closest following value.
  mdData_smr <- spp_dataPrep_OK(crwOut_ALL[[1]], spp_telem_covs[[1]])      
  mdData_wtr <- spp_dataPrep_OK(crwOut_ALL[[2]], spp_telem_covs[[2]])      
  cougData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[7]], spp_telem_covs[[7]]) 
  cougData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[8]], spp_telem_covs[[8]])
  wolfData_smr_OK <- spp_dataPrep_OK(crwOut_ALL[[11]], spp_telem_covs[[11]])
  wolfData_wtr_OK <- spp_dataPrep_OK(crwOut_ALL[[12]], spp_telem_covs[[12]])
  
  #'  NORTHEAST data sets
  spp_dataPrep_NE <- function(crwOut, telem_covs){
    #'  Merge crawlOut data with extracted covariate data
    crwlMerge <- crawlMerge(crwOut, telem_covs, Time.name = "time")
    #'  Make categorical variables factors
    crwlMerge$crwPredict$StudyArea <- as.factor(crwlMerge$crwPredict$StudyArea)
    crwlMerge$crwPredict$Sex <- as.factor(crwlMerge$crwPredict$Sex)
    crwlMerge$crwPredict$Season <- as.factor(crwlMerge$crwPredict$Season)
    crwlMerge$crwPredict$SnowCover <- as.factor(crwlMerge$crwPredict$SnowCover)
    crwlMerge$crwPredict$daytime <- as.factor(crwlMerge$crwPredict$daytime)
    #'  Standardize continuous variables
    crwlMerge$crwPredict$Dist2Road <- scale(crwlMerge$crwPredict$Dist2Road)
    crwlMerge$crwPredict$PercOpen <- scale(crwlMerge$crwPredict$PercOpen)
    crwlMerge$crwPredict$TRI <- scale(crwlMerge$crwPredict$TRI)
    crwlMerge$crwPredict$ELK_RSF <- scale(crwlMerge$crwPredict$ELK_RSF)
    crwlMerge$crwPredict$WTD_RSF <- scale(crwlMerge$crwPredict$WTD_RSF)
    crwlMerge$crwPredict$COUG_RSF <- scale(crwlMerge$crwPredict$COUG_RSF)
    crwlMerge$crwPredict$WOLF_RSF <- scale(crwlMerge$crwPredict$WOLF_RSF)
    crwlMerge$crwPredict$hour <- as.integer(crwlMerge$crwPredict$hour)
    crwlMerge$crwPredict$hour_fix <- as.integer(crwlMerge$crwPredict$hour_fix)
    crwlMerge$crwPredict$hour3 <- as.integer(crwlMerge$crwPredict$hour3)
    #'  Prep crwlMerge data for fitHMM function
    Data <- prepData(data = crwlMerge, covNames = c("Dist2Road", "PercOpen",  
                                                    "SnowCover", "TRI", "ELK_RSF", 
                                                    "WTD_RSF", "COUG_RSF", "WOLF_RSF", 
                                                    "hour", "hour_fix", "hour3", "daytime",
                                                    "Sex", "StudyArea", "Season")) 
    return(Data)
  }
  #'  Run season & species-specific data from the Northeast through prep function
  #'  Warnings are due to missing data for interpolated locations. prepData 
  #'  command automatically fills in values with closest following value.
  elkData_smr <- spp_dataPrep_NE(crwOut_ALL[[3]], spp_telem_covs[[3]])       
  elkData_wtr <- spp_dataPrep_NE(crwOut_ALL[[4]], spp_telem_covs[[4]])
  wtdData_smr <- spp_dataPrep_NE(crwOut_ALL[[5]], spp_telem_covs[[5]])
  wtdData_wtr <- spp_dataPrep_NE(crwOut_ALL[[6]], spp_telem_covs[[6]])
  cougData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[9]], spp_telem_covs[[9]])
  cougData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[10]], spp_telem_covs[[10]])
  wolfData_smr_NE <- spp_dataPrep_NE(crwOut_ALL[[13]], spp_telem_covs[[13]])
  wolfData_wtr_NE <- spp_dataPrep_NE(crwOut_ALL[[14]], spp_telem_covs[[14]])
  
  #'  Save data prepped for HMMs
  hmm_data <- list(mdData_smr, mdData_wtr, elkData_smr, elkData_wtr, wtdData_smr, 
                   wtdData_wtr, cougData_smr_OK, cougData_wtr_OK, cougData_smr_NE, 
                   cougData_wtr_NE, wolfData_smr_OK, wolfData_wtr_OK, wolfData_smr_NE, 
                   wolfData_wtr_NE)
  
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_for_pubs.RData")
  
  names(hmm_data) <- c("mdData_smr", "mdData_wtr", "elkData_smr", "elkData_wtr", "wtdData_smr", "wtdData_wtr",
                  "cougData_smr_OK", "cougData_wtr_OK", "cougData_smr_NE", "cougData_wtr_NE",
                  "wolfData_smr_OK", "wolfData_wtr_OK", "wolfData_smr_NE", "wolfData_wtr_NE")
  
  #'  Correlation Matrix
  #'  ==================
  #'  Function to create correlation matrix for all continuous covariates at once
  cov_correlation_OK <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "MD_RSF", "COUG_RSF", #"NDVI", 
                      "WOLF_RSF")) 
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (md_smr_corr <- cov_correlation_OK(mdData_smr)) 
  (md_wtr_corr <- cov_correlation_OK(mdData_wtr)) 
  (coug_smr_OK_corr <- cov_correlation_OK(cougData_smr_OK))
  (coug_wtr_OK_corr <- cov_correlation_OK(cougData_wtr_OK)) 
  (wolf_smr_OK_corr <- cov_correlation_OK(wolfData_smr_OK))
  (wolf_wtr_OK_corr <- cov_correlation_OK(wolfData_wtr_OK)) 
  
  cov_correlation_NE <- function(dat) {
    covs <- dat %>%
      dplyr::select(c("Dist2Road", "PercOpen", "TRI", "ELK_RSF", "WTD_RSF",  
                      "COUG_RSF", "WOLF_RSF"))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  #'  Generate correlation matrix for each species and season
  (elk_smr_corr <- cov_correlation_NE(elkData_smr)) 
  (elk_wtr_corr <- cov_correlation_NE(elkData_wtr))
  (wtd_smr_corr <- cov_correlation_NE(wtdData_smr)) 
  (wtd_wtr_corr <- cov_correlation_NE(wtdData_wtr)) 
  (coug_smr_NE_corr <- cov_correlation_NE(cougData_smr_NE)) 
  (coug_wtr_NE_corr <- cov_correlation_NE(cougData_wtr_NE))
  (wolf_smr_NE_corr <- cov_correlation_NE(wolfData_smr_NE))  
  (wolf_wtr_NE_corr <- cov_correlation_NE(wolfData_wtr_NE))
  
  #' #'  Visualize data to inform initial parameter specifications
  mean(mdData_smr$step, na.rm = T); sd(mdData_smr$step, na.rm = T)
  mean(mdData_wtr$step, na.rm = T); sd(mdData_wtr$step, na.rm = T)
  mean(elkData_smr$step, na.rm = T); sd(elkData_smr$step, na.rm = T)
  mean(elkData_wtr$step, na.rm = T); sd(elkData_wtr$step, na.rm = T)
  mean(wtdData_smr$step, na.rm = T); sd(wtdData_smr$step, na.rm = T)
  mean(wtdData_wtr$step, na.rm = T); sd(wtdData_wtr$step, na.rm = T)
  mean(cougData_smr_OK$step, na.rm = T); sd(cougData_smr_OK$step, na.rm = T)
  mean(cougData_wtr_OK$step, na.rm = T); sd(cougData_wtr_OK$step, na.rm = T)
  mean(cougData_smr_NE$step, na.rm = T); sd(cougData_smr_NE$step, na.rm = T)
  mean(cougData_wtr_NE$step, na.rm = T); sd(cougData_wtr_NE$step, na.rm = T)
  mean(wolfData_smr_OK$step, na.rm = T); sd(wolfData_smr_OK$step, na.rm = T)
  mean(wolfData_wtr_OK$step, na.rm = T); sd(wolfData_wtr_OK$step, na.rm = T)
  mean(wolfData_smr_NE$step, na.rm = T); sd(wolfData_smr_NE$step, na.rm = T)
  mean(wolfData_wtr_NE$step, na.rm = T); sd(wolfData_wtr_NE$step, na.rm = T)
  
  #'  Visualize data to identify potential temporal autocorrelation
  #'  lag.max is measured in hours
  acf(mdData_smr$step[!is.na(mdData_smr$step)],lag.max=100)
  acf(mdData_wtr$step[!is.na(mdData_wtr$step)],lag.max=100)
  acf(elkData_smr$step[!is.na(elkData_smr$step)],lag.max=100)
  acf(elkData_wtr$step[!is.na(elkData_wtr$step)],lag.max=100)
  acf(wtdData_smr$step[!is.na(wtdData_smr$step)],lag.max=100)
  acf(wtdData_wtr$step[!is.na(wtdData_wtr$step)],lag.max=100)
  acf(cougData_smr_OK$step[!is.na(cougData_smr_OK$step)],lag.max=100)
  acf(cougData_wtr_OK$step[!is.na(cougData_wtr_OK$step)],lag.max=100)
  acf(cougData_smr_NE$step[!is.na(cougData_smr_NE$step)],lag.max=100)
  acf(cougData_wtr_NE$step[!is.na(cougData_wtr_NE$step)],lag.max=100)
  acf(wolfData_smr_OK$step[!is.na(wolfData_smr_OK$step)],lag.max=100)
  acf(wolfData_wtr_OK$step[!is.na(wolfData_wtr_OK$step)],lag.max=100)
  acf(wolfData_smr_NE$step[!is.na(wolfData_smr_NE$step)],lag.max=100)
  acf(wolfData_wtr_NE$step[!is.na(wolfData_wtr_NE$step)],lag.max=100)
 
  #'  What's up with the ACF? Plot step lengths against hour to look for patterns
  mdData_smr <- hmm_data[[1]] 
  mdData_wtr <- hmm_data[[2]] 
  elkData_smr <- hmm_data[[3]]
  elkData_wtr <- hmm_data[[4]]
  wtdData_smr <- hmm_data[[5]] 
  wtdData_wtr <- hmm_data[[6]] 
  cougData_smr_OK <- hmm_data[[7]] 
  cougData_wtr_OK <- hmm_data[[8]] 
  cougData_smr_NE <- hmm_data[[9]]
  cougData_wtr_NE <- hmm_data[[10]]
  wolfData_smr_OK <- hmm_data[[11]]
  wolfData_wtr_OK <- hmm_data[[12]]
  wolfData_smr_NE <- hmm_data[[13]] 
  wolfData_wtr_NE <- hmm_data[[14]] 
  
  # write.csv(mdData_smr, "./Outputs/Telemetry_crwOut/mdData_smr_crwOut.csv")
  # write.csv(mdData_wtr, "./Outputs/Telemetry_crwOut/mdData_wtr_crwOut.csv")
  # write.csv(wtdData_smr, "./Outputs/Telemetry_crwOut/wtdData_smr_crwOut.csv")
  # write.csv(wtdData_wtr, "./Outputs/Telemetry_crwOut/wtdData_wtr_crwOut.csv")
  # write.csv(cougData_smr_OK, "./Outputs/Telemetry_crwOut/cougData_smr_OK_crwOut.csv")
  # write.csv(cougData_wtr_OK, "./Outputs/Telemetry_crwOut/cougData_wtr_OK_crwOut.csv")
  # write.csv(wolfData_smr_NE, "./Outputs/Telemetry_crwOut/wolfData_smr_NE_crwOut.csv")
  # write.csv(wolfData_wtr_NE, "./Outputs/Telemetry_crwOut/wolfData_wtr_NE_crwOut.csv")
  
  
  ####  Initial model set up  ####
  #'  ============================
  #'  Define initial parameters associated with each distribution & each state
  #'  Species-specific parameters based on viewing plotted data and mean step lengths
  #'  Providing value close to mean step length as "exploratory" mean & SD
  Par0_m1_md <- list(step = c(100, 250, 100, 250, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_elk <- list(step = c(100, 450, 100, 450, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_wtd <- list(step = c(100, 260, 100, 260, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_coug <- list(step = c(100, 650, 100, 650, 0.01, 0.005), angle = c(0.1, 0.5))  
  Par0_m1_wolf <- list(step = c(100, 1600, 100, 1600), angle = c(0.1, 0.5))  
  #'  Step arguments: report 2 means then the 2 SD for the two different states
  #'  Gamma distribution: mean & standard deviation of step lengths for each state
  #'  Michelot & Langrock 2019 recommend using same value for mean and SD per state
  #'  Wrapped Cauchy distribution: concentration of turning angles for each state
  #'  Include zero-mass parameters when there are 0s in the data w/gamma, Weibull, 
  #'  etc. distributions, e.g., zeromass0 <- c(0.1,0.05) # step zero-mass
  #'  Applies to mule deer, elk, white-tailed deer, and cougars
  
  #'  Label states
  stateNames <- c("encamped", "exploratory")
  
  ####  Models describing State-Dependent Distributions  ####
  #' Distributions for observation processes
  #' Step length: gamma or Weibull; Turning angle: von Mises or wrapped Cauchy
  #' State dwell time: geometric distribution
  #' Weibull = "weibull"; von Mises = "vm"
  dists_wc <- list(step = "gamma", angle = "wrpcauchy")  
  dists_vm <- list(step = "gamma", angle = "vm")
  
  #'  Define formula(s) to be applied to state-dependent distributions
  #'  Covariates that help describe movement patterns of a given state
  #'  Add zeromass = formula for species that need zeromass parameters above
  DM_formula_null <- ~1

  #'  Create pseudo-design matrices for state-dependent distributions
  DM_null <- list(step = list(mean = ~1, sd = ~1), angle = list(concentration = ~1))
  DM_null_ZeroMass <- list(step = list(mean = ~1, sd = ~1, zeromass = ~1), angle = list(concentration = ~1)) # includes zeromass parameters
  DM_time <- list(step = list(mean = ~daytime + cosinor(hour_fix, period = 12), sd = ~daytime + cosinor(hour_fix, period = 12)), angle = list(concentration = ~1))
  DM_Zerotime <- list(step = list(mean = ~daytime + cosinor(hour_fix, period = 12), sd = ~daytime + cosinor(hour_fix, period = 12), zeromass = ~1), angle = list(concentration = ~1)) # includes zeromass parameters
  
  ####  Models describing Transition Probabilities  ####
  #'  Define formula(s) to be applied to transition probabilities
  #'  Covariates affecting probability of transitioning from one state to another
  #'  and associated with behavioral states
  trans_formula_null <- ~1
  #'  For prey species
  trans_formula_smr_all <- ~TRI + PercOpen + Dist2Road + COUG_RSF + WOLF_RSF 
  trans_formula_wtr_all <- ~TRI + PercOpen + Dist2Road + SnowCover + COUG_RSF + WOLF_RSF  
  #'  For predator species
  trans_formula_smr_OK <- ~TRI + PercOpen + Dist2Road + MD_RSF 
  trans_formula_wtr_OK <- ~TRI + PercOpen + Dist2Road + SnowCover + MD_RSF 
  trans_formula_wtr_OK_noMD <- ~TRI + PercOpen + Dist2Road + SnowCover
  trans_formula_wtr_OK_noTRI <- ~PercOpen + Dist2Road + SnowCover + MD_RSF
  trans_formula_smr_NE <- ~TRI + PercOpen + Dist2Road + ELK_RSF + WTD_RSF 
  trans_formula_wtr_NE <- ~TRI + PercOpen + Dist2Road + SnowCover + ELK_RSF + WTD_RSF 
  
  
  ####  It's H[a]MM[er] Time!  ####
  #'  =============================
  #'  Keep in mind I can fit covariates on the state transition probabilities, 
  #'  meaning the variables that influence whether an animal will transition from
  #'  one state to the other, or on the state-dependent observation distributions,
  #'  meaning variables that influence step length and/or turning angle for each
  #'  of the states. 

  #'  Use retryFits argument to specify the number of attempts to minimize the 
  #'  negative log-likelihood based on random perturbations of the parameter 
  #'  estimates at the current minimum- helps ensure convergence
  
  #'  Function to run data through null and global HMM for each species
  HMM_fit <- function(Data, dists, Par0_m1, dm, tformula, fits) { 
    
    #' Fit basic model with no covariates
    m1 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m1,
                 estAngleMean = list(angle = FALSE), stateNames = stateNames,
                 retryFits = fits)
    
    #'  Get new initial parameter values for global model based on nested m1 model
    Par0_m2 <- getPar0(model = m1, DM = dm, formula = tformula)   
    
    #'  Fit model with sex covariate on transition probability
    m2 <- fitHMM(data = Data, nbStates = 2, dist = dists, Par0 = Par0_m2$Par,
                 stateNames = stateNames, DM = dm, beta0 = Par0_m2$beta, formula = tformula) 
    
    #'  What proportion of the locations fall within each state?
    states <- viterbi(m2)
    print(table(states)/nrow(Data))
    
    #'  Model selection with AIC
    print(AIC(m1,m2))
    
    #'  Model summary and covariate effects
    print(m2)
    
    global_est <- CIbeta(m2, alpha = 0.95)
    print(global_est[[3]])
    
    return(m2)
    
  }
  ####  MULE DEER HMMS  ####     
  #'  Summer
  md_HMM_smr <- HMM_fit(mdData_smr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(md_HMM_smr, lag.max = 100, ncores = 4)
  pr_md_HMM_smr <- pseudoRes(md_HMM_smr)
  acf(pr_md_HMM_smr$stepRes[is.finite(pr_md_HMM_smr$stepRes)], lag.max = 100)
  plot(md_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Winter
  md_HMM_wtr <- HMM_fit(mdData_wtr, dists_vm, Par0_m1_md, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(md_HMM_wtr, lag.max = 100, ncores = 4)
  pr_md_HMM_wtr <- pseudoRes(md_HMM_wtr)
  acf(pr_md_HMM_wtr$stepRes[is.finite(pr_md_HMM_wtr$stepRes)], lag.max = 100)
  plot(md_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  ####  ELK HMMS  ####
  #'  Summer
  elk_HMM_smr <- HMM_fit(elkData_smr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(elk_HMM_smr, lag.max = 100, ncores = 4)
  pr_elk_HMM_smr <- pseudoRes(elk_HMM_smr)
  acf(pr_elk_HMM_smr$stepRes[!is.na(pr_elk_HMM_smr$stepRes)],lag.max = 100)
  plot(elk_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE) 
  
  #'  Winter
  elk_HMM_wtr <- HMM_fit(elkData_wtr, dists_vm, Par0_m1_elk, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(elk_HMM_wtr, lag.max = 100, ncores = 4)
  pr_elk_HMM_wtr <- pseudoRes(elk_HMM_wtr)
  acf(pr_elk_HMM_wtr$stepRes[!is.na(pr_elk_HMM_wtr$stepRes)],lag.max = 100)
  plot(elk_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  ####  WHITE-TAILED DEER HMMS  ####
  #'  Summer
  wtd_HMM_smr <- HMM_fit(wtdData_smr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_smr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(wtd_HMM_smr, lag.max = 100, ncores = 4)
  pr_wtd_HMM_smr <- pseudoRes(wtd_HMM_smr)
  acf(pr_wtd_HMM_smr$stepRes[!is.na(pr_wtd_HMM_smr$stepRes)],lag.max = 100)
  plot(wtd_HMM_smr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Winter
  wtd_HMM_wtr <- HMM_fit(wtdData_wtr, dists_vm, Par0_m1_wtd, DM_Zerotime, trans_formula_wtr_all, fits = 1)
  #'  Inspect residuals and plot
  plotPR(wtd_HMM_wtr, lag.max = 100, ncores = 4)
  pr_wtd_HMM_wtr <- pseudoRes(wtd_HMM_wtr)
  acf(pr_wtd_HMM_wtr$stepRes[!is.na(pr_wtd_HMM_wtr$stepRes)],lag.max = 100)
  plot(wtd_HMM_wtr, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  ####  COUGAR HMMS  ####       
  #'  Okanogan Summer
  coug_HMM_smr_OK <- HMM_fit(cougData_smr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_smr_OK, lag.max = 100, ncores = 4)
  pr_coug_HMM_smr_OK <- pseudoRes(coug_HMM_smr_OK)
  acf(pr_coug_HMM_smr_OK$stepRes[!is.na(pr_coug_HMM_smr_OK$stepRes)],lag.max = 100)
  plot(coug_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Okanogan Winter
  coug_HMM_wtr_OK <- HMM_fit(cougData_wtr_OK, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_OK_noMD, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_wtr_OK, lag.max = 100, ncores = 4)
  pr_coug_HMM_wtr_OK <- pseudoRes(coug_HMM_wtr_OK)
  acf(pr_coug_HMM_wtr_OK$stepRes[!is.na(pr_coug_HMM_wtr_OK$stepRes)],lag.max = 100)
  plot(coug_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Northeast Summer
  coug_HMM_smr_NE <- HMM_fit(cougData_smr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_smr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_smr_NE, lag.max = 100, ncores = 4)
  pr_coug_HMM_smr_NE <- pseudoRes(coug_HMM_smr_NE)
  acf(pr_coug_HMM_smr_NE$stepRes[!is.na(pr_coug_HMM_smr_NE$stepRes)],lag.max = 100)
  plot(coug_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Northeast Winter
  coug_HMM_wtr_NE <- HMM_fit(cougData_wtr_NE, dists_vm, Par0_m1_coug, DM_Zerotime, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(coug_HMM_wtr_NE, lag.max = 100, ncores = 4)
  pr_coug_HMM_wtr_NE <- pseudoRes(coug_HMM_wtr_NE)
  acf(pr_coug_HMM_wtr_NE$stepRes[!is.na(pr_coug_HMM_wtr_NE$stepRes)],lag.max = 100)
  plot(coug_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  ####  WOLF HMMS  ####
  wolf_HMM_smr_OK <- HMM_fit(wolfData_smr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_smr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_smr_OK, lag.max = 100, ncores = 4)
  pr_wolf_HMM_smr_OK <- pseudoRes(wolf_HMM_smr_OK)
  acf(pr_wolf_HMM_smr_OK$stepRes[!is.na(pr_wolf_HMM_smr_OK$stepRes)],lag.max = 100)
  plot(wolf_HMM_smr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Okanogan Winter
  wolf_HMM_wtr_OK <- HMM_fit(wolfData_wtr_OK, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtr_OK, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_wtr_OK, lag.max = NULL, ncores = 4)
  pr_wolf_HMM_wtr_OK <- pseudoRes(wolf_HMM_wtr_OK)
  acf(pr_wolf_HMM_wtr_OK$stepRes[!is.na(pr_wolf_HMM_wtr_OK$stepRes)],lag.max = 100)
  plot(wolf_HMM_wtr_OK, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Northeast Summer
  wolf_HMM_smr_NE <- HMM_fit(wolfData_smr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_smr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_smr_NE, lag.max = NULL, ncores = 4)
  pr_wolf_HMM_smr_NE <- pseudoRes(wolf_HMM_smr_NE)
  acf(pr_wolf_HMM_smr_NE$stepRes[!is.na(pr_wolf_HMM_smr_NE$stepRes)],lag.max = 100)
  plot(wolf_HMM_smr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  #'  Northeast Winter
  wolf_HMM_wtr_NE <- HMM_fit(wolfData_wtr_NE, dists_vm, Par0_m1_wolf, DM_time, trans_formula_wtr_NE, fits = 1)
  #'  QQplot of residuals
  plotPR(wolf_HMM_wtr_NE, lag.max = NULL, ncores = 4)
  pr_wolf_HMM_wtr_NE <- pseudoRes(wolf_HMM_wtr_NE)
  acf(pr_wolf_HMM_wtr_NE$stepRes[!is.na(pr_wolf_HMM_wtr_NE$stepRes)],lag.max = 100)
  plot(wolf_HMM_wtr_NE, ask = TRUE, animals = 1, breaks = 20, plotCI = TRUE)
  
  
  #'  Save model results
  spp_HMM_output <- list(md_HMM_smr, md_HMM_wtr, elk_HMM_smr, elk_HMM_wtr, wtd_HMM_smr, 
                         wtd_HMM_wtr, coug_HMM_smr_OK, coug_HMM_wtr_OK, coug_HMM_smr_NE, 
                         coug_HMM_wtr_NE, wolf_HMM_smr_OK, wolf_HMM_wtr_OK, 
                         wolf_HMM_smr_NE, wolf_HMM_wtr_NE)
  
  save(spp_HMM_output, file = "./Outputs/HMM_output/spp_HMM_output.RData")
  

  ####  Summarize Results  ####
  #'  Review model output
  print(spp_HMM_output[[1]]) # md_HMM_smr
  print(spp_HMM_output[[2]]) # md_HMM_wtr
  print(spp_HMM_output[[3]]) # elk_HMM_smr
  print(spp_HMM_output[[4]]) # elk_HMM_wtr
  print(spp_HMM_output[[5]]) # wtd_HMM_smr
  print(spp_HMM_output[[6]]) # wtr_HMM_wtr
  print(spp_HMM_output[[7]]) # coug_HMM_smr_OK
  print(spp_HMM_output[[8]]) # coug_HMM_wtr_OK
  print(spp_HMM_output[[9]]) # coug_HMM_smr_NE
  print(spp_HMM_output[[10]]) # coug_HMM_wtr_NE
  print(spp_HMM_output[[11]]) # wolf_HMM_smr_OK
  print(spp_HMM_output[[12]]) # wolf_HMM_wtr_OK
  print(spp_HMM_output[[13]]) # wolf_HMM_smr_NE
  print(spp_HMM_output[[14]]) # wolf_HMM_wtr_NE
  

  ####  State-Dependent Distributions  ####
  #'  Function to report state-dependent distribution parameters, including zero-mass parameters
  step_turn_parms_zmass <- function(mod, spp, season, area){ 
    #'  Pull out turning angle parameters
    step_out <- as.data.frame(mod$mle[[1]])
    step_out$Species <- spp
    step_out$Season <- season
    step_out$StudyArea <- area    
    colnames(step_out) <- c("State1 Intercept_mu", "State1 Daylight_mu", "State1 Cos_mu", "State1 Sin_mu", 
                            "State2 Intercept_mu", "State2 Daylight_mu", "State2 Cos_mu", "State2 Sin_mu",
                            "State1 Intercept_sd", "State1 Daylight_sd", "State1 Cos_sd", "State1 Sin_sd", 
                            "State2 Intercept_sd", "State2 Daylight_sd", "State2 Cos_sd", "State2 Sin_sd",
                            "State1 Intercept_zmass", "State2 Intercept_zmass",  
                            "Species", "Season", "StudyArea")
    #'  Wrangle parameters into an interpret-able table
    step_table <- step_out %>%
      pivot_longer(!c(Species, Season, StudyArea), names_to = "Parameter", values_to = "Estimate") %>%
      separate(Parameter, c("State", "Parameter"), sep = " ") %>%
      pivot_wider(names_from = "State", values_from = "Estimate") %>%
      separate(Parameter, c("Coefficient", "Parameter"), sep = "_") %>%
      pivot_wider(names_from = "Parameter", values_from = c("State1", "State2"))
    #'  Create separate tables for state 1 & 2 parameters
    state1 <- step_table[,1:7]
    state1$State <- "Encamped"
    colnames(state1) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "Zeromass", "State")
    state2 <- step_table[,c(1:4,8:10)]
    state2$State <- "Exploratory"
    colnames(state2) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "Zeromass", "State")
    #'  Merge into one single table of step length parameters
    step_out_tbl <- rbind(state1, state2) %>%
      relocate(State, .before = "Coefficient")

    #'  Turning angles parameters
    turn_out <- as.data.frame(mod$mle[[2]])
    turn_out$Species <- spp
    turn_out$Season <- season
    turn_out$StudyArea <- area
    turn_out_tbl <- turn_out %>%
      relocate(Species, .before = "encamped") %>%
      relocate(Season, .after = "Species") %>%
      relocate(StudyArea, .after = "Season") %>%
      rownames_to_column(var = "Parameter") %>%
      relocate(Parameter, .after = "StudyArea") %>%
      mutate(Parameter = ifelse(Parameter == "mean", "Mean", "Concentration"))
    colnames(turn_out_tbl) <- c("Species", "Season", "Study Area", "Parameter", "Encamped", "Exploratory")
      
    #'  List parameter tables together
    params_out <- list(step_out_tbl, turn_out_tbl)
    return(params_out)

  }
  #'  Create parameter tables for species that included zero-mass parameters
  md_smr_params <- step_turn_parms_zmass(spp_HMM_output[[1]], spp = "Mule Deer", season = "Summer", area = "Okanogan")
  md_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[2]], spp = "Mule Deer", season = "Winter", area = "Okanogan")
  elk_smr_params <- step_turn_parms_zmass(spp_HMM_output[[3]], spp = "Elk", season = "Summer", area = "Northeast")
  elk_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[4]], spp = "Elk", season = "Winter", area = "Northeast")
  wtd_smr_params <- step_turn_parms_zmass(spp_HMM_output[[5]], spp = "White-tailed Deer", season = "Summer", area = "Northeast")
  wtd_wtr_params <- step_turn_parms_zmass(spp_HMM_output[[6]], spp = "White-tailed Deer", season = "Winter", area = "Northeast")
  coug_smr_params_OK <- step_turn_parms_zmass(spp_HMM_output[[7]], spp = "Cougar", season = "Summer", area = "Okanogan")
  coug_wtr_params_OK <- step_turn_parms_zmass(spp_HMM_output[[8]], spp = "Cougar", season = "Winter", area = "Okanogan")
  coug_smr_params_NE <- step_turn_parms_zmass(spp_HMM_output[[9]], spp = "Cougar", season = "Summer", area = "Northeast")
  coug_wtr_params_NE <- step_turn_parms_zmass(spp_HMM_output[[10]], spp = "Cougar", season = "Winter", area = "Northeast")
  

  #'  Function to report state-dependent distribution parameters, excluding zero-mass parameters
  step_turn_parms <- function(mod, spp, season, area){ 
    #'  Pull out turning angle parameters
    step_out <- as.data.frame(mod$mle[[1]])
    step_out$Species <- spp
    step_out$Season <- season
    step_out$StudyArea <- area
    colnames(step_out) <- c("State1 Intercept_mu", "State1 Daylight_mu", "State1 Cos_mu", "State1 Sin_mu", 
                            "State2 Intercept_mu", "State2 Daylight_mu", "State2 Cos_mu", "State2 Sin_mu",
                            "State1 Intercept_sd", "State1 Daylight_sd", "State1 Cos_sd", "State1 Sin_sd", 
                            "State2 Intercept_sd", "State2 Daylight_sd", "State2 Cos_sd", "State2 Sin_sd",
                            "Species", "Season", "StudyArea")
    #'  Wrangle parameters into an interpret-able table
    step_table <- step_out %>%
      pivot_longer(!c(Species, Season, StudyArea), names_to = "Parameter", values_to = "Estimate") %>%
      separate(Parameter, c("State", "Parameter"), sep = " ") %>%
      pivot_wider(names_from = "State", values_from = "Estimate") %>%
      separate(Parameter, c("Coefficient", "Parameter"), sep = "_") %>%
      pivot_wider(names_from = "Parameter", values_from = c("State1", "State2"))
    #'  Create separate tables for state 1 & 2 parameters
    state1 <- step_table[,1:6]
    state1$State <- "Encamped"
    colnames(state1) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "State")
    state2 <- step_table[,c(1:4,7:8)]
    state2$State <- "Exploratory"
    colnames(state2) <- c("Species", "Season", "Study Area", "Coefficient", "Mean", "SD", "State")
    #'  Merge into one single table of step length parameters
    step_out_tbl <- rbind(state1, state2) %>%
      relocate(State, .before = "Coefficient") %>%
      mutate(zmass = NA)
    colnames(step_out_tbl) <- c("Species", "Season", "Study Area", "State", "Coefficient", 
                                "Mean", "SD", "Zeromass")
    
    #'  Turning angles parameters
    turn_out <- as.data.frame(mod$mle[[2]])
    turn_out$Species <- spp
    turn_out$Season <- season
    turn_out$StudyArea <- area
    turn_out_tbl <- turn_out %>%
      relocate(Species, .before = "encamped") %>%
      relocate(Season, .after = "Species") %>%
      relocate(StudyArea, .after = "Season") %>%
      rownames_to_column(var = "Parameter") %>%
      relocate(Parameter, .after = "StudyArea") %>%
      mutate(Parameter = ifelse(Parameter == "mean", "Mean", "Concentration"))
    colnames(turn_out_tbl) <- c("Species", "Season", "Study Area", "Parameter", "Encamped", "Exploratory")
    
    #'  List parameter tables together
    params_out <- list(step_out_tbl, turn_out_tbl)
    return(params_out)
    
  }
  #'  Create parameter tables for species that don't include zero-mass parameters
  wolf_smr_params_OK <- step_turn_parms(spp_HMM_output[[11]], spp = "Wolf", season = "Summer", area = "Okanogan")
  wolf_wtr_params_OK <- step_turn_parms(spp_HMM_output[[12]], spp = "Wolf", season = "Winter", area = "Okanogan")
  wolf_smr_params_NE <- step_turn_parms(spp_HMM_output[[13]], spp = "Wolf", season = "Summer", area = "Northeast")
  wolf_wtr_params_NE <- step_turn_parms(spp_HMM_output[[14]], spp = "Wolf", season = "Winter", area = "Northeast")
  
  #'  Make single giant table of all step length parameters
  all_steps <- bind_rows(md_smr_params[[1]], md_wtr_params[[1]], elk_smr_params[[1]], 
                         elk_wtr_params[[1]], wtd_smr_params[[1]], wtd_wtr_params[[1]],
                         coug_smr_params_OK[[1]], coug_wtr_params_OK[[1]], 
                         coug_smr_params_NE[[1]], coug_wtr_params_NE[[1]], 
                         wolf_smr_params_OK[[1]], wolf_wtr_params_OK[[1]], 
                         wolf_smr_params_NE[[1]], wolf_wtr_params_NE[[1]]) %>%
    mutate(Mean = round(Mean, 2),
           SD = round(SD, 2),
           Zeromass = round(Zeromass, 2)) %>%
    arrange(Species)
  
  #'  Make single giant table of all turning angles parameters
  all_turns <- bind_rows(md_smr_params[[2]], md_wtr_params[[2]], elk_smr_params[[2]], 
                         elk_wtr_params[[2]], wtd_smr_params[[2]], wtd_wtr_params[[2]],
                         coug_smr_params_OK[[2]], coug_wtr_params_OK[[2]], 
                         coug_smr_params_NE[[2]], coug_wtr_params_NE[[2]], 
                         wolf_smr_params_OK[[2]], wolf_wtr_params_OK[[2]], 
                         wolf_smr_params_NE[[2]], wolf_wtr_params_NE[[2]]) %>%
    mutate(Encamped = round(Encamped, 2),
           Exploratory = round(Exploratory, 2)) %>%
    arrange(Species)

  
  ####  Transition Probabilities  ####
  #'  Function to report transition probability coefficients in a table
  rounddig <- 2
  hmm_out <- function(mod, spp, season, area) {
    #'  Extract estimates, standard error, and 95% Confidence Intervals for effect
    #'  of each covariate on transition probabilities
    est_out <- CIbeta(mod, alpha = 0.95)
    beta1.2 <- formatC(round(est_out[[3]]$est[,1], rounddig), rounddig, format="f")
    beta2.1 <- formatC(round(est_out[[3]]$est[,2], rounddig), rounddig, format="f")
    se1.2 <- formatC(round(est_out[[3]]$se[,1], rounddig), rounddig, format="f")
    se2.1 <- formatC(round(est_out[[3]]$se[,2], rounddig), rounddig, format="f")
    lci1.2 <- formatC(round(est_out[[3]]$lower[,1], rounddig), rounddig, format="f")
    lci2.1 <- formatC(round(est_out[[3]]$lower[,2], rounddig), rounddig, format="f")
    uci1.2 <- formatC(round(est_out[[3]]$upper[,1], rounddig), rounddig, format="f")
    uci2.1 <- formatC(round(est_out[[3]]$upper[,2], rounddig), rounddig, format="f")
    #'  Merge into a data frame and organize
    out1.2 <- as.data.frame(cbind(beta1.2, se1.2, lci1.2, uci1.2)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        StudyArea = rep(area, nrow(.)),
        Transition = rep("Trans.1->2", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta1.2) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(StudyArea, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out1.2) <- c("Species", "Season", "Study Area", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out2.1 <- as.data.frame(cbind(beta2.1, se2.1, lci2.1, uci2.1)) %>%
      mutate(
        Parameter = row.names(est_out[[3]]$est),
        Species = rep(spp, nrow(.)),
        Season = rep(season, nrow(.)),
        StudyArea = rep(area, nrow(.)),
        Transition = rep("Trans.2->1", nrow(.))
      ) %>%
      relocate(Parameter, .before = beta2.1) %>%
      relocate(Species, .before = Parameter) %>%
      relocate(Season, .before = Parameter) %>%
      relocate(StudyArea, .before = Parameter) %>%
      relocate(Transition, .before = Parameter)
    colnames(out2.1) <- c("Species", "Season", "Study Area", "Transition", "Parameter", "Estimate", "SE", "Lower", "Upper")
    out <- as.data.frame(rbind(out1.2, out2.1))
    return(out)
  }
  #'  Run each season and species-specific model through function
  md_smr_hmm <- hmm_out(spp_HMM_output[[1]], "Mule Deer", "Summer", "Okanogan") #md_HMM_smr
  md_wtr_hmm <- hmm_out(spp_HMM_output[[2]], "Mule Deer", "Winter", "Okanogan") #md_HMM_wtr
  elk_smr_hmm <- hmm_out(spp_HMM_output[[3]], "Elk", "Summer", "Northeast") #elk_HMM_smr
  elk_wtr_hmm <- hmm_out(spp_HMM_output[[4]], "Elk", "Winter", "Northeast") #elk_HMM_wtr
  wtd_smr_hmm <- hmm_out(spp_HMM_output[[5]], "White-tailed Deer", "Summer", "Northeast") #wtd_HMM_smr
  wtd_wtr_hmm <- hmm_out(spp_HMM_output[[6]], "White-tailed Deer", "Winter", "Northeast") #wtd_HMM_wtr
  coug_smr_hmm_OK <- hmm_out(spp_HMM_output[[7]], "Cougar", "Summer", "Okanogan") #coug_HMM_smr_OK
  coug_wtr_hmm_OK <- hmm_out(spp_HMM_output[[8]], "Cougar", "Winter", "Okanogan") #coug_HMM_wtr_OK
  coug_smr_hmm_NE <- hmm_out(spp_HMM_output[[9]], "Cougar", "Summer", "Northeast") #coug_HMM_smr_NE
  coug_wtr_hmm_NE <- hmm_out(spp_HMM_output[[10]], "Cougar", "Winter", "Northeast") #coug_HMM_wtr_NE
  wolf_smr_hmm_OK <- hmm_out(spp_HMM_output[[11]], "Wolf", "Summer", "Okanogan") #wolf_HMM_smr_OK
  wolf_wtr_hmm_OK <- hmm_out(spp_HMM_output[[12]], "Wolf", "Winter", "Okanogan") #wolf_HMM_wtr_OK
  wolf_smr_hmm_NE <- hmm_out(spp_HMM_output[[13]], "Wolf", "Summer", "Northeast") #wolf_HMM_smr_NE
  wolf_wtr_hmm_NE <- hmm_out(spp_HMM_output[[14]], "Wolf", "Winter", "Northeast") #wolf_HMM_wtr_NE
  
  #'  Gather prey and predator results to put into a single results table
  results_hmm_TransPr <- rbind(md_smr_hmm, md_wtr_hmm, elk_smr_hmm, elk_wtr_hmm, 
                               wtd_smr_hmm, wtd_wtr_hmm, coug_smr_hmm_OK, coug_wtr_hmm_OK, 
                               coug_smr_hmm_NE, coug_wtr_hmm_NE, wolf_smr_hmm_OK, 
                               wolf_wtr_hmm_OK, wolf_smr_hmm_NE, wolf_wtr_hmm_NE)
  
  results_hmm_TransPr_prey <- rbind(md_smr_hmm, md_wtr_hmm, elk_smr_hmm, elk_wtr_hmm, 
                                    wtd_smr_hmm, wtd_wtr_hmm) %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    mutate(
      Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
      Parameter = ifelse(Parameter == "TRI", "Terrain Ruggedness", Parameter),
      Parameter = ifelse(Parameter == "PercOpen", "Percent Open", Parameter),
      Parameter = ifelse(Parameter == "Dist2Road", "Nearest Road", Parameter),
      Parameter = ifelse(Parameter == "SnowCover1", "Snow Cover (Y)", Parameter),
      Parameter = ifelse(Parameter == "COUG_RSF", "Pr(Cougar)", Parameter),
      Parameter = ifelse(Parameter == "WOLF_RSF", "Pr(Wolf)", Parameter)
    ) 
  colnames(results_hmm_TransPr_prey) <- c("Species", "Season", "Study Area", 
                                          "Transition", "Parameter", "Estimate", 
                                          "SE", "CI95")
  results_hmm_TransPr_pred <- rbind(coug_smr_hmm_OK, coug_wtr_hmm_OK, coug_smr_hmm_NE, 
                                    coug_wtr_hmm_NE, wolf_smr_hmm_OK, wolf_wtr_hmm_OK, 
                                    wolf_smr_hmm_NE, wolf_wtr_hmm_NE) %>%
    unite(CI95, Lower, Upper, sep = ", ") %>%
    mutate(
      Parameter = ifelse(Parameter == "(Intercept)", "Intercept", Parameter),
      Parameter = ifelse(Parameter == "TRI", "Terrain Ruggedness", Parameter),
      Parameter = ifelse(Parameter == "PercOpen", "Percent Open", Parameter),
      Parameter = ifelse(Parameter == "Dist2Road", "Nearest Road", Parameter),
      Parameter = ifelse(Parameter == "SnowCover1", "Snow Cover (Y)", Parameter),
      Parameter = ifelse(Parameter == "MD_RSF", "Pr(Mule Deer)", Parameter),
      Parameter = ifelse(Parameter == "ELK_RSF", "Pr(Elk)", Parameter),
      Parameter = ifelse(Parameter == "WTD_RSF", "Pr(White-tailed Deer)", Parameter)
    ) 
  colnames(results_hmm_TransPr_pred) <- c("Species", "Season", "Study Area", 
                                          "Transition", "Parameter", "Estimate",
                                          "SE", "CI95")
  
  #'  Spread results so the coefficient effects are easier to compare between 
  #'  transition probabilities and across species
  #'  Prey HMM results
  results_hmm_wide_TransPr_prey <- results_hmm_TransPr_prey %>%  
    mutate(
      SE = paste0("(", SE, ")"),
    ) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("Intercept", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Terrain Ruggedness", c("Terrain Ruggedness (SE)", "Terrain Ruggedness 95% CI"), sep = "_") %>%
    separate("Percent Open", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
    separate("Nearest Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
    separate("Snow Cover (Y)", c("Snow Cover (Y) (SE)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
    separate("Pr(Cougar)", c("Pr(Cougar) (SE)", "Pr(Cougar) 95% CI"), sep = "_") %>%
    separate("Pr(Wolf)", c("Pr(Wolf) (SE)", "Pr(Wolf) 95% CI"), sep = "_") %>%
    arrange(match(Species, c("Mule Deer", "Elk", "White-tailed Deer")))

  #'  Predators HMM results
  results_hmm_wide_TransPr_pred <- results_hmm_TransPr_pred %>% 
    mutate(
      SE = paste0("(", SE, ")"),
    ) %>%
    unite(Est_SE, Estimate, SE, sep = " ") %>%
    unite(Est_SE_CI, Est_SE, CI95, sep = "_") %>%
    spread(Parameter, Est_SE_CI) %>%
    separate("Intercept", c("Intercept (SE)", "Intercept 95% CI"), sep = "_") %>%
    separate("Terrain Ruggedness", c("Terrain Ruggedness (SE)", "Terrain Ruggedness 95% CI"), sep = "_") %>%
    separate("Percent Open", c("Percent Open (SE)", "Percent Open 95% CI"), sep = "_") %>%
    separate("Nearest Road", c("Nearest Road (SE)", "Nearest Road 95% CI"), sep = "_") %>%
    separate("Snow Cover (Y)", c("Snow Cover (Y) (SE)", "Snow Cover (Y) 95% CI"), sep = "_") %>%
    separate("Pr(Mule Deer)", c("Pr(Mule Deer) (SE)", "Pr(Mule Deer) 95% CI"), sep = "_") %>%
    separate("Pr(Elk)", c("Pr(Elk) (SE)", "Pr(Elk) 95% CI"), sep = "_") %>%
    separate("Pr(White-tailed Deer)", c("Pr(White-tailed Deer) (SE)", "Pr(White-tailed Deer) 95% CI"), sep = "_") %>%
    group_by(Species) %>%
    arrange(match(`Study Area`, c("Okanogan", "Northeast")), .by_group = TRUE) %>%
    ungroup()

 
  ####  Back-transformed Results  ####
  #'  Back-transform HMM results to the real (natural) scale of the data
  #'  Extract parameter means, SE, and 95% CI on natural scale when all covariates
  #'  are held at their mean value (i.e., 0 since covariates are scaled)
  backtrans_params <- function(mod, spp, season, area) {
    
    #'  CIreal has 4 lists: [[1]] step length params, [[2]] turning angle concentration,
    #'  [[3]] transition probabilities, and [[4]] initial state for each track.
    #'  Step length includes 4-6 lists depending on if zeromass parameter is needed
    #'  Lists 1:3 are State1 mean, sd, zeromass, 4:6 are State2 mean, sd, zeromass
    #'  Transition probability included 4 lists: [[1]] staying in State1, [[2]] 
    #'  transition from State1 to State2, [[3]] transitioning from State2 to State1,
    #'  and [[4]] staying in State2
    ci_nat <- CIreal(mod)
    
    #'  Table of step lengths (in meters) and 95% CI
    steps_state1 <- c(ci_nat[[1]]$est[[1]], ci_nat[[1]]$lower[[1]], ci_nat[[1]]$upper[[1]])
    steps_state2 <- c(ci_nat[[1]]$est[[4]], ci_nat[[1]]$lower[[4]], ci_nat[[1]]$upper[[4]])
    steps_real <- as.data.frame(rbind(steps_state1, steps_state2))
    colnames(steps_real) <- c("Mean", "Lower", "Upper")
    steps_real <- rownames_to_column(steps_real, var = "State") %>%
      mutate(Species = spp,
             Season = season,
             StudyArea = area,
             State = ifelse(State == "steps_state1", "Encamped", "Exploratory"),
             Mean = round(Mean, 2),
             Lower = round(Lower, 2),
             Upper = round(Upper, 2)) %>%
      unite("95%CI", Lower:Upper, sep = " - ") %>%
      relocate(Species, .before = "State") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea")
    
    #'  Table of turning angles and concentrations
    turn_matrix <- matrix(c(0, 0, ci_nat[[2]]$est[[1]], ci_nat[[2]]$est[[2]]),nrow=2,ncol=2,byrow=TRUE)
    colnames(turn_matrix) <- c("Encamped", "Exploratory")
    rownames(turn_matrix) <- c("Mean", "Concentration")
    turn_real <- rownames_to_column(as.data.frame(turn_matrix), var = "Parameter") %>%
      mutate(Species = spp,
             Season = season,
             StudyArea = area,
             Encamped = round(Encamped, 2),
             Exploratory = round(Exploratory, 2)) %>%
      relocate(Species, .before = "Parameter") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea")
      
    #'  Table of transition probabilties and 95% CI
    trans_probs <- ci_nat[[3]]$est
    trans_real <- rownames_to_column(as.data.frame(trans_probs), var = "States") %>%
      mutate(States = ifelse(States == "encamped", "Encamped", "Exploratory"),
             Species = spp,
             Season = season,
             StudyArea = area,
             encamped = round(encamped, 2),
             exploratory = round(exploratory, 2)) %>%
      relocate(Species, .before = "States") %>%
      relocate(StudyArea, .after = "Species") %>%
      relocate(Season, .after = "StudyArea") 
    colnames(trans_real) <- c("Species", "Study Area", "Season", "Start State", "Pr(To Encamped)", "Pr(To Exploratory)")
 
    print(round(ci_nat[[1]]$est, 2))
    print(round(ci_nat[[2]]$est, 2))
    print(round(ci_nat[[3]]$est, 2))
    
    table_list <- list(steps_real, turn_real, trans_real)
    
    return(table_list)
  }
  md_smr_backtrans <- backtrans_params(spp_HMM_output[[1]], spp = "Mule Deer", season = "Summer", area = "Okanogan")
  md_wtr_backtrans <- backtrans_params(spp_HMM_output[[2]], spp = "Mule Deer", season = "Winter", area = "Okanogan")
  elk_smr_backtrans <- backtrans_params(spp_HMM_output[[3]], spp = "Elk", season = "Summer", area = "Northeast")
  elk_wtr_backtrans <- backtrans_params(spp_HMM_output[[4]], spp = "Elk", season = "Winter", area = "Northeast")
  wtd_smr_backtrans <- backtrans_params(spp_HMM_output[[5]], spp = "White-tailed Deer", season = "Summer", area = "Northeast")
  wtd_wtr_backtrans <- backtrans_params(spp_HMM_output[[6]], spp = "White-tailed Deer", season = "Winter", area = "Northeast")
  coug_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[7]], spp = "Cougar", season = "Summer", area = "Okanogan")
  coug_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[8]], spp = "Cougar", season = "Winter", area = "Okanogan")
  coug_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[9]], spp = "Cougar", season = "Summer", area = "Northeast")
  coug_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[10]], spp = "Cougar", season = "Winter", area = "Northeast")
  wolf_smr_backtrans_OK <- backtrans_params(spp_HMM_output[[11]], spp = "Wolf", season = "Summer", area = "Okanogan")
  wolf_wtr_backtrans_OK <- backtrans_params(spp_HMM_output[[12]], spp = "Wolf", season = "Winter", area = "Okanogan")
  wolf_smr_backtrans_NE <- backtrans_params(spp_HMM_output[[13]], spp = "Wolf", season = "Summer", area = "Northeast")
  wolf_wtr_backtrans_NE <- backtrans_params(spp_HMM_output[[14]], spp = "Wolf", season = "Winter", area = "Northeast")
  
  #'  Table for back-transformed step lengths
  all_steps_backtrans <- bind_rows(md_smr_backtrans[[1]], md_wtr_backtrans[[1]], 
                                   elk_smr_backtrans[[1]], elk_wtr_backtrans[[1]], 
                                   wtd_smr_backtrans[[1]], wtd_wtr_backtrans[[1]],
                                   coug_smr_backtrans_OK[[1]], coug_wtr_backtrans_OK[[1]], 
                                   coug_smr_backtrans_NE[[1]], coug_wtr_backtrans_NE[[1]], 
                                   wolf_smr_backtrans_OK[[1]], wolf_wtr_backtrans_OK[[1]], 
                                   wolf_smr_backtrans_NE[[1]], wolf_wtr_backtrans_NE[[1]]) %>%
    arrange(Species) %>%
    pivot_wider(names_from = "State", values_from = c("Mean", "95%CI")) %>%
    relocate('95%CI_Encamped', .after = "Mean_Encamped") 
  colnames(all_steps_backtrans) <- c("Species", "Study Area", "Season", "Mean Encamped", "95% CI Encamped", "Mean Exploratory", "95% CI Exploratory")
  
  #'  Table for back-transformed turning angles
  all_turns_backtrans <- bind_rows(md_smr_backtrans[[2]], md_wtr_backtrans[[2]], 
                                   elk_smr_backtrans[[2]], elk_wtr_backtrans[[2]], 
                                   wtd_smr_backtrans[[2]], wtd_wtr_backtrans[[2]],
                                   coug_smr_backtrans_OK[[2]], coug_wtr_backtrans_OK[[2]], 
                                   coug_smr_backtrans_NE[[2]], coug_wtr_backtrans_NE[[2]], 
                                   wolf_smr_backtrans_OK[[2]], wolf_wtr_backtrans_OK[[2]], 
                                   wolf_smr_backtrans_NE[[2]], wolf_wtr_backtrans_NE[[2]]) %>%
    arrange(Species) 
  colnames(all_turns_backtrans) <- c("Species", "Study Area", "Season", "Parameter", "Encamped", "Exploratory")
  
  #'  Table for back-transformed transition probabilities
  all_TransPr_backtrans <- bind_rows(md_smr_backtrans[[3]], md_wtr_backtrans[[3]], 
                                   elk_smr_backtrans[[3]], elk_wtr_backtrans[[3]], 
                                   wtd_smr_backtrans[[3]], wtd_wtr_backtrans[[3]],
                                   coug_smr_backtrans_OK[[3]], coug_wtr_backtrans_OK[[3]], 
                                   coug_smr_backtrans_NE[[3]], coug_wtr_backtrans_NE[[3]], 
                                   wolf_smr_backtrans_OK[[3]], wolf_wtr_backtrans_OK[[3]], 
                                   wolf_smr_backtrans_NE[[3]], wolf_wtr_backtrans_NE[[3]]) %>%
    arrange(Species) 

  
  
  
  #'  Next up: Stationary_State_Plots.R to visualize movement results

  