  #'  ======================================================
  #'  Code from Bassing et al. 2024 "Predator-prey space-use  
  #'  and landscape features influence animal movement  
  #'  behaviors in a large-mammal community". Ecology.
  #'  
  #'  Plotting Plot Stationary-State Probabilities
  #'  Washington Predator-Prey Project
  #'  ======================================================
  #'  Predict and plot stationary-state probabilities across range of covariate 
  #'  values. Create figures for publication.
  #'  ======================================================

  #'  Load libraries
  library(momentuHMM)
  library(ggplot2)
  library(tidyverse)  
  library(purrr)
  library(patchwork)
  
  #'  Load raw data with standardized covariates
  load("./Outputs/Telemetry_crwOut/crwOut_ALL_wCovs_for_pubs.RData") 

  #'  Load HMM results from Hidden_Markov_Movement_Models.R
  load("./Outputs/HMM_output/spp_HMM_output.RData") 
  
  
  ####  Stationary State Probs with MEAN Covariate Values  ####
  #'  -----------------------------------------------------
  #'  Functions to extract stationary state probabilities & plot predicted responses
  stay_probs_prey <- function(hmmm, snow) {
    #'  Calculate stationary state probs. for each state based on covariate data
    #'  for each time step
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    #'  Calculate stationary state probs. for each state when covariate data are
    #'  held at their mean value (0 b/c data are centered and scaled)
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0,
                                                   COUG_RSF = 0, WOLF_RSF = 0))
    print(stay_mu0)
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0,
                                                  COUG_RSF = 0, WOLF_RSF = 0),
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities for deer and elk
  stay_md_smr <- stay_probs_prey(spp_HMM_output[[1]], snow = 0)
  stay_md_wtr <- stay_probs_prey(spp_HMM_output[[2]], snow = 1)
  stay_elk_smr <- stay_probs_prey(spp_HMM_output[[3]], snow = 0)
  stay_elk_wtr <- stay_probs_prey(spp_HMM_output[[4]], snow = 1)
  stay_wtd_smr <- stay_probs_prey(spp_HMM_output[[5]], snow = 0)
  stay_wtd_wtr <- stay_probs_prey(spp_HMM_output[[6]], snow = 1)
  
  #'  Stationary probabilities for predators in the Okanogan
  stay_probs_pred_OK <- function(hmmm, snow) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0, MD_RSF = 0)) 
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0, MD_RSF = 0),  
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE) 
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[7]], snow = 0)
  stay_coug_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[8]], snow = 1)
  stay_wolf_smr_OK <- stay_probs_pred_OK(spp_HMM_output[[11]], snow = 0)
  stay_wolf_wtr_OK <- stay_probs_pred_OK(spp_HMM_output[[12]], snow = 1)
  
  #'  Stationary state probabilities for predators in the Northeast
  stay_probs_pred_NE <- function(hmmm, snow) {
    stay_pr <- stationary(hmmm)
    stay_pr <- stay_pr[[1]]
    stay_mu0 <- stationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                   SnowCover = snow, TRI = 0, 
                                                   ELK_RSF = 0, WTD_RSF = 0))
    print(stay_mu0) 
    #'  Plot stationary state probabilities and extract predicted estimates
    fig <- plotStationary(hmmm, covs = data.frame(Dist2Road = 0, PercOpen = 0,
                                                  SnowCover = snow, TRI = 0, 
                                                  ELK_RSF = 0, WTD_RSF = 0),    
                          col = c("red", "blue"), plotCI = TRUE, alpha = 0.95, return =  TRUE)
    stationary_probs <- list(stay_pr, fig)
    
    return(stationary_probs)
  }
  #'  Extract stationary state probabilities
  stay_coug_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[9]], snow = 0)
  stay_coug_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[10]], snow = 1)
  stay_wolf_smr_NE <- stay_probs_pred_NE(spp_HMM_output[[13]], snow = 0)
  stay_wolf_wtr_NE <- stay_probs_pred_NE(spp_HMM_output[[14]], snow = 1)
  
  #'  Function to extract stationary state probabilities for prettier plotting
  stay_covs <- function(stay, season, spp, area) {
    #'  Extract list of calculated stationary states for range of covariate values 
    #'  from HMM stationary output
    stay_covs <- stay[[2]]
    #'  Create empty list
    covs_out <- list()
    #'  Loop through all list elements (results for each covariate)
    for(l in 1:length(stay_covs)){
      #'  Hold list of interest
      cov <- stay_covs[[l]]
      #'  Add column indicating which behavioral state values belong to
      cov[[1]]$State <- "State 1"
      cov[[2]]$State <- "State 2"
      #'  Convert to data frame instead of list
      cov <- rbind(as.data.frame(cov[[1]]), as.data.frame(cov[[2]]))
      cov$State <- as.factor(cov$State)
      cov$Species <- spp
      cov$Season <- season
      cov$StudyArea <- area
      #'  Append to new list of data frames
      covs_out[[l]] <- cov
    }
    #'  Rename list elements based on covariate
    names(covs_out) <- names(stay_covs)
    return(covs_out)
  }
  md_smr_PrStay <- stay_covs(stay_md_smr, season = "Summer", spp = "Mule Deer", area = "Okanogan")
  md_wtr_PrStay <- stay_covs(stay_md_wtr, season = "Winter", spp = "Mule Deer", area = "Okanogan")
  elk_smr_PrStay <- stay_covs(stay_elk_smr, season = "Summer", spp = "Elk", area = "Northeast")
  elk_wtr_PrStay <- stay_covs(stay_elk_wtr, season = "Winter", spp = "Elk", area = "Northeast")
  wtd_smr_PrStay <- stay_covs(stay_wtd_smr, season = "Summer", spp = "White-tailed Deer", area = "Northeast")
  wtd_wtr_PrStay <- stay_covs(stay_wtd_wtr, season = "Winter", spp = "White-tailed Deer", area = "Northeast")
  coug_smr_OK_PrStay <- stay_covs(stay_coug_smr_OK, season = "Summer", spp = "Cougar", area = "Okanogan")
  coug_wtr_OK_PrStay <- stay_covs(stay_coug_wtr_OK, season = "Winter", spp = "Cougar", area = "Okanogan")
  coug_smr_NE_PrStay <- stay_covs(stay_coug_smr_NE, season = "Summer", spp = "Cougar", area = "Northeast")
  coug_wtr_NE_PrStay <- stay_covs(stay_coug_wtr_NE, season = "Winter", spp = "Cougar", area = "Northeast")
  wolf_smr_OK_PrStay <- stay_covs(stay_wolf_smr_OK, season = "Summer", spp = "Wolf", area = "Okanogan")
  wolf_wtr_OK_PrStay <- stay_covs(stay_wolf_wtr_OK, season = "Winter", spp = "Wolf", area = "Okanogan")
  wolf_smr_NE_PrStay <- stay_covs(stay_wolf_smr_NE, season = "Summer", spp = "Wolf", area = "Northeast")
  wolf_wtr_NE_PrStay <- stay_covs(stay_wolf_wtr_NE, season = "Winter", spp = "Wolf", area = "Northeast")
  
  
  ####  Predator-Prey Stationary State Plots  ####
  #'  ----------------------------------------
  #'  Color schemes
  #'  Mule deer, elk, white-tailed deer ["#40B0A6", "#E66100", "#5D3A9B"]
  #'  Cougar, wolf ["#D41159", "#0C7BDC"]
  
  #'  Ungulate stationary state ~ Cougar RSF
  coug_effects <- rbind(md_smr_PrStay$COUG_RSF, md_wtr_PrStay$COUG_RSF,
                        elk_smr_PrStay$COUG_RSF, elk_wtr_PrStay$COUG_RSF,
                        wtd_smr_PrStay$COUG_RSF, wtd_wtr_PrStay$COUG_RSF) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_coug_plot <- ggplot(coug_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(text = element_text(size = 14)) +
    #theme(legend.position="bottom") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled cougar RSF value") +
    ylab("Probability of prey being in faster state") 
  
  #'  Ungulate stationary state ~ Wolf RSF
  wolf_effects <- rbind(md_smr_PrStay$WOLF_RSF, 
                        elk_smr_PrStay$WOLF_RSF, elk_wtr_PrStay$WOLF_RSF,
                        wtd_smr_PrStay$WOLF_RSF, wtd_wtr_PrStay$WOLF_RSF) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_wolf_plot <- ggplot(wolf_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 14)) +
    xlim(-1.5, 2) + ylim(0, 1.0) +
    xlab("Scaled wolf RSF value") 
  
  
  #'  Predator stationary state ~ Mule Deer RSF
  md_effects <- rbind(coug_smr_OK_PrStay$MD_RSF, 
                      wolf_smr_OK_PrStay$MD_RSF) %>%  
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_md_plot <- ggplot(md_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    # theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 14)) +
    theme(legend.position = "none") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled mule deer RSF value") +
    ylab("Probability of predator being in faster state") 
  
  #'  Predator stationary state ~ Elk RSF
  elk_effects <- rbind(coug_smr_NE_PrStay$ELK_RSF, coug_wtr_NE_PrStay$ELK_RSF,
                       wolf_smr_NE_PrStay$ELK_RSF) %>% 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_elk_plot <- ggplot(elk_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +    
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 14)) +
    theme(legend.position = "none") +
    xlim(-2, 2.5) + ylim(0, 1.0) +
    xlab("Scaled elk RSF value") 
  
  #'  Predator stationary state ~ White-tailed Deer RSF
  wtd_effects <- rbind(coug_smr_NE_PrStay$WTD_RSF, 
                       wolf_wtr_NE_PrStay$WTD_RSF) %>%  
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_wtd_plot <- ggplot(wtd_effects, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +    
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(text = element_text(size = 14)) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2.5, 2) + ylim(0, 1.0) +
    xlab("Scaled white-tailed deer RSF value") 

  
  #'  Patchwork figures together in panels
  PredEffect_onPrey_fig <- prey_coug_plot + prey_wolf_plot + 
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 2) +
      plot_annotation(tag_levels = 'a', 
                      title = 'Ungulate stationary state probabilities',
                      theme = theme(plot.title = element_text(size = 18))) & 
    theme(plot.tag = element_text(size = 18)) 
  PreyEffect_onPred_fig <- pred_md_plot + pred_elk_plot + pred_wtd_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'a', 
                    title = 'Predator stationary state probabilities',
                    theme = theme(plot.title = element_text(size = 18))) & 
    theme(plot.tag = element_text(size = 18)) 
  
  PredPrey_patchwork <- PreyEffect_onPred_fig / PredEffect_onPrey_fig +
    plot_annotation(tag_levels = 'a',
                    theme = theme(plot.title = element_text(size = 18))) &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16)) &
    theme(plot.tag = element_text(size = 18)) 
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PredEffect_onPrey_StationaryProb_plot.tiff", PredEffect_onPrey_fig, width = 11, height = 7, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PreyEffect_onPred_StationaryProb_plot.tiff", PreyEffect_onPred_fig, width = 11, height = 7, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/PredPrey_StationaryProb_plot.tiff", PredPrey_patchwork, width = 13, height = 14, dpi = 600, units = "in", device = 'tiff', compression = 'lzw')
  
  
  
  ####  Landscape Effect Stationary State Plots  ####
  #'  ----------------------------------------
  #'  Color schemes
  #'  Mule deer, elk, white-tailed deer ["#40B0A6", "#E66100", "#5D3A9B"]
  #'  Cougar, wolf ["#D41159", "#0C7BDC"]
  
  #'  Ungulate stationary state ~ TRI (only significant relationships)
  tri_effects_prey <- rbind(md_smr_PrStay$TRI, md_wtr_PrStay$TRI,
                            elk_smr_PrStay$TRI, elk_wtr_PrStay$TRI,
                            wtd_smr_PrStay$TRI) %>% # wtd_wtr_PrStay$TRI
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_tri_plot <- ggplot(tri_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Ungulate movement in response to habitat complexity")
  
  #'  Ungulate stationary state ~ Distance to Road (only significant relationships)
  dist2rd_effects_prey <- rbind(md_smr_PrStay$Dist2Road, md_wtr_PrStay$Dist2Road,
                                elk_smr_PrStay$Dist2Road, #elk_wtr_PrStay$Dist2Road,
                                wtd_smr_PrStay$Dist2Road) %>% #, wtd_wtr_PrStay$Dist2Road
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_dist2rd_plot <- ggplot(dist2rd_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) + 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100", "#5D3A9B")) +
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    #theme(legend.position="bottom") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Prey Species"), color = guide_legend(title = "Prey Species")) +
    xlab("Scaled distance to road") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Ungulate stationary state probabilities in response to distance to nearest road")
  
  #'  Ungulate stationary state ~ Habitat Openness (only significant relationships)
  percopen_effects_prey <- rbind(md_smr_PrStay$PercOpen, md_wtr_PrStay$PercOpen,  
                                 elk_smr_PrStay$PercOpen) %>% #, elk_wtr_PrStay$PercOpen,
                                 #wtd_smr_PrStay$PercOpen, wtd_wtr_PrStay$PercOpen) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  prey_open_plot <- ggplot(percopen_effects_prey, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#40B0A6", "#E66100")) +  #, "#5D3A9B"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#40B0A6", "#E66100")) + #, "#5D3A9B"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(legend.position="none") +
    #theme(legend.position="bottom") +
    xlim(-1.5, 2.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Prey Species"), color = guide_legend(title = "Prey Species")) +
    xlab("Scaled percent open habitat") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Ungulate stationary state probabilities in response to percentage of open habitat")
  
  #'  Predator stationary state ~ TRI (only significant relationships)
  tri_effects_pred_OK <- rbind(coug_smr_OK_PrStay$TRI, coug_wtr_OK_PrStay$TRI,
                               wolf_smr_OK_PrStay$TRI) %>% #, wolf_wtr_OK_PrStay$TRI
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_tri_plot <- ggplot(tri_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +   #, "#FFC20A"  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI, Okanogan") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to habitat complexity, Okanogan")
  
  tri_effects_pred_NE <- rbind(coug_smr_NE_PrStay$TRI, coug_wtr_NE_PrStay$TRI,
                               wolf_smr_NE_PrStay$TRI, wolf_wtr_NE_PrStay$TRI) %>%
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_tri_plot <- ggplot(tri_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +   #, "#FFC20A"  
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-2, 4) + ylim(0, 1.0) +
    xlab("Scaled TRI, Northeast") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to habitat complexity, Northeast")
  
  #'  Predator stationary state ~ Distance to Road (only significant relationshp)
  dist2rd_effects_pred_OK <- rbind(coug_smr_OK_PrStay$Dist2Road, coug_wtr_OK_PrStay$Dist2Road,
                                   wolf_wtr_OK_PrStay$Dist2Road) %>% #wolf_smr_OK_PrStay$Dist2Road
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_dist2rd_plot <- ggplot(dist2rd_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled distance to road, OK") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to distance to nearest road, Okanogan")
  
  dist2rd_effects_pred_NE <- rbind(coug_smr_NE_PrStay$Dist2Road, coug_wtr_NE_PrStay$Dist2Road,
                                   wolf_smr_NE_PrStay$Dist2Road) %>% #, wolf_wtr_NE_PrStay$Dist2Road
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_dist2rd_plot <- ggplot(dist2rd_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A" 
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"     
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    #theme(legend.position="bottom") +
    xlim(-1, 4.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled distance to road, NE") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to distance to nearest road, Northeast")
  
  #'  Predator stationary state ~ Open Habitat (only significant relationships)
  percopen_effects_pred_OK <- rbind(coug_wtr_OK_PrStay$PercOpen, #coug_smr_OK_PrStay$PercOpen, 
                                    wolf_smr_OK_PrStay$PercOpen, wolf_wtr_OK_PrStay$PercOpen) %>% 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_OK_open_plot <- ggplot(percopen_effects_pred_OK, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    #theme(legend.position="bottom") +
    xlim(-1.5, 1.5) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled percent open, OK") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to percentage of open habitat, Okanogan")
  
  percopen_effects_pred_NE <- rbind(coug_smr_NE_PrStay$PercOpen, coug_wtr_NE_PrStay$PercOpen) %>% 
                                    #wolf_smr_NE_PrStay$PercOpen, wolf_wtr_NE_PrStay$PercOpen 
    filter(!State == "State 1") %>%
    mutate(Species = as.factor(Species),
           Season = as.factor(Season),
           StudyArea = as.factor(StudyArea)) %>% 
    dplyr::select(-c(StudyArea, State))
  pred_NE_open_plot <- ggplot(percopen_effects_pred_NE, aes(x = cov, y = est, colour = Species, linetype = Season)) + 
    geom_line(size = 0.75) + 
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c("#D41159", "#0C7BDC")) +  #, "#FFC20A"
    #'  Add confidence intervals
    geom_ribbon(aes(ymin = lci, ymax = uci, fill = Species), alpha = 0.3, colour = NA) +
    scale_fill_manual(values=c("#D41159", "#0C7BDC")) +     #, "#FFC20A"
    #'  Get rid of lines and gray background
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(color = 'black')) +
    theme(axis.title.y = element_blank()) +
    theme(legend.position="none") +
    xlim(-1, 3) + ylim(0, 1.0) +
    guides(fill = guide_legend(title = "Predator Species"), color = guide_legend(title = "Predator Species")) +
    xlab("Scaled percent open, NE") +
    ylab("Probability of being in faster state") #+
    #labs(title = "Predator stationary state probabilities in response to percentage of open habitat, Northeast")
  
  #'  Save
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_tri_plot.tiff", prey_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_dsit2rd_plot.tiff", prey_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/prey_percopen_plot.tiff", prey_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_tri_plot.tiff", pred_OK_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_dsit2rd_plot.tiff", pred_OK_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_OK_percopen_plot.tiff", pred_OK_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_tri_plot.tiff", pred_NE_tri_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_dsit2rd_plot.tiff", pred_NE_dist2rd_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/pred_NE_percopen_plot.tiff", pred_NE_open_plot, width = 7, height = 7, dpi = 800, units = "in", device = 'tiff')
  
  
  ##'  Patchwork figures together in panels
  tri_fig <- prey_tri_plot + pred_OK_tri_plot + pred_NE_tri_plot +
      plot_layout(guides = 'collect') + 
      plot_layout(ncol = 3) +
      plot_annotation(tag_levels = 'a') &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16)) &
    theme(plot.tag = element_text(size = 18)) 
  
  open_road_fig <- prey_dist2rd_plot + pred_OK_dist2rd_plot + pred_NE_dist2rd_plot +
    prey_open_plot + pred_OK_open_plot + pred_NE_open_plot +
    plot_layout(guides = 'collect') + 
    plot_layout(ncol = 3) +
    plot_annotation(tag_levels = 'a',
                    title = 'Effect of distance to nearst road and open habitat on stationary state probabilities') & 
    theme(title = element_text(size = 18)) &
    theme(axis.title = element_text(size = 16)) &
    theme(axis.text = element_text(size = 16)) &
    theme(legend.text = element_text(size = 16)) &
    theme(legend.title = element_text(size = 16)) &
    theme(plot.tag = element_text(size = 18)) 
  
  ggsave("./Outputs/Figures for ms/HMM Stationary States/TRI_StationaryProb_plot.tiff", tri_fig, width = 13, height = 7, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  ggsave("./Outputs/Figures for ms/HMM Stationary States/OpenHabitat-RoadDist_StationaryProb_plot.tiff", open_road_fig, width = 13, height = 13, dpi = 800, units = "in", device = 'tiff', compression = 'lzw')
  
  
  
  
 
  