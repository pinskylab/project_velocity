  
dom_centroid <- function(years, dat){
  time=years
  test_dat=dat
  dom_table <- data.frame(niche=as.numeric(), modelrun=as.numeric(), rcp=as.numeric())
  mod_full <- lm(centroid~nicheModel + modelrun_fact + rcp_fact + modelrun_fact:rcp_fact:nicheModel, data=test_dat[test_dat$years == time,])

  mod_niche <- lm(centroid~nicheModel, data=test_dat[test_dat$years == time,])
  mod_modelrun <- lm(centroid~modelrun_fact, data=test_dat[test_dat$years == time,])
  mod_rcp <- lm(centroid~rcp_fact, data=test_dat[test_dat$years == time,])

    dom_table[1,] <-  data.frame(niche=summary(mod_niche)$adj.r.squared, modelrun=summary(mod_modelrun)$adj.r.squared, rcp=summary(mod_rcp)$adj.r.squared)

  mod_niche_modelrun <- lm(centroid~nicheModel + modelrun_fact + modelrun_fact:nicheModel, data=test_dat[test_dat$years == time,])
  mod_niche_rcp <- lm(centroid~nicheModel + rcp_fact + rcp_fact:nicheModel, data=test_dat[test_dat$years == time,])
  mod_modelrun_rcp <- lm(centroid~modelrun_fact + rcp_fact + modelrun_fact:rcp_fact, data=test_dat[test_dat$years == time,])

  dom_table[2,] <- data.frame(niche=((summary(mod_niche_modelrun)$adj.r.squared - summary(mod_modelrun)$adj.r.squared) + (summary(mod_niche_rcp)$adj.r.squared - summary(mod_rcp)$adj.r.squared))/2, 
                    modelrun=((summary(mod_niche_modelrun)$adj.r.squared - summary(mod_niche)$adj.r.squared) + (summary(mod_modelrun_rcp)$adj.r.squared - summary(mod_rcp)$adj.r.squared))/2,  
                    rcp=((summary(mod_niche_rcp)$adj.r.squared - summary(mod_niche)$adj.r.squared) + (summary(mod_modelrun_rcp)$adj.r.squared - summary(mod_modelrun)$adj.r.squared))/2)  

  dom_table[3,] <- data.frame(niche=(summary(mod_full)$adj.r.squared - summary(mod_modelrun_rcp)$adj.r.squared), 
                            modelrun=(summary(mod_full)$adj.r.squared - summary(mod_niche_rcp)$adj.r.squared),  
                            rcp=(summary(mod_full)$adj.r.squared - summary(mod_niche_modelrun)$adj.r.squared))  

  totalSS <- sum(anova(mod_full)[2])
  param_SS <- data.matrix(anova(mod_full)[2])[5]

  dom_table[4,] <- data.frame(niche=mean(dom_table$niche), modelrun=mean(dom_table$modelrun), rcp=mean(dom_table$rcp))
  dom_table[5,] <- data.frame(niche=dom_table$niche[4]/sum(dom_table[4,]), modelrun=dom_table$modelrun[4]/sum(dom_table[4,]), rcp=dom_table$rcp[4]/sum(dom_table[4,]))
  final_tab <- data.frame(years = time, nicheModel=dom_table$niche[5]*(totalSS-param_SS), modelrun_fact=dom_table$modelrun[5]*(totalSS-param_SS), rcp_fact=dom_table$rcp[5]*(totalSS-param_SS), parameter=param_SS, stringsAsFactors = F)
  return(final_tab)
}


# SIMILAR FUNCTION BELOW, BUT FOR USE WITH MIXED EFFECTS MODELS FOR HABITAT CALCULATION
 
dom_habitat <- function(years, dat){
  time=years
  test_dat=dat
  dom_table <- data.frame(niche=as.numeric(), modelrun=as.numeric(), rcp=as.numeric())
  mod_full <- lmer(habitat~nicheModel + modelrun_fact + rcp_fact + modelrun_fact:rcp_fact:nicheModel + (1|method), data=test_dat[test_dat$years == time,])
  
  mod_niche <- lmer(habitat~nicheModel + (1|method), data=test_dat[test_dat$years == time,])
  mod_modelrun <- lmer(habitat~modelrun_fact + (1|method), data=test_dat[test_dat$years == time,])
  mod_rcp <- lmer(habitat~rcp_fact + (1|method), data=test_dat[test_dat$years == time,])
  
  dom_table[1,] <-  data.frame(niche=r.squaredGLMM(mod_niche)[1], modelrun=r.squaredGLMM(mod_modelrun)[1], rcp=r.squaredGLMM(mod_rcp)[1])
  
  mod_niche_modelrun <- lmer(habitat~nicheModel + modelrun_fact + modelrun_fact:nicheModel + (1|method), data=test_dat[test_dat$years == time,])
  mod_niche_rcp <- lmer(habitat~nicheModel + rcp_fact + rcp_fact:nicheModel + (1|method), data=test_dat[test_dat$years == time,])
  mod_modelrun_rcp <- lmer(habitat~modelrun_fact + rcp_fact + modelrun_fact:rcp_fact + (1|method), data=test_dat[test_dat$years == time,])
  
  dom_table[2,] <- data.frame(niche=((r.squaredGLMM(mod_niche_modelrun)[1] - r.squaredGLMM(mod_modelrun)[1]) + (r.squaredGLMM(mod_niche_rcp)[1] - r.squaredGLMM(mod_rcp)[1]))/2, 
                              modelrun=((r.squaredGLMM(mod_niche_modelrun)[1] - r.squaredGLMM(mod_niche)[1]) + (r.squaredGLMM(mod_modelrun_rcp)[1] - r.squaredGLMM(mod_rcp)[1]))/2,  
                              rcp=((r.squaredGLMM(mod_niche_rcp)[1] - r.squaredGLMM(mod_niche)[1]) + (r.squaredGLMM(mod_modelrun_rcp)[1] - r.squaredGLMM(mod_modelrun)[1]))/2)  
  
  dom_table[3,] <- data.frame(niche=(r.squaredGLMM(mod_full)[1] - r.squaredGLMM(mod_modelrun_rcp)[1]), 
                              modelrun=(r.squaredGLMM(mod_full)[1] - r.squaredGLMM(mod_niche_rcp)[1]),  
                              rcp=(r.squaredGLMM(mod_full)[1] - r.squaredGLMM(mod_niche_modelrun)[1]))  
  
  totalSS <- sum(anova(mod_full)[2]) + (sum(resid(mod_full)^2))
  param_SS <- sum(resid(mod_full)^2)
  
  dom_table[4,] <- data.frame(niche=mean(dom_table$niche), modelrun=mean(dom_table$modelrun), rcp=mean(dom_table$rcp))
  dom_table[5,] <- data.frame(niche=dom_table$niche[4]/sum(dom_table[4,]), modelrun=dom_table$modelrun[4]/sum(dom_table[4,]), rcp=dom_table$rcp[4]/sum(dom_table[4,]))
  final_tab <- data.frame(years = time, nicheModel=dom_table$niche[5]*(totalSS-param_SS), modelrun_fact=dom_table$modelrun[5]*(totalSS-param_SS), rcp_fact=dom_table$rcp[5]*(totalSS-param_SS), parameter=param_SS, stringsAsFactors = F)
  return(final_tab)
}



# SIMILAR FUNCTION BUT FOR PERCENTAGE CHANGE OF HABITAT
dom_percChan <- function(years, dat){
  time=years
  test_dat=dat
  dom_table <- data.frame(niche=as.numeric(), modelrun=as.numeric(), rcp=as.numeric())
  mod_full <- lm(habitat~nicheModel + modelrun_fact + rcp_fact + modelrun_fact:rcp_fact:nicheModel, data=test_dat[test_dat$years == time,])
  
  mod_niche <- lm(habitat~nicheModel, data=test_dat[test_dat$years == time,])
  mod_modelrun <- lm(habitat~modelrun_fact, data=test_dat[test_dat$years == time,])
  mod_rcp <- lm(habitat~rcp_fact, data=test_dat[test_dat$years == time,])
  
  dom_table[1,] <- data.frame(niche=summary(mod_niche)$adj.r.squared, modelrun=summary(mod_modelrun)$adj.r.squared, rcp=summary(mod_rcp)$adj.r.squared)
  
  mod_niche_modelrun <- lm(habitat~nicheModel + modelrun_fact + modelrun_fact:nicheModel, data=test_dat[test_dat$years == time,])
  mod_niche_rcp <- lm(habitat~nicheModel + rcp_fact + rcp_fact:nicheModel, data=test_dat[test_dat$years == time,])
  mod_modelrun_rcp <- lm(habitat~modelrun_fact + rcp_fact + modelrun_fact:rcp_fact, data=test_dat[test_dat$years == time,])
  
  dom_table[2,] <- data.frame(niche=((summary(mod_niche_modelrun)$adj.r.squared - summary(mod_modelrun)$adj.r.squared) + (summary(mod_niche_rcp)$adj.r.squared - summary(mod_rcp)$adj.r.squared))/2, 
                              modelrun=((summary(mod_niche_modelrun)$adj.r.squared - summary(mod_niche)$adj.r.squared) + (summary(mod_modelrun_rcp)$adj.r.squared - summary(mod_rcp)$adj.r.squared))/2,  
                              rcp=((summary(mod_niche_rcp)$adj.r.squared - summary(mod_niche)$adj.r.squared) + (summary(mod_modelrun_rcp)$adj.r.squared - summary(mod_modelrun)$adj.r.squared))/2)  
  
  dom_table[3,] <- data.frame(niche=(summary(mod_full)$adj.r.squared - summary(mod_modelrun_rcp)$adj.r.squared), 
                              modelrun=(summary(mod_full)$adj.r.squared - summary(mod_niche_rcp)$adj.r.squared),  
                              rcp=(summary(mod_full)$adj.r.squared - summary(mod_niche_modelrun)$adj.r.squared))  
  
  totalSS <- sum(anova(mod_full)[2])
  param_SS <- data.matrix(anova(mod_full)[2])[5]
  
  dom_table[4,] <- data.frame(niche=mean(dom_table$niche), modelrun=mean(dom_table$modelrun), rcp=mean(dom_table$rcp))
  dom_table[5,] <- data.frame(niche=dom_table$niche[4]/sum(dom_table[4,]), modelrun=dom_table$modelrun[4]/sum(dom_table[4,]), rcp=dom_table$rcp[4]/sum(dom_table[4,]))
  final_tab <- data.frame(years = time, nicheModel=dom_table$niche[5]*(totalSS-param_SS), modelrun_fact=dom_table$modelrun[5]*(totalSS-param_SS), rcp_fact=dom_table$rcp[5]*(totalSS-param_SS), parameter=param_SS, stringsAsFactors = F)
  return(final_tab)
}

 
