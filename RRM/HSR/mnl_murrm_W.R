# ################################################################# #
#### LOAD LIBRARY AND DEFINE CORE SETTINGS                       ####
# ################################################################# #
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\6_Mixed model - SP Rail\\Model")

### Clear memory
rm(list = ls())

### Load Apollo library
library(apollo)
library(tidyverse)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database <-  read.csv('database.hsr_calib.csv',header = T)
#database <-  read.csv('database.hsr_valid.csv',header = T)

database <- database %>% 
  filter(PURP == 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'mnl_murrm 1',
  modelDescr   = "MNL CRRM model",
  indivID      = "ID",
  outputDirectory = "C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\6_Mixed model - SP Rail/Model/Second Paper/rrm/HSR/W"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car     = 0,
                asc_bus     = 0,
                asc_hsr     = 0,
                asc_air     = 0,
                b_tt        = 0,
                b_fr        = 0,
                b_co        = 0,
                b_age2     = 0,
                b_age2_car     = 0,
                b_age3_car     = 0,
                b_age3_bus     = 0,
                b_age3_air     = 0,
                b_inc2         = 0,
                b_inc2_air     = 0,
                b_inc3         = 0,
                b_inc3_car     = 0,
                b_mu           = 7
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_hsr")

# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  
  ### List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  R = list()
  
  R[['car']]  =  asc_car + b_mu*(- 
            log(1 + exp(b_tt/b_mu*(TT_BUS - TT_CAR)))*BUS_AV -
            log(1 + exp(b_tt/b_mu*(TT_HSR - TT_CAR)))*HSR_AV -
            log(1 + exp(b_tt/b_mu*(TT_AIR - TT_CAR)))*AIR_AV -
            log(1 + exp(b_co/b_mu*(FA_BUS - CO_CAR)))*BUS_AV -
            log(1 + exp(b_co/b_mu*(FA_HSR - CO_CAR)))*HSR_AV -
            log(1 + exp(b_co/b_mu*(FA_AIR - CO_CAR)))*AIR_AV) + b_age2_car*AGE_2 + b_age3_car*AGE_3 +
    b_inc2*INC_2 + b_inc3_car*INC_3  
            
  R[['bus']]  = asc_bus + b_mu*(- 
    log(1 + exp(b_tt/b_mu*(TT_CAR - TT_BUS)))*DRIV_LIC -
    log(1 + exp(b_tt/b_mu*(TT_HSR - TT_BUS)))*HSR_AV -
    log(1 + exp(b_tt/b_mu*(TT_AIR - TT_BUS)))*AIR_AV -
    log(1 + exp(b_co/b_mu*(CO_CAR - FA_BUS)))*DRIV_LIC -
    log(1 + exp(b_co/b_mu*(FA_HSR - FA_BUS)))*HSR_AV -
    log(1 + exp(b_co/b_mu*(FA_AIR - FA_BUS)))*AIR_AV -
    log(1 + exp(b_fr/b_mu*(FR_HSR - FR_BUS)))*HSR_AV -
    log(1 + exp(b_fr/b_mu*(FR_AIR - FR_BUS)))*AIR_AV) + b_age2*AGE_2 + b_age3_bus*AGE_3 +
      b_inc2*INC_2 + b_inc3*INC_3  
  
  R[['hsr']] = asc_hsr + b_mu*(- 
    log(1 + exp(b_tt/b_mu*(TT_CAR - TT_HSR)))*DRIV_LIC -
    log(1 + exp(b_tt/b_mu*(TT_BUS - TT_HSR)))*BUS_AV -
    log(1 + exp(b_tt/b_mu*(TT_AIR - TT_HSR)))*AIR_AV -
    log(1 + exp(b_co/b_mu*(CO_CAR - FA_HSR)))*DRIV_LIC -
    log(1 + exp(b_co/b_mu*(FA_BUS - FA_HSR)))*BUS_AV -
    log(1 + exp(b_co/b_mu*(FA_AIR - FA_HSR)))*AIR_AV -
    log(1 + exp(b_fr/b_mu*(FR_BUS - FR_HSR)))*BUS_AV -
    log(1 + exp(b_fr/b_mu*(FR_AIR - FR_HSR)))*AIR_AV)
  
  
  R[['air']] = asc_air + b_mu*(- 
    log(1 + exp(b_tt/b_mu*(TT_CAR - TT_AIR)))*DRIV_LIC -
    log(1 + exp(b_tt/b_mu*(TT_BUS - TT_AIR)))*BUS_AV -
    log(1 + exp(b_tt/b_mu*(TT_HSR - TT_AIR)))*HSR_AV -
    log(1 + exp(b_co/b_mu*(CO_CAR - FA_AIR)))*DRIV_LIC -
    log(1 + exp(b_co/b_mu*(FA_BUS - FA_AIR)))*BUS_AV -
    log(1 + exp(b_co/b_mu*(FA_HSR - FA_AIR)))*HSR_AV -
    log(1 + exp(b_fr/b_mu*(FR_BUS - FR_AIR)))*BUS_AV -
    log(1 + exp(b_fr/b_mu*(FR_HSR - FR_AIR)))*HSR_AV) + b_age2*AGE_2 + b_age3_air*AGE_3 +
      b_inc2_air*INC_2 + b_inc3*INC_3 
  
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(car=1,bus=2,hsr=4,air=5),
    avail        = list(car=DRIV_LIC,bus=BUS_AV,hsr=HSR_AV,air=AIR_AV),
    choiceVar    = CHOICE,
    V            = R
  )
  
  ### Compute probabilities using MNL model
  P[['model']] = apollo_mnl(mnl_settings, functionality)
  
  P = apollo_panelProd(P,apollo_inputs,functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### MODEL ESTIMATION                                            ####
# ################################################################# #

model_crrm  <-  apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                              estimate_settings = list(printLevel = 3))

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

apollo_modelOutput(model_crrm, modelOutput_settings = list(printClassical = TRUE,printPVal = TRUE))

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #

apollo_saveOutput(model_crrm,saveOutput_settings = list(printPVal=T ,printT1=T,
                                                         printDiagnostics = T)) #ver saveoutput list para mais configuracoes

#### VTTS ####

#VTTS car
dtt_car <- (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_BUS - database$TT_CAR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_HSR - database$TT_CAR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_AIR - database$TT_CAR))+1)

dco_car <- (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_BUS - database$CO_CAR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_HSR - database$CO_CAR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_AIR - database$CO_CAR))+1)

vtts_car <-  dtt_car/dco_car*60

#VTTS bus
dtt_bus <- (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_CAR - database$TT_BUS))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_HSR - database$TT_BUS))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_AIR - database$TT_BUS))+1)

dco_bus <- (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$CO_CAR - database$FA_BUS))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_HSR - database$FA_BUS))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_AIR - database$FA_BUS))+1)

vtts_bus <- dtt_bus/dco_bus*60

#VTTS air
dtt_air <- (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_CAR - database$TT_AIR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_BUS - database$TT_AIR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_HSR - database$TT_AIR))+1)

dco_air <- (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$CO_CAR - database$FA_AIR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_BUS - database$FA_AIR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_HSR - database$FA_AIR))+1)

vtts_air <- dtt_air/dco_air*60

#VTTS rail
dtt_rail <- (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_CAR - database$TT_HSR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_BUS - database$TT_HSR))+1) +
  (-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_tt']]/model_crrm$estimate[['b_mu']])*(database$TT_AIR - database$TT_HSR))+1)

dco_rail <- (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$CO_CAR - database$FA_HSR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_BUS - database$FA_HSR))+1) +
  (-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])/(exp((-model_crrm$estimate[['b_co']]/model_crrm$estimate[['b_mu']])*(database$FA_AIR - database$FA_HSR))+1)

vtts_rail <- dtt_rail/dco_rail*60

mean(vtts_car)
sd(vtts_car)

mean(vtts_bus)
sd(vtts_bus)

mean(vtts_air)
sd(vtts_air)

mean(vtts_rail)
sd(vtts_rail)

vtts_models <- cbind(vtts_muRRM,vtts_muRRM_paper)
write_csv(vtts_models,'muRRM_RP.csv')

#### Elasticities ####
predictions_base = apollo_prediction(model_crrm, apollo_probabilities, apollo_inputs)

#tt car
elast_ttcar <-       (model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_car))+1) +
                     model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_air - database$time_car))+1) +
                     model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_car))+1) +                  
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_car))+1))*predictions_base[,'bus'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_car))+1))*predictions_base[,'bus'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_car))+1))*predictions_base[,'bus'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_car))+1))*predictions_base[,'air'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_car))+1))*predictions_base[,'air'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_car))+1))*predictions_base[,'air'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_car))+1))*predictions_base[,'rail'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_car))+1))*predictions_base[,'rail'] +
                     (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_car))+1))*predictions_base[,'rail'])*database$time_car

agg_ttcar <- sum(elast_ttcar*(predictions_base[,'car']/sum(predictions_base[,'car'])))

#bus
elast_ttbus <-          (model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_car - database$time_bus))+1) +
                        model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_air - database$time_bus))+1) +
                        model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_bus))+1) +                  
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_bus))+1))*predictions_base[,'car'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_bus))+1))*predictions_base[,'car'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_bus))+1))*predictions_base[,'car'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_bus))+1))*predictions_base[,'air'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_bus))+1))*predictions_base[,'air'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_bus))+1))*predictions_base[,'air'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_bus))+1))*predictions_base[,'rail'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_bus))+1))*predictions_base[,'rail'] +
                        (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_bus))+1))*predictions_base[,'rail'])*database$time_bus

agg_ttbus <- sum(elast_ttbus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

#metro
elast_ttair <-         (model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_car - database$time_air))+1) +
                       model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_air))+1) +
                       model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_air))+1) +                  
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_air))+1))*predictions_base[,'car'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_air))+1))*predictions_base[,'car'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_air))+1))*predictions_base[,'car'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_air))+1))*predictions_base[,'bus'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_air))+1))*predictions_base[,'bus'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_air))+1))*predictions_base[,'bus'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_air))+1))*predictions_base[,'rail'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_air))+1))*predictions_base[,'rail'] +
                       (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_rail - database$time_air))+1))*predictions_base[,'rail'])*database$time_air

agg_ttair <- sum(elast_ttair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

#metro
elast_ttrail <-           (model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_car - database$time_rail))+1) +
                          model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_rail))+1) +
                          model_crrm$estimate[['b_tt']]/(exp((-model_crrm$estimate[['b_tt']])*(database$time_air - database$time_rail))+1) +                  
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_rail))+1))*predictions_base[,'car'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_rail))+1))*predictions_base[,'car'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_rail))+1))*predictions_base[,'car'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_rail))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_rail))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_rail))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_car - database$time_rail))+1))*predictions_base[,'rail'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_bus - database$time_rail))+1))*predictions_base[,'rail'] +
                          (model_crrm$estimate[['b_tt']]/(exp((model_crrm$estimate[['b_tt']])*(database$time_air - database$time_rail))+1))*predictions_base[,'rail'])*database$time_rail

agg_ttrail <- sum(elast_ttrail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))

#co auto
elast_cocar <-            (model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_car))+1) +
                          model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_car))+1) +
                          model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_car))+1) +                  
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_car))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_car))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_car))+1))*predictions_base[,'bus'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_car))+1))*predictions_base[,'air'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_car))+1))*predictions_base[,'air'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_car))+1))*predictions_base[,'air'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_car))+1))*predictions_base[,'rail'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_car))+1))*predictions_base[,'rail'] +
                          (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_car))+1))*predictions_base[,'rail'])*database$cost_car

agg_cocar <- sum(elast_cocar*(predictions_base[,'car']/sum(predictions_base[,'car'])))

#co bus
elast_cobus <-               (model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_bus))+1) +
                              model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_bus))+1) +
                              model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_bus))+1) +                  
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_bus))+1))*predictions_base[,'car'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_bus))+1))*predictions_base[,'car'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_bus))+1))*predictions_base[,'car'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_bus))+1))*predictions_base[,'air'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_bus))+1))*predictions_base[,'air'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_bus))+1))*predictions_base[,'air'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_bus))+1))*predictions_base[,'rail'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_bus))+1))*predictions_base[,'rail'] +
                             (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_bus))+1))*predictions_base[,'rail'])*database$cost_bus

agg_cobus <- sum(elast_cobus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

#co air
elast_coair <-                  (model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_air))+1) +
                                model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_air))+1) +
                                model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_air))+1) +                  
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_air))+1))*predictions_base[,'car'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_air))+1))*predictions_base[,'car'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_air))+1))*predictions_base[,'car'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_air))+1))*predictions_base[,'bus'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_air))+1))*predictions_base[,'bus'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_air))+1))*predictions_base[,'bus'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_air))+1))*predictions_base[,'rail'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_air))+1))*predictions_base[,'rail'] +
                                (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_rail - database$cost_air))+1))*predictions_base[,'rail'])*database$cost_air

agg_coair <- sum(elast_coair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

#co rail
elast_corail <-                    (model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_rail))+1) +
                                   model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_rail))+1) +
                                   model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_rail))+1) +                  
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_rail))+1))*predictions_base[,'car'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_rail))+1))*predictions_base[,'car'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_rail))+1))*predictions_base[,'car'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_rail))+1))*predictions_base[,'bus'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_rail))+1))*predictions_base[,'bus'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_rail))+1))*predictions_base[,'bus'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_car - database$cost_rail))+1))*predictions_base[,'rail'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_bus - database$cost_rail))+1))*predictions_base[,'rail'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$cost_air - database$cost_rail))+1))*predictions_base[,'rail'])*database$cost_rail

agg_corail <- sum(elast_corail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))

#acc bus
elast_accbus <-              (model_crrm$estimate[['b_acc']]/(exp((-model_crrm$estimate[['b_acc']])*(database$access_air - database$access_bus))+1) +
                              model_crrm$estimate[['b_acc']]/(exp((-model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_bus))+1) +                  
                             (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_air - database$access_bus))+1))*predictions_base[,'air'] +
                             (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_bus))+1))*predictions_base[,'air'] +
                             (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_air - database$access_bus))+1))*predictions_base[,'rail'] +
                             (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_bus))+1))*predictions_base[,'rail'])*database$access_bus

agg_accbus <- sum(elast_accbus*(predictions_base[,'bus']/sum(predictions_base[,'bus'])))

#acc air
elast_accair <-                 (model_crrm$estimate[['b_acc']]/(exp((-model_crrm$estimate[['b_acc']])*(database$access_bus - database$access_air))+1) +
                                 model_crrm$estimate[['b_acc']]/(exp((-model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_air))+1) +                  
                                (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_bus - database$access_air))+1))*predictions_base[,'bus'] +
                                (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_air))+1))*predictions_base[,'bus'] +
                                (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_bus - database$access_air))+1))*predictions_base[,'rail'] +
                                (model_crrm$estimate[['b_acc']]/(exp((model_crrm$estimate[['b_acc']])*(database$access_rail - database$access_air))+1))*predictions_base[,'rail'])*database$access_air

agg_accair <- sum(elast_accair*(predictions_base[,'air']/sum(predictions_base[,'air'])))

#acc rail
elast_accrail <-                   (model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$access_bus - database$access_rail))+1) +
                                    model_crrm$estimate[['b_co']]/(exp((-model_crrm$estimate[['b_co']])*(database$access_air - database$access_rail))+1) +                  
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$access_bus - database$access_rail))+1))*predictions_base[,'bus'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$access_air - database$access_rail))+1))*predictions_base[,'bus'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$access_bus - database$access_rail))+1))*predictions_base[,'rail'] +
                                   (model_crrm$estimate[['b_co']]/(exp((model_crrm$estimate[['b_co']])*(database$access_air - database$access_rail))+1))*predictions_base[,'rail'])*database$access_rail

agg_accrail <- sum(elast_accrail*(predictions_base[,'rail']/sum(predictions_base[,'rail'])))

#enumerado
agg_ttcar
agg_ttbus
agg_ttair
agg_ttrail
agg_cocar
agg_cobus
agg_coair
agg_corail
agg_accbus
agg_accair
agg_accrail
