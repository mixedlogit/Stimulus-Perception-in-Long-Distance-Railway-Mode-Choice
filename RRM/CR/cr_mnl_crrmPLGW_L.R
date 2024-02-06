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

database <-  read.csv('database.cr_calib.csv',header = T)
#database <-  read.csv('database.cr_valid.csv',header = T)

database <- database %>% 
  filter(PURP != 1)

### Initialise code
apollo_initialise()

### Set core controls
apollo_control = list(
  modelName    = 'mnl_crrmPLGW 2',
  modelDescr   = "MNL CRRM model",
  indivID      = "ID",
  outputDirectory = "C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\6_Mixed model - SP Rail/Model/Second Paper/rrm/CR/L"
)

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #

### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(asc_car    = 0,
                asc_bus    = 0,
                asc_cr     = 0,
                asc_air    = 0,
                b_tt       = 0,
                b_fr             = 0,
                b_co             = 0,
                b_cocar          = 0,
                b_coair          = 0,
                b_age2           = 0,
                b_age2_air       = 0,
                b_age3           = 0,
                b_age3_air        = 0,
                b_inc2       = 0,
                b_inc3       = 0,
                b_inc3_air       = 0,
                b_fi             = 1
)

### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c("asc_cr")

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
  
  R[['car']]  =  asc_car - 
    log(1 + exp(b_tt*((TT_BUS - TT_CAR)/(TT_CAR**b_fi))))*BUS_AV -
    log(1 + exp(b_tt*((TT_CR - TT_CAR)/(TT_CAR**b_fi))))*CR_AV -
    log(1 + exp(b_tt*((TT_AIR - TT_CAR)/(TT_CAR**b_fi))))*AIR_AV -
    log(1 + exp((b_co*FA_BUS - b_cocar*CO_CAR)/(CO_CAR**b_fi)))*BUS_AV -
    log(1 + exp(b_co*FA_CR - b_cocar*CO_CAR)/(CO_CAR**b_fi))*CR_AV -
    log(1 + exp(b_coair*FA_AIR - b_cocar*CO_CAR)/(CO_CAR**b_fi))*AIR_AV +
    b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  R[['bus']]  = asc_bus -
    log(1 + exp(b_tt*((TT_CAR - TT_BUS)/(TT_BUS**b_fi))))*DRIV_LIC -
    log(1 + exp(b_tt*((TT_CR - TT_BUS)/(TT_BUS**b_fi))))*CR_AV -
    log(1 + exp(b_tt*((TT_AIR - TT_BUS)/(TT_BUS**b_fi))))*AIR_AV -
    log(1 + exp((b_cocar*CO_CAR - b_co*FA_BUS)/(FA_BUS**b_fi)))*DRIV_LIC -
    log(1 + exp(b_co*FA_CR - b_co*FA_BUS)/(FA_BUS**b_fi))*CR_AV -
    log(1 + exp(b_coair*FA_AIR - b_co*FA_BUS)/(FA_BUS**b_fi))*AIR_AV -
    log(1 + exp(b_fr*(FR_CR - FR_BUS)/(FR_BUS**b_fi)))*CR_AV -
    log(1 + exp(b_fr*(FR_AIR - FR_BUS)/(FR_BUS**b_fi)))*AIR_AV + 
    b_age2*AGE_2 + b_age3*AGE_3 +
    b_inc2*INC_2 + b_inc3*INC_3
  
  R[['cr']] = asc_cr - 
    log(1 + exp(b_tt*((TT_CAR - TT_CR)/(TT_CR**b_fi))))*DRIV_LIC -
    log(1 + exp(b_tt*((TT_BUS - TT_CR)/(TT_CR**b_fi))))*BUS_AV -
    log(1 + exp(b_tt*((TT_AIR - TT_CR)/(TT_CR**b_fi))))*AIR_AV -
    log(1 + exp(b_cocar*CO_CAR - b_co*FA_CR)/(FA_CR**b_fi))*DRIV_LIC -
    log(1 + exp(b_co*FA_BUS - b_co*FA_CR))*BUS_AV -
    log(1 + exp(b_coair*FA_AIR - b_co*FA_CR)/(FA_CR**b_fi))*AIR_AV -
    log(1 + exp(b_fr*(FR_BUS - FR_CR)/(FR_CR**b_fi)))*BUS_AV -
    log(1 + exp(b_fr*(FR_AIR - FR_CR)/(FR_CR**b_fi)))*AIR_AV
  
  
  R[['air']] = asc_air - 
    log(1 + exp(b_tt*((TT_CAR - TT_AIR)/(TT_AIR**b_fi))))*DRIV_LIC -
    log(1 + exp(b_tt*((TT_BUS - TT_AIR)/(TT_AIR**b_fi))))*BUS_AV -
    log(1 + exp(b_tt*((TT_CR - TT_AIR)/(TT_AIR**b_fi))))*HSR_AV -
    log(1 + exp(b_cocar*CO_CAR - b_coair*FA_AIR)/(FA_AIR**b_fi))*DRIV_LIC -
    log(1 + exp(b_co*FA_BUS - b_coair*FA_AIR)/(FA_AIR**b_fi))*BUS_AV -
    log(1 + exp(b_co*FA_CR - b_coair*FA_AIR)/(FA_AIR**b_fi))*HSR_AV -
    log(1 + exp(b_fr*(FR_BUS - FR_AIR)/(FR_AIR**b_fi)))*BUS_AV -
    log(1 + exp(b_fr*(FR_CR - FR_AIR)/(FR_AIR**b_fi)))*HSR_AV + 
    b_age2_air*AGE_2 + b_age3_air*AGE_3 +
    b_inc2*INC_2 + b_inc3_air*INC_3
  
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(car=1,bus=2,cr=3,air=5),
    avail        = list(car=DRIV_LIC,bus=BUS_AV,cr=CR_AV,air=AIR_AV),
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
dtt_car <- (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_BUS - database$TT_CAR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CR - database$TT_CAR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_AIR - database$TT_CAR))+1)

dco_car <- (-model_crrm$estimate[['b_co']])/(exp((-model_crrm$estimate[['b_co']])*(database$FA_BUS - database$CO_CAR))+1) +
  (-model_crrm$estimate[['b_co']])/(exp(-(model_crrm$estimate[['b_cocr']]*database$FA_CR - model_crrm$estimate[['b_co']]*database$CO_CAR))+1) +
  (-model_crrm$estimate[['b_co']])/(exp(-(model_crrm$estimate[['b_coair']]*database$FA_AIR - model_crrm$estimate[['b_co']]*database$CO_CAR))+1)

vtts_car <-  dtt_car/dco_car*60

#VTTS bus
dtt_bus <- (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CAR - database$TT_BUS))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CR - database$TT_BUS))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_AIR - database$TT_BUS))+1)

dco_bus <- (-model_crrm$estimate[['b_co']])/(exp((-model_crrm$estimate[['b_co']])*(database$CO_CAR - database$FA_BUS))+1) +
  (-model_crrm$estimate[['b_co']])/(exp(-(model_crrm$estimate[['b_cocr']]*database$FA_CR - model_crrm$estimate[['b_co']]*database$FA_BUS))+1) +
  (-model_crrm$estimate[['b_co']])/(exp(-(model_crrm$estimate[['b_coair']]*database$FA_AIR - model_crrm$estimate[['b_co']]*database$FA_BUS))+1)

vtts_bus <- dtt_bus/dco_bus*60

#VTTS air
dtt_air <- (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CAR - database$TT_AIR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_BUS - database$TT_AIR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CR - database$TT_AIR))+1)

dco_air <- (-model_crrm$estimate[['b_coair']])/(exp(-(model_crrm$estimate[['b_co']]*database$CO_CAR - model_crrm$estimate[['b_coair']]*database$FA_AIR))+1) +
  (-model_crrm$estimate[['b_coair']])/(exp(-(model_crrm$estimate[['b_co']]*database$FA_BUS - model_crrm$estimate[['b_coair']]*database$FA_AIR))+1) +
  (-model_crrm$estimate[['b_coair']])/(exp(-(model_crrm$estimate[['b_cocr']]*database$FA_CR - model_crrm$estimate[['b_coair']]*database$FA_AIR))+1)

vtts_air <- dtt_air/dco_air*60

#VTTS CR
dtt_cr <- (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_CAR - database$TT_CR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_BUS - database$TT_CR))+1) +
  (-model_crrm$estimate[['b_tt']])/(exp((-model_crrm$estimate[['b_tt']])*(database$TT_AIR - database$TT_CR))+1)

dco_cr <- (-model_crrm$estimate[['b_cocr']])/(exp(-(model_crrm$estimate[['b_co']]*database$CO_CAR - model_crrm$estimate[['b_cocr']]*database$FA_CR))+1) +
  (-model_crrm$estimate[['b_cocr']])/(exp(-(model_crrm$estimate[['b_co']]*database$FA_BUS - model_crrm$estimate[['b_cocr']]*database$FA_CR))+1) +
  (-model_crrm$estimate[['b_cocr']])/(exp(-(model_crrm$estimate[['b_coair']]*database$FA_AIR - model_crrm$estimate[['b_cocr']]*database$FA_CR))+1)

vtts_cr <- dtt_cr/dco_cr*60

mean(vtts_car)
mean(vtts_bus)
mean(vtts_air)
mean(vtts_cr)
sd(vtts_car)
sd(vtts_bus)
sd(vtts_air)
sd(vtts_cr)
  
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
