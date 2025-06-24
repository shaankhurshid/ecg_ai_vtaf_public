# Depends
library(data.table)
library(stringr)
library(plyr)
library(sktools)
library(lubridate)
library(APtools)
library(nricens)

# Load baseline
baseline <- fread(file='baseline.txt')

# Load outcomes
outcomes <- fread(file='AF_incidence.txt')

# Load ECG-AI inferences
vtaf_inf <- fread(file='ecg2af_mgh_survival_curve_af_inference_v2025_04_16.tsv')

# Classes
baseline[,':='(baseline_date = mdy(fstvisit_date),
               last_fu_date = mdy(end_date))]
outcomes[,':='(af_dx_date = mdy(af_dx_date))]

# Remove prevalent AF
baseline_no_af <- baseline[PrevalentAF==0]

# Merge in inferences and outcomes
setkey(vtaf_inf,mrn); setkey(baseline_no_af,mrn); setkey(outcomes,mrn)
baseline_no_af[vtaf_inf, ':='(ecg_ai_pred = i.survival_curve_af_prediction,
                              ecg_ai_datetime = i.datetime)]
baseline_no_af[outcomes,':='(af_dx_date = i.af_dx_date)]
baseline_no_af[,':='(incd_af = ifelse(!is.na(af_dx_date) & (af_dx_date >= baseline_date),1,0),
                     prev_af = ifelse(!is.na(af_dx_date) & (af_dx_date <= baseline_date),1,0))]
baseline_no_af <- baseline_no_af[prev_af==0]

# Remove missing ECG or missing CHARGE
analysis <- baseline_no_af[!is.na(ecg_ai_pred)]
analysis <- analysis[!is.na(CHARGE_AF)]

# Save out exclusion set
excluded <- baseline_no_af[is.na(ecg_ai_pred) | is.na(CHARGE_AF)]
write.csv(excluded,file='excluded_set.csv',row.names=F)

# Define CH-AI
analysis[,':='(ecg_ai_logit = log((ecg_ai_pred)/(1-ecg_ai_pred)))]
analysis[,':='(ch_ai = 0.35655*CHARGE_AF+0.44266*ecg_ai_logit)]
analysis[,':='(ecg_ai_lp = 0.61270*ecg_ai_logit)]

# Quantilize
analysis[,':='(ecg_ai_quartile = quantilize(ecg_ai_logit,4),
               charge_quartile = quantilize(CHARGE_AF,4),
               ch_ai_quartile = quantilize(ch_ai,4))]

analysis[,':='(ecg_ai_decile = quantilize(ecg_ai_logit,10),
               charge_decile = quantilize(CHARGE_AF,10),
               ch_ai_decile = quantilize(ch_ai,10))]

# Time var
analysis[,':='(time_to_af = ifelse(incd_af==1,(as.Date(af_dx_date) - as.Date(baseline_date))/365.25,
                                   (as.Date(last_fu_date) - as.Date(baseline_date) + 365.25)/365.25))]

# Discrimination
charge_ap <- APSurv(stime=analysis$time_to_af,status=analysis$incd_af,
                    marker=analysis$CHARGE_AF,t0.list=c(1,1.999),method='bootstrap',B=500,Plot=F)
ecgai_ap <- APSurv(stime=analysis$time_to_af,status=analysis$incd_af,
                   marker=analysis$ecg_ai_logit,t0.list=c(1,1.999),method='bootstrap',B=500,Plot=F)
chai_ap <- APSurv(stime=analysis$time_to_af,status=analysis$incd_af,
                  marker=analysis$ch_ai,t0.list=c(1,1.999),method='bootstrap',B=500,Plot=F)

# CompareAP
chai_v_charge <- CompareAP(stime=analysis$time_to_af,status=analysis$incd_af,
                           marker1=analysis$CHARGE_AF,marker2=analysis$ch_ai,
                           t0.list=c(1,1.999),method='bootstrap',B=500,Plot=F)

# Save out
write.csv(analysis,file='analysis.csv',row.names=F)

# Generate quartile estimates
output <- data.table(score=NA,quartile=NA,n_af_control=NA,n_af_intervention=NA,
                     n_total_control=NA,n_total_intervention=NA,
                     af_rate_control=NA,af_rate_control_lower=NA,af_rate_control_upper=NA,
                     af_rate_intervention=NA,af_rate_intervention_lower=NA,af_rate_intervention_upper=NA,
                     af_rate_difference=NA,af_rate_difference_lower=NA,af_rate_difference_upper=NA)
for (i in c('ecg_ai_quartile','charge_quartile','ch_ai_quartile')){
  for (j in 1:4){
    subset <- analysis[get(i)==j]
    control <- subset[AsTreated=='C']
    intervention <- subset[AsTreated=='I']
    af_control <- nrow(control[incd_af==1]); af_intervention <- nrow(intervention[incd_af==1])
    total_control <- nrow(control); total_intervention <- nrow(intervention)
    control_rate = nrow(control[incd_af==1])/nrow(control)
    intervention_rate = nrow(intervention[incd_af==1])/nrow(intervention)
    control_prop <- prop.test(x=af_control,n=total_control)
    intervention_prop <- prop.test(x=af_intervention,n=total_intervention)
    difference = intervention_rate - control_rate
    test <- prop.test(x=c(af_intervention,af_control),n=c(total_intervention,total_control))
    out <- data.table(score=i,quartile=j,n_af_control=af_control,n_af_intervention=af_intervention,
                      n_total_control=total_control,n_total_intervention=total_intervention,
                      af_rate_control=control_rate,af_rate_control_lower=control_prop$conf.int[1],af_rate_control_upper=control_prop$conf.int[2],
                      af_rate_intervention=intervention_rate,af_rate_intervention_lower=intervention_prop$conf.int[1],af_rate_intervention_upper=intervention_prop$conf.int[2],
                      af_rate_difference=difference,af_rate_difference_lower=test$conf.int[1],af_rate_difference_upper=test$conf.int[2])
    output <- rbind(output,out)
  }
}

# Generate decile estimates
output_decile <- data.table(score=NA,quartile=NA,n_af_control=NA,n_af_intervention=NA,
                            n_total_control=NA,n_total_intervention=NA,
                            af_rate_control=NA,af_rate_control_lower=NA,af_rate_control_upper=NA,
                            af_rate_intervention=NA,af_rate_intervention_lower=NA,af_rate_intervention_upper=NA,
                            af_rate_difference=NA,af_rate_difference_lower=NA,af_rate_difference_upper=NA)
for (i in c('ecg_ai_decile','charge_decile','ch_ai_decile')){
  for (j in 1:10){
    subset <- analysis[get(i)==j]
    control <- subset[AsTreated=='C']
    intervention <- subset[AsTreated=='I']
    af_control <- nrow(control[incd_af==1]); af_intervention <- nrow(intervention[incd_af==1])
    total_control <- nrow(control); total_intervention <- nrow(intervention)
    control_rate = nrow(control[incd_af==1])/nrow(control)
    intervention_rate = nrow(intervention[incd_af==1])/nrow(intervention)
    control_prop <- prop.test(x=af_control,n=total_control)
    intervention_prop <- prop.test(x=af_intervention,n=total_intervention)
    difference = intervention_rate - control_rate
    test <- prop.test(x=c(af_intervention,af_control),n=c(total_intervention,total_control))
    out <- data.table(score=i,quartile=j,n_af_control=af_control,n_af_intervention=af_intervention,
                      n_total_control=total_control,n_total_intervention=total_intervention,
                      af_rate_control=control_rate,af_rate_control_lower=control_prop$conf.int[1],af_rate_control_upper=control_prop$conf.int[2],
                      af_rate_intervention=intervention_rate,af_rate_intervention_lower=intervention_prop$conf.int[1],af_rate_intervention_upper=intervention_prop$conf.int[2],
                      af_rate_difference=difference,af_rate_difference_lower=test$conf.int[1],af_rate_difference_upper=test$conf.int[2])
    output_decile <- rbind(output_decile,out)
  }
}

# Generate quartile estimates (ITT)
output_itt <- data.table(score=NA,quartile=NA,n_af_control=NA,n_af_intervention=NA,
                         n_total_control=NA,n_total_intervention=NA,
                         af_rate_control=NA,af_rate_control_lower=NA,af_rate_control_upper=NA,
                         af_rate_intervention=NA,af_rate_intervention_lower=NA,af_rate_intervention_upper=NA,
                         af_rate_difference=NA,af_rate_difference_lower=NA,af_rate_difference_upper=NA)
for (i in c('ecg_ai_quartile','charge_quartile','ch_ai_quartile')){
  for (j in 1:4){
    subset <- analysis[get(i)==j]
    control <- subset[group=='C']
    intervention <- subset[group=='I']
    af_control <- nrow(control[incd_af==1]); af_intervention <- nrow(intervention[incd_af==1])
    total_control <- nrow(control); total_intervention <- nrow(intervention)
    control_rate = nrow(control[incd_af==1])/nrow(control)
    intervention_rate = nrow(intervention[incd_af==1])/nrow(intervention)
    control_prop <- prop.test(x=af_control,n=total_control)
    intervention_prop <- prop.test(x=af_intervention,n=total_intervention)
    difference = intervention_rate - control_rate
    test <- prop.test(x=c(af_intervention,af_control),n=c(total_intervention,total_control))
    out <- data.table(score=i,quartile=j,n_af_control=af_control,n_af_intervention=af_intervention,
                      n_total_control=total_control,n_total_intervention=total_intervention,
                      af_rate_control=control_rate,af_rate_control_lower=control_prop$conf.int[1],af_rate_control_upper=control_prop$conf.int[2],
                      af_rate_intervention=intervention_rate,af_rate_intervention_lower=intervention_prop$conf.int[1],af_rate_intervention_upper=intervention_prop$conf.int[2],
                      af_rate_difference=difference,af_rate_difference_lower=test$conf.int[1],af_rate_difference_upper=test$conf.int[2])
    output_itt <- rbind(output_itt,out)
  }
}

# Generate decile estimates (ITT)
output_decile_itt <- data.table(score=NA,quartile=NA,n_af_control=NA,n_af_intervention=NA,
                                n_total_control=NA,n_total_intervention=NA,
                                af_rate_control=NA,af_rate_control_lower=NA,af_rate_control_upper=NA,
                                af_rate_intervention=NA,af_rate_intervention_lower=NA,af_rate_intervention_upper=NA,
                                af_rate_difference=NA,af_rate_difference_lower=NA,af_rate_difference_upper=NA)
for (i in c('ecg_ai_decile','charge_decile','ch_ai_decile')){
  for (j in 1:10){
    subset <- analysis[get(i)==j]
    control <- subset[group=='C']
    intervention <- subset[group=='I']
    af_control <- nrow(control[incd_af==1]); af_intervention <- nrow(intervention[incd_af==1])
    total_control <- nrow(control); total_intervention <- nrow(intervention)
    control_rate = nrow(control[incd_af==1])/nrow(control)
    intervention_rate = nrow(intervention[incd_af==1])/nrow(intervention)
    control_prop <- prop.test(x=af_control,n=total_control)
    intervention_prop <- prop.test(x=af_intervention,n=total_intervention)
    difference = intervention_rate - control_rate
    test <- prop.test(x=c(af_intervention,af_control),n=c(total_intervention,total_control))
    out <- data.table(score=i,quartile=j,n_af_control=af_control,n_af_intervention=af_intervention,
                      n_total_control=total_control,n_total_intervention=total_intervention,
                      af_rate_control=control_rate,af_rate_control_lower=control_prop$conf.int[1],af_rate_control_upper=control_prop$conf.int[2],
                      af_rate_intervention=intervention_rate,af_rate_intervention_lower=intervention_prop$conf.int[1],af_rate_intervention_upper=intervention_prop$conf.int[2],
                      af_rate_difference=difference,af_rate_difference_lower=test$conf.int[1],af_rate_difference_upper=test$conf.int[2])
    output_decile_itt <- rbind(output_decile_itt,out)
  }
}

# Define missing ECG subset
missing_ecg <- baseline_no_af[is.na(ecg_ai_pred)]
control <- missing_ecg[group=='C']
intervention <- missing_ecg[group=='I']
af_control <- nrow(control[incd_af==1]); af_intervention <- nrow(intervention[incd_af==1])
total_control <- nrow(control); total_intervention <- nrow(intervention)
control_rate = nrow(control[incd_af==1])/nrow(control)
intervention_rate = nrow(intervention[incd_af==1])/nrow(intervention)
difference = intervention_rate - control_rate
test <- prop.test(x=c(af_intervention,af_control),n=c(total_intervention,total_control))

# Cox PH and interaction testing
mod_chai <- glm(incd_af ~ AsTreated + ch_ai + AsTreated*ch_ai_decile,family='binomial',data=analysis)
mod_chai_bin <- glm(incd_af ~ AsTreated*ifelse(ch_ai_decile==10,1,0),family='binomial',data=analysis)

# Event rate
km <- survfit(Surv(time_to_af,incd_af) ~ 1,data=analysis)
est <- 1-stepfun(km$time, c(1, km$surv))(2)
upper <- 1-stepfun(km$time, c(1, km$lower))(2)
lower <- 1-stepfun(km$time, c(1, km$upper))(2)
c(est,upper,lower)

# NRI
analysis[,dummy := 1]
analysis[,chai_top_decile := ifelse(ch_ai_decile==10,1,0)]
analysis[,charge_top_decile := ifelse(charge_decile==10,1,0)]

nri_chai_65 <- nricens(p.std=analysis$dummy,
                       p.new=analysis$chai_top_decile,
                       time=analysis$time_to_af, event=analysis$incd_af,
                       updown='category',cut=0.5,niter=500,t0=1.999)

nri_chai_charge <- nricens(p.std=analysis$charge_top_decile,
                           p.new=analysis$chai_top_decile,
                           time=analysis$time_to_af, event=analysis$incd_af,
                           updown='category',cut=0.5,niter=500,t0=1.999)

# Write outputs
write.csv(output,file='subgroup_as_treated.csv',row.names=F)
write.csv(output_itt,file='subgroup_itt.csv',row.names=F)
write.csv(output_decile,file='output_decile.csv',row.names=F)
write.csv(output_decile_itt,file='output_decile_itt.csv',row.names=F)