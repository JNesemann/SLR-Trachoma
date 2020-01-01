#Final SLR Data analysis

#To do:
  #do bootstrapping for TI - JN Done
  #likelihood ration or PPV/NPV of test - JN Done w/ bca
  #change village level prevalence CIs to binomial and do for TI - JN done
  #Wilcoxon sign rank test to compare paired prevalence estimates for each village ?and overall? - JN Done
  #stratify by age and compare SLR vs smartphone for TI and TF prevalence (SLR probs better for 0-2yo) - 
  #comparing SLR/field vs smartphone/field kappas--wrtie function in R that calculates difference in ks then bootstrap the difference. 
    #is it 0 or not? - JN attempted
  #figure out an LCA as gold standard - JN attempted
  #use LCA gold standard to calculate sens and spec - JN Done
  #bootstrap CIs around sens/spec calculations using "yardstick" package - JN Done with bcas

library(tidyverse)
library(skimr)
library(readr)
library(readxl)
library(DescTools)
library(irr)
library(purrr)

################################
##     IMPORT AND CLEAN      ##
################################
master_key_import <- read_excel("Masked Number Master Key (1).xlsx")
TIRET3_PCR_import <- read_excel("TIRET3_PCR_for John.xlsx")
SLR_import <- read_csv("SLRTrachoma_DATA_2019-11-27_1325.csv")
TIRET3_PhotoExamData13Villages <- read_excel("TIRET3-PhotoExamDatain13Villages.xlsx")

#cleaning master key
master_key <- master_key_import %>%
  dplyr::select(number, camera, `masked number`, `repeat`) %>%
  rename(mask = `masked number`,
         repeat.instance=`repeat`) %>%
  mutate(repeat.instance=if_else(is.na(repeat.instance), 1, repeat.instance))

#cleaning and renaming TIRET 3 PCR
PCR <- TIRET3_PCR_import %>%
  dplyr::select(Group, `Random#`, PoolID, Pool_PCR_Results, Individual_PCR_Results, `Exam(0/1)`, Examiner) %>%
  rename(group = Group,
         number = `Random#`,
         poolid=PoolID,
         pool_pcr = Pool_PCR_Results,
         individual_pcr = Individual_PCR_Results,
         exam_yn = `Exam(0/1)`,  
         examiner = Examiner) %>%
  #creating new pcr variable
  mutate(newindpcr=if_else(pool_pcr=="Negative", 0,
                           if_else(pool_pcr=="Positive" & individual_pcr %in% c("Negative","equivocal"), 0,
                                   if_else(pool_pcr=="Positive" & individual_pcr=="Positive", 1, NA_real_))))

#cleaning and renaming TIRET3_PhotoExamData13Villages
TIRET3_villages <- TIRET3_PhotoExamData13Villages %>%
  dplyr::select(Clinicalexam, Photo, `Random#`, Sex, `Stateteam-code`, Tuber, Age, ID) %>%
  rename(number = `Random#`,
         id=ID,
         clinical_exam = Clinicalexam,
         photo = Photo,
         sex = Sex,
         state_code = `Stateteam-code`,
         age = Age,
         tuber = Tuber,
         photo = Photo) %>%
  #creating numeric fieldgrade
  mutate(clinic_ti_yn = if_else(!is.na(clinical_exam), as.numeric(grepl("TI", clinical_exam)), NA_real_),
                  clinic_tf_yn = if_else(!is.na(clinical_exam), as.numeric(grepl("TF", clinical_exam)), NA_real_),
                  clinic_tfti_yn = if_else(!is.na(clinical_exam), as.numeric(grepl("TF/TI", clinical_exam)), NA_real_),
                  number=as.numeric(number))

#merging PCR data with complete data in TIRET3_PhotoExamData13Villages, this will provide the stateteam-code which I am using as a proxy for village ID
PCR_exam <- full_join(PCR, TIRET3_villages, by = "number")

#creating one large dataset 
master1 <- SLR_import %>%
  #clean and separate graders 1 and 2
  rename(id = record_id,
         complete = trachoma_grading_complete) %>%
  # JK: This mutate is new...
  mutate(tf_di=if_else(tf %in% c(1, 2),1,if_else(tf %in% c(3, 4),0,NA_real_)),
         ti_di=if_else(ti %in% c(1, 2),1,if_else(ti %in% c(3, 4),0,NA_real_))) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  # JK: Get rid of the 2 "test" id's, then convert to numeric
  filter(!(grepl("est", id))) %>%
  # JK: This group_by is new; 
  mutate(id=as.numeric(id)) %>%
  group_by(id) %>%
  mutate(sumtf=sum(tf_di, na.rm=TRUE),
         sumti=sum(ti_di, na.rm=TRUE),
         totalid=n(), # xtabs(data=xyz, ~sumtf+totalid, addNA=TRUE)
         tf_yn=case_when(totalid==3 & sumtf>=2 ~ 1,
                         totalid==3 & sumtf<2 ~ 0,
                         totalid==2 & sumtf==2 ~ 1,
                         totalid==2 & sumtf==0 ~ 0,
                         TRUE ~ NA_real_),
         ti_yn=case_when(totalid==3 & sumti>=2 ~ 1,
                         totalid==3 & sumti<2 ~ 0,
                         totalid==2 & sumti==2 ~ 1,
                         totalid==2 & sumti==0 ~ 0,
                         TRUE ~ NA_real_)) %>%
  dplyr::select(-complete, -sumtf, -sumti, -totalid) %>%
  dplyr::select(id, tf_yn, ti_yn, dde, tf_di, ti_di, everything()) %>% 
  #creating a wide dataset
  gather(field, value, tf_di:notes) %>%
  mutate(field_dde = paste(field, dde, sep = "_")) %>%
  dplyr::select(-field, -dde) %>%
  spread(field_dde, value, convert = TRUE) %>%
  # JK: Note that when you enter the dot below it means to use the current data, so you're merging into the current dataframe
  right_join(., master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  dplyr::select(number, camera.instance, id, tf_yn, ti_yn, everything()) %>%
  dplyr::select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_NA) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  dplyr::select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE) %>%
  left_join(., PCR_exam, by="number")
# JK: Now this is one big dataset. The variable names follow the pattern:
# SLR vs smartphone, then...
# whether it's a repeat (1 is the first instance; 2 is the second instance; so primary analyses should only be with 1)
# variable (tf /ti etc)
# grader (1=Blake, 2=John, NA=Jeremy)

################################
##          ANALYSES          ##
################################

#ARVO abstract cross tabulations
xtabs(data = master1, ~SLR_1_tf_yn + clinic_tf_yn + smartphone_1_tf_yn)
xtabs(data = master1, ~SLR_1_ti_yn + smartphone_1_ti_yn + clinic_ti_yn)

#xtabs determining which villages had positive PCR
xtabs(data=master1, ~state_code + newindpcr)

#summarizing demographic information
demoMF <- master1 %>%
  group_by(sex) %>%
  dplyr::summarize(n=n(),
            mean.age=mean(age))

################################
##         PREVALENCE         ##
################################
# The bootstrapping really only needed to account for inter-village clustering
# This creates a nested data frame (D), where all data from the same village is put on the same line
# So if we resample this data frame D, we will choose a random village from D, and then take all kids from that village
D <- master1 %>% 
  dplyr::select(number, age, sex, state_code, age, SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) %>%
  filter(!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn) & !is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn)) %>% 
  nest(-state_code)
head(D)
library(rsample)
set.seed(154234)
# The bs object is the boostrap object; we are creating separate populations with resampling
# You could alter the "times" option; usually use small number of replications as testing code because faster
# But then change to a larger number (9999?) for the final analysis
bs <- bootstraps(D, times = 9)
library(coxed) # This is for the bca confidence intervals, but it masks the "summarize" command so that means you have to use dplyr:: beforehand
library(purrr)

bs_mean <- map(bs$splits, ~as_tibble(.) %>% unnest %>%  # Note that the unnest here because you need to un-nest the data in each bootstrap resample (ie bs$splits), 
                 # and then take the mean. Then this essentially gives you the mean from each bootstrap resample in each of the rows
                 # From this distribution we can come up with the 95%CI. We can get it different ways, a simple way is to just
                 # take the lower 2.5% and upper 97.5% percentiles
                 # But what I have been working on the last 2 days is getting a bias corrected and accelerated CI, which is supposed to be better
                 dplyr::summarize(mean_age=mean(age, na.rm=TRUE),
                                  prop_slr_tf=mean(SLR_1_tf_yn, na.rm=TRUE),
                                  prop_smartphone_tf=mean(smartphone_1_tf_yn, na.rm=TRUE),
                                  prop_clinic_tf=mean(clinic_tf_yn, na.rm=TRUE),
                                  prop_slr_ti=mean(SLR_1_ti_yn, na.rm = TRUE),
                                  prop_smartphone_ti=mean(smartphone_1_ti_yn, na.rm = TRUE),
                                  prop_clinic_ti=mean(clinic_ti_yn, na.rm = TRUE))) %>% 
  bind_rows(.id = 'boots') 
# This bca gives the approximate bootstrap confidence interval --> analytic version of BCa applying to smoothly defined parameters in exponential families
# They are second order accurate and second order correct
# See DiCiccio and Efron, Statistical Science 1996; 11:189
ci_age <- bca(bs_mean$mean_age)
ci_slr.tf <- bca(bs_mean$prop_slr_tf)
ci_smart.tf <- bca(bs_mean$prop_smartphone_tf)
ci_clin.tf <- bca(bs_mean$prop_clinic_tf)
ci_slr.ti <- bca(bs_mean$prop_slr_ti)
ci_smart.ti <- bca(bs_mean$prop_smartphone_ti)
ci_clin.ti <- bca(bs_mean$prop_clinic_ti)
lowerupper <- c("low", "up")
means <- master1 %>%
  dplyr::summarize(mean_age=mean(age, na.rm=TRUE),
                   mean_slr.tf=mean(SLR_1_tf_yn, na.rm=TRUE),
                   mean_smart.tf=mean(smartphone_1_tf_yn, na.rm=TRUE),
                   mean_clin.tf=mean(clinic_tf_yn, na.rm=TRUE),
                   mean_slr.ti=mean(SLR_1_ti_yn, na.rm = TRUE),
                   mean_smart.ti=mean(smartphone_1_ti_yn, na.rm = TRUE),
                   mean_clin.ti=mean(clinic_ti_yn, na.rm = TRUE),
                   num_age=sum(!is.na(age)),
                   num_slr.tf=sum(!is.na(SLR_1_tf_yn)),
                   num_smart.tf=sum(!is.na(smartphone_1_tf_yn)),
                   num_clin.tf=sum(!is.na(clinic_tf_yn)),
                   num_slr.ti=sum(!is.na(SLR_1_ti_yn)),
                   num_smart.ti=sum(!is.na(smartphone_1_ti_yn)),
                   num_clin.ti=sum(!is.na(clinic_ti_yn))) %>%
  gather(stat_field, value, mean_age:num_clin.ti) %>%
  separate(stat_field, into=c("stat", "field"), sep="_") %>%
  spread(stat, value, convert=TRUE)

bs_mean_table <- data.frame(lowerupper, ci_age, ci_slr.tf, ci_smart.tf, ci_clin.tf, ci_slr.ti, ci_smart.ti, ci_clin.ti, stringsAsFactors=FALSE) %>%
  gather(ci_field, value, ci_age:ci_clin.ti) %>%
  separate(ci_field, into=c("ci", "field"), sep="_") %>%
  dplyr::select(-ci) %>%
  spread(lowerupper, value, convert=TRUE) %>%
  full_join(., means, by="field") %>%
  dplyr::select(field, num, mean, low, up)

# To check
confint(lm(data=filter(master1,!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)), age ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)), SLR_1_tf_yn ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)), smartphone_1_tf_yn ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)), clinic_tf_yn ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn)), age ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn)), SLR_1_ti_yn ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn)), smartphone_1_tf_yn ~ 1))
confint(lm(data=filter(master1,!is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn)), clinic_ti_yn ~ 1))

#graphing overall prevalence estimates with bootstrapped CIs
overall_means_plot <- bs_mean_table %>%
  filter(field != "age") %>%
  separate(field, into = c("method", "sign")) %>%
  ggplot(aes(x = sign, y = mean, fill= method)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin=low, ymax=up), width=0.2, position=position_dodge(0.9))

#wilcoxon signed rank test for overall prevalences
master1_filtered <- master1 %>%
  dplyr::select(number, age, sex, state_code, age, SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) %>%
  filter(!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn) & !is.na(SLR_1_ti_yn) & !is.na(smartphone_1_ti_yn) & !is.na(clinic_ti_yn))

overallwil_slr_smart_tf <- wilcox.test(master1_filtered$SLR_1_tf_yn, master1_filtered$smartphone_1_tf_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95) 
  #significant positive difference
overallwil_slr_clinic_tf <- wilcox.test(master1_filtered$SLR_1_tf_yn, master1_filtered$clinic_tf_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95)
  #no significant difference
overallwil_smart_clinic_tf <- wilcox.test(master1_filtered$smartphone_1_tf_yn, master1_filtered$clinic_tf_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95)
  #significant negative difference
overallwil_slr_smart_ti <- wilcox.test(master1_filtered$SLR_1_ti_yn, master1_filtered$smartphone_1_ti_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95)
overallwil_slr_clinic_ti <- wilcox.test(master1_filtered$SLR_1_ti_yn, master1_filtered$clinic_ti_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95)
overallwil_smart_clinic_ti <- wilcox.test(master1_filtered$smartphone_1_ti_yn, master1_filtered$clinic_ti_yn, alternative = "two.sided", mu=0, paired = T, conf.int = T, conf.level = 0.95)

#now for village level prevalences and binom confidence intervals
#JN comment: I added TI into this table, not sure if you wanted me to do them separately...
villagenumsjk <- master1 %>%
  group_by(state_code) %>%
  dplyr::summarize(slr.tf_count=sum(SLR_1_tf_yn==1),
                   slr.tf_total=sum(SLR_1_tf_yn==1 | SLR_1_tf_yn==0),
                   smart.tf_count=sum(smartphone_1_tf_yn==1),
                   smart.tf_total=sum(smartphone_1_tf_yn==1 | smartphone_1_tf_yn==0),
                   clin.tf_count=sum(clinic_tf_yn==1),
                   clin.tf_total=sum(clinic_tf_yn==1 | clinic_tf_yn==0),
                   slr.ti_count=sum(SLR_1_ti_yn==1),
                   slr.ti_total=sum(SLR_1_ti_yn==1 | SLR_1_ti_yn==0),
                   smart.ti_count=sum(smartphone_1_ti_yn==1),
                   smart.ti_total=sum(smartphone_1_ti_yn==1 | smartphone_1_ti_yn==0),
                   clin.ti_count=sum(clinic_ti_yn==1),
                   clin.ti_total=sum(clinic_ti_yn==1 | clinic_ti_yn==0)) %>%
  gather(methodfield, value, slr.tf_count:clin.ti_total) %>%
  separate(methodfield, into=c("method", "field"), sep="_") %>%
  spread(field, value, convert=TRUE)

villagemeansjk <- as.tibble(cbind(villagenumsjk, binconf(villagenumsjk$count, villagenumsjk$total)))
villagemeansjk %>% group_by(method) %>% dplyr::summarize(sum_total=sum(total)) # Confirming everyone accounted for...

#wilcoxon rank sum to compare prevalence estimates within villages
  #JN : actually now that I think about it a wilcoxon signed rank test is probably better as the data are paired (different photos of the same individual)
#attempting to do this for results nested by community in D, instead of writing it out for evey vilage
village_wilcox <- D %>%
  mutate(slr_smart_tf=map(data, ~wilcox.test(.$SLR_1_tf_yn, .$smartphone_1_tf_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)),
         slr_clinic_tf=map(data, ~wilcox.test(.$SLR_1_tf_yn, .$clinic_tf_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)),
         smart_clinic_tf=map(data, ~wilcox.test(.$smartphone_1_tf_yn, .$clinic_tf_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)),
         slr_smart_ti=map(data, ~wilcox.test(.$SLR_1_ti_yn, .$smartphone_1_ti_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)),
         slr_clinic_ti=map(data, ~wilcox.test(.$SLR_1_ti_yn, .$clinic_ti_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)),
         smart_clinic_ti=map(data, ~wilcox.test(.$smartphone_1_ti_yn, .$clinic_ti_yn, alternative = "two.sided", mu = 0, paired = TRUE, exact = F, 
                                             conf.int = F)))

#I can call each new column to see the results
village_wilcox$slr_smart_tf #significant difference in 2, 4, 6, 12, 
village_wilcox$slr_clinic_tf #significant difference in 8
village_wilcox$smart_clinic_tf #significant difference in 5, 8, 9, 12
village_wilcox$slr_smart_ti #5 & 9 p=NA
village_wilcox$slr_clinic_ti #significant difference in 1, 3, 4, 9=NA, 12
village_wilcox$smart_clinic_ti #significant difference in 1, 3, 4, 9=NA, 10=NA, 13=NA
#JN: I see what you meant when you said the resulst would be difficult to interpret there are a far number of p=NA or p=1
  #I cannot think of a better test to compare, a paired t-test would be a stretch since the distribution is not normal

#graphing village level prev estimates with CIs
village_means_plot <- villagemeansjk %>%
  separate(method, into = c("method", "sign")) %>%
  ggplot(aes(x = state_code, y = PointEst, fill = method)) +
  geom_bar(position="dodge", stat = "identity", alpha=0.6) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~sign, nrow = 1) +
  coord_flip()

################################
##         SENS/SPEC          ##
################################
#I made four tables for two index tests (smartphone, SLR) grading TI and TF, using the yardstick package and bca
#currently using clinical grades as glod standard but eventally will use an LCA per TL
library(yardstick)
options(yardstick.event_first = FALSE)
library(rsample)

addmargins(xtabs(data = master1_filtered, ~SLR_1_tf_yn + clinic_tf_yn))
master1_filtered %>% sens(., truth = as.factor(clinic_tf_yn), estimate = as.factor(SLR_1_tf_yn))
master1_filtered %>% spec(., truth = as.factor(clinic_tf_yn), estimate = as.factor(SLR_1_tf_yn))

#defining all the metrics I want to calculate
class_metrics <- metric_set(sens, spec, ppv, npv)

#now to get bootstrapped CIs, following same method as for prevalence/means above
#JN: will use the D object that is nested by village to account for clustering
set.seed(3323)
bs2 <- bootstraps(D, times = 9) #object containing the seperate resamples

#SLR TF
bs_slr_tf_metrics <- map(bs2$splits, ~as.tibble(.) %>% unnest %>%
  class_metrics(., truth = as.factor(clinic_tf_yn), estimate = as.factor(SLR_1_tf_yn))) %>%
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=T)
  #now calculating bca confidence intervals, using the same workaround
ci_sens <- bca(bs_slr_tf_metrics$sens)
ci_spec <- bca(bs_slr_tf_metrics$spec)
ci_npv <- bca(bs_slr_tf_metrics$npv)
ci_ppv <- bca(bs_slr_tf_metrics$ppv)
#calculating estimates 
slr_tf_estimates <- class_metrics(master1_filtered, truth = as.factor(clinic_tf_yn), estimate = as.factor(SLR_1_tf_yn)) %>%
  dplyr::select(-.estimator) %>%
  mutate(field=.metric) %>%
  dplyr::select(-.metric)
  #creating one data table with the metrics for SLR_tf
slr_tf_metrics <- data.frame(lowerupper, ci_sens, ci_spec, ci_npv, ci_ppv, stringsAsFactors = F) %>%
  gather(ci_field, value, ci_sens:ci_ppv) %>%
  separate(ci_field, into=c("ci", "field"), sep="_") %>%
  dplyr::select(-ci) %>%
  spread(lowerupper, value, convert = T) %>%
  full_join(., slr_tf_estimates, by = "field") 
#checking
Conf(x = master1$SLR_1_tf_yn, ref = master1$clinic_tf_yn) #also gives NPV/PPV which compare nicely
BinomCI(260, 260+41) #sens CI   est    lwr.ci    upr.ci
                           # 0.8637874 0.8204256 0.8979806
BinomCI(144, 155+43) #spec CI   est    lwr.ci    upr.ci
                           #  0.7272727 0.6613545 0.78454

#smartphone TF
bs_smartphone_tf_metrics <- map(bs2$splits, ~as.tibble(.) %>% unnest %>%
                           class_metrics(., truth = as.factor(clinic_tf_yn), estimate = as.factor(smartphone_1_tf_yn))) %>%
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=T)
#JN : note these CIs are named the same as those above. I did this to simplify the separating and merging process to get the final table
ci_sens <- bca(bs_smartphone_tf_metrics$sens)
ci_spec <- bca(bs_smartphone_tf_metrics$spec)
ci_npv <- bca(bs_smartphone_tf_metrics$npv)
ci_ppv <- bca(bs_smartphone_tf_metrics$ppv)
#calculating estimates 
smart_tf_estimates <- class_metrics(master1_filtered, truth = as.factor(clinic_tf_yn), estimate = as.factor(smartphone_1_tf_yn)) %>%
  dplyr::select(-.estimator) %>%
  mutate(field=.metric) %>%
  dplyr::select(-.metric)
#creating one data table with the metrics for SLR_tf
smart_tf_metrics <- data.frame(lowerupper, ci_sens, ci_spec, ci_npv, ci_ppv, stringsAsFactors = F) %>%
  gather(ci_field, value, ci_sens:ci_ppv) %>%
  separate(ci_field, into=c("ci", "field"), sep="_") %>%
  dplyr::select(-ci) %>%
  spread(lowerupper, value, convert = T) %>%
  full_join(., smart_tf_estimates, by = "field")
#checking
Conf(x = master1$smartphone_1_tf_yn, ref = master1$clinic_tf_yn) #more sensitive for TF than SLR
BinomCI(279, 279+22) #sens CI   est    lwr.ci    upr.ci
                            # 0.9269103 0.891821 0.9512402
BinomCI(127, 127+71) #spec CI   est    lwr.ci    upr.ci
                          # 0.6414141 0.572506 0.7049395

#SLR TI
bs_slr_ti_metrics <- map(bs2$splits, ~as.tibble(.) %>% unnest %>%
                           class_metrics(., truth = as.factor(clinic_ti_yn), estimate = as.factor(SLR_1_ti_yn))) %>%
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=T)
#now calculating bca confidence intervals, using the same workaround
ci_sens <- bca(bs_slr_ti_metrics$sens)
ci_spec <- bca(bs_slr_ti_metrics$spec)
ci_npv <- bca(bs_slr_ti_metrics$npv)
ci_ppv <- bca(bs_slr_ti_metrics$ppv)
#calculating estimates 
slr_ti_estimates <- class_metrics(master1_filtered, truth = as.factor(clinic_ti_yn), estimate = as.factor(SLR_1_ti_yn)) %>%
  dplyr::select(-.estimator) %>%
  mutate(field=.metric) %>%
  dplyr::select(-.metric)
#creating one data table with the metrics for SLR_tf
slr_ti_metrics <- data.frame(lowerupper, ci_sens, ci_spec, ci_npv, ci_ppv, stringsAsFactors = F) %>%
  gather(ci_field, value, ci_sens:ci_ppv) %>%
  separate(ci_field, into=c("ci", "field"), sep="_") %>%
  dplyr::select(-ci) %>%
  spread(lowerupper, value, convert = T) %>%
  full_join(., slr_ti_estimates, by = "field")
#checking
Conf(x = master1$SLR_1_ti_yn, ref = master1$clinic_ti_yn) 
BinomCI(417, 417+41) #sens CI   est    lwr.ci    upr.ci
                          #  0.9104803 0.8808059 0.9333263
BinomCI(31, 31+10) #spec CI   est    lwr.ci    upr.ci
                          # 0.7560976 0.6065666 0.86175

#Smartphone TI
bs_smartphone_ti_metrics <- map(bs2$splits, ~as.tibble(.) %>% unnest %>%
                                  class_metrics(., truth = as.factor(clinic_ti_yn), estimate = as.factor(smartphone_1_ti_yn))) %>%
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=T)
#JN : note these CIs are named the same as those above. I did this to simplify the separating and merging process to get the final table
ci_sens <- bca(bs_smartphone_ti_metrics$sens)
ci_spec <- bca(bs_smartphone_ti_metrics$spec)
ci_npv <- bca(bs_smartphone_ti_metrics$npv)
ci_ppv <- bca(bs_smartphone_ti_metrics$ppv)
#calculating estimates 
smart_ti_estimates <- class_metrics(master1_filtered, truth = as.factor(clinic_ti_yn), estimate = as.factor(smartphone_1_ti_yn)) %>%
  dplyr::select(-.estimator) %>%
  mutate(field=.metric) %>%
  dplyr::select(-.metric)
#creating one data table with the metrics for SLR_tf
smart_ti_metrics <- data.frame(lowerupper, ci_sens, ci_spec, ci_npv, ci_ppv, stringsAsFactors = F) %>%
  gather(ci_field, value, ci_sens:ci_ppv) %>%
  separate(ci_field, into=c("ci", "field"), sep="_") %>%
  dplyr::select(-ci) %>%
  spread(lowerupper, value, convert = T) %>%
  full_join(., smart_ti_estimates, by = "field")
#checking
Conf(x = master1$smartphone_1_ti_yn, ref = master1$clinic_ti_yn)
BinomCI(427, 427+31) #sens CI   est    lwr.ci    upr.ci
                          #  0.9323144 0.9055278 0.9519093
BinomCI(29, 29+12) #spec CI   est    lwr.ci    upr.ci
                          # 0.7073171 0.5552053 0.8239081

################################
##             LCA            ##
################################
#example of poLCA https://www.r-bloggers.com/example-8-21-latent-class-analysis/
# another example https://statistics.ohlsen-web.de/latent-class-analysis-polca/
library(poLCA)

#latent groups: +/- TF (will build separate model later for TI)
#Response variables: SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, newindvpcr -- if + all would indicate latent/true trachoma infection
#JN note: maybe incorporate covariates: state_code, sex, age. We may want to do a clustered analysis using state_code (mixed effects linear regression?), 
  #but not quite sure how to do this

#new df with relevant variables
LCAtf_df <- master1_filtered %>%
  dplyr::select(SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn) %>%  #consider adding in later: state_code, age, sex
  #and recoding as positive integers as required for poLCA
  mutate(SLR_tf_yn=case_when(SLR_1_tf_yn == 0 ~ as.integer(1),
                             SLR_1_tf_yn == 1 ~ as.integer(2)),
         smartphone_tf_yn=case_when(smartphone_1_tf_yn == 0 ~ as.integer(1),
                                    smartphone_1_tf_yn == 1 ~ as.integer(2)),
         clinic_tf=case_when(clinic_tf_yn == 0 ~ as.integer(1),
                             clinic_tf_yn == 1 ~ as.integer(2))) %>%
  #now 1=no and 2-yes
  dplyr::select(-SLR_1_tf_yn, -smartphone_1_tf_yn, -clinic_tf_yn)
#running the LCA, can change the iterations and reps for the final analysis
LCA_tf_model1 = poLCA(cbind(SLR_tf_yn, smartphone_tf_yn, clinic_tf) ~ 1, maxiter = 5000, nclass = 2, nrep = 10, data = LCAtf_df)
#JN : It seems class 1 would be +TF (0.374) and class 2 would be -TF (0.626). 
#reordering latent classes so 2 = +TF and 1 = -TF
probs.start <- LCA_tf_model1$probs.start
new.probs.start <- poLCA.reorder(probs.start, c(2,1)) #here class 1 (+TF) becomes 2 and class 2 (-TF) becomes 1
LCA_tf_model_reordered = poLCA(cbind(SLR_tf_yn, smartphone_tf_yn, clinic_tf) ~ 1, probs.start = new.probs.start, 
                               maxiter = 5000, nclass = 2, nrep = 10, data = LCAtf_df)
#Not sure how I would get the sens/spec/npv/ppv from the results...
#inspecting the latent class memberships:
table(LCA_tf_model1$predclass) #could calc the sens/spec etc from this but bootstrapping would be an issue
#JN: I would like to have a column of predicted class memberships I can bind to the df used in the LCA model
LCA_tf_model_reordered$predcell #this give predicted and observed for each combination of response variables
LCA_tf_model_reordered$posterior #class membership probabilities for each individual
LCA_tf_model_reordered$predclass  #here it is!!
skim(LCA_tf_model_reordered$predclass) #499, everyone accounted for

tf_with_latentclass <- LCAtf_df %>%
  cbind(LCA_tf_model_reordered$predclass) %>%
  rename(predclass = 'LCA_tf_model_reordered$predclass') %>%
  #reformating to binary
  mutate(SLR_tf_yn=SLR_tf_yn-1,
         smartphone_tf_yn=smartphone_tf_yn-1,
         clinic_tf=clinic_tf-1,
         predclass=predclass-1) #with this I should be able to calculate sens/spec with bca using the yardstick backage
#just trying it out now. will do for all and calculate bcas once we have final latent class
class_metrics(tf_with_latentclass, truth = as.factor(predclass), estimate = as.factor(SLR_tf_yn))
class_metrics(tf_with_latentclass, truth = as.factor(predclass), estimate = as.factor(smartphone_tf_yn))
class_metrics(tf_with_latentclass, truth = as.factor(predclass), estimate = as.factor(clinic_tf))
#JN : I remember Tom saying something about including each of the villages as separate inputs. Does this mean I should create a separate model 
  #with the SLR/smart/clinic grades for each village as separate response variables? 
  #This would probably help our degrees of freedom
  #I would have to merge a lot of tables to get predicted class for each individuals 
  #but maybe I could use the posterior class membership probabilities from the model to determine the class membership for the new data
  #https://datascience.stackexchange.com/questions/40632/how-do-i-use-the-model-generated-by-the-r-package-polca-to-classify-new-data-as

################################
##           KAPPAS           ##
################################
# JK: So the kappas can be done with the same dataset, which is nice for internal consistency:
# Using 4-level, but without weighting (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
agreement_matrix <- tibble(
  "1" = c(1, 0.9, 0.2, 0.1),
  "2" = c(0.9, 1, 0.5, 0.2),
  "3" = c(0.2, 0.5, 1, 0.9),
  "4" = c(0.1, 0.2, 0.9, 1))
#Repeats
# Blake
xtabs(data=master1, ~ SLR_1_tf_1+SLR_2_tf_1, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_1, master1$SLR_2_tf_1, conf.level=0.95)
# Using dichotomous variable
xtabs(data=master1, ~ SLR_1_tf_di_1 +SLR_2_tf_di_1 , addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_1 , master1$SLR_2_tf_di_1 , conf.level=0.95)
# John
# Using 4-level, but without weighting 
xtabs(data=master1, ~ SLR_1_tf_2+SLR_2_tf_2, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_2, master1$SLR_2_tf_2, conf.level=0.95)
# Using dichotomous variable
xtabs(data=master1, ~ SLR_1_tf_di_2 +SLR_2_tf_di_2 , addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_2 , master1$SLR_2_tf_di_2 , conf.level=0.95)
# Jeremy
# Using 4-level, but without weighting 
xtabs(data=master1, ~ SLR_1_tf_NA+SLR_2_tf_NA, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_NA, master1$SLR_2_tf_NA, conf.level=0.95)
# Using dichotomous variable
xtabs(data=master1, ~ SLR_1_tf_di_NA +SLR_2_tf_di_NA , addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_NA , master1$SLR_2_tf_di_NA , conf.level=0.95)
#TRYING TO REPRODUCE WHAT YOU DID BELOW:
#among slr photos john vs blake
xtabs(data=master1, ~ SLR_1_tf_di_1+SLR_1_tf_di_2, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_1, master1$SLR_1_tf_di_2, conf.level=0.95)
# among smart photos john vs blake
xtabs(data = master1, ~smartphone_1_tf_di_1+smartphone_1_tf_di_2, addNA=TRUE)
CohenKappa(master1$smartphone_1_tf_di_1, master1$smartphone_1_tf_di_2, conf.level=0.95)
# Blake SLR vs field
xtabs(data=master1, ~ SLR_1_tf_di_1+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_1, master1$clinic_tf_yn, conf.level=0.95)
# John SLR vs field
xtabs(data=master1, ~ SLR_1_tf_di_2+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_2, master1$clinic_tf_yn, conf.level=0.95)
# Blake smart vs field
xtabs(data=master1, ~smartphone_1_tf_di_1+clinic_tf_yn, addNA=TRUE)
CohenKappa(master1$smartphone_1_tf_di_1, master1$clinic_tf_yn, conf.level=0.95)
# John smart vs field
xtabs(data=master1, ~smartphone_1_tf_di_2+clinic_tf_yn, addNA=TRUE)
CohenKappa(master1$smartphone_1_tf_di_2, master1$clinic_tf_yn, conf.level=0.95)
# Consensus SLR vs field
xtabs(data=master1, ~ SLR_1_tf_yn+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# Consensus smart vs field
xtabs(data=master1, ~ smartphone_1_tf_yn+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$smartphone_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# Consensus SLR vs Consensus smartphone
xtabs(data=master1, ~ smartphone_1_tf_yn+SLR_1_tf_yn, addNA=TRUE) 
CohenKappa(master1$smartphone_1_tf_yn, master1$SLR_1_tf_yn, conf.level=0.95)
#for TI
#among slr photos john vs blake
xtabs(data=master1, ~ SLR_1_ti_di_1+SLR_1_ti_di_2, addNA=TRUE) 
CohenKappa(master1$SLR_1_ti_di_1, master1$SLR_1_ti_di_2, conf.level=0.95)
# among smart photos john vs blake
xtabs(data = master1, ~smartphone_1_ti_di_1+smartphone_1_ti_di_2, addNA=TRUE)
CohenKappa(master1$smartphone_1_ti_di_1, master1$smartphone_1_ti_di_2, conf.level=0.95)
# Consensus SLR vs field
xtabs(data=master1, ~ SLR_1_ti_yn+clinic_ti_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_ti_yn, master1$clinic_ti_yn, conf.level=0.95)
# Consensus smart vs field
xtabs(data=master1, ~ smartphone_1_ti_yn+clinic_ti_yn, addNA=TRUE) 
CohenKappa(master1$smartphone_1_ti_yn, master1$clinic_ti_yn, conf.level=0.95)

################################
##      KAPPA DIFFERENCE     ##
################################
#Note to write function to determine difference between kappas
#write function in R that calculates difference in ks then bootstrap the difference. check if it is 0
#JN : My thought here is to calculate kappas for bootstrapped samples then calculate the difference between the two (e.g., clinic and slr)
  #and determine the 95%CI from this distribution of differences. 

set.seed(2020) 
#JN note: can use D (nested by community) if we need to account for nesting in the kappa analysis
bs3 <- bootstraps(master1_filtered, times = 999) #object containing the seperate resamples

#Difference function
dummy_df <- as.tibble(list(x=5:8, y=1:4))
dummy_dif <- function(x,y) {
  x - y }
dummy_dif(5,3) #works
dummy_dif(dummy_df$x, dummy_df$y) #also works

#two functions that take intput from split and return kappas for smart-clinic and slr-clinic
  #starting for slr-clinic and smart-clinic tf, reiterate for other combinations if it works.
bootstrapkappa_smart <- function(splits) {
  x <- analysis(splits)
  CohenKappa(x$smartphone_1_tf_yn, x$clinic_tf_yn, conf.level=NA)
}
bootstrapkappa_slr <- function(splits) {
  x <- analysis(splits)
  CohenKappa(x$SLR_1_tf_yn, x$clinic_tf_yn, conf.level=NA)
}
#creating bootstrapped kappas from bs3 in order to calculate differences
bootstrap_kappa <- bs3 %>% 
  mutate(
    kappa_slr = map(splits, bootstrapkappa_slr), 
    kappa_smart = map(splits, bootstrapkappa_smart))
bootstrap_kappa
#now applying the dummy_dif function
bootstrap_kappa1 <- mapply(dummy_dif, bootstrap_kappa$kappa_smart, bootstrap_kappa$kappa_slr)

#alternatively I can build a function that calculates the difference in ks then bootstrap the difference
bootstrapkappa_dif <- function(splits) {
  x <- analysis(splits)
  CohenKappa(x$smartphone_1_tf_yn, x$clinic_tf_yn, conf.level=NA) - CohenKappa(x$SLR_1_tf_yn, x$clinic_tf_yn, conf.level=NA)
}
#now to bootstrap it
bootstrap_kappa2 <- bs3 %>%
  mutate(kappa_dif = map(splits, bootstrapkappa_dif))
bootstrap_kappa2$kappa_dif #the same as the results of mapply(dummy_dif). This way seems simpler and faster to me.

#I think the next step would be to take the 95% CI of this distribution and see if it contains 0 to determine if there is a difference
ci_kappa_dif1 <- bca(bootstrap_kappa1, conf.level = 0.95)
ci_kappa_dif2 <- bca(unlist(bootstrap_kappa2$kappa_dif), conf.level = 0.95)
#same results and both contain 0

################################
##             ICC            ##
################################
# JK: This looks fine. I like the Stata help file for icc to help me understand the options, 
# just google Stata help icc and then open up the pdf and read the "Introduction" section (note, not the html; there is more info in the pdf)
# Not sure what we do here. Most people will want to see kappas. Tom is probably the only person who would want to see an ICC. 
# Not sure it really matters
#JN notes from Stata help ICC: 
  #one way random effects = each target is rated by a different set of independent raters
  #two way random effects model=each target is rated by the same set of independent raters who are randonly drawn from the population of raters. 
    #the random effects in this model are target and rater and possibly their interaction (cannot determine without repeated measurements for each rater on each target)
  #two way mixed effects model = each target is rated by the same set of independent raters, thus rater is a fixed effect.
    #random effects are target and possible target-rater interaction
  #there is also consistency of agreement vs absolute agreement. I think we are more interested in absolute agreement since we have a binary outcome/no ranking
#JN: seems like we want a two-way random-effects model as we are using the the combined data from three judges to compare the groups (SLR/smart) and want to generalize this to
  #the population of raters (allow our study to be generalizable)
  #we would also want an average ICC as a team of raters (3) are used to rate the target in each group (SLR/smart)

  #blake repeats
  #blake slr rep
master1 %>%
  dplyr::select(SLR_1_tf_1, SLR_2_tf_1) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #blake smart rep
master1 %>%
  dplyr::select(smartphone_1_tf_1, smartphone_2_tf_1) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #john repeats
  #johnslrrep
master1 %>%
  dplyr::select(SLR_1_tf_2, SLR_2_tf_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #johnsmartrep 
master1 %>%
  dplyr::select(smartphone_1_tf_2, smartphone_2_tf_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #JK repeats
  #JK SLR reps
master1 %>%
  dplyr::select(SLR_1_tf_NA, SLR_2_tf_NA) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #JK smart rep
master1 %>%
  dplyr::select(smartphone_1_tf_NA, smartphone_2_tf_NA) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #among slr photos john vs blake
master1 %>%
  dplyr::select(SLR_1_tf_di_1, SLR_1_tf_di_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # among smart photos john vs blake
master1 %>%
  dplyr::select(smartphone_1_tf_di_1, smartphone_1_tf_di_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Blake SLR vs field
master1 %>%
  dplyr::select(SLR_1_tf_di_1, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # John SLR vs field
master1 %>%
  dplyr::select(SLR_1_tf_di_2, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Blake smart vs field
master1 %>%
  dplyr::select(smartphone_1_tf_di_1, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # John smart vs field
master1 %>%
  dplyr::select(smartphone_1_tf_di_2, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus SLR vs field
master1 %>%
  dplyr::select(SLR_1_tf_yn, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus smart vs field
master1 %>%
  dplyr::select(smartphone_1_tf_yn, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus SLR vs Consensus smartphone
master1 %>%
  dplyr::select(SLR_1_tf_yn, smartphone_1_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
#Did we want to calculate Kappas and ICCs compared to each field grader (there are 9 of them)? 
# JK: I don't think I would bother