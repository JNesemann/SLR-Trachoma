#Final SLR Data analysis

#Analyses of interest:
  #Kappas:
    #among slr photos john vs blake--JN: Done
    #among smartphone photos john vs blake--JN: Done
    #among slr repeat photos john vs john, blake vs blake (we can eventually combine in an ICC maybe)--JN: Done
    #among smartphone repeat photos john vs john, blake vs blake (we can eventually combine in an ICC maybe)--JN: Done
    #consensus/trump grade smartphone vs consensus SLR--JN: Done
    #JN: I also compared all of the above to field grades
  #Village-level
    #Prevalence of TF, TI, TF+/-TI per village.
      #SLR: JN Done
      #Smartphone  JN done
      #PCR  JN done
    #Then plot different combinations of prevalences JN done
    #Regression/correlation coefficient to assess correlation between prevalences
    #Specifically looking at 5% threshold for TF, 10% threshold for TF
  #Sensitivity/Specificity--JN Done
      #Reference standard: TF by field grade--JN Done
      #Reference standard: PCR--JN Done
      #Do separately for 2 index tests (smartphone, SLR)
#Remember also: let’s try 2 different formulas for consensus photo grades: one in which I trump you two--JN: done
          # and the other where it’s a consensus. I’d be curious how often these are different--JN: done

library(tidyverse)
library(skimr)
library(readr)
library(readxl)
library(DescTools)
library(irr)
library(purrr)

#importing
master_key_import <- read_excel("Masked Number Master Key (1).xlsx")
TIRET3_PCR_import <- read_excel("TIRET3_PCR_for John.xlsx")
SLR_import <- read_csv("SLRTrachoma_DATA_2019-11-27_1325.csv")
TIRET3_PhotoExamData13Villages <- read_excel("TIRET3-PhotoExamDatain13Villages.xlsx")

#cleaning master key
master_key <- master_key_import %>%
  select(number, camera, `masked number`, `repeat`) %>%
  rename(mask = `masked number`,
         repeat.instance=`repeat`) %>%
  mutate(repeat.instance=if_else(is.na(repeat.instance), 1, repeat.instance))

#cleaning and renaming TIRET 3 PCR
PCR <- TIRET3_PCR_import %>%
  select(Group, `Random#`, PoolID, Pool_PCR_Results, Individual_PCR_Results, `Exam(0/1)`, Examiner) %>%
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
  select(Clinicalexam, Photo, `Random#`, Sex, `Stateteam-code`, Tuber, Age, ID) %>%
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
#JN note--need to filter out smartphone and SlR images that were ungradeable/unavailable (quality = 3 or 4)
xtabs(data=master1, ~SLR_1_quality_1 + SLR_1_quality_2)
xtabs(data=master1, ~SLR_1_quality_1 + SLR_1_quality_2 + SLR_1_quality_NA) #3 JK ruled as not available (qualit = 4) none quality = 3. will exclude
xtabs(data=master1, ~smartphone_1_quality_1 + smartphone_1_quality_2) #1 q=4, 6 q=3, and 6 borderline q=2/3 cases
xtabs(data=master1, ~smartphone_1_quality_1 + smartphone_1_quality_2 + smartphone_1_quality_NA)
xtabs(data=master1, ~examiner, addNA = TRUE)

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
         #creating sum score for quality
         totalid=n(), 
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
  #JN: new binary quality variable
  mutate(qual_di=if_else(quality %in% c(1, 2),1,if_else(quality %in% c(3, 4),0,NA_real_))) %>%
  mutate(sumqual=sum(qual_di, na.rm = TRUE),
         #JN: new consensus score for quality. 0 = ungradeable or no image, 1 = gradeable
         qual_yn=case_when(totalid == 3 & sumqual >= 2 ~ 1,
                           totalid == 3 & sumqual < 2 ~ 0,
                           totalid == 2 & sumqual == 2 ~ 1,
                           totalid == 2 & sumqual == 0 ~ 0)) %>%
  select(-complete, -sumtf, -sumti, -sumqual, -totalid) %>%
  select(id, tf_yn, ti_yn, qual_yn, dde, tf_di, ti_di, qual_di, everything()) %>% 
  #creating a wide dataset
  gather(field, value, tf_di:notes) %>%
  mutate(field_dde = paste(field, dde, sep = "_")) %>%
  select(-field, -dde) %>%
  spread(field_dde, value, convert = TRUE) %>%
  # JK: Note that when you enter the dot below it means to use the current data, so you're merging into the current dataframe
  right_join(., master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  #JN: added qual_yn
  select(number, camera.instance, id, tf_yn, ti_yn, qual_yn, everything()) %>%
  select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_NA) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE) %>%
  left_join(., PCR_exam, by="number") %>%
  #JN filtering out missing and ungradeable images. xtabs(data = master1, ~SLR_1_qual_yn + smartphone_1_qual_yn) 7 smartphone and 6 SLR were ungradeable
  filter(SLR_1_qual_yn != 0 & smartphone_1_qual_yn != 0)

# JK: Now this is one big dataset. The variable names follow the pattern:
# SLR vs smartphone, then...
# whether it's a repeat (1 is the first instance; 2 is the second instance; so primary analyses should only be with 1)
# variable (tf /ti etc)
# grader (1=Blake, 2=John, NA=Jeremy)

#ARVO abstract cross tabulations
xtabs(data = master1, ~SLR_1_tf_yn + clinic_tf_yn + smartphone_1_tf_yn)
xtabs(data = master1, ~SLR_1_ti_yn + smartphone_1_ti_yn + clinic_ti_yn)

#Paper xtabs
xtabs(data=master1, ~state_code + newindpcr)

#summarizing demographic information
demoMF <- master1 %>%
  group_by(sex) %>%
  summarize(n=n())

set.seed(23) #Temp in lima 12/9/19
demonest <- master1 %>%
  nest(-sex) %>%
  mutate(mean_age = map(.x = data, ~mean(x = .x$age)),
         sd_age = map(.x = data, ~sd(x = .x$age)),
         #running bootstrapped mean and CIs 
         slr_tf = map(data, ~DescTools::BootCI(.$SLR_1_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         smart_tf = map(.x = data, ~BootCI(.x$smartphone_1_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         clinic_tf = map(data, ~DescTools::BootCI(.x$clinic_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         slr_ti = map(.x = data, ~BootCI(.x$SLR_1_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         smart_ti = map(.x = data, ~BootCI(.x$smartphone_1_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         clinic_ti = map(data, ~DescTools::BootCI(.x$clinic_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000))) %>%
  #mutating bootstrapped results into their own columns
  mutate(slr_tf_prev=map(slr_tf, ~.x[1]),
         slr_tf_low=map(slr_tf, ~.x[2]),
         slr_tf_up=map(slr_tf, ~.x[3]),
         smart_tf_prev=map(.x = smart_tf, ~.x[1]),
         smart_tf_low=map(.x = smart_tf, ~.x[2]),
         smart_tf_up=map(.x=smart_tf, ~.x[3]),
         clinic_tf_prev=map(.x = clinic_tf, ~.x[1]),
         clinic_tf_low=map(.x = clinic_tf, ~.x[2]),
         clinic_tf_up=map(.x=clinic_tf, ~.x[3]),
         slr_ti_prev=map(.x = slr_ti, ~.x[1]),
         slr_ti_low=map(.x = slr_ti, ~.x[2]),
         slr_ti_up=map(.x=slr_ti, ~.x[3]),
         smart_ti_prev=map(.x = smart_ti, ~.x[1]),
         smart_ti_low=map(.x = smart_ti, ~.x[2]),
         smart_ti_up=map(.x=smart_ti, ~.x[3]),
         clinic_ti_prev=map(.x = clinic_ti, ~.x[1]),
         clinic_ti_low=map(.x = clinic_ti, ~.x[2]),
         clinic_ti_up=map(.x=clinic_ti, ~.x[3])) %>%
  #removing unecessary columns
  select(-data, -(slr_tf:clinic_ti)) %>%
  unnest()

#JN: do chi-squared to determine differences

# JK: So the kappas can be done with the same dataset, which is nice for internal consistency:
# Using 4-level, but without weighting (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
agreement_matrix <- tibble(
  "1" = c(1, 0.9, 0.2, 0.1),
  "2" = c(0.9, 1, 0.5, 0.2),
  "3" = c(0.2, 0.5, 1, 0.9),
  "4" = c(0.1, 0.2, 0.9, 1))
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
#JN comment--doing for TI
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

#JN comment--there are 9 field graders, going to calculate a Fleiss' kappa for them
  #Fleiss' kappa specifically allows that although there are a fixed number of raters (e.g., 9), different items may be rated by different individuals
  #another option is Krippendorff's alpha
field_exam <- master1 %>%
  select(id, examiner, clinic_ti_yn, clinic_tf_yn) %>%
  group_by(examiner) %>%
  gather(field, value, examiner:clinic_tf_yn) %>%
  spread(field, value)

#JN: Doing an ICC. Best to do an ICC2 (two-way random-effects model), look at absolute agreement
  #extent to which different graders assign the same score to the same subject (could be wrong on this)
  #redoing the comparisons above
  #blake repeats
  #blake slr rep
master1 %>%
  select(SLR_1_tf_1, SLR_2_tf_1) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #blake smart rep
master1 %>%
  select(smartphone_1_tf_1, smartphone_2_tf_1) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #john repeats
  #johnslrrep
master1 %>%
  select(SLR_1_tf_2, SLR_2_tf_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #johnsmartrep 
master1 %>%
  select(smartphone_1_tf_2, smartphone_2_tf_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #JK repeats
  #JK SLR reps
master1 %>%
  select(SLR_1_tf_NA, SLR_2_tf_NA) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #JK smart rep
master1 %>%
  select(smartphone_1_tf_NA, smartphone_2_tf_NA) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  #among slr photos john vs blake
master1 %>%
  select(SLR_1_tf_di_1, SLR_1_tf_di_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # among smart photos john vs blake
master1 %>%
  select(smartphone_1_tf_di_1, smartphone_1_tf_di_2) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Blake SLR vs field
master1 %>%
  select(SLR_1_tf_di_1, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # John SLR vs field
master1 %>%
  select(SLR_1_tf_di_2, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Blake smart vs field
master1 %>%
  select(smartphone_1_tf_di_1, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # John smart vs field
master1 %>%
  select(smartphone_1_tf_di_2, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus SLR vs field
master1 %>%
  select(SLR_1_tf_yn, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus smart vs field
master1 %>%
  select(smartphone_1_tf_yn, clinic_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
  # Consensus SLR vs Consensus smartphone
master1 %>%
  select(SLR_1_tf_yn, smartphone_1_tf_yn) %>%
  icc(model = "twoway", type = "agreement", conf.level = 0.95)
#Did we want to calculate Kappas and ICCs compared to each field grader (there are 9 of them)? 

#now for village level prevalences and boostraped confidence intervals
#first I have to figure out how to nest all individuals from the same village together
library(boot)
set.seed(22) #Temp in Lima 12/6/19
dummy_nest <- master1 %>%
  select(number, id, SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn, 
         examiner, age, sex, state_code, newindpcr) %>%
  group_by(state_code) %>%
  nest()
#I found this way of booting ci's online so this is method 1
#making a boot function
boot_mean <- function(d, i) {
  mean(d[i])
}

dummy_boot1 <- dummy_nest %>%
  #doing just tf for now
  mutate(booted_slr_tf=map(.x=data,
                         ~boot::boot(data = .x$SLR_1_tf_yn,
                               statistic = boot_mean,
                               R = 10000,
                               stype = "i")),
         booted_slr_tf_ci=map(.x=booted_slr_tf,  #this is the list-column containing bootstraped samples from which I will derive my confidence intervals
                       ~boot::boot.ci(.x, 
                               conf = 0.95,
                               type = "basic")), #not sure if basic is the right type to use
         booted_smart_tf=map(.x = data, 
                            ~boot::boot(data = .x$smartphone_1_tf_yn,
                               statistic = boot_mean,
                               R = 10000,
                               stype = "i")),
         booted_smart_tf_ci=map(.x = booted_smart_tf,
                            ~boot::boot.ci(.x,
                                conf = 0.95,
                                type = "basic")),
         booted_clinic_tf=map(.x = data,
                            ~boot::boot(data = .x$clinic_tf_yn,
                                statistic = boot_mean,
                                R = 10000,
                                stype = "i")),
         booted_clinic_tf_ci=map(.x = booted_clinic_tf,
                             ~boot::boot.ci(.x,
                                conf = 0.95,
                                type = "basic"))) %>%
  #extracting the relative data from the mutated boot_ci
  mutate(slr_tf_prev = map(.x = booted_slr_tf_ci,
                          ~.x$t0),        #extracting the point estimate 
         slr_tf_lower_ci = map(.x = booted_slr_tf_ci,
                        ~ .x$basic[[4]]),  #extracting the lower 2.5% limit
         slr_tf_upper_ci = map(.x = booted_slr_tf_ci,
                        ~ .x$basic[[5]]),
         smart_tf_prev = map(.x = booted_smart_tf_ci,
                           ~.x$t0),        
         smart_tf_lower_ci = map(.x = booted_smart_tf_ci,
                               ~ .x$basic[[4]]), 
         smart_tf_upper_ci = map(.x = booted_smart_tf_ci,
                               ~ .x$basic[[5]]),
         clinic_tf_prev = map(.x = booted_clinic_tf_ci,
                           ~.x$t0),        
         clinic_tf_lower_ci = map(.x = booted_clinic_tf_ci,
                               ~ .x$basic[[4]]), 
         clinic_tf_upper_ci = map(.x = booted_clinic_tf_ci,
                               ~ .x$basic[[5]])) %>%
  #drop list columns as they are no longer needed
  select(-data, -booted_slr_tf, -booted_slr_tf_ci, -booted_smart_tf, -booted_smart_tf_ci, -booted_clinic_tf, -booted_clinic_tf_ci) %>%
  unnest()

  #method 2
dummy_boot2 <- dummy_nest %>%
  #not quite sure why but the only way I was able to get it to work was specify the fucntion ~DescTools::BootCI
  mutate(slr_tf_boot = map(data, ~DescTools::BootCI(.x$SLR_1_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         smart_tf_boot = map(data, ~DescTools::BootCI(.x$smartphone_1_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         clinic_tf_boot = map(data, ~DescTools::BootCI(.x$clinic_tf_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         slr_ti_boot = map(data, ~DescTools::BootCI(.x$SLR_1_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         smart_ti_boot = map(data, ~DescTools::BootCI(.x$smartphone_1_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         clinic_ti_boot = map(data, ~DescTools::BootCI(.x$clinic_ti_yn, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000))) %>%
  #extracting the mean [1], lower ci [2], and upper ci [3] from the boot strapped TF means
  mutate(slr_tf_prev=map(.x = slr_tf_boot, ~.x[1]),
         slr_tf_lowci=map(.x = slr_tf_boot, ~.x[2]),
         slr_tf_upci=map(.x = slr_tf_boot, ~.x[3]),
         smart_tf_prev=map(.x = smart_tf_boot, ~.x[1]),
         smart_tf_lowci=map(.x = smart_tf_boot, ~.x[2]),
         smart_tf_upci=map(.x = smart_tf_boot, ~.x[3]),
         clinic_tf_prev=map(.x = clinic_tf_boot, ~.x[1]),
         clinic_tf_lowci=map(.x = clinic_tf_boot, ~.x[2]),
         clinic_tf_upci=map(.x = clinic_tf_boot, ~.x[3]),
         #same thing for the bootstrapped ti
         slr_ti_prev=map(.x = slr_ti_boot, ~.x[1]),
         slr_ti_lowci=map(.x = slr_ti_boot, ~.x[2]),
         slr_ti_upci=map(.x = slr_ti_boot, ~.x[3]),
         smart_ti_prev=map(.x = smart_ti_boot, ~.x[1]),
         smart_ti_lowci=map(.x = smart_ti_boot, ~.x[2]),
         smart_ti_upci=map(.x = smart_ti_boot, ~.x[3]),
         clinic_ti_prev=map(.x = clinic_ti_boot, ~.x[1]),
         clinic_ti_lowci=map(.x = clinic_ti_boot, ~.x[2]),
         clinic_ti_upci=map(.x = clinic_ti_boot, ~.x[3])) %>%
  #removing unecessary columns
  select(-(data:clinic_ti_boot)) %>%
  unnest(cols = c(slr_tf_prev, slr_tf_lowci, slr_tf_upci, smart_tf_prev, smart_tf_lowci, 
                  smart_tf_upci, clinic_tf_prev, clinic_tf_lowci, clinic_tf_upci, 
                  slr_ti_prev, slr_ti_lowci, slr_ti_upci, smart_ti_prev, smart_ti_lowci, 
                  smart_ti_upci, clinic_ti_prev, clinic_ti_lowci, clinic_ti_upci))
glimpse(dummy_boot2) #it works!

cor.test(dummy_boot2$slr_tf_prev, dummy_boot2$smart_tf_prev, alternative = "two.sided", method = "spearman", conf.level = 0.95)
    
#restructuring the prevalences to get the graph I want
village_prev <- dummy_boot2 %>%
  select(state_code:clinic_ti_upci) %>%
  transmute(slr_tf_prev=slr_tf_prev, slr_tf_lowci=slr_tf_lowci, slr_tf_upci=slr_tf_upci,
            smart_tf_prev=smart_tf_prev, smart_tf_lowci=smart_tf_lowci, smart_tf_upci=smart_tf_upci,
            clinic_tf_prev=clinic_tf_prev, clinic_tf_lowci=clinic_tf_lowci, clinic_tf_upci=clinic_tf_upci,
            slr_ti_prev=slr_ti_prev, slr_ti_lowci=slr_ti_lowci, slr_ti_upci=slr_ti_upci,
            smart_ti_prev=smart_ti_prev, smart_ti_lowci=smart_ti_lowci, smart_ti_upci=smart_ti_upci,
            clinic_ti_prev=clinic_ti_prev, clinic_ti_lowci=clinic_ti_lowci, clinic_ti_upci=clinic_ti_upci) %>%
  gather(field, value, slr_tf_prev:clinic_ti_upci) %>%
  separate(field, into = c("method", "sign", "measure"), sep = "_", convert = TRUE) %>%
  spread(measure, value)

  #Graph showing village prev + 95%CIs for TF and TI seperately
  ggplot(data = village_prev, mapping = aes(x=state_code, y=prev, fill=method)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin = lowci, ymax = upci), width=0.2, position = position_dodge(0.9)) +
  facet_wrap(~sign, nrow = 1) +
  coord_flip()

xtabs(data = master1, ~state_code + sex)

#making separate data table in order to graph overall prevalence estimates
overall_prev <- master1 %>%
  select(id, SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) %>%
  transmute(id=id, SLR_tf=SLR_1_tf_yn, smartphone_tf=smartphone_1_tf_yn, clinic_tf=clinic_tf_yn,
            SLR_ti=SLR_1_ti_yn, smartphone_ti=smartphone_1_ti_yn, clinic_ti=clinic_ti_yn) %>%
  gather(field, value, SLR_tf:clinic_ti) %>%
  separate(field, into = c("method", "sign"), sep = "_", convert = TRUE) %>%
  spread(method, value) %>%
  group_by(sign) %>%
  nest() %>%
  mutate(slr_boot=map(data, ~DescTools::BootCI(.x$SLR, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         smart_boot=map(data, ~DescTools::BootCI(.x$smartphone, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000)),
         clinic_boot=map(data, ~DescTools::BootCI(.x$clinic, FUN=mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 1000))) %>%
  mutate(slr_prev=map(.x=slr_boot, ~.x[1]),
         slr_low=map(.x=slr_boot, ~.x[2]),
         slr_up=map(.x=slr_boot, ~.x[3]),
         smart_prev=map(.x=smart_boot, ~.x[1]),
         smart_low=map(.x=smart_boot, ~.x[2]),
         smart_up=map(.x=smart_boot, ~.x[3]),
         clinic_prev=map(.x=clinic_boot, ~.x[1]),
         clinic_low=map(.x=clinic_boot, ~.x[2]),
         clinic_up=map(.x=clinic_boot, ~.x[3])) %>%
  select(-(data:clinic_boot)) %>%
  unnest() %>%
  gather(field, value, slr_prev:clinic_up) %>%
  separate(field, into = c("method", "measure"), sep = "_") %>%
  spread(measure, value)

  ggplot(data = overall_prev, mapping = aes(x = sign, y = prev, fill = method)) +
    geom_bar(position="dodge", stat = "identity") +
    geom_errorbar(aes(ymin = low, ymax = up), width=0.2, position = position_dodge(0.9))
    
#Sensitivity/Specificity, Do separately for 2 index tests (smartphone, SLR), alternatively we could do an LCA per TL--JN Tried below but was unable to
#unsure how to calculate bootstrapped confidence intervals for these so I just used BinomCI
#Reference standard: TF by field grade
#SLR
Conf(x = master1$SLR_1_tf_yn, ref = master1$clinic_tf_yn)
BinomCI(260, 260+41) #sens CI   est    lwr.ci    upr.ci
                          # 0.8637874 0.8204256 0.8979806
BinomCI(144, 155+43) #spec CI   est    lwr.ci    upr.ci
                          # 0.7272727 0.6613545 0.78454
                      #Pos Pred Value : 0.8581
                      # Neg Pred Value : 0.7908
#smartphone
Conf(x = master1$smartphone_1_tf_yn, ref = master1$clinic_tf_yn) #more sensitive for TF than SLR
BinomCI(279, 279+22) #sens CI   est    lwr.ci    upr.ci
                          # 0.9269103 0.891821 0.9512402
BinomCI(127, 127+71) #spec CI   est    lwr.ci    upr.ci
                          # 0.6414141 0.572506 0.7049395
                      #Pos Pred Value : 0.7971
                      #Neg Pred Value : 0.8523
#Reference standard: TI by field grade
#SLR
Conf(x = master1$SLR_1_ti_yn, ref = master1$clinic_ti_yn)
BinomCI(417, 417+41) #sens CI   est    lwr.ci    upr.ci
                          #  0.9104803 0.8808059 0.9333263
BinomCI(31, 31+10) #spec CI   est    lwr.ci    upr.ci
                          # 0.7560976 0.6065666 0.86175
                    #Pos Pred Value : 0.9766
                    #Neg Pred Value : 0.4306
#smartphone
Conf(x = master1$smartphone_1_ti_yn, ref = master1$clinic_ti_yn) #slightly more sensitive than SLR
BinomCI(427, 427+31) #sens CI   est    lwr.ci    upr.ci
                          #  0.9323144 0.9055278 0.9519093
BinomCI(29, 29+12) #spec CI   est    lwr.ci    upr.ci
                          # 0.7073171 0.5552053 0.8239081
                    #Pos Pred Value : 0.9727
                    #Neg Pred Value : 0.4833

#Attempting an LCA
#install.packages("poLCA")
library(poLCA)

#latent groups: +/- TF (will build separate model later for TI)
#Response variables: SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, newindvpcr -- if + all would indicate latent/true trachoma infection
#covariates: state_code, sex, age. We may want to do a clustered analysis using state_code (mixed effects linear regression?), but not quite sure how to do this
#new df with relevant variables
LCAtf <- master1 %>%
  dplyr::select(SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn) %>%  #consider adding in later: state_code, age, sex
  #and recoding as positive integers as required for poLCA
  mutate(SLR_tf_yn=case_when(SLR_1_tf_yn == 0 ~ as.integer(1),
                               SLR_1_tf_yn == 1 ~ as.integer(2)),
         smartphone_tf_yn=case_when(smartphone_1_tf_yn == 0 ~ as.integer(1),
                                    smartphone_1_tf_yn == 1 ~ as.integer(2)),
         clinic_tf=case_when(clinic_tf_yn == 0 ~ as.integer(1),
                             clinic_tf_yn == 1 ~ as.integer(2))) %>%
  dplyr::select(-SLR_1_tf_yn, -smartphone_1_tf_yn, -clinic_tf_yn)
#building a function reflecting this
f <- cbind(LCAtf$SLR_tf_yn, LCAtf$smartphone_tf_yn, LCAtf$clinic_tf) ~ 1
#model 1
M1 <- poLCA(f, LCAtf, nclass = 2) #still not working for me

Cor(x = master1$SLR_1_tf_yn, y=master1$smartphone_1_tf_yn, method = c("pearson"), use = "complete.obs")

#JN STOPPED HERE--OLD CODE
#Regression/correlation coefficient to assess correlation between prevalences
#I think we talked about taking into account how field grades will be clustered by examiner and village, so we could use a mixed effects linear regression
#but I am not sure this is what you meant... https://m-clark.github.io/mixed-models-with-R/random_intercepts.html

    #Correlations & regressions
ggplot(data = village_prevalence) +
  geom_point(mapping = aes(x = clinic_tf, y = pcr)) 
qqplot(x = village_prevalence$clinic_tf, y = village_prevalence$pcr) #I repeated the plot for each of the combos listed below. 
      #SLR tf cons vs clinic. Roughly normal
      Cor(x = village_prevalence$clinic_tf, y = village_prevalence$slr_tf_cons, method = c("pearson"), use = "complete.obs") #0.75, comparing prevalences
      Cor(x = photo_pcr_master$clinic_tf_yn, y = photo_pcr_master$SLR_tf_yn_cons, method = c("pearson"), use = "complete.obs") #0.67, comparing grades used to calculate prevalences
      lm1 = lm(clinic_tf ~ slr_tf_cons, data = village_prevalence) %>%
        summary()  #0.50, comparing prevalences
      lm(clinic_tf_yn~SLR_tf_yn_cons, data = photo_pcr_master) %>%
        summary()  #0.67, comparing grades used to calculate prevalences
      ggplot(data=village_prevalence, mapping = aes(x=clinic_tf, y = slr_tf_cons)) +
        geom_point(color = 'red') +
        geom_smooth(method = 'lm', se = FALSE) 
      #SLR tf trump vs clinic. Roughly normal
      Cor(x = village_prevalence$clinic_tf, y = village_prevalence$slr_tf_trump, method = c("pearson"), use = "complete.obs")
      lm(clinic_tf ~ slr_tf_trump, data = village_prevalence) %>%
        summary() 
      #smart tf cons vs clinic.
      Cor(x = village_prevalence$clinic_tf, y = village_prevalence$smart_tf_cons, method = c("pearson"), use = "complete.obs")
      lm(clinic_tf ~ smart_tf_cons, data = village_prevalence) %>%
        summary()
      #smart tf trump vs clinic. 
      Cor(x = village_prevalence$clinic_tf, y = village_prevalence$smart_tf_trump, method = c("pearson"), use = "complete.obs")
      lm(clinic_tf ~ smart_tf_trump, data = village_prevalence) %>%
        summary()
      #slr tf cons vs smart tf cons
      Cor(x = village_prevalence$slr_tf_cons, y = village_prevalence$smart_tf_cons, method = c("pearson"), use = "complete.obs")
      lm(slr_tf_cons ~ smart_tf_cons, data = village_prevalence) %>%
        summary()
      #slr tf trump vs smart tf trump
      Cor(x = village_prevalence$slr_tf_trump, y = village_prevalence$smart_tf_trump, method = c("pearson"), use = "complete.obs")
      lm(slr_tf_trump ~ smart_tf_trump, data = village_prevalence) %>%
        summary()
      #slr tf cons vs pcr. does not appear to be linear
      Cor(x = village_prevalence$slr_tf_cons, y = village_prevalence$pcr, method = c("spearman"), use = "complete.obs")
      lm(slr_tf_cons ~ pcr, data = village_prevalence) %>%
        summary()
      #slr tf trump vs pcr. does not appear to be linear  
      Cor(x = village_prevalence$slr_tf_trump, y = village_prevalence$pcr, method = c("spearman"), use = "complete.obs")
      lm(slr_tf_trump ~ pcr, data = village_prevalence) %>%
        summary()
      #smart tf cons vs pcr
      Cor(x = village_prevalence$smart_tf_cons, y = village_prevalence$pcr, method = c("spearman"), use = "complete.obs")  #smartphone seems to correlate poorly with pcr
      lm(smart_tf_cons ~ pcr, data = village_prevalence) %>%
        summary()
      #smart tf trump vs pcr
      Cor(x = village_prevalence$smart_tf_trump, y = village_prevalence$pcr, method = c("spearman"), use = "complete.obs")
      lm(smart_tf_trump ~ pcr, data = village_prevalence) %>%
        summary()
      #clinic tf vs pcr
      Cor(x = village_prevalence$clinic_tf, y = village_prevalence$pcr, method = c("spearman"), use = "complete.obs")
      lm(clinic_tf ~ pcr, data = village_prevalence) %>%
        summary()