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

#importing
master_key_import <- read_excel("Masked Number Master Key (1).xlsx")
TIRET3_PCR_import <- read_excel("TIRET3_PCR_for John.xlsx")
SLR_import <- read_csv("SLRTrachoma_DATA_2019-11-27_1325.csv")
TIRET3_PhotoExamData13Villages <- read_excel("TIRET3-PhotoExamDatain13Villages.xlsx")
skim(TIRET3_PhotoExamData13Villages)  #this dataset has 13 unique stateteam-codes, corresponds to the 13 villages

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
skim(TIRET3_villages)    #JN: even when using this data set to calculate clinical grades we are still missing grades for 80 individuals...

#merging PCR data with complete data in TIRET3_PhotoExamData13Villages, this will provide the stateteam-code which I am using as a proxy for village ID
PCR_exam <- full_join(PCR, TIRET3_villages, by = "number")
skim(PCR)
skim(TIRET3_villages)
skim(PCR_exam)

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
  select(-complete, -sumtf, -sumti, -totalid) %>%
  select(id, tf_yn, ti_yn, dde, tf_di, ti_di, everything()) %>% 
  #creating a wide dataset
  gather(field, value, tf_di:notes) %>%
  mutate(field_dde = paste(field, dde, sep = "_")) %>%
  select(-field, -dde) %>%
  spread(field_dde, value, convert = TRUE) %>%
  # JK: Note that when you enter the dot below it means to use the current data, so you're merging into the current dataframe
  right_join(., master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  select(number, camera.instance, id, tf_yn, ti_yn, everything()) %>%
  select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_NA) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE) %>%
  left_join(., PCR_exam, by="number")
# JK: Now this is one big dataset. The variable names follow the pattern:
# SLR vs smartphone, then...
# whether it's a repeat (1 is the first instance; 2 is the second instance; so primary analyses should only be with 1)
# variable (tf /ti etc)
# grader (1=Blake, 2=John, NA=Jeremy)

#ARVO abstract cross tabulations
xtabs(data = master1, ~SLR_1_tf_yn + clinic_tf_yn + smartphone_1_tf_yn)
xtabs(data = master1, ~SLR_1_ti_yn + smartphone_1_ti_yn + clinic_ti_yn)

# JK: So the kappas can be done with the same dataset, which is nice for internal consistency:
# Using 4-level, but without weighting (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
    #JN comment--make agreement matrix for 1-4 grading
    #JN comment--?nest by grader to account for grades clustering by grader (included examiner from TIRET data sets)
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
#JN comment--do we want to do TI?
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

#JN: Doing an ICC. Best to do an ICC2 (two-way random-effects model), look at absolute agreement
  #extent to which different graders assign the same score to the same subject (could be wrong on this)
library(irr)
  #redoing the CohenK comparisons above
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
#Did we want to calculate Kappas and ICCs compared to each field grader? 

#now for village level prevalences and boostraped confidence intervals
#first I have to figure out how to nest all individuals from the same village together
library(boot)
set.seed(22) #setting seed to make results replicable
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
                         ~boot(data = .x$SLR_1_tf_yn,
                               statistic = boot_mean,
                               R = 10000,
                               stype = "i")),
         booted_slr_tf_ci=map(.x=booted_slr_tf,  #this is the list-column containing bootstraped samples from which I will derive my confidence intervals
                       ~ boot.ci(.x, 
                               conf = 0.95,
                               type = "basic")), #not sure if basic is the right type to use but it sounds right
         booted_smart_tf=map(.x = data, 
                            ~boot(data = .x$smartphone_1_tf_yn,
                               statistic = boot_mean,
                               R = 10000,
                               stype = "i")),
         booted_smart_tf_ci=map(.x = booted_smart_tf,
                            ~ boot.ci(.x,
                                conf = 0.95,
                                type = "basic")),
         booted_clinic_tf=map(.x = data,
                            ~boot(data = .x$clinic_tf_yn,
                                statistic = boot_mean,
                                R = 10000,
                                stype = "i")),
         booted_clinic_tf_ci=map(.x = booted_clinic_tf,
                             ~ boot.ci(.x,
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
    
#restructuring the prevalences to get the map I want
dummy_map <- dummy_boot2 %>%
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

#Then plot different combinations of prevalences
  #Graph showing village prev + 95%CIs for TF and TI seperately
ggplot(data = dummy_map, mapping = aes(x=state_code, y=prev, fill=method)) +
  geom_bar(position="dodge", stat = "identity") +
  geom_errorbar(aes(ymin = lowci, ymax = upci), width=0.2, position = position_dodge(0.9)) +
  facet_wrap(~sign, nrow = 1) +
  coord_flip()

#Sensitivity/Specificity, Do separately for 2 index tests (smartphone, SLR), alternatively we could do an LCA per TL
#unsure how to calculate bootstrapped confidence intervals for these so I just used BinomCI
#Reference standard: TF by field grade
#SLR
Conf(x = master1$SLR_1_tf_yn, ref = master1$clinic_tf_yn)
BinomCI(260, 260+41) #sens CI   est    lwr.ci    upr.ci
                          # 0.8637874 0.8204256 0.8979806
BinomCI(144, 155+43) #spec CI   est    lwr.ci    upr.ci
                          # 0.7272727 0.6613545 0.78454
#smartphone
Conf(x = master1$smartphone_1_tf_yn, ref = master1$clinic_tf_yn) #more sensitive for TF than SLR
BinomCI(279, 279+22) #sens CI   est    lwr.ci    upr.ci
                          # 0.9269103 0.891821 0.9512402
BinomCI(127, 127+71) #spec CI   est    lwr.ci    upr.ci
                          # 0.6414141 0.572506 0.7049395
#Reference standard: TI by field grade
#SLR
Conf(x = master1$SLR_1_ti_yn, ref = master1$clinic_ti_yn)
BinomCI(417, 417+41) #sens CI   est    lwr.ci    upr.ci
                          #  0.9104803 0.8808059 0.9333263
BinomCI(31, 31+10) #spec CI   est    lwr.ci    upr.ci
                          # 0.7560976 0.6065666 0.86175
#smartphone
Conf(x = master1$smartphone_1_ti_yn, ref = master1$clinic_ti_yn) #slightly more sensitive than SLR
BinomCI(427, 427+31) #sens CI   est    lwr.ci    upr.ci
                          #  0.9323144 0.9055278 0.9519093
BinomCI(29, 29+12) #spec CI   est    lwr.ci    upr.ci
                          # 0.7073171 0.5552053 0.8239081

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
M1 <- poLCA(f, LCAtf, nclass = 2)




#JN STOPPED HERE
#Regression/correlation coefficient to assess correlation between prevalences
#I think we talked about taking into account how field grades will be clustered by examiner and village, so we could use a mixed effects linear regression
#but I am not sure this is what you meant... https://m-clark.github.io/mixed-models-with-R/random_intercepts.html
gradernest <- master1 %>%
  select(id, SLR_1_tf_1:SLR_1_ti_yn, SLR_2_tf_1:SLR_2_ti_yn, smartphone_1_tf_1:smartphone_1_ti_yn, smartphone_2_tf_1:smartphone_2_ti_yn, newindpcr, examiner) %>%
  nest(-examiner)
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
      