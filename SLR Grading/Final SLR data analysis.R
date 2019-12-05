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
#per TL can consider an ICC? JN will redo once final image grades are in
# Using 4-level, but without weighting 
    #JN comment--make agreement matrix for 1-4 grading
    #JN comment--?nest by grader to account for grades clustering by grader (included examiner from TIRET data sets)
# (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
# Blake
xtabs(data=master1, ~ SLR_1_tf_1+SLR_2_tf_1, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_1, master1$SLR_2_tf_1, conf.level=0.95)
# Using dichotomous variable
xtabs(data=master1, ~ SLR_1_tf_di_1 +SLR_2_tf_di_1 , addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_1 , master1$SLR_2_tf_di_1 , conf.level=0.95)
# John
# Using 4-level, but without weighting 
# (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
xtabs(data=master1, ~ SLR_1_tf_2+SLR_2_tf_2, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_2, master1$SLR_2_tf_2, conf.level=0.95)
# Using dichotomous variable
xtabs(data=master1, ~ SLR_1_tf_di_2 +SLR_2_tf_di_2 , addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_2 , master1$SLR_2_tf_di_2 , conf.level=0.95)
# Jeremy
# Using 4-level, but without weighting 
# (in reality we'd probably want to weight 1 and 2 as being more similar, and 3 and 4 more similar)
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
# Blake SLR vs field
xtabs(data=master1, ~ SLR_1_tf_di_1+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_1, master1$clinic_tf_yn, conf.level=0.95)
# John SLR vs field
xtabs(data=master1, ~ SLR_1_tf_di_2+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_di_2, master1$clinic_tf_yn, conf.level=0.95)
# Blake smart vs field
# John smart vs field
# Consensus SLR vs field
xtabs(data=master1, ~ SLR_1_tf_yn+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$SLR_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# Consensus smart vs field
xtabs(data=master1, ~ smartphone_1_tf_yn+clinic_tf_yn, addNA=TRUE) 
CohenKappa(master1$smartphone_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# Consensus SLR vs Consensus smartphone
xtabs(data=master1, ~ smartphone_1_tf_yn+SLR_1_tf_yn, addNA=TRUE) 
CohenKappa(master1$smartphone_1_tf_yn, master1$SLR_1_tf_yn, conf.level=0.95)


#now for village level prevalences and boostraped confidence intervals
#first I have to figure out how to nest all individuals from the same village together
library(boot)
set.seed(22) #setting seed to make results replicable

#I found this way of booting ci's online so this is method 1
#making a boot function
boot_mean <- function(d, i) {
  mean(d[i])
}

dummy_boot1 <- master1 %>%
  #selecting out unused columns--can undo this later
  select(number, id, SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn, 
         examiner, age, sex, state_code, newindpcr) %>%
  group_by(state_code) %>%
  nest()  %>%
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

  #method 2--START HERE TOMORROW
dummy_boot2

dummy_boot <- dummy_nest %>%
  mutate(booted_prev=map(.x=data,
                         ~boot(data = .x$SLR_1_tf_yn,
                               statistic = boot_mean,
                               R = 10000,
                               stype = "i")),
         booted_ci=map(.x=booted_prev,  #this is the list-column containing bootstraped samples from which I will derive my confidence intervals
                       ~ boot.ci(.x, 
                                 conf = 0.95,
                                 type = "basic"))) %>%  #not sure if basic is the right type to use but it sounds right
  mutate(prevalence = map(.x = booted_ci,
                         ~.x$t0),        #extracting the point estimate 
         lower_ci = map(.x = booted_ci,
                        ~ .x$basic[[4]]),  #extracting the lower 2.5% limit
         upper_ci = map(.x=booted_ci,
                        ~ .x$basic[[5]])) %>%
  #drop list columns as they are no longer needed
  select(-data, -booted_prev, -booted_ci) %>%
  unnest()
         
#extracting and tidying the data
str(dummy_boot$booted_ci[[1]])

#inspecting the results
boot.plots <- map(.x = dummy_boot$booted.prev,
                  ~ plot(.x))
prints <- map(.x = dummy_boot$booted_ci,
              ~ print(.x))


  mutate(slr_tf_prev.CI = map(,mean(SLR_1_ti_yn)))

  
dummy_map <- dummy_nest %>%
  mutate(mean = map(data, "SLR_1_tf_yn") %>% map_dbl(BootCI(mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 10000)))

map(dummy_nest$results, ~.x mean(SLR_1_tf_yn))

BootCI(SLR_1_tf_yn, mean, bci.method = "norm", conf.level = 0.95, sides = "two.sided", R = 10000)

skim(dummy_nest)
glimpse(dummy_boot)
glimpse(dummy_boot$booted.prev)
glimpse(dummy_boot$data[[1]]) #by doing this it seems each row (village) in dummy_prev has a associated list of all the variables and their values in the data column
#SLR, smartphone, clinical, and PCR prevalences per village
village_prev1 <- master1 %>%
  group_by(state_code) %>%
    summarise(slr_tf_prev=mean(SLR_1_tf_yn, na.rm = TRUE),
            smart_tf_prev=mean(smartphone_1_tf_yn, na.rm = TRUE),
            clinic_tf_prev=mean(clinic_tf_yn, na.rm = TRUE),
            slr_ti_prev=mean(SLR_1_tf_yn, na.rm = TRUE),
            smart_ti_prev=mean(smartphone_1_tf_yn, na.rm = TRUE),
            clinic_ti_prev=mean(clinic_ti_yn, na.rm = TRUE),
            pcr_prev=mean(newindpcr, na.rm = TRUE))

#alternatively aclculating village level exam prevalence--might be useful to compare to that calculated by the sample of 499 kids

#making overall prev
overall_prev <- xyz %>%
  summarise(slr_tf_prev=mean(SLR_1_tf_yn, na.rm = TRUE),
            smart_tf_prev=mean(smartphone_1_tf_yn, na.rm = TRUE),
            clinic_tf_prev=mean(clinic_tf_yn, na.rm = TRUE),
            slr_ti_prev=mean(SLR_1_ti_yn, na.rm = TRUE),
            smart_ti_prev=mean(smartphone_1_ti_yn, na.rm = TRUE),
            clinic_ti_prev=mean(clinic_ti_yn, na.rm = TRUE),
            pcr_prev=mean(newindpcr, na.rm = TRUE))
#restructuring for overall prevalence graph
dummy_overall <- overall_prev %>%
  select(-pcr_prev) %>%
  transmute(slr.tf_prev=slr_tf_prev,
         smart.tf_prev=smart_tf_prev,
         clinic.tf_prev=clinic_tf_prev,
         slr.ti_prev=slr_ti_prev,
         smart.ti_prev=smart_ti_prev,
         clinic.ti_prev=clinic_ti_prev) %>%
  gather(field, prev, slr.tf_prev:clinic.ti_prev) %>%
  separate(field, into =c("method", "sign")) %>%
  cbind(lowerCI, upperCI)
  
  ggplot(dummy_overall, aes(x=method, y=prev)) +
    geom_bar(stat= "identity") + 
    facet_wrap(sign, nrow = 1)
  
    geom_errorbar(mapping = aes(x = sign, ymin = lowerCI, ymax = upperCI, position_dodge()))

#Then plot different combinations of prevalences
#restructuring the tibble so I can make the graphs I want
dummy_village <- village_prev %>%
  select(-pcr_prev) %>%  #I chose not to include pcr as only three villages had non-zero prevalences
  mutate(slr.cons_tf=slr_tf_cons,
         slr.trump_tf=slr_tf_trump,
         smart.cons_tf=smart_tf_cons,
         smart.trump_tf=smart_tf_trump,
         slr.cons_ti=slr_ti_cons,
         slr.trump_ti=slr_ti_trump,
         smart.cons_ti=smart_ti_cons,
         smart.trump_ti=smart_ti_trump) %>%
  select(-(slr_tf_cons:smart_tf_trump), -(slr_ti_cons:smart_ti_trump)) %>%
  gather(field, prevalence, clinic_tf:smart.trump_ti) %>%
  separate(field, into=c("method", "sign"), sep = "_", remove = TRUE, convert = TRUE)

#overall TF and TI prevalence estimated by different methods
ggplot(data = dummy_village) +  
  geom_point(mapping = aes(x = method, y = prevalence, color = sign)) +
  facet_wrap(~state_code, nrow=2) +
  coord_flip()
#showing TF prevalence estimated via different methods for each village
dummy_village %>%  
  filter(sign == "tf") %>%
  ggplot() +
  geom_bar(mapping = aes(x = method, y = prevalence), stat = "identity") +
  facet_wrap(~state_code, nrow=2) +
  coord_flip()
#showing TI prevalence estimated via different methods for each village
dummy_village %>%  
  filter(sign == "ti") %>%
  ggplot() +
  geom_bar(mapping = aes(x = method, y = prevalence), stat = "identity") +
  facet_wrap(~state_code, nrow=2) +
  coord_flip()
#overall TF prevalence by different methods
dummy_village %>%  
  filter(sign == "tf") %>%
  ggplot() +
  geom_bar(mapping = aes(x = method, y = prevalence), stat = "identity") #doesn't seem to be a big difference between trump and consensus
#overall TI prevalence by different methods
dummy_village %>%  
  filter(sign == "tf") %>%
  ggplot() +
  geom_bar(mapping = aes(x = method, y = prevalence), stat = "identity") #doesn't seem to be a big difference between trump and consensus

#Regression/correlation coefficient to assess correlation between prevalences
  #Specifically looking at 5% threshold for TF, 10% threshold for TF. JN comment: not sure how to set a threshold
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
      
#Sensitivity/Specificity, Do separately for 2 index tests (smartphone, SLR)
  #Reference standard: TF by field grade
    #SLR
Conf(x = photo_pcr_master$SLR_tf_yn_cons, ref = photo_pcr_master$clinic_tf_yn)
BinomCI(225, 225+29) #sens CI   est    lwr.ci    upr.ci
                          # 0.8858268 0.8408377 0.9193194
BinomCI(123, 123+36) #spec CI   est    lwr.ci    upr.ci
                          # 0.7735849 0.7025283 0.8317336
Conf(x = photo_pcr_master$SLR_tf_yn_trump, ref = photo_pcr_master$clinic_tf_yn) #quite similar
    #smartphone
Conf(x = photo_pcr_master$smart_tf_yn_cons, ref = photo_pcr_master$clinic_tf_yn) #more sensitive for TF than SLR
BinomCI(241, 241+15) #sens CI   est    lwr.ci    upr.ci
                          # 0.9414062 0.9055877 0.9641734
BinomCI(103, 103 + 54) #spec CI   est    lwr.ci    upr.ci
                            # 0.656051 0.5788177 0.7258301
Conf(x = photo_pcr_master$smart_tf_yn_trump, ref = photo_pcr_master$clinic_tf_yn)
  #Reference standard: PCR
    #SLR
Conf(x = photo_pcr_master$SLR_tf_yn_cons, ref = photo_pcr_master$newindpcr)
Conf(x = photo_pcr_master$SLR_tf_yn_trump, ref = photo_pcr_master$newindpcr)
    #smartphone
Conf(x = photo_pcr_master$smart_tf_yn_cons, ref = photo_pcr_master$newindpcr)  #pretty specific with good PPV
Conf(x = photo_pcr_master$smart_tf_yn_trump, ref = photo_pcr_master$newindpcr)

#JN note: I know TF is the most important outcome for public health purposes but I repeated the analysis for TI too
  #Reference standard: TF by field grade
    #SLR
Conf(x = photo_pcr_master$SLR_ti_yn_cons, ref = photo_pcr_master$clinic_ti_yn)
BinomCI(343, 343+35) #sens CI   est    lwr.ci    upr.ci
                          #  0.9074074 0.8739479 0.9326696
BinomCI(27, 27+8) #spec CI   est    lwr.ci    upr.ci
                          # 0.7714286 0.6098268 0.8793412
Conf(x = photo_pcr_master$SLR_ti_yn_trump, ref = photo_pcr_master$clinic_ti_yn) #quite similar
    #smartphone
Conf(x = photo_pcr_master$smart_ti_yn_cons, ref = photo_pcr_master$clinic_ti_yn) #slightly more sensitive than SLR
BinomCI(352, 352+26) #sens CI   est    lwr.ci    upr.ci
                          #  0.9312169 0.901126 0.9526315
BinomCI(25, 25+10) #spec CI   est    lwr.ci    upr.ci
                          # 0.7142857 0.5494507 0.8367346
Conf(x = photo_pcr_master$smart_ti_yn_trump, ref = photo_pcr_master$clinic_ti_yn)
  #Reference standard: PCR
    #SLR
Conf(x = photo_pcr_master$SLR_ti_yn_cons, ref = photo_pcr_master$newindpcr)
Conf(x = photo_pcr_master$SLR_ti_yn_trump, ref = photo_pcr_master$newindpcr) #same, spec of 1... is this correct?
    #smartphone
Conf(x = photo_pcr_master$smart_ti_yn_cons, ref = photo_pcr_master$newindpcr)  #slightly better sens than SLR but worse spec
Conf(x = photo_pcr_master$smart_ti_yn_trump, ref = photo_pcr_master$newindpcr)

pospcr <- master1 %>%
  select(newindpcr) %>%
  filter(newindpcr == 1)