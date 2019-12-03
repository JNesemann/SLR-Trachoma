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
    #Prevalence of TF, TI, TF+/-TI per village. JN comment: how to divide into villages, look at other TIRET 3 datasets
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
skim(TIRET3_PhotoExamData13Villages)  #this dataset has 13 unique stateteam-codes, I am assuming this corresponds to the 13 villages

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
  mutate(clinic_ti_yn = as.numeric(grepl("TI", clinical_exam)),
         clinic_tf_yn = as.numeric(grepl("TF", clinical_exam)),
         clinic_tfti_yn = as.numeric(grepl("TF/TI", clinical_exam)),
         number=as.numeric(number))
skim(TIRET3_villages)    #JN: even when using this data set to calculate clinical grades we are still missing grades for 80 individuals...

#merging PCR data with complete data in TIRET3_PhotoExamData13Villages, this will provide the stateteam-code which I am using as a proxy for village ID
PCR_exam <- left_join(PCR, TIRET3_villages, by = "number")
skim(PCR)
skim(TIRET3_villages)
skim(PCR_exam)

#selecting blakes images to calculate SLR and smartphone ICC
slr_blake <- SLR_import %>%
  rename(id = record_id, complete = trachoma_grading_complete) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  #filtering out dde 2, to give just blakes images
  filter(dde == 1) %>%
  select(-dde) %>%
  filter(!(grepl("est", id))) %>%
  mutate(id=as.numeric(id)) %>%
  #mutating grades
  mutate(tf_yn=case_when(tf %in% c(1, 2) ~ as.integer(1),
                         tf %in% c(3, 4) ~ as.integer(0)),
         ti_yn=case_when(ti %in% c(1, 2) ~ as.integer(1),
                         ti %in% c(3, 4)  ~ as.integer(0)),
         tfti_yn=case_when(tf %in% c(1, 2) & ti %in% c(1, 2) ~ as.integer(1),
                           tf %in% c(3, 4) | ti %in% c(3, 4) ~ as.integer(1)))
  #checking for good quality images without grades
  dummy_blake <- slr_blake %>%
    filter(quality %in% c(1,2) & (is.na(tf) | is.na(ti)))
  skim(dummy_blake) #4 images
  #checking for duplicates
  blake_dups <- slr_blake %>%
    group_by(id) %>%
    mutate(dup=n()) #xtabs(data=blake_dups, ~dup, addNA=TRUE) none
  
  #Merging blakes image grades and master key and putting SLR and iphone grades in the same row (widening data set)
slr_blake_wide <- right_join(slr_blake, master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  select(number, camera.instance, id, everything()) %>%
  select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_yn) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE)
skim(slr_blake_wide)
  #Blakes intra-rater reliability
    #SLR TF intra-rater reliability
xtabs(data=slr_blake_wide, ~ SLR_1_tf_yn+SLR_2_tf_yn, addNA=TRUE) 
CohenKappa(slr_blake_wide$SLR_1_tf_yn, slr_blake_wide$SLR_2_tf_yn, conf.level=0.95)
    #SLR TI intra-rater reliability
xtabs(data=slr_blake_wide, ~ SLR_1_ti_yn+SLR_2_ti_yn, addNA=TRUE) 
CohenKappa(slr_blake_wide$SLR_1_ti_yn, slr_blake_wide$SLR_2_ti_yn, conf.level=0.95)
    #smartphone TF intra-rater reliability
xtabs(data=slr_blake_wide, ~ smartphone_1_tf_yn+smartphone_2_tf_yn, addNA=TRUE) 
CohenKappa(slr_blake_wide$smartphone_1_tf_yn, slr_blake_wide$smartphone_2_tf_yn, conf.level=0.95)
    #smartphone TI intra-rater reliability
xtabs(data=slr_blake_wide, ~ smartphone_1_ti_yn+smartphone_2_ti_yn, addNA=TRUE) 
CohenKappa(slr_blake_wide$smartphone_1_ti_yn, slr_blake_wide$smartphone_2_ti_yn, conf.level=0.95)

#Repeating the process for John
slr_john <- SLR_import %>%
  rename(id = record_id, complete = trachoma_grading_complete) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  #filtering out dde 2, to give just johns images
  filter(dde == 2) %>%
  select(-dde) %>%
  filter(!(grepl("est", id))) %>%
  mutate(id=as.numeric(id)) %>%
  #mutating grades
  mutate(tf_yn=case_when(tf %in% c(1, 2) ~ as.integer(1),
                         tf %in% c(3, 4) ~ as.integer(0)),
         ti_yn=case_when(ti %in% c(1, 2) ~ as.integer(1),
                         ti %in% c(3, 4)  ~ as.integer(0)),
         tfti_yn=case_when(tf %in% c(1, 2) & ti %in% c(1, 2) ~ as.integer(1),
                           tf %in% c(3, 4) | ti %in% c(3, 4) ~ as.integer(1)))
  #checking for good quality images without grades
dummy_john <- slr_john %>%
  filter(quality %in% c(1,2) & (is.na(tf) | is.na(ti)))
skim(dummy_john) #0 images
  #checking for duplicates
john_dups <- slr_john %>%
  group_by(id) %>%
  mutate(dup=n()) #xtabs(data=john_dups, ~dup, addNA=TRUE) none

  #Merging johns image grades and master key and putting SLR and iphone grades in the same row (widening data set)
slr_john_wide <- right_join(slr_john, master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  select(number, camera.instance, id, everything()) %>%
  select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_yn) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE)
skim(slr_john_wide)
  #Johns intra-rater reliability
    #SLR TF intra-rater reliability
xtabs(data=slr_john_wide, ~ SLR_1_tf_yn+SLR_2_tf_yn, addNA=TRUE) 
CohenKappa(slr_john_wide$SLR_1_tf_yn, slr_john_wide$SLR_2_tf_yn, conf.level=0.95)
    #SLR TI intra-rater reliability
xtabs(data=slr_john_wide, ~ SLR_1_ti_yn+SLR_2_ti_yn, addNA=TRUE) 
CohenKappa(slr_john_wide$SLR_1_ti_yn, slr_john_wide$SLR_2_ti_yn, conf.level=0.95)
    #smartphone TF intra-rater reliability
xtabs(data=slr_john_wide, ~ smartphone_1_tf_yn+smartphone_2_tf_yn, addNA=TRUE) 
CohenKappa(slr_john_wide$smartphone_1_tf_yn, slr_john_wide$smartphone_2_tf_yn, conf.level=0.95)
    #smartphone TI intra-rater reliability
xtabs(data=slr_john_wide, ~ smartphone_1_ti_yn+smartphone_2_ti_yn, addNA=TRUE) 
CohenKappa(slr_john_wide$smartphone_1_ti_yn, slr_john_wide$smartphone_2_ti_yn, conf.level=0.95)
    #not that great...

#cleaning SLR import
SLR <- SLR_import %>%
  #clean and separate graders 1 and 2
  rename(
    id = record_id,
    complete = trachoma_grading_complete) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  #creating a wide dataset
  gather(field, value, quality:complete) %>%
  mutate(field_dde = paste(field, dde, sep = "_")) %>%
  select(-field, -dde) %>%
  spread(field_dde, value, convert = TRUE) %>%
  #Get rid of the 2 "test" id's, then convert to numeric
  filter(!(grepl("est", id))) %>%
  mutate(id=as.numeric(id))
skim(SLR)

#merging datasets
slr_master_wide <- right_join(SLR, master_key, by=c("id" = "mask")) %>%
  #JK - get number as the first field, camera as second, then doesn't matter
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  select(number, camera.instance, id, everything()) %>%
  select(-camera, -repeat.instance) %>%
  gather(field, value, id:ti_NA) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE)
  
slr_master_wide_grades <- slr_master_wide %>%
  #selecting columns relevant to analysis and removing grades for the duplicate images (e.g., "SLR/smartphone_2_")
  select(number, SLR_1_tf_1:SLR_1_ti_NA, smartphone_1_tf_1:smartphone_1_ti_NA) %>%
  #creating new SLR trump grades
  mutate(SLR_tf_yn_trump=case_when(SLR_1_tf_NA %in% c(1,2) ~ as.integer(1),
                                   SLR_1_tf_NA %in% c(3,4) ~ as.integer(0),
                                   SLR_1_tf_1 %in% c(1, 2) & SLR_1_tf_2 %in% c(1, 2) ~ as.integer(1),
                                   SLR_1_tf_1 %in% c(3, 4) & SLR_1_tf_2 %in% c(3, 4) ~ as.integer(0),
                               TRUE ~ NA_integer_),
         SLR_ti_yn_trump=case_when(SLR_1_ti_NA %in% c(1,2) ~ as.integer(1),
                                   SLR_1_ti_NA %in% c(3,4) ~ as.integer(0),
                                   SLR_1_ti_1 %in% c(1, 2) & SLR_1_ti_2 %in% c(1, 2) ~ as.integer(1),
                                   SLR_1_ti_1 %in% c(3, 4) & SLR_1_ti_2 %in% c(3, 4) ~ as.integer(0),
                               TRUE ~ NA_integer_),
         SLR_tfti_yn_trump=case_when(SLR_1_ti_NA %in% c(1,2) & SLR_1_tf_NA %in% c(1,2) ~ as.integer(1),
                                     SLR_1_ti_NA %in% c(3,4) & SLR_1_tf_NA %in% c(3,4) ~ as.integer(0),
                                     SLR_1_ti_1 %in% c(1,2) & SLR_1_tf_1 %in% c(1,2) & SLR_1_tf_2 %in% c(1,2) & SLR_1_ti_2 %in% c(1,2) ~ as.integer(1),
                                     SLR_1_ti_1 %in% c(3, 4) | SLR_1_ti_2 %in% c(3, 4) | SLR_1_tf_1 %in% c(3, 4) | SLR_1_tf_2 %in% c(3, 4)  ~ as.integer(0),
                                     TRUE ~ NA_integer_)) %>%
  #creating new smartphone trump grades
  mutate(smart_tf_yn_trump=case_when(smartphone_1_tf_NA %in% c(1,2) ~ as.integer(1),
                                     smartphone_1_tf_NA %in% c(3,4) ~ as.integer(0),
                                     smartphone_1_tf_1 %in% c(1, 2) & smartphone_1_tf_2 %in% c(1, 2) ~ as.integer(1),
                                     smartphone_1_tf_1 %in% c(3, 4) & smartphone_1_tf_2 %in% c(3, 4) ~ as.integer(0),
                                   TRUE ~ NA_integer_),
         smart_ti_yn_trump=case_when(smartphone_1_ti_NA %in% c(1,2) ~ as.integer(1),
                                     smartphone_1_ti_NA %in% c(3,4) ~ as.integer(0),
                                     smartphone_1_ti_1 %in% c(1, 2) & smartphone_1_ti_2 %in% c(1, 2) ~ as.integer(1),
                                     smartphone_1_ti_1 %in% c(3, 4) & smartphone_1_ti_2 %in% c(3, 4) ~ as.integer(0),
                                   TRUE ~ NA_integer_),
         smart_tfti_yn_trump=case_when(smartphone_1_ti_NA %in% c(1,2) & smartphone_1_tf_NA %in% c(1,2) ~ as.integer(1),
                                       smartphone_1_ti_NA %in% c(3,4) & smartphone_1_tf_NA %in% c(3,4) ~ as.integer(0),
                                       smartphone_1_ti_1 %in% c(1,2) & smartphone_1_tf_1 %in% c(1,2) & smartphone_1_tf_2 %in% c(1,2) & smartphone_1_ti_2 %in% c(1,2) ~ as.integer(1),
                                       smartphone_1_ti_1 %in% c(3,4) | smartphone_1_ti_2 %in% c(3,4) | smartphone_1_tf_1 %in% c(3,4) | smartphone_1_tf_2 %in% c(3,4)  ~ as.integer(0),
                                     TRUE ~ NA_integer_)) %>%
  #creating new SLR consensus grades
  mutate(SLR_tf_yn_cons=case_when(SLR_1_tf_NA %in% c(1,2) & SLR_1_tf_1 %in% c(1,2) | SLR_1_tf_NA %in% c(1,2) & SLR_1_tf_2 %in% c(1,2) | SLR_1_tf_1 %in% c(1,2) & SLR_1_tf_2 %in% c(1,2) ~ as.integer(1),
                                  SLR_1_tf_NA %in% c(3,4) & SLR_1_tf_1 %in% c(3,4) | SLR_1_tf_NA %in% c(3,4) & SLR_1_tf_2 %in% c(3,4) | SLR_1_tf_1 %in% c(3,4) & SLR_1_tf_2 %in% c(3,4) ~ as.integer(0),
                                  TRUE ~ NA_integer_),
         SLR_ti_yn_cons=case_when(SLR_1_ti_NA %in% c(1,2) & SLR_1_ti_1 %in% c(1,2) | SLR_1_ti_NA %in% c(1,2) & SLR_1_ti_2 %in% c(1,2) | SLR_1_ti_1 %in% c(1,2) & SLR_1_ti_2 %in% c(1,2) ~ as.integer(1),
                                  SLR_1_ti_NA %in% c(3,4) & SLR_1_ti_1 %in% c(3,4) | SLR_1_ti_NA %in% c(3,4) & SLR_1_ti_2 %in% c(3,4) | SLR_1_ti_1 %in% c(3,4) & SLR_1_ti_2 %in% c(3,4) ~ as.integer(0),
                                  TRUE ~ NA_integer_)) %>%
  mutate(SLR_tfti_yn_cons=if_else(SLR_tf_yn_cons == 1 & SLR_ti_yn_cons == 1, as.integer(1), as.integer(0))) %>%
  #creating new smartphone consensus grades
  mutate(smart_tf_yn_cons=case_when(smartphone_1_tf_NA %in% c(1,2) & smartphone_1_tf_1 %in% c(1,2) | smartphone_1_tf_NA %in% c(1,2) & smartphone_1_tf_2 %in% c(1,2) | smartphone_1_tf_1 %in% c(1,2) & smartphone_1_tf_2 %in% c(1,2) ~ as.integer(1),
                                    smartphone_1_tf_NA %in% c(3,4) & smartphone_1_tf_1 %in% c(3,4) | smartphone_1_tf_NA %in% c(3,4) & smartphone_1_tf_2 %in% c(3,4) | smartphone_1_tf_1 %in% c(3,4) & smartphone_1_tf_2 %in% c(3,4) ~ as.integer(0),
                                  TRUE ~ NA_integer_),
         smart_ti_yn_cons=case_when(smartphone_1_ti_NA %in% c(1,2) & smartphone_1_ti_1 %in% c(1,2) | smartphone_1_ti_NA %in% c(1,2) & smartphone_1_ti_2 %in% c(1,2) | smartphone_1_ti_1 %in% c(1,2) & smartphone_1_ti_2 %in% c(1,2) ~ as.integer(1),
                                    smartphone_1_ti_NA %in% c(3,4) & smartphone_1_ti_1 %in% c(3,4) | smartphone_1_ti_NA %in% c(3,4) & smartphone_1_ti_2 %in% c(3,4) | smartphone_1_ti_1 %in% c(3,4) & smartphone_1_ti_2 %in% c(3,4) ~ as.integer(0),
                                  TRUE ~ NA_integer_)) %>%
  mutate(smart_tfti_yn_cons=if_else(smart_tf_yn_cons == 1 & smart_ti_yn_cons == 1, as.integer(1), as.integer(0))) %>%
  #converting individual grades from 1-4 to binary
  #first SLRs
  mutate(slr_tf_1=case_when(SLR_1_tf_1 %in% c(1,2) ~ as.integer(1),
                            SLR_1_tf_1 %in% c(3,4) ~ as.integer(0)),
         slr_tf_2=case_when(SLR_1_tf_2 %in% c(1,2) ~ as.integer(1),
                            SLR_1_tf_2 %in% c(3,4) ~ as.integer(0)),
         slr_tf_NA=case_when(SLR_1_tf_NA %in% c(1,2) ~ as.integer(1),
                            SLR_1_tf_NA %in% c(3,4) ~ as.integer(0)),
         slr_ti_1=case_when(SLR_1_ti_1 %in% c(1,2) ~ as.integer(1),
                            SLR_1_ti_1 %in% c(3,4) ~ as.integer(0)),
         slr_ti_2=case_when(SLR_1_ti_2 %in% c(1,2) ~ as.integer(1),
                            SLR_1_ti_2 %in% c(3,4) ~ as.integer(0)),
         slr_ti_NA=case_when(SLR_1_ti_NA %in% c(1,2) ~ as.integer(1),
                             SLR_1_ti_NA %in% c(3,4) ~ as.integer(0)),
         #and now smartphones
         smart_tf_1=case_when(smartphone_1_tf_1 %in% c(1,2) ~ as.integer(1),
                              smartphone_1_tf_1 %in% c(3,4) ~ as.integer(0)),
         smart_tf_2=case_when(smartphone_1_tf_2 %in% c(1,2) ~ as.integer(1),
                              smartphone_1_tf_2 %in% c(3,4) ~ as.integer(0)),
         smart_tf_NA=case_when(smartphone_1_tf_NA %in% c(1,2) ~ as.integer(1),
                               smartphone_1_tf_NA %in% c(3,4) ~ as.integer(0)),
         smart_ti_1=case_when(smartphone_1_ti_1 %in% c(1,2) ~ as.integer(1),
                              smartphone_1_ti_1 %in% c(3,4) ~ as.integer(0)),
         smart_ti_2=case_when(smartphone_1_ti_2 %in% c(1,2) ~ as.integer(1),
                              smartphone_1_ti_2 %in% c(3,4) ~ as.integer(0)),
         smart_ti_NA=case_when(smartphone_1_ti_NA %in% c(1,2) ~ as.integer(1),
                               smartphone_1_ti_NA %in% c(3,4) ~ as.integer(0)))
skim(slr_master_wide_grades)
head(slr_master_wide_grades)

#Kappas:
  #among slr photos john vs blake
    #SLR TF, 0.654
    xtabs(data=slr_master_wide_grades, ~ slr_tf_1+slr_tf_2, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$slr_tf_1, slr_master_wide_grades$slr_tf_2, conf.level=0.95)
    #SLR TI, 0.711
    xtabs(data=slr_master_wide_grades, ~ slr_ti_1+slr_ti_2, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$slr_ti_1, slr_master_wide_grades$slr_ti_2, conf.level=0.95)
  #among smartphone photos john vs blake
    #smartphone TF 0.55...
    xtabs(data=slr_master_wide_grades, ~ smart_tf_1+smart_tf_2, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_tf_1, slr_master_wide_grades$smart_tf_2, conf.level=0.95)
    #smartphone TI 0.74
    xtabs(data=slr_master_wide_grades, ~ smart_ti_1+smart_ti_2, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_ti_1, slr_master_wide_grades$smart_ti_2, conf.level=0.95)
  #consensus grade smartphone vs consensus SLR
    #consensus TF  0.72
    xtabs(data=slr_master_wide_grades, ~ smart_tf_yn_cons + SLR_tf_yn_cons, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_tf_yn_cons, slr_master_wide_grades$SLR_tf_yn_cons, conf.level=0.95)
    #consensus TI 0.809
    xtabs(data=slr_master_wide_grades, ~ smart_ti_yn_cons + SLR_ti_yn_cons, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_ti_yn_cons, slr_master_wide_grades$SLR_ti_yn_cons, conf.level=0.95)
  #trump grade smartphone vs trump SLR
    #consensus TF
    xtabs(data=slr_master_wide_grades, ~ smart_tf_yn_trump + SLR_tf_yn_trump, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_tf_yn_trump, slr_master_wide_grades$SLR_tf_yn_trump, conf.level=0.95)
    #consensus TI
    xtabs(data=slr_master_wide_grades, ~ smart_ti_yn_trump + SLR_ti_yn_trump, addNA=TRUE) 
    CohenKappa(slr_master_wide_grades$smart_ti_yn_trump, slr_master_wide_grades$SLR_ti_yn_trump, conf.level=0.95)
    
#merging in PCR data
photo_pcr_master <- left_join(slr_master_wide_grades, PCR_exam, by = "number")
skim(photo_pcr_master)

#kappas comparing image grades with field grades
  #blake images vs fieldgrades
    #SLR TF vs in field TF
  xtabs(data = photo_pcr_master, ~slr_tf_1 + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$slr_tf_1, photo_pcr_master$clinic_tf_yn, conf.level=0.95)
    #smartphone TF vs in field TF
  xtabs(data = photo_pcr_master, ~smart_tf_1 + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_tf_1, photo_pcr_master$clinic_tf_yn, conf.level=0.95)    
    #SLR TI vs in field TI    
  xtabs(data = photo_pcr_master, ~slr_ti_1 + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$slr_ti_1, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
    #smartphone TI vs in field TI
  xtabs(data = photo_pcr_master, ~smart_ti_1 + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_ti_1, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
  #johns images vs fieldgrades
    #SLR TF vs in field TF
  xtabs(data = photo_pcr_master, ~slr_tf_2 + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$slr_tf_2, photo_pcr_master$clinic_tf_yn, conf.level=0.95)
    #smartphone TF vs in field TF
  xtabs(data = photo_pcr_master, ~smart_tf_2 + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_tf_2, photo_pcr_master$clinic_tf_yn, conf.level=0.95)    
    #SLR TI vs in field TI    
  xtabs(data = photo_pcr_master, ~slr_ti_2 + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$slr_ti_2, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
    #smartphone TI vs in field TI
  xtabs(data = photo_pcr_master, ~smart_ti_2 + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_ti_2, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
  #consensus grades vs infield
    #SLR TF vs in field TF
  xtabs(data = photo_pcr_master, ~SLR_tf_yn_cons + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$SLR_tf_yn_cons, photo_pcr_master$clinic_tf_yn, conf.level=0.95)
    #smartphone TF vs in field TF
  xtabs(data = photo_pcr_master, ~smart_tf_yn_cons + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_tf_yn_cons, photo_pcr_master$clinic_tf_yn, conf.level=0.95)    
    #SLR TI vs in field TI    
  xtabs(data = photo_pcr_master, ~SLR_ti_yn_cons + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$SLR_ti_yn_cons, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
    #smartphone TI vs in field TI
  xtabs(data = photo_pcr_master, ~smart_ti_yn_cons + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_ti_yn_cons, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
  #trump grades vs infield
    #SLR TF vs in field TF
  xtabs(data = photo_pcr_master, ~SLR_tf_yn_trump + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$SLR_tf_yn_trump, photo_pcr_master$clinic_tf_yn, conf.level=0.95)
    #smartphone TF vs in field TF
  xtabs(data = photo_pcr_master, ~smart_tf_yn_trump + clinic_tf_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_tf_yn_trump, photo_pcr_master$clinic_tf_yn, conf.level=0.95)    
    #SLR TI vs in field TI    
  xtabs(data = photo_pcr_master, ~SLR_ti_yn_trump + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$SLR_ti_yn_trump, photo_pcr_master$clinic_ti_yn, conf.level=0.95)
    #smartphone TI vs in field TI
  xtabs(data = photo_pcr_master, ~smart_ti_yn_trump + clinic_ti_yn, addNA = TRUE)
  CohenKappa(photo_pcr_master$smart_ti_yn_trump, photo_pcr_master$clinic_ti_yn, conf.level=0.95)

#SLR, smartphone, and PCR prevalences per village
village_image_prev <- photo_pcr_master %>%
  group_by(state_code) %>%
  summarise(slr_tf_cons=mean(SLR_tf_yn_cons, na.rm = TRUE),
            slr_tf_trump=mean(SLR_tf_yn_trump, na.rm = TRUE),
            smart_tf_cons=mean(smart_tf_yn_cons, na.rm = TRUE),
            smart_tf_trump=mean(smart_tf_yn_trump, na.rm = TRUE),
            slr_ti_cons=mean(SLR_ti_yn_cons, na.rm = TRUE),
            slr_ti_trump=mean(SLR_ti_yn_trump, na.rm = TRUE),
            smart_ti_cons=mean(smart_ti_yn_cons, na.rm = TRUE),
            smart_ti_trump=mean(smart_ti_yn_trump, na.rm = TRUE),
            pcr=mean(newindpcr, na.rm = TRUE))  #these are all binary so the mean should reflect the prevalence

slr_tf_consCI <- MeanCI(photo_pcr_master$SLR_tf_yn_cons, sd = NULL, na.rm = TRUE, conf.level = 0.95, sides = "two.sided")
slr_tf_trumpCI <- MeanCI(photo_pcr_master$SLR_tf_yn_trump, sd = NULL, na.rm = TRUE, conf.level = 0.95, sides = "two.sided")

slrtf_pos <- sum(photo_pcr_master$SLR_tf_yn_cons == 1, na.rm = TRUE)
slrtf_tot <- length(photo_pcr_master$SLR_tf_yn_cons)
prop.test(slrtf_pos, slrtf_tot)

skim(photo_pcr_master)

#calculating village level prevalence from TIRET 3, clinical exams had a larger n (~1300 kids) so calculating clinical prevalence seperately
clinical_prev <- TIRET3_villages%>%
  group_by(state_code) %>%
  summarise(clinic_tf=mean(clinic_tf_yn),
            clinic_ti=mean(clinic_ti_yn))

#merging clinical and image prevalence described above.
village_prevalence <- left_join(village_image_prev, clinical_prev, by = "state_code")
#Alternatively I can use the Conf() function from the DescTools package which should give the prevalence w/ conf intervals and also sens and spec. See sens/spec section below

#Then plot different combinations of prevalences
#restructuring the tibble so I can make the graphs I want
dummy_village <- village_prevalence %>%
  select(-pcr) %>%  #I chose not to include pcr as only three villages had non-zero prevalences
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
