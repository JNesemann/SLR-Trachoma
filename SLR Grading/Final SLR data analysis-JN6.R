#Final SLR Data analysis

#JN : major addition is the calculation of ppv/npv with BCA CIs in lines ***

#  library(venneuler) # Note that on a Mac you need to also install the Java SE Development Kit
                   # And look at the error message but for this version it required the 11.0.1 dmg file
                   # Google "Java SE development kit archive" to find it
library(eulerr)
library(tidyverse)
library(skimr)
library(readr)
library(readxl)
library(DescTools)
library(irr)
library(purrr)
library(ICCbin)

################################
##     IMPORT AND CLEAN      ##
################################
master_key_import <- read_excel("Masked Number Master Key (1).xlsx")
TIRET3_PCR_import <- read_excel("TIRET3_PCR_for John.xlsx")
SLR_import <- read_csv("SLRTrachoma_DATA_2020-01-03_1209.csv")
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
         clinic_tfti_yn = if_else(!is.na(clinical_exam) & (clinic_tf_yn==1 | clinic_ti_yn==1), 1, 0),
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
         ti_di=if_else(ti %in% c(1, 2),1,if_else(ti %in% c(3, 4),0,NA_real_)),
         tfti_di=if_else(tf %in% c(1, 2) | ti %in% c(1, 2),1,if_else(tf %in% c(3, 4) & ti %in% c(3, 4),0,NA_real_)),
         ungradable=if_else(quality==3, 1, 0),
         nophotos=if_else(quality==4, 1, 0)) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  # JK: Get rid of the 2 "test" id's, then convert to numeric
  filter(!(grepl("est", id))) %>%
  # JK: This group_by is new; 
  mutate(id=as.numeric(id)) %>%
  group_by(id) %>%
  mutate(sumtf=sum(tf_di, na.rm=TRUE),
         sumti=sum(ti_di, na.rm=TRUE),
         sumtfti=sum(tfti_di, na.rm=TRUE),
         sumungrad=sum(ungradable, na.rm=TRUE),
         sumnophotos=sum(nophotos, na.rm=TRUE),
         totalid=n(), # xtabs(data=master1, ~sumtf+totalid, addNA=TRUE)
         tf_yn=case_when(totalid==3 & sumtf>=2 ~ 1,
                         totalid==3 & sumtf<2 ~ 0,
                         totalid==2 & sumtf==2 ~ 1,
                         totalid==2 & sumtf==0 ~ 0,
                         TRUE ~ NA_real_),
         ti_yn=case_when(totalid==3 & sumti>=2 ~ 1,
                         totalid==3 & sumti<2 ~ 0,
                         totalid==2 & sumti==2 ~ 1,
                         totalid==2 & sumti==0 ~ 0,
                         TRUE ~ NA_real_),
         tfti_yn=case_when(totalid==3 & sumtfti>=2 ~ 1,
                           totalid==3 & sumtfti<2 ~ 0,
                           totalid==2 & sumtfti==2 ~ 1,
                           totalid==2 & sumtfti==0 ~ 0,
                           TRUE ~ NA_real_),
         ungrad_yn=case_when(totalid==3 & sumungrad>=2 ~ 1,
                             totalid==3 & sumungrad<2 ~ 0,
                             totalid==2 & sumungrad==2 ~ 1,
                             totalid==2 & sumungrad==0 ~ 0,
                             TRUE ~ NA_real_),
         nophoto_yn=case_when(totalid==3 & sumnophotos>=2 ~ 1,
                              totalid==3 & sumnophotos<2 ~ 0,
                              totalid==2 & sumnophotos==2 ~ 1,
                              totalid==2 & sumnophotos==0 ~ 0,
                              TRUE ~ NA_real_)) %>%
  dplyr::select(-complete, -sumtf, -sumti, -sumtfti, -sumungrad, -sumnophotos, -totalid) %>%
  dplyr::select(id, tf_yn, ti_yn, tfti_yn, ungrad_yn, nophoto_yn, dde, tf_di, ti_di, tfti_di, ungradable, nophotos, everything()) %>% 
  #creating a wide dataset
  gather(field, value, tf_di:notes) %>%
  mutate(field_dde = paste(field, dde, sep = "_")) %>%
  dplyr::select(-field, -dde) %>%
  spread(field_dde, value, convert = TRUE) %>%
  # JK: Note that when you enter the dot below it means to use the current data, so you're merging into the current dataframe
  right_join(., master_key, by=c("id" = "mask")) %>%
  mutate(camera.instance=paste(camera, repeat.instance, sep="_")) %>%
  dplyr::select(number, camera.instance, id, tf_yn, ti_yn, tfti_yn, ungrad_yn, nophoto_yn, everything()) %>%
  dplyr::select(-camera, -repeat.instance) %>%
  gather(field, value, id:ungradable_NA) %>%
  mutate(camera_field=paste(camera.instance, field, sep="_")) %>%
  dplyr::select(-camera.instance, -field) %>%
  spread(camera_field, value, convert=TRUE) %>%
  left_join(., PCR_exam, by="number") %>%
  # JK: Let's restrict to only those with PCR for simplicity
  # Also, can run without "& smartphone_1_ungrad_yn!=1" to get the 412 children. Exclude to get final sample of 408
  filter(!is.na(newindpcr) & SLR_1_nophoto_yn!=1 & smartphone_1_nophoto_yn!=1 & smartphone_1_ungrad_yn!=1)
# JK: Now this is one big dataset. The variable names follow the pattern:
# SLR vs smartphone, then...
# whether it's a repeat (1 is the first instance; 2 is the second instance; so primary analyses should only be with 1)
# variable (tf /ti etc)
# grader (1=Blake, 2=John, NA=Jeremy)

master1long <- master1 %>%
  gather(field, value, SLR_1_angle_1:clinic_tfti_yn)
master1longdups <- master1 %>%
  filter(!is.na(smartphone_2_id) | !is.na(SLR_2_id)) %>%
  gather(field, value, SLR_1_angle_1:clinic_tfti_yn)

# Finding missing TF grades (2=John, 1=Blake)
# ts2.sm <- master1 %>%
#   filter((is.na(smartphone_2_tf_di_2) & !is.na(smartphone_2_id)) | is.na(smartphone_1_tf_di_2) & !is.na(smartphone_1_id)) %>%
#   select(smartphone_1_id, smartphone_2_id, smartphone_1_tf_di_1, smartphone_1_tf_di_2, smartphone_2_tf_di_1, smartphone_2_tf_di_2)
# ts1.sm <- master1 %>%
#   filter((is.na(smartphone_2_tf_di_1) & !is.na(smartphone_2_id)) | is.na(smartphone_1_tf_di_1) & !is.na(smartphone_1_id)) %>%
#   select(smartphone_1_id, smartphone_2_id, smartphone_1_tf_di_1, smartphone_1_tf_di_2, smartphone_2_tf_di_1, smartphone_2_tf_di_2)
# ts2.sl <- master1 %>%
#   filter((is.na(SLR_2_tf_di_2) & !is.na(SLR_2_id)) | is.na(SLR_1_tf_di_2) & !is.na(SLR_1_id)) %>%
#   select(SLR_1_id, SLR_2_id, SLR_1_tf_di_1, SLR_1_tf_di_2, SLR_2_tf_di_1, SLR_2_tf_di_2)
# ts1.sl <- master1 %>%
#   filter((is.na(SLR_2_tf_di_1) & !is.na(SLR_2_id)) | is.na(SLR_1_tf_di_1) & !is.na(SLR_1_id)) %>%
#   select(SLR_1_id, SLR_2_id, SLR_1_tf_di_1, SLR_1_tf_di_2, SLR_2_tf_di_1, SLR_2_tf_di_2)
# Finding missing TI grades (2=John, 1=Blake)
# ts2.sm <- master1 %>%
#   filter((is.na(smartphone_2_ti_di_2) & !is.na(smartphone_2_id)) | is.na(smartphone_1_ti_di_2) & !is.na(smartphone_1_id)) %>%
#   select(smartphone_1_id, smartphone_2_id, smartphone_1_ti_di_1, smartphone_1_ti_di_2, smartphone_2_ti_di_1, smartphone_2_ti_di_2)
# ts1.sm <- master1 %>%
#   filter((is.na(smartphone_2_ti_di_1) & !is.na(smartphone_2_id)) | is.na(smartphone_1_ti_di_1) & !is.na(smartphone_1_id)) %>%
#   select(smartphone_1_id, smartphone_2_id, smartphone_1_ti_di_1, smartphone_1_ti_di_2, smartphone_2_ti_di_1, smartphone_2_ti_di_2)
# ts2.sl <- master1 %>%
#   filter((is.na(SLR_2_ti_di_2) & !is.na(SLR_2_id)) | is.na(SLR_1_ti_di_2) & !is.na(SLR_1_id)) %>%
#   select(SLR_1_id, SLR_2_id, SLR_1_ti_di_1, SLR_1_ti_di_2, SLR_2_ti_di_1, SLR_2_ti_di_2)
# ts1.sl <- master1 %>%
#   filter((is.na(SLR_2_ti_di_1) & !is.na(SLR_2_id)) | is.na(SLR_1_ti_di_1) & !is.na(SLR_1_id)) %>%
#   select(SLR_1_id, SLR_2_id, SLR_1_ti_di_1, SLR_1_ti_di_2, SLR_2_ti_di_1, SLR_2_ti_di_2)

#TI totals
addmargins(xtabs(data = master1, ~SLR_1_ti_yn + smartphone_1_ti_yn))
addmargins(xtabs(data = master1, ~SLR_1_ti_yn + clinic_ti_yn))
addmargins(xtabs(data = master1, ~smartphone_1_ti_yn + clinic_ti_yn))
#TF totals
addmargins(xtabs(data = master1, ~SLR_1_tf_yn + smartphone_1_tf_yn))
addmargins(xtabs(data = master1, ~SLR_1_tf_yn + clinic_tf_yn))
addmargins(xtabs(data = master1, ~smartphone_1_tf_yn + clinic_tf_yn))
# Examiner complete data?
addmargins(xtabs(data = master1, ~examiner+smartphone_1_tf_yn, addNA=TRUE))

################################
##          ANALYSES          ##
################################
# PAPER TEXT
# First paragraph
#summarizing demographic information
addmargins(xtabs(data=master1, ~sex,addNA=TRUE))
master1 %>% dplyr::summarize(meanage=mean(age), sdage=sd(age))
# Note that the below is after excluding the 4 kids with ungradable smartphone photos (above in master1 code)
addmargins(xtabs(data=master1, ~SLR_1_ungrad_yn+smartphone_1_ungrad_yn,addNA=TRUE))
addmargins(xtabs(data=master1, ~SLR_1_nophoto_yn+smartphone_1_nophoto_yn,addNA=TRUE))

################################
##  FIGURE 1: VENN DIAGRAM    ##
################################
# ETABLE 1: NUMBERS (NOT SURE WE NEED)
  #JN: I would not include
#ARVO abstract cross tabulations
xtabs(data = master1, ~SLR_1_tf_yn + smartphone_1_tf_yn + clinic_tf_yn)
xtabs(data = master1, ~SLR_1_ti_yn + smartphone_1_ti_yn + clinic_ti_yn)
# TF
# , , clinic_tf_yn = 0
# 
# smartphone_1_tf_yn
# SLR_1_tf_yn   0   1   C-SL-Sm
#           0 220   3   0-0-0 = 220= D         
#           1  16  12   0-1-0 = 16 = B --> SLR
#                       0-0-1 = 3  = C --> smartphone
#                       0-1-1 = 12 = B&C
# , , clinic_tf_yn = 1
# 
# smartphone_1_tf_yn
# SLR_1_tf_yn   0   1   C-SL-Sm
#           0  32   4   1-0-0 = 32 = A --> field exam
#           1  23  98   1-1-0 = 23 = A&B
#                       1-0-1 = 4  = A&C 
#                       1-1-1 = 98 = A&B&C

tfvenn <- euler(c(A=32, B=16, C=3, D=220, "A&B"=23, "A&C"=4, "B&C"=12, "A&B&C"=98, "D&A"=0, "D&B"=0, "D&C"=0), shape="ellipse")
tfvennplot <- plot(tfvenn, quantities=TRUE)
tfvennplot

tfvenn2 <- euler(c("D&A"=32, "D&B"=16, "D&C"=3, D=220, "D&A&B"=23, "D&A&C"=4, "D&B&C"=12, "D&A&B&C"=98), shape="ellipse")
tfvennplot2 <- plot(tfvenn2, quantities=TRUE)
tfvennplot2

# , , clinic_ti_yn = 0
# 
# smartphone_1_ti_yn
# SLR_1_ti_yn   0   1   C-SL-Sm
#           0 335   6   0-0-0 = 335= D      
#           1  13  20   0-1-0 = 13 = B --> SLR
#                       0-0-1 = 6  = C --> smartphone
#                       0-1-1 = 20 = B&C
# , , clinic_ti_yn = 1
# 
# smartphone_1_ti_yn
# SLR_1_ti_yn   0   1   C-SL-Sm
#           0   8   0   1-0-0 = 8  = A --> field exam
#           1   2  24   1-1-0 = 2  = A&B
#                       1-0-1 = 0  = A&C 
#                       1-1-1 = 24 = A&B&C
tivenn <- euler(c(A=8, B=13, C=6, D=335, "A&B"=2, "A&C"=0, "B&C"=20, "A&B&C"=24, "D&A"=0, "D&B"=0, "D&C"=0), shape="ellipse")
tivennplot <- plot(tivenn, quantities=TRUE)
tivennplot
tivenn2 <- euler(c(D=335, "D&A"=8, "D&B"=13, "D&C"=6, "A&B&D"=2, "A&C&D"=0, "B&C&D"=20, "A&B&C&D"=24), shape="ellipse")
tivennplot2 <- plot(tivenn2, quantities=TRUE)
tivennplot2

################################
##        TABLE 1: ICC        ##
################################
# TF
# Consensus SLR vs Consensus smartphone
addmargins(xtabs(data=master1, ~ smartphone_1_tf_yn+SLR_1_tf_yn, addNA=TRUE))
# CohenKappa(master1$smartphone_1_tf_yn, master1$SLR_1_tf_yn, conf.level=0.95)
# # Note that the following line is the original ICC in the code, but not specifically for binary data
# master1 %>% select(SLR_1_tf_yn, smartphone_1_tf_yn) %>% icc(model = "twoway", type = "agreement", conf.level = 0.95)
# ICC for binary data
master1long.slsm.tf <- master1long %>%
  filter(field %in% c("smartphone_1_tf_yn","SLR_1_tf_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.slsm.tf, method = "rm", ci.type = "rm")

# Consensus smart vs field
addmargins(xtabs(data=master1, ~ smartphone_1_tf_yn+clinic_tf_yn, addNA=TRUE))
# CohenKappa(master1$smartphone_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# ICC for binary data
master1long.smf.tf <- master1long %>%
  filter(field %in% c("smartphone_1_tf_yn","clinic_tf_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.smf.tf, method = "rm", ci.type = "rm")

# Consensus SLR vs field
addmargins(xtabs(data=master1, ~ SLR_1_tf_yn+clinic_tf_yn, addNA=TRUE))
# CohenKappa(master1$SLR_1_tf_yn, master1$clinic_tf_yn, conf.level=0.95)
# master1 %>% select(SLR_1_tf_yn, clinic_tf_yn) %>% icc(model = "twoway", type = "agreement", conf.level = 0.95)
master1long.slf.tf <- master1long %>%
  filter(field %in% c("SLR_1_tf_yn","clinic_tf_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.slf.tf, method = "rm", ci.type = "rm")

# Grader 1 vs grader 2 (smartphone)
addmargins(xtabs(data=master1, ~ smartphone_1_tf_di_1+smartphone_1_tf_di_2, addNA=TRUE)) # FYI: 1 missing value for grader 1, maybe we ask him to grade it
# CohenKappa(master1$smartphone_1_tf_di_1, master1$smartphone_1_tf_di_2, conf.level=0.95)
master1long.g1g2sm.tf <- master1long %>%
  filter(field %in% c("smartphone_1_tf_di_1","smartphone_1_tf_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g2sm.tf, method = "rm", ci.type = "rm")

# Grader 1 vs grader 2 (SLR)
addmargins(xtabs(data=master1, ~ SLR_1_tf_di_1+SLR_1_tf_di_2, addNA=TRUE)) # FYI: 1 missing value for grader 1, maybe we ask him to grade it
# CohenKappa(master1$SLR_1_tf_di_1, master1$SLR_1_tf_di_2, conf.level=0.95)
master1long.g1g2sl.tf <- master1long %>%
  filter(field %in% c("SLR_1_tf_di_1","SLR_1_tf_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g2sl.tf, method = "rm", ci.type = "rm")

# Grader 1 vs grader 1 (smartphone)
addmargins(xtabs(data=filter(master1, !is.na(smartphone_2_id)), ~ smartphone_1_tf_di_1+smartphone_2_tf_di_1, addNA=TRUE)) 
master1long.g1g1sm.tf <- master1longdups %>%
  filter(field %in% c("smartphone_1_tf_di_1","smartphone_2_tf_di_1")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g1sm.tf, method = "rm", ci.type = "rm")

# Grader 1 vs grader 1 (SLR)
addmargins(xtabs(data=filter(master1, !is.na(SLR_2_id)), ~ SLR_1_tf_di_1+SLR_2_tf_di_1, addNA=TRUE)) 
master1long.g1g1slr.tf <- master1longdups %>%
  filter(field %in% c("SLR_1_tf_di_1","SLR_2_tf_di_1")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g1slr.tf, method = "rm", ci.type = "rm")

# Grader 2 vs grader 2 (smartphone)
addmargins(xtabs(data=filter(master1, !is.na(smartphone_2_id)), ~ smartphone_1_tf_di_2+smartphone_2_tf_di_2, addNA=TRUE)) 
master1long.g2g2sm.tf <- master1longdups %>%
  filter(field %in% c("smartphone_1_tf_di_2","smartphone_2_tf_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g2g2sm.tf, method = "rm", ci.type = "rm")

# Grader 2 vs grader 2 (SLR)
addmargins(xtabs(data=filter(master1, !is.na(SLR_2_id)), ~ SLR_1_tf_di_2+SLR_2_tf_di_2, addNA=TRUE)) 
master1long.g2g2slr.tf <- master1longdups %>%
  filter(field %in% c("SLR_1_tf_di_2","SLR_2_tf_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g2g2slr.tf, method = "rm", ci.type = "rm")

#creating dataset for ICCs over time
master2wide <- SLR_import %>%
  #clean and separate graders 1 and 2
  rename(id = record_id,
         complete = trachoma_grading_complete) %>%
  # JK: This mutate is new...
  mutate(tf_di=if_else(tf %in% c(1, 2),1,if_else(tf %in% c(3, 4),0,NA_real_)),
         ti_di=if_else(ti %in% c(1, 2),1,if_else(ti %in% c(3, 4),0,NA_real_)),
         tfti_di=if_else(tf %in% c(1, 2) | ti %in% c(1, 2),1,if_else(tf %in% c(3, 4) & ti %in% c(3, 4),0,NA_real_)),
         ungradable=if_else(quality==3, 1, 0),
         nophotos=if_else(quality==4, 1, 0)) %>%
  separate(id, into = c("id", "dde"), sep = "--", remove = TRUE, convert = TRUE) %>%
  # JK: Get rid of the 2 "test" id's, then convert to numeric
  filter(!(grepl("est", id))) %>%
  mutate(id=as.factor(id)) %>%
  filter(!is.na(tf_di) & !is.na(dde)) %>%
  dplyr::select(id, dde, tf_di) %>%
  group_by(id) %>%
  mutate(count=n()) %>%
  filter(count>1) %>%
  dplyr::select(-count) %>%
  gather(field, value, tf_di) %>%
  mutate(fielddde=paste(field, dde, sep="_")) %>%
  dplyr::select(-field, -dde) %>%
  spread(fielddde, value, convert=TRUE)
# Can't get the iccbin to work for some reason...
master2wide %>% ungroup() %>% filter(as.numeric(id)<200) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
master2wide %>% ungroup() %>% filter(as.numeric(id)>=200 & as.numeric(id)<400) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
master2wide %>% ungroup() %>% filter(as.numeric(id)>=400 & as.numeric(id)<600) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
master2wide %>% ungroup() %>% filter(as.numeric(id)>=600 & as.numeric(id)<800) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
master2wide %>% ungroup() %>% filter(as.numeric(id)>=800 & as.numeric(id)<1000) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
master2wide %>% ungroup() %>% filter(as.numeric(id)>=1000) %>% dplyr::select(tf_di_1, tf_di_2) %>% irr::icc(model = "twoway", type = "agreement", conf.level = 0.95)
# This suggests less agreement earlier (first 400)

# TI
# Consensus SLR vs Consensus smartphone
addmargins(xtabs(data=master1, ~ smartphone_1_ti_yn+SLR_1_ti_yn, addNA=TRUE))
# CohenKappa(master1$smartphone_1_ti_yn, master1$SLR_1_ti_yn, conf.level=0.95)
master1long.slsm.ti <- master1long %>%
  filter(field %in% c("smartphone_1_ti_yn","SLR_1_ti_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.slsm.ti, method = "rm", ci.type = "rm")

# Consensus smart vs field
addmargins(xtabs(data=master1, ~ smartphone_1_ti_yn+clinic_ti_yn, addNA=TRUE))
# CohenKappa(master1$smartphone_1_ti_yn, master1$clinic_ti_yn, conf.level=0.95)
# ICC for binary data
master1long.smf.ti <- master1long %>%
  filter(field %in% c("smartphone_1_ti_yn","clinic_ti_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.smf.ti, method = "rm", ci.type = "rm")

# Consensus SLR vs field
addmargins(xtabs(data=master1, ~ SLR_1_ti_yn+clinic_ti_yn, addNA=TRUE))
# CohenKappa(master1$SLR_1_ti_yn, master1$clinic_ti_yn, conf.level=0.95)
master1long.slf.ti <- master1long %>%
  filter(field %in% c("SLR_1_ti_yn","clinic_ti_yn")) %>% mutate(value=as.numeric(value)) %>% dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.slf.ti, method = "rm", ci.type = "rm")

# Grader 1 vs grader 2 (smartphone)
addmargins(xtabs(data=master1, ~ smartphone_1_ti_di_1+smartphone_1_ti_di_2, addNA=TRUE)) # FYI: 1 missing value for grader 1, maybe we ask him to grade it
# CohenKappa(master1$smartphone_1_ti_di_1, master1$smartphone_1_ti_di_2, conf.level=0.95)
master1long.g1g2sm.ti <- master1long %>%
  filter(field %in% c("smartphone_1_ti_di_1","smartphone_1_ti_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g2sm.ti, method = "rm", ci.type = "rm")

# Grader 1 vs grader 2 (SLR)
addmargins(xtabs(data=master1, ~ SLR_1_ti_di_1+SLR_1_ti_di_2, addNA=TRUE)) # FYI: 1 missing value for grader 1, maybe we ask him to grade it
# CohenKappa(master1$SLR_1_ti_di_1, master1$SLR_1_ti_di_2, conf.level=0.95)
master1long.g1g2sl.ti <- master1long %>%
  filter(field %in% c("SLR_1_ti_di_1","SLR_1_ti_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g2sl.ti, method = "rm", ci.type = "rm")

# Grader 1 vs grader 1 (smartphone)
addmargins(xtabs(data=filter(master1, !is.na(smartphone_2_id)), ~ smartphone_1_ti_di_1+smartphone_2_ti_di_1, addNA=TRUE)) 
master1long.g1g1sm.ti <- master1longdups %>%
  filter(field %in% c("smartphone_1_ti_di_1","smartphone_2_ti_di_1")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g1sm.ti, method = "rm", ci.type = "rm")

# Grader 1 vs grader 1 (SLR)
addmargins(xtabs(data=filter(master1, !is.na(SLR_2_id)), ~ SLR_1_ti_di_1+SLR_2_ti_di_1, addNA=TRUE)) 
master1long.g1g1slr.ti <- master1longdups %>%
  filter(field %in% c("SLR_1_ti_di_1","SLR_2_ti_di_1")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g1g1slr.ti, method = "rm", ci.type = "rm")

# Grader 2 vs grader 2 (smartphone)
addmargins(xtabs(data=filter(master1, !is.na(smartphone_2_id)), ~ smartphone_1_ti_di_2+smartphone_2_ti_di_2, addNA=TRUE)) 
master1long.g2g2sm.ti <- master1longdups %>%
  filter(field %in% c("smartphone_1_ti_di_2","smartphone_2_ti_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g2g2sm.ti, method = "rm", ci.type = "rm")

# Grader 2 vs grader 2 (SLR)
addmargins(xtabs(data=filter(master1, !is.na(SLR_2_id)), ~ SLR_1_ti_di_2+SLR_2_ti_di_2, addNA=TRUE)) 
master1long.g2g2slr.ti <- master1longdups %>%
  filter(field %in% c("SLR_1_ti_di_2","SLR_2_ti_di_2")) %>% 
  mutate(value=as.numeric(value)) %>% 
  dplyr::select(-field)
iccbin(cid=number, y=value, data=master1long.g2g2slr.ti, method = "rm", ci.type = "rm")

###################################
##  SENS/SPEC w/ field grades   ##
###################################
# ETABLE 2, SENS/SPEC RELATIVE TO FIELD GRADES
#JN comment: I don't think we will include this in the final manuscript since we are calculating all of these metrics relative to the LCA
library(yardstick)
options(yardstick.event_first = FALSE)
class_metrics <- metric_set(sens, spec, ppv, npv) # JN added ppv and npv
# Basic yardstick command, SLR compared to field grades reference standard
dx_estimates_tf <- master1 %>%
  class_metrics(., truth = factor(clinic_tf_yn), estimate = factor(SLR_1_tf_yn))

# BOOTSTRAPPING THE CONFIDENCE INTERVALS
# First, set up the bootstraps
library(rsample)
library(coxed) # This is for the bca confidence intervals, but it masks the "summarize" command so that means you have to use dplyr:: beforehand
library(purrr)
D <- master1 %>% 
  dplyr::select(number, age, sex, state_code, age, SLR_1_tf_yn,smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn,smartphone_1_ti_yn, clinic_ti_yn, examiner) %>%
  filter(!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)) %>% 
  nest(-examiner)
set.seed(154234)
# The bs object is the boostrap object; we are creating separate populations with resampling
bs_nested <- bootstraps(D, times = 9)
bs_unnested <- bootstraps(master1, times = 9)

# SLR TF
bs_sensspec_sl.tf <- map(bs_nested$splits, ~as_tibble(.) %>% unnest %>% 
                           # group_by(camera) %>% 
                           class_metrics(., truth = factor(clinic_tf_yn), estimate = factor(SLR_1_tf_yn))) %>% 
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=TRUE)

sensspec_estimates_tf <- master1 %>%
  class_metrics(., truth = factor(clinic_tf_yn), estimate = factor(SLR_1_tf_yn))
sensspec_estimates_tf
coxed::bca(bs_sensspec_sl.tf$sens)
coxed::bca(bs_sensspec_sl.tf$spec)
coxed::bca(bs_sensspec_sl.tf$ppv)
coxed::bca(bs_sensspec_sl.tf$npv)

# SLR TI
bs_sensspec_sl.ti <- map(bs_nested$splits, ~as_tibble(.) %>% unnest %>% 
                           # group_by(camera) %>% 
                           class_metrics(., truth = factor(clinic_ti_yn), estimate = factor(SLR_1_ti_yn))) %>% 
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=TRUE)
sensspec_estimates_ti <- master1 %>%
  class_metrics(., truth = factor(clinic_ti_yn), estimate = factor(SLR_1_ti_yn))
sensspec_estimates_ti
coxed::bca(bs_sensspec_sl.ti$sens)
coxed::bca(bs_sensspec_sl.ti$spec)

# Smartphone TF
bs_sensspec_sm.tf <- map(bs_nested$splits, ~as_tibble(.) %>% unnest %>% 
                           # group_by(camera) %>% 
                           class_metrics(., truth = factor(clinic_tf_yn), estimate = factor(smartphone_1_tf_yn))) %>% 
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=TRUE)
sensspec_estimates_tf <- master1 %>%
  class_metrics(., truth = factor(clinic_tf_yn), estimate = factor(smartphone_1_tf_yn))
sensspec_estimates_tf
coxed::bca(bs_sensspec_sm.tf$sens)
coxed::bca(bs_sensspec_sm.tf$spec)

# Smartphone TI
bs_sensspec_sm.ti <- map(bs_nested$splits, ~as_tibble(.) %>% unnest %>% 
                           # group_by(camera) %>% 
                           class_metrics(., truth = factor(clinic_ti_yn), estimate = factor(smartphone_1_ti_yn))) %>% 
  bind_rows(.id = 'boots') %>%
  dplyr::select(-.estimator) %>%
  spread(.metric, .estimate, convert=TRUE)
sensspec_estimates_ti <- master1 %>%
  class_metrics(., truth = factor(clinic_ti_yn), estimate = factor(smartphone_1_ti_yn))
sensspec_estimates_ti
coxed::bca(bs_sensspec_sm.ti$sens)
coxed::bca(bs_sensspec_sm.ti$spec)

################################
##             LCA            ##
################################
# There are at least 2 packages for LCA; experimenting with both
library(poLCA)
library(randomLCA) # This one allows random effect for subject, good for diagnostic tests with multiple graders
# Note we could have used that here actually but we did a consensus grade so not necessary
# The following is an example for learning...
# # FIRST, an example without random effects (myocardial)
# # Note below that myocardial dataframe has 4 dichotomous variables,
# # and the fifth variable is freq=the frequency of each from a rxc contingency table
# # nclass is how many classes to fit
# myocardial.lca1 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq, nclass = 1)
# myocardial.lca2 <- randomLCA(myocardial[, 1:4], freq = myocardial$freq, nclass = 2)
# myocardial.bic <- data.frame(classes = 1:2, bic = c(BIC(myocardial.lca1), BIC(myocardial.lca2)))
# print(myocardial.bic, row.names = FALSE)
# # BIC lower for 2-class model, (so if use BIC as a selection method, then selects the 2 class model, 
# # indicating a breakdown into diseased and non-diseased
# # Alternatively simulate with parametric bootstrap to figure out number of classes
# nsims <- 999
# obslrt <- 2 * (logLik(myocardial.lca2) - logLik(myocardial.lca1))
# thesims <- simulate(myocardial.lca1, nsim = nsims)
# simlrt <- as.vector(lapply(thesims, function(x) {
#   submodel <- refit(myocardial.lca1, newpatterns = x)
#   fullmodel <- refit(myocardial.lca2, newpatterns = x)
#   return(2 * (logLik(fullmodel) - logLik(submodel)))
#   }))
# print((sum(simlrt >= obslrt) + 1)/(nsims + 1)) # Gives a p-value, not entirely sure but he says significant p-value means you should use a 2-class model
# summary(myocardial.lca2) # Shows that class 2 is disease and class 1 is not diseased
# outcomeProbs(myocardial.lca2) # For individual outcome probabilities
# # Note that Q.wave in class 1 has bad CIs. The parametric bootstrap with boot = TRUE will produce improved results 
# # or alternatively the value of the penalty argument could be increased.
# probs <- outcomeProbs(myocardial.lca2, boot = TRUE)
# # For sensitivity if class 1 is positive disease
# probs[[1]]
# # For specificity if class 2 is negative disease
# 1-probs[[2]]
# # His code for sens/spec; this is the same as output above
# diseased <- ifelse(probs[[1]]$Outcome[1] < probs[[2]]$Outcome[1], 2, 1)
# notdiseased <- 3 - diseased
# sens <- apply(probs[[diseased]], 1, function(x) sprintf("%3.2f (%3.2f, %3.2f)", x[1], x[2], x[3]))
# spec <- apply(probs[[notdiseased]], 1, function(x) sprintf("%3.2f (%3.2f, %3.2f)", 1 - x[1], 1 - x[3], 1 - x[2]))
# stable <- data.frame(sens, spec)
# names(stable) <- c("Sensitivity", "Specificity")
# print(stable, row.names = TRUE)
# # If want the probability of disease/not disease for each combination of tests:
# print(postClassProbs(myocardial.lca2), row.names = FALSE)
# plot(myocardial.lca2, type = "b", pch = 1:2, xlab = "Test",
#      ylab = "Outcome Probability",
#      scales = list(x = list(at = 1:4, labels = names(myocardial)[1:4])),
#      key = list(corner = c(0.05, .95), border = TRUE, cex = 1.2,
#                   text = list(c("Class 1", "Class 2")),
#                   col = trellis.par.get()$superpose.symbol$col[1:2],
#                   points = list(pch = 1:2)))
# # SECOND, an example with random effects (dentistry)
# # Here, 5 different dentists have graded presence/abence of caries
# dentistry.lca1       <- randomLCA(dentistry[, 1:5], freq = dentistry$freq, nclass = 1)
# dentistry.lca1       <- randomLCA(dentistry[,"V1":"V5"], freq = dentistry$freq, nclass = 1)
# 
# dentistry.lca1random <- randomLCA(dentistry[, 1:5], freq = dentistry$freq, nclass = 1, random = TRUE, probit = TRUE)
# dentistry.lca1random2 <- randomLCA(dentistry[, 1:5],freq = dentistry$freq, nclass = 1, random = TRUE, probit = TRUE, constload = FALSE)
# dentistry.lca2random <- randomLCA(dentistry[, 1:5], freq = dentistry$freq, nclass = 2, random = TRUE, quadpoints = 71, probit = TRUE)
# probs <- outcomeProbs(dentistry.lca1random2, boot = TRUE)
# probs

# FOR TRACHOMA DATA
tfxtab <- as.data.frame(xtabs(data=master1, ~clinic_tf_yn+SLR_1_tf_yn+smartphone_1_tf_yn))
tf.lca1 <- randomLCA(tfxtab[, 1:3], freq = tfxtab$Freq, nclass = 1)                   
tf.lca2 <- randomLCA(tfxtab[, 1:3], freq = tfxtab$Freq, nclass = 2)
tf.bic <- data.frame(classes = 1:2, bic = c(BIC(tf.lca1), BIC(tf.lca2)))
print(tf.bic, row.names = FALSE)
summary(tf.lca2) # Suggests class 2 is disease (TF)
# However I later learned that if you repeat the randomLCA command, sometimes class 1 is disease and sometimes clas 2 is disease
# But this randomLCA nice because gives confidence intervals for sens/spec
print(postClassProbs(tf.lca2), row.names = FALSE)
plot(tf.lca2, type = "b", pch = 1:2, xlab = "Test",
     ylab = "Outcome Probability",
     scales = list(x = list(at = 1:4, labels = names(tfxtab)[1:4])),
     key = list(corner = c(0.05, .95), border = TRUE, cex = 1.2,
                text = list(c("Class 1", "Class 2")),
                col = trellis.par.get()$superpose.symbol$col[1:2],
                points = list(pch = 1:2)))
tfprobs <- outcomeProbs(tf.lca2, boot = TRUE, R = 9999)
tfprobs
# For sensitivity if class 2 is positive disease
tfprobs[[2]]
# For specificity if class 1 is negative disease
1-tfprobs[[1]]
#JN: maybe these change with increasing bootstrap samples but none of these results match with the sens/spec reported in the paper

# TRY REPRODUCING WITH poLCA (DOWNSIDES: HARDER TO TELL WHICH CLASS IS DISEASE AND NO CIs; UPSIDE: EASIER FOR BOOTSTRAPPING THE DIFFERENCES)
tfdata <- master1 %>%
  dplyr::select(clinic_tf_yn, SLR_1_tf_yn, smartphone_1_tf_yn) %>%
  mutate(clinic_tf_yn=clinic_tf_yn+1,
         SLR_1_tf_yn=SLR_1_tf_yn+1,
         smartphone_1_tf_yn=smartphone_1_tf_yn+1) # This is because poLCA needs 1/2 variable, not 0/1
tf.form <- cbind(clinic_tf_yn,SLR_1_tf_yn,smartphone_1_tf_yn)~1
tf.polca2 <- poLCA(tf.form,tfdata,nclass=2) # This gives the same result; look at the Pr(2) column
tf.polca2.probs.start.new <- poLCA.reorder(tf.polca2$probs.start,order(tf.polca2$P,decreasing=TRUE))
tf.polca2 <- poLCA(tf.form, tfdata, nclass=2, verbose=FALSE, probs.start=tf.polca2.probs.start.new)
tf.polca2 # The "tf.polca2.probs.start.new" makes reproducible results; in this case that 
# class 2 is disease (hard to know from this input, easier to tell by matching the numbers from the randomLCA output)
# these match with what I find in the paper
D <- master1 %>% 
  dplyr::select(number, age, sex, state_code, age, SLR_1_tf_yn,smartphone_1_tf_yn, clinic_tf_yn, SLR_1_ti_yn,smartphone_1_ti_yn, clinic_ti_yn, examiner) %>%
  filter(!is.na(SLR_1_tf_yn) & !is.na(smartphone_1_tf_yn) & !is.na(clinic_tf_yn)) %>% 
  nest(-examiner)
set.seed(12346)
# The bs object is the boostrap object; we are creating separate populations with resampling
bs_nested <- bootstraps(D, times = 9999)
bs_unnested <- bootstraps(master1, times = 9999)

# BOOTSTRAP CIs FOR SENS/SPEC DIFFERENCES BASED ON LCA 
# First make a formula
lca.fx <- function(d) {
  tf.form <- cbind(clinic_tf_yn,SLR_1_tf_yn,smartphone_1_tf_yn)~1
  d2 <- d %>% 
    dplyr::select(clinic_tf_yn,SLR_1_tf_yn,smartphone_1_tf_yn) %>%
    dplyr::mutate(clinic_tf_yn=(clinic_tf_yn+1),
                  SLR_1_tf_yn=(SLR_1_tf_yn+1),
                  smartphone_1_tf_yn=(smartphone_1_tf_yn+1))
  tf.polca <- poLCA(tf.form, d2, nclass=2, verbose=FALSE)
  # Note that the classes are arbitrary, assigned by initial values of EM algorithm.
  # So sometimes class 1=disease and sometimes class 2 = disease
  # So this next line defines starting values so we get reproducible results
  # Especially important for the bootstrapping
  tf.polca.probs.start.new <- poLCA.reorder(tf.polca$probs.start,order(tf.polca$P,decreasing=TRUE))
  tf.polca <- poLCA(tf.form, d2, nclass=2, verbose=FALSE, probs.start=tf.polca.probs.start.new)
  tf.polca.sensspec <- as.data.frame(tf.polca$probs) %>%
    rownames_to_column() %>%
    dplyr::mutate(rowname=if_else(grepl("class 2:", rowname), "sens", "spec"),
                  sens.slr.minus.clin=if_else(rowname=="sens", (SLR_1_tf_yn.Pr.2. - clinic_tf_yn.Pr.2.), NA_real_ ),
                  sens.slr.minus.smart=if_else(rowname=="sens", (SLR_1_tf_yn.Pr.2. - smartphone_1_tf_yn.Pr.2.), NA_real_ ),
                  sens.clin.minus.smart=if_else(rowname=="sens", (clinic_tf_yn.Pr.2. - smartphone_1_tf_yn.Pr.2.), NA_real_ ),
                  spec.slr.minus.clin=if_else(rowname=="spec", (SLR_1_tf_yn.Pr.1. - clinic_tf_yn.Pr.1.), NA_real_ ),
                  spec.slr.minus.smart=if_else(rowname=="spec", (SLR_1_tf_yn.Pr.1. - smartphone_1_tf_yn.Pr.1.), NA_real_ ),
                  spec.clin.minus.smart=if_else(rowname=="spec", (clinic_tf_yn.Pr.1. - smartphone_1_tf_yn.Pr.1.), NA_real_ )) %>%
    dplyr::summarize(sens.slr.minus.clin=min(sens.slr.minus.clin, na.rm=TRUE),
                     sens.slr.minus.smart=min(sens.slr.minus.smart, na.rm=TRUE),
                     sens.clin.minus.smart=min(sens.clin.minus.smart, na.rm=TRUE),
                     spec.slr.minus.clin=min(spec.slr.minus.clin, na.rm=TRUE),
                     spec.slr.minus.smart=min(spec.slr.minus.smart, na.rm=TRUE),
                     spec.clin.minus.smart=min(spec.clin.minus.smart, na.rm=TRUE))
  tf.polca.sensspec
} 

lca.fx(master1) # function works
# Making sure it works within a pipe...
master1 %>% lca.fx(.)

bs_lca <- map(bs_nested$splits, ~as_tibble(.) %>% 
                unnest %>%  
                lca.fx(.)) %>%
  # Then bind together the results (the single row from the summarize above)
  bind_rows(.id = 'boots') 
# So finally the mean difference and CIs for TABLE 1:
master1 %>% lca.fx(.)
bca(bs_lca$sens.slr.minus.clin)
bca(bs_lca$sens.slr.minus.smart)
bca(bs_lca$sens.clin.minus.smart)
bca(bs_lca$spec.slr.minus.clin)
bca(bs_lca$spec.slr.minus.smart)
bca(bs_lca$spec.clin.minus.smart)

#JN: Here is what I did to calculate pdlr and ndlr and their BCas using the poLCA and DTComPair packages
library(DTComPair)
#JN this package also allows for comparison between sens/spec/plr/nlr

lca_tf_df <- master1 %>%
  dplyr::select(SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn) %>%
  mutate(SLR_1_tf_yn=SLR_1_tf_yn+1,
         smartphone_1_tf_yn=smartphone_1_tf_yn+1,
         clinic_tf_yn=clinic_tf_yn+1)
lca_tf_model1 = poLCA(cbind(SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn) ~ 1, nclass = 2, data = lca_tf_df)
probs.start <- lca_tf_model1$probs.start
new.probs.start <- poLCA.reorder(lca_tf_model1$probs.start,order(lca_tf_model1$P,decreasing=TRUE)) #here class 1 (+TF) becomes 2 and class 2 (-TF) becomes 1
lca_tf_model1 = poLCA(cbind(SLR_1_tf_yn, smartphone_1_tf_yn, clinic_tf_yn) ~ 1, nclass = 2, data = lca_tf_df, probs.start = new.probs.start)
#this gives the predicted class
lca_tf_model1$predclass

tf_with_lc <- lca_tf_df %>%
  cbind(lca_tf_model1$predclass) %>%
  rename(predclass = 'lca_tf_model1$predclass') %>%
  #reformating to binary
  mutate(SLR_1_tf_yn=SLR_1_tf_yn-1,
         smartphone_1_tf_yn=smartphone_1_tf_yn-1,
         clinic_tf_yn=clinic_tf_yn-1,
         predclass=predclass-1)

#Calculating PLR/NLR
#PLR sens/1-spec
#NLR 1-sens/spec
# SLR TF vs LC
slr_tf_tab <- tab.1test(predclass, SLR_1_tf_yn, data=tf_with_lc)     
acc.1test(slr_tf_tab, alpha = 0.05)
#smart TF vs LC
smart_tf_tab <- tab.1test(predclass, smartphone_1_tf_yn, data=tf_with_lc)
acc.1test(smart_tf_tab, alpha = 0.05)
#clinic TF vs LC
clinic_tf_tab <- tab.1test(predclass, clinic_tf_yn, data=tf_with_lc)
acc.1test(clinic_tf_tab, alpha = 0.05)
#function to get PDLR and NDLR estimates

tf.lr <- function(d) {
tf.data <- d %>%
  dplyr::select(clinic_tf_yn, SLR_1_tf_yn, smartphone_1_tf_yn)
d2 <- d %>% 
  dplyr::select(clinic_tf_yn, SLR_1_tf_yn, smartphone_1_tf_yn) %>%
  dplyr::mutate(clinic_tf_yn=(clinic_tf_yn + 1),
                SLR_1_tf_yn=(SLR_1_tf_yn + 1),
                smartphone_1_tf_yn=(smartphone_1_tf_yn + 1))
tf.polca.model <- poLCA(cbind(clinic_tf_yn, SLR_1_tf_yn, smartphone_1_tf_yn) ~ 1, nclass = 2, data = d2, verbose = F)
new.probs.start <- poLCA.reorder(tf.polca.model$probs.start, order(tf.polca.model$P, decreasing = T))
tf.polca.model <- poLCA(cbind(clinic_tf_yn, SLR_1_tf_yn, smartphone_1_tf_yn) ~ 1, nclass = 2, data = d2, verbose = F, probs.start=new.probs.start)
tf.lc <- d2 %>%
  cbind(tf.polca.model$predclass) %>%
  rename(predclass = 'tf.polca.model$predclass') %>%
  #reformating to binary
  mutate(clinic_tf_yn=clinic_tf_yn-1,
         SLR_1_tf_yn=SLR_1_tf_yn-1,
         smartphone_1_tf_yn=smartphone_1_tf_yn-1,
         predclass=predclass-1)
#calculating estimates for each
slr_tf_tab <- tab.1test(predclass, SLR_1_tf_yn, data=tf.lc)     
slr <- acc.1test(slr_tf_tab, alpha = 0.05)
smart_tf_tab <- tab.1test(predclass, smartphone_1_tf_yn, data=tf.lc)
smart <- acc.1test(smart_tf_tab, alpha = 0.05)
clinic_tf_tab <- tab.1test(predclass, clinic_tf_yn, data=tf.lc)
clinic <- acc.1test(clinic_tf_tab, alpha = 0.05)
#subsetting the PDLR estimates
slr.pdlr <- slr$pdlr[[1]]
smart.pdlr <- smart$pdlr[[1]]
clinic.pdlr <- clinic$pdlr[[1]]
#subsetting the NDLR estimates
slr.ndlr <- slr$ndlr[[1]]
smart.ndlr <- smart$ndlr[[1]]
clinic.ndlr <- clinic$ndlr[[1]]
#merging the estimates into one data table
pdlr.ndlr.estimates <- as_tibble(cbind(slr.pdlr, smart.pdlr, clinic.pdlr, slr.ndlr, smart.ndlr, clinic.ndlr)) %>%
  mutate(pdlr.slr.minus.clinic=slr.pdlr-clinic.pdlr,
         pdlr.slr.minus.smart=slr.pdlr-smart.pdlr,
         pdlr.clinic.minus.smart=clinic.pdlr-smart.pdlr,
         ndlr.slr.minus.clinic=slr.ndlr-clinic.ndlr,
         ndlr.slr.minus.smart=slr.ndlr-smart.ndlr,
         ndlr.clinic.minus.smart=clinic.ndlr-smart.ndlr)
pdlr.ndlr.estimates 
}
view(master1 %>%tf.lr(.))
#bootstrapping
bs_lca_lr <- map(bs_nested$splits, ~as_tibble(.) %>% 
                        unnest %>%  
                        tf.lr(.)) %>%
  bind_rows(.id = 'boots') 
#BCas
bca(bs_lca_lr$slr.pdlr) # output is NaN NaN
bca(is.finite(bs_lca_lr$slr.pdlr)) #now 0 to 0...
#for now will use the CIs given by acc.1test
bca(bs_lca_lr$smart.pdlr)
bca(bs_lca_lr$clinic.pdlr)
bca(bs_lca_lr$slr.ndlr)
bca(bs_lca_lr$smart.ndlr)
bca(bs_lca_lr$clinic.ndlr)
#differeces
bca(bs_lca_lr$pdlr.slr.minus.clinic)
bs_plr_slr_min_clinic <- bs_lca_lr %>%
  dplyr::select(pdlr.slr.minus.clinic) %>%
  filter(is.finite(pdlr.slr.minus.clinic)) #gets rid of 242 vallues
bca(bs_plr_slr_min_clinic$pdlr.slr.minus.clinic) 
#alternative to filtering these out is write a function that redraws a boot until all 9999 are finite
bca(bs_lca_lr$pdlr.slr.minus.smart) #lots of NaN and Inf numbers
bs_plr_slr_min_smart <- bs_lca_lr %>%
  dplyr::select(pdlr.slr.minus.smart) %>%
  filter(is.finite(pdlr.slr.minus.smart)) #gets rid of 640 values
bca(bs_plr_slr_min_smart$pdlr.slr.minus.smart)

bca(bs_lca_lr$pdlr.clinic.minus.smart)
bs_plr_clin_min_smart <- bs_lca_lr %>%
  dplyr::select(pdlr.clinic.minus.smart) %>%
  filter(is.finite(pdlr.clinic.minus.smart)) #removes 630
bca(bs_plr_clin_min_smart$pdlr.clinic.minus.smart)

bca(bs_lca_lr$ndlr.slr.minus.clinic)
bca(bs_lca_lr$ndlr.slr.minus.smart)
bca(bs_lca_lr$ndlr.clinic.minus.smart)

#doing for TI
LCAti_df <- master1 %>%
  dplyr::select(SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) %>%
  mutate(SLR_1_ti_yn=SLR_1_ti_yn+1,
         smartphone_1_ti_yn=smartphone_1_ti_yn+1,
         clinic_ti_yn=clinic_ti_yn+1)
LCA_ti_model1 = poLCA(cbind(SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) ~ 1, nclass = 2, data = LCAti_df)
probs.start <- LCA_ti_model1$probs.start
new.probs.start <- poLCA.reorder(LCA_ti_model1$probs.start,order(LCA_ti_model1$P,decreasing=TRUE)) #here class 1 (+TF) becomes 2 and class 2 (-TF) becomes 1
LCA_ti_model1 = poLCA(cbind(SLR_1_ti_yn, smartphone_1_ti_yn, clinic_ti_yn) ~ 1, nclass = 2, data = LCAti_df, probs.start = new.probs.start)
#this gives the predicted class
LCA_ti_model1$predclass

ti_with_lc <- LCAti_df %>%
  cbind(LCA_ti_model1$predclass) %>%
  rename(predclass = 'LCA_ti_model1$predclass') %>%
  #reformating to binary
  mutate(SLR_1_ti_yn=SLR_1_ti_yn-1,
         smartphone_1_ti_yn=smartphone_1_ti_yn-1,
         clinic_ti_yn=clinic_ti_yn-1,
         predclass=predclass-1) #with this I should be able to calculate sens/spec with bca using the yardstick backage
#SLR TI vs LC
slr_ti_tab <- tab.1test(predclass, SLR_1_ti_yn, data=ti_with_lc)     
acc.1test(slr_ti_tab, alpha = 0.05) #sens and NPV = 1 and NDLR =0...
#smart TI vs LC
smart_ti_tab <- tab.1test(predclass, smartphone_1_ti_yn, data=ti_with_lc)
acc.1test(smart_ti_tab, alpha = 0.05)
#clinic TI vs LC
clinic_ti_tab <- tab.1test(predclass, clinic_ti_yn, data=ti_with_lc)
acc.1test(clinic_ti_tab, alpha = 0.05)

#creating a function for point estimates and differences between plr and nlr
ti.lr <- function(d) {
  ti.data <- d %>%
    dplyr::select(clinic_ti_yn, SLR_1_ti_yn, smartphone_1_ti_yn)
  d2 <- d %>% 
    dplyr::select(clinic_ti_yn, SLR_1_ti_yn, smartphone_1_ti_yn) %>%
    dplyr::mutate(clinic_ti_yn=(clinic_ti_yn + 1),
                  SLR_1_ti_yn=(SLR_1_ti_yn + 1),
                  smartphone_1_ti_yn=(smartphone_1_ti_yn + 1))
  ti.polca.model <- poLCA(cbind(clinic_ti_yn, SLR_1_ti_yn, smartphone_1_ti_yn) ~ 1, nclass = 2, data = d2, verbose = F)
  new.probs.start <- poLCA.reorder(ti.polca.model$probs.start, order(ti.polca.model$P, decreasing = T))
  ti.polca.model <- poLCA(cbind(clinic_ti_yn, SLR_1_ti_yn, smartphone_1_ti_yn) ~ 1, nclass = 2, data = d2, verbose = F, probs.start=new.probs.start)
  ti.lc <- d2 %>%
    cbind(ti.polca.model$predclass) %>%
    rename(predclass = 'ti.polca.model$predclass') %>%
    #reformating to binary
    mutate(clinic_ti_yn=clinic_ti_yn-1,
           SLR_1_ti_yn=SLR_1_ti_yn-1,
           smartphone_1_ti_yn=smartphone_1_ti_yn-1,
           predclass=predclass-1)
  #calculating estimates for each
  slr_ti_tab <- tab.1test(predclass, SLR_1_ti_yn, data=ti.lc)     
  slr <- acc.1test(slr_ti_tab, alpha = 0.05)
  smart_ti_tab <- tab.1test(predclass, smartphone_1_ti_yn, data=ti.lc)
  smart <- acc.1test(smart_ti_tab, alpha = 0.05)
  clinic_ti_tab <- tab.1test(predclass, clinic_ti_yn, data=ti.lc)
  clinic <- acc.1test(clinic_ti_tab, alpha = 0.05)
  #subsetting the PDLR estimates
  slr.pdlr <- slr$pdlr[[1]]
  smart.pdlr <- smart$pdlr[[1]]
  clinic.pdlr <- clinic$pdlr[[1]]
  #subsetting the NDLR estimates
  slr.ndlr <- slr$ndlr[[1]]
  smart.ndlr <- smart$ndlr[[1]]
  clinic.ndlr <- clinic$ndlr[[1]]
  #merging the estimates into one data table
  pdlr.ndlr.estimates <- as_tibble(cbind(slr.pdlr, smart.pdlr, clinic.pdlr, slr.ndlr, smart.ndlr, clinic.ndlr)) %>%
    mutate(pdlr.slr.minus.clinic=slr.pdlr-clinic.pdlr,
           pdlr.slr.minus.smart=slr.pdlr-smart.pdlr,
           pdlr.clinic.minus.smart=clinic.pdlr-smart.pdlr,
           ndlr.slr.minus.clinic=slr.ndlr-clinic.ndlr,
           ndlr.slr.minus.smart=slr.ndlr-smart.ndlr,
           ndlr.clinic.minus.smart=clinic.ndlr-smart.ndlr)
  pdlr.ndlr.estimates 
}
view(master1 %>%ti.lr(.))
#bootstrapping
bs_lca_lr_ti <- map(bs_nested$splits, ~as_tibble(.) %>% 
                   unnest %>%  
                   ti.lr(.)) %>%
  bind_rows(.id = 'boots') 
#BCas
bca(bs_lca_lr_ti$slr.pdlr)
bca(bs_lca_lr_ti$smart.pdlr)
bca(bs_lca_lr_ti$clinic.pdlr)
bca(bs_lca_lr_ti$slr.ndlr)
bca(bs_lca_lr_ti$smart.ndlr)
bca(bs_lca_lr_ti$clinic.ndlr)

bca(bs_lca_lr_ti$pdlr.slr.minus.clinic)
pdlr_slr_min_clinic <- bs_lca_lr_ti %>%
  dplyr::select(pdlr.slr.minus.clinic) %>%
  filter(is.finite(pdlr.slr.minus.clinic)) #removes 1097
bca(pdlr_slr_min_clinic$pdlr.slr.minus.clinic)

bca(bs_lca_lr_ti$pdlr.slr.minus.smart)
pdlr_slr_min_smart <- bs_lca_lr_ti %>%
  dplyr::select(pdlr.slr.minus.smart) %>%
  filter(is.finite(pdlr.slr.minus.smart)) #removes 857
bca(pdlr_slr_min_smart$pdlr.slr.minus.smart)

bca(bs_lca_lr_ti$pdlr.clinic.minus.smart)
pdlr_clin_min_smart <- bs_lca_lr_ti %>%
  dplyr::select(pdlr.clinic.minus.smart) %>%
  filter(is.finite(pdlr.clinic.minus.smart)) #removes 323
bca(pdlr_clin_min_smart$pdlr.clinic.minus.smart)

bca(bs_lca_lr_ti$ndlr.slr.minus.clinic)
bca(bs_lca_lr_ti$ndlr.slr.minus.smart)
bca(bs_lca_lr_ti$ndlr.clinic.minus.smart)

#Resume JK's code
# And for TI:
# FOR TRACHOMA DATA
tixtab <- as.data.frame(xtabs(data=master1, ~clinic_ti_yn+SLR_1_ti_yn+smartphone_1_ti_yn))
ti.lca1 <- randomLCA(tixtab[, 1:3], freq = tixtab$Freq, nclass = 1)                   
ti.lca2 <- randomLCA(tixtab[, 1:3], freq = tixtab$Freq, nclass = 2)
ti.bic <- data.frame(classes = 1:2, bic = c(BIC(ti.lca1), BIC(ti.lca2)))
print(ti.bic, row.names = FALSE)
summary(ti.lca2) # Suggests class 2 is disease (ti)
# However I later learned that if you repeat the randomLCA command, sometimes class 1 is disease and sometimes clas 2 is disease
# But this randomLCA nice because gives confidence intervals for sens/spec
tiprobs <- outcomeProbs(ti.lca2, boot = TRUE, R = 9999)
tiprobs
# For sensitivity if class 2 is positive disease
tiprobs[[2]]
# For specificity if class 1 is negative disease
1-tiprobs[[1]]

# BOOTSTRAP CIs FOR SENS/SPEC DIFFERENCES BASED ON LCA 
# First make a formula
lca.fx.ti <- function(d) {
  ti.form <- cbind(clinic_ti_yn,SLR_1_ti_yn,smartphone_1_ti_yn)~1
  d2 <- d %>% 
    dplyr::select(clinic_ti_yn,SLR_1_ti_yn,smartphone_1_ti_yn) %>%
    dplyr::mutate(clinic_ti_yn=(clinic_ti_yn+1),
                  SLR_1_ti_yn=(SLR_1_ti_yn+1),
                  smartphone_1_ti_yn=(smartphone_1_ti_yn+1))
  ti.polca <- poLCA(ti.form, d2, nclass=2, verbose=FALSE)
  # Note that the classes are arbitrary, assigned by initial values of EM algorithm.
  # So sometimes class 1=disease and sometimes class 2 = disease
  # So this next line defines starting values so we get reproducible results
  # Especially important for the bootstrapping
  ti.polca.probs.start.new <- poLCA.reorder(ti.polca$probs.start,order(ti.polca$P,decreasing=TRUE))
  ti.polca <- poLCA(ti.form, d2, nclass=2, verbose=FALSE, probs.start=ti.polca.probs.start.new)
  ti.polca.sensspec <- as.data.frame(ti.polca$probs) %>%
    rownames_to_column() %>%
    dplyr::mutate(rowname=if_else(grepl("class 2:", rowname), "sens", "spec"),
                  sens.slr.minus.clin=if_else(rowname=="sens", (SLR_1_ti_yn.Pr.2. - clinic_ti_yn.Pr.2.), NA_real_ ),
                  sens.slr.minus.smart=if_else(rowname=="sens", (SLR_1_ti_yn.Pr.2. - smartphone_1_ti_yn.Pr.2.), NA_real_ ),
                  sens.clin.minus.smart=if_else(rowname=="sens", (clinic_ti_yn.Pr.2. - smartphone_1_ti_yn.Pr.2.), NA_real_ ),
                  spec.slr.minus.clin=if_else(rowname=="spec", (SLR_1_ti_yn.Pr.1. - clinic_ti_yn.Pr.1.), NA_real_ ),
                  spec.slr.minus.smart=if_else(rowname=="spec", (SLR_1_ti_yn.Pr.1. - smartphone_1_ti_yn.Pr.1.), NA_real_ ),
                  spec.clin.minus.smart=if_else(rowname=="spec", (clinic_ti_yn.Pr.1. - smartphone_1_ti_yn.Pr.1.), NA_real_ )) %>%
    dplyr::summarize(sens.slr.minus.clin=min(sens.slr.minus.clin, na.rm=TRUE),
                     sens.slr.minus.smart=min(sens.slr.minus.smart, na.rm=TRUE),
                     sens.clin.minus.smart=min(sens.clin.minus.smart, na.rm=TRUE),
                     spec.slr.minus.clin=min(spec.slr.minus.clin, na.rm=TRUE),
                     spec.slr.minus.smart=min(spec.slr.minus.smart, na.rm=TRUE),
                     spec.clin.minus.smart=min(spec.clin.minus.smart, na.rm=TRUE))
  ti.polca.sensspec
} 
lca.fx.ti(master1) # function works
bs_lca_ti <- map(bs_nested$splits, ~as_tibble(.) %>% 
                   unnest %>%  
                   lca.fx.ti(.)) %>%
  # Then bind together the results (the single row from the summarize above)
  bind_rows(.id = 'boots') 
# So finally the mean difference and CIs for TABLE 1:
master1 %>% lca.fx.ti(.)
bca(bs_lca_ti$sens.slr.minus.clin)
bca(bs_lca_ti$sens.slr.minus.smart)
bca(bs_lca_ti$sens.clin.minus.smart)
bca(bs_lca_ti$spec.slr.minus.clin)
bca(bs_lca_ti$spec.slr.minus.smart)
bca(bs_lca_ti$spec.clin.minus.smart)

#slr smart tf
slr.smart.tf.tab <- tab.paired(predclass, SLR_1_tf_yn, smartphone_1_tf_yn, data=tf_with_lc) #test1 = slr, test2= smart
#comparing sens/spec
sesp.mcnemar(slr.smart.tf.tab) #Mcnemar chi-squared comparison
sesp.exactbinom(slr.smart.tf.tab) #Exact binom method for comparison
#comparing dlrp/dlrn
dlr.regtest(slr.smart.tf.tab, alpha = 0.05) #uses regression model approach proposed by Gu and Pepe (2009)

#slr clinic tf
slr.clinic.tf.tab <- tab.paired(predclass, SLR_1_tf_yn, clinic_tf_yn, data=tf_with_lc) #test1 = slr, test2= clinic
#comparing sens/spec
sesp.mcnemar(slr.clinic.tf.tab) 
sesp.exactbinom(slr.clinic.tf.tab)
#comparing dlrp/dlrn
dlr.regtest(slr.clinic.tf.tab, alpha = 0.05)

#smart clinic tf
smart.clinic.tf.tab <- tab.paired(predclass, smartphone_1_tf_yn, clinic_tf_yn, data=tf_with_lc) #test1 = smart, test2= clinic
#comparing sens/spec
sesp.mcnemar(smart.clinic.tf.tab)
sesp.exactbinom(smart.clinic.tf.tab)
#comparing dlrp/dlrn
dlr.regtest(smart.clinic.tf.tab, alpha = 0.05)

#slr smart ti
slr.smart.ti.tab <- tab.paired(predclass, SLR_1_ti_yn, smartphone_1_ti_yn, data=ti_with_lc) #test1 = slr, test2= smart
#comparing sens/spec
sesp.mcnemar(slr.smart.ti.tab)
sesp.exactbinom(slr.smart.ti.tab)
#comparing dlrp/dlrn
dlr.regtest(slr.smart.ti.tab)

#slr clinic ti
slr.clinic.ti.tab <- tab.paired(predclass, SLR_1_ti_yn, clinic_ti_yn, data=ti_with_lc) #test1 = slr, test2= clinic
#comparing sens/spec
sesp.mcnemar(slr.clinic.ti.tab) 
sesp.exactbinom(slr.clinic.ti.tab)
#comparing dlrp/dlrn
dlr.regtest(slr.clinic.ti.tab, alpha = 0.05)

#smart clinic ti
smart.clinic.ti.tab <- tab.paired(predclass, smartphone_1_ti_yn, clinic_ti_yn, data=ti_with_lc) #test1 = smart, test2= clinic
#comparing sens/spec
sesp.mcnemar(smart.clinic.ti.tab) 
sesp.exactbinom(smart.clinic.ti.tab)
#comparing dlrp/dlrn
dlr.regtest(smart.clinic.ti.tab, alpha = 0.05)

# # THIS WHOLE SECTION IS NOT REALLY BEING USED IF WE ONLY DO DIFFERENCE FROM LCA
# # BUT LEAVING HERE JUST AS A RECORD
# # THIS WAS TO CALCULATE DIFFERENCE (95%CI) OF SENS/SPEC ASSUMING FIELD GRADE GOLD STANDARD
# # FUNCTION A: Here, I don't have the bang-bang (!!), so when running the function you need to specify the object and column
# spec.diff.fx.original <- function(d,a,b,c) {
#   spec.b <- spec(data=d, truth = factor(a), estimate = factor(b))
#   spec.c <- spec(data=d, truth = factor(a), estimate = factor(c))
#   spec.diff = spec.b$.estimate - spec.c$.estimate
#   spec.diff
# } 
# spec.diff.fx.original(master1,master1$clinic_tf_yn, master1$smartphone_1_tf_yn, master1$SLR_1_tf_yn )
# 
# # FUNCTION B: But here, I use the bang-bang (ie, unquote), so I first have to quote the arguments; 
# # can do this by using the quo() when running the function
# spec.diff.fx <- function(d,a,b,c) {
#   spec.b <- spec(data=d, truth = factor(!!a), estimate = factor(!!b))
#   spec.c <- spec(data=d, truth = factor(!!a), estimate = factor(!!c))
#   spec.diff = spec.b$.estimate - spec.c$.estimate
#   spec.diff
# } 
# spec.diff.fx(master1,quo(clinic_tf_yn), quo(smartphone_1_tf_yn), quo(SLR_1_tf_yn))
# 
# # FUNCTION C: Here I put the quo() inside the function to make it nicer
# spec.diff.fx <- function(d,ref,test1,test2) {
#   ref <- enquo(ref)
#   test1 <- enquo(test1)
#   test2 <- enquo(test2)
#   spec.test1 <- spec(d, truth = factor(!!ref), estimate = factor(!!test1))
#   spec.test2 <- spec(d, truth = factor(!!ref), estimate = factor(!!test2))
#   spec.diff = spec.test1$.estimate - spec.test2$.estimate
#   names(spec.diff)="Diff"
#   spec.diff
# } 
# master1 %>% spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# 
# # JK: OK now feed the function into the command. (It works!)
# bs_diff <- map(bs_nested$splits, ~as_tibble(.) %>% # the map takes each split and turns it into a tibble (data frame)
#                  unnest %>%  # Then because bs_nested was nested above we need to unnest.
#                              # Note though that this is still inside the parentheses of the map
#                  # Now each of the unnested splits is a dataframe and you can do whatever statistic you want
#                  # We will just do our function, which calculates the specificity
#                  # Make sure to use summarize so that the end result is a single row
#                  dplyr::summarize(specdiff=spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn))) %>%
#   # Then bind together the results (the single row from the summarize above)
#   bind_rows(.id = 'boots') 
# bca(bs_diff$specdiff)
# 
# # FUNCTION D: For sensitivity
# # Finally here I put the quo() inside the function to make it nicer
# sens.diff.fx <- function(d,ref,test1,test2) {
#   ref <- enquo(ref)
#   test1 <- enquo(test1)
#   test2 <- enquo(test2)
#   sens.test1 <- sens(d, truth = factor(!!ref), estimate = factor(!!test1))
#   sens.test2 <- sens(d, truth = factor(!!ref), estimate = factor(!!test2))
#   sens.diff = sens.test1$.estimate - sens.test2$.estimate
#   names(sens.diff)="Diff"
#   sens.diff
# } 
# master1 %>% sens.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# 
# 
# # JK: OK now do sensitivity and specificity with same bootstrap resamples
# bs_diff <- map(bs_nested$splits, ~as_tibble(.) %>% 
#                  unnest %>%  
#                  dplyr::summarize(specdiff.tf=spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn),
#                                   sensdiff.tf=sens.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn),
#                                   specdiff.ti=spec.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn),
#                                   sensdiff.ti=sens.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn))) %>%
#   # Then bind together the results (the single row from the summarize above)
#   bind_rows(.id = 'boots') 
# # So finally the mean difference and CIs for TABLE 1:
# master1 %>% spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# bca(bs_diff$specdiff.tf)
# master1 %>% sens.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# bca(bs_diff$sensdiff.tf)
# master1 %>% spec.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn)
# bca(bs_diff$specdiff.ti)
# master1 %>% sens.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn)
# bca(bs_diff$sensdiff.ti)
# 
# # Visualize the distributions
# ggplot(data=bs_diff, aes(specdiff.tf)) +
#   geom_density()
# ggplot(data=bs_diff, aes(sensdiff.tf)) +
#   geom_density()
# ggplot(data=bs_diff, aes(specdiff.ti)) +
#   geom_density()
# ggplot(data=bs_diff, aes(sensdiff.ti)) +
#   geom_density()
# 
# # Let's try it but not accounting for clustering by examiner.
# # We would expect the CIs to be narrower here
# # Note that it's the same code except no "unnest"
# bs_diff2 <- map(bs_unnested$splits, ~as_tibble(.) %>% 
#                  # unnest %>% 
#                  dplyr::summarize(specdiff.tf=spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn),
#                                   sensdiff.tf=sens.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn),
#                                   specdiff.ti=spec.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn),
#                                   sensdiff.ti=sens.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn))) %>%
#   # Then bind together the results (the single row from the summarize above)
#   bind_rows(.id = 'boots') 
# # Mean difference and CIs without accounting for clustering by examiner:
# # Generally narrower but not for specificity of TF
# master1 %>% spec.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# bca(bs_diff2$specdiff.tf)
# master1 %>% sens.diff.fx(., ref=clinic_tf_yn,test1=smartphone_1_tf_yn, test2=SLR_1_tf_yn)
# bca(bs_diff2$sensdiff.tf)
# master1 %>% spec.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn)
# bca(bs_diff2$specdiff.ti)
# master1 %>% sens.diff.fx(., ref=clinic_ti_yn,test1=smartphone_1_ti_yn, test2=SLR_1_ti_yn)
# bca(bs_diff2$sensdiff.ti)

# FIGURE 2
#now for village level prevalences and boostraped confidence intervals
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
                   clin.ti_total=sum(clinic_ti_yn==1 | clinic_ti_yn==0),
                   slr.tfti_count=sum(SLR_1_tfti_yn==1),
                   slr.tfti_total=sum(SLR_1_tfti_yn==1 | SLR_1_tfti_yn==0),
                   smart.tfti_count=sum(smartphone_1_tfti_yn==1),
                   smart.tfti_total=sum(smartphone_1_tfti_yn==1 | smartphone_1_tfti_yn==0),
                   clin.tfti_count=sum(clinic_tfti_yn==1),
                   clin.tfti_total=sum(clinic_tfti_yn==1 | clinic_tfti_yn==0)) %>%
  gather(methodfield, value, slr.tf_count:clin.tfti_total) %>%
  separate(methodfield, into=c("method", "field"), sep="_") %>%
  spread(field, value, convert=TRUE)
villagemeansjk <- as_tibble(cbind(villagenumsjk, binconf(villagenumsjk$count, villagenumsjk$total))) %>% 
  # Note that the above line creates the binomial exact confidence intervals for each village
  separate(method, into=c("method", "grade")) %>%
  mutate(method=case_when(method=="smart" ~ "Smartphone",
                          method=="slr" ~ "SLR",
                          method=="clin" ~ "Field grading",
                          TRUE ~ NA_character_),
         grade=case_when(grade=="tf" ~ "TF",
                         grade=="ti" ~ "TI",
                         grade=="tfti" ~ "TF±TI",
                         TRUE ~ NA_character_))
villagemeansjk %>% group_by(method) %>% dplyr::summarize(sum_total=sum(total)) # Confirming everyone accounted for...
villagepcr <- master1 %>%
  group_by(state_code) %>%
  dplyr::summarize(pcr_count=sum(newindpcr==1),
                   pcr_total=sum(newindpcr==1 | newindpcr==0),
                   pcr_prev=pcr_count/pcr_total) %>%
  dplyr::select(state_code, pcr_count, pcr_total, pcr_prev)
villagemeansjk_wide <- villagemeansjk %>%
  dplyr::select(state_code, method, grade, count, total, PointEst) %>%
  gather(field, value, count:PointEst) %>%
  mutate(methodfield=paste(method, field, sep="_")) %>%
  dplyr::select(-method, -field) %>%
  spread(methodfield, value, convert=TRUE) %>%
  full_join(., villagepcr, by="state_code")
villagemeansjk_cl.wide <- villagemeansjk %>%
  group_by(state_code, grade) %>%
  mutate(PointEstField=if_else(method=="Field grading", PointEst, NA_real_),
         PointEstField=max(PointEstField, na.rm=TRUE)) %>%
  filter(method!="Field grading") %>%
  full_join(., villagepcr, by="state_code")
# Wilcoxon sign rank test (paired) comparing prevalences
# For TF
villagemeansjk_wide_tf <- villagemeansjk_wide %>% filter(grade=="TF")
villagemeansjk_wide_ti <- villagemeansjk_wide %>% filter(grade=="TI")
wilcox.test(villagemeansjk_wide_tf$SLR_PointEst, villagemeansjk_wide_tf$`Field grading_PointEst`, paired=TRUE) 
wilcox.test(villagemeansjk_wide_tf$Smartphone_PointEst, villagemeansjk_wide_tf$`Field grading_PointEst`, paired=TRUE) 
wilcox.test(villagemeansjk_wide_tf$SLR_PointEst, villagemeansjk_wide_tf$Smartphone_PointEst, paired=TRUE) 
# TI
wilcox.test(villagemeansjk_wide_ti$SLR_PointEst, villagemeansjk_wide_ti$`Field grading_PointEst`, paired=TRUE) 
wilcox.test(villagemeansjk_wide_ti$Smartphone_PointEst, villagemeansjk_wide_ti$`Field grading_PointEst`, paired=TRUE) 
wilcox.test(villagemeansjk_wide_ti$SLR_PointEst, villagemeansjk_wide_ti$Smartphone_PointEst, paired=TRUE) 

# Bootstrap the difference
# Make a function
prevdiff.fx <- function(d) {
  d %>% 
    mutate(clin.minus.slr=(`Field grading_PointEst`-SLR_PointEst),
           clin.minus.smart=(`Field grading_PointEst`-Smartphone_PointEst),
           slr.minus.smart=(SLR_PointEst-Smartphone_PointEst)) %>%
    dplyr::summarize(mean.clin.minus.slr=mean(clin.minus.slr),
                     mean.clin.minus.smart=mean(clin.minus.smart),
                     mean.slr.minus.smart=mean(slr.minus.smart))
} 
# Make sure function works
prevdiff.fx(filter(villagemeansjk_wide, grade=="TF"))
prevdiff.fx(filter(villagemeansjk_wide, grade=="TI"))
# Make bootstrap splits
bs_villageprev <- bootstraps(villagemeansjk_wide, times = 9999)

bs_prevdiff <- map(bs_villageprev$splits, ~as_tibble(.) %>% 
                     group_by(grade) %>%
                     prevdiff.fx(.)) %>%
  # Then bind together the results (the single row from the summarize above)
  bind_rows(.id = 'boots') %>%
  gather(field,value, mean.clin.minus.slr:mean.slr.minus.smart ) %>%
  mutate(fieldgrade=paste(field, grade, sep="_")) %>%
  dplyr::select(-field, -grade) %>%
  spread(fieldgrade, value, convert=TRUE)
# Mean difference in prevalence between 2 methods, with 95%CI
# TF
# Clin minus SLR
bs_prevdiff %>% dplyr::summarize(mean.clin.minus.slr.tf=mean(mean.clin.minus.slr_TF))
bca(bs_prevdiff$mean.clin.minus.slr_TF)
# Clin minus Smart
bs_prevdiff %>% dplyr::summarize(mean.clin.minus.smart.tf=mean(mean.clin.minus.smart_TF))
bca(bs_prevdiff$mean.clin.minus.smart_TF)
# SLR minus Smart
bs_prevdiff %>% dplyr::summarize(mean.slr.minus.smart.tf=mean(mean.slr.minus.smart_TF))
bca(bs_prevdiff$mean.slr.minus.smart_TF)
# TI
# Clin minus SLR
bs_prevdiff %>% dplyr::summarize(mean.clin.minus.slr.ti=mean(mean.clin.minus.slr_TI))
bca(bs_prevdiff$mean.clin.minus.slr_TI)
# Clin minus Smart
bs_prevdiff %>% dplyr::summarize(mean.clin.minus.smart.ti=mean(mean.clin.minus.smart_TI))
bca(bs_prevdiff$mean.clin.minus.smart_TI)
# SLR minus Smart
bs_prevdiff %>% dplyr::summarize(mean.slr.minus.smart.ti=mean(mean.slr.minus.smart_TI))
bca(bs_prevdiff$mean.slr.minus.smart_TI)

# # Plot with CIs: This is not finished but leaving here if need starting code in future for a bar graph with confidence intervals that uses points instead of dots
# village_prevplotjk <- ggplot(data=villagemeansjk, aes(x=state_code, y=PointEst, fill=method)) +
#   geom_point(aes(color=method), position=position_dodge(width=0.5)) +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper, color=method), width=0.1, position = position_dodge(0.5), alpha=0.6) +
#   scale_y_continuous(labels = scales::percent) +
#   xlab("Community") + ylab("Prevalence")
#   # facet_wrap(~sign, nrow = 1) +
#   # coord_flip()
# village_prevplotjk

# Scatter plot
village_scatterjk <- ggplot(data=filter(villagemeansjk_cl.wide, grade!="TF±TI"), aes(x=PointEstField, y=PointEst, size=total, group=method)) +
  geom_abline(slope=1, intercept = 0, color="gray", linetype="dotted") +
  # geom_abline(slope=0, intercept = 0.05, color="gray", linetype="dashed") +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color=method), size=0.5, alpha=0.3, position=position_dodge(width=0.0025)) +
  geom_point(aes(shape=method, fill=method, color=method), position=position_dodge(width=0.0025)) +
  # geom_smooth(aes(fill=method, color=method), size=0.5, method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  # geom_text(aes(label=pcr_prev),hjust=0, vjust=0) + # Just to identify which villages have positive PCR
  geom_point(data=filter(villagemeansjk_cl.wide, grade!="TF±TI" & pcr_prev!=0), aes(x=PointEstField, y=pcr_prev, size=total, group=method), shape=8, color="gray") +
  scale_size(range=c(0.5, 3)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 5L)) +   
  scale_y_continuous(labels = scales::percent_format(accuracy = 5L)) +
  scale_shape_manual(values=c(21,21)) +
  scale_color_manual(values=c("blue", "red")) + 
  scale_fill_manual(values=c("white", "white")) + 
  xlab("Field prevalence") + ylab("Photo-prevalence") +
  labs(shape="Photo method", fill="Photo method", color="Photo method" ) +
  guides(size=FALSE,
         shape=guide_legend(keywidth=0.1, keyheight=0.1, default.unit="cm"),
         color=guide_legend(keywidth=0.1, keyheight=0.1, default.unit="cm"),
         fill=guide_legend(keywidth=0.1, keyheight=0.1, default.unit="cm")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border= element_rect(colour="black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color="black", size=8), 
        legend.title=element_text(size=8), 
        legend.text=element_text(size=8),
        legend.position=c(0, 1),
        legend.justification=c(0,1), 
        legend.background = element_rect(fill="white", size=0.25, linetype="solid", colour ="black"),
        legend.direction="vertical",
        legend.key=element_blank(),
        legend.spacing.y = unit(0.1, 'cm')) +
  # coord_fixed(ratio=1) +
  facet_wrap(~grade, scales = "free", dir="v")
  # facet_grid(grade ~ .)
village_scatterjk
ggsave("village_scatterjk_pcrprev.svg", height=5, width=3.3, units='in')
# Note you still made some changes in Illustrator
