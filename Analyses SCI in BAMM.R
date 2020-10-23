
#Code for analyses of: "Striatal functional connectivity in psychosis relapse: A comparison between antipsychotic adherent and non-adherent patients at the time of relapse"
#Completed October 2020
#Jose M Rubio 'jrubio13@northwell.edu' 

library(tidyverse)
#Read from various data files and create master data file to work with
bamm.clin.dat<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/57bcaa6cb52d828c60cc7b9fc16fad3ca8f56441/BAMM%20clinical.csv")
bamm.clin.dat<-bamm.clin.dat %>%
  mutate(LQ=(TotRight - TotLeft)/(TotRight + TotLeft) * 100)%>%
  mutate(CD_RISC_Total=CD_RISC1+CD_RISC2+CD_RISC3+CD_RISC4+CD_RISC5+CD_RISC6+CD_RISC7+CD_RISC8+CD_RISC9+CD_RISC10+CD_RISC11+CD_RISC12+CD_RISC13+CD_RISC14+CD_RISC15+CD_RISC16+CD_RISC17+CD_RISC18+CD_RISC19+CD_RISC20+CD_RISC21+CD_RISC22+CD_RISC23+CD_RISC24+CD_RISC25)  %>%
  mutate(CDSS_Total=CDSS_Depression+CDSS_Hopelessness+CDSS_SelfDepreciation+CDSS_GuiltyIdeas+CDSS_PathologicalGuilt+CDSS_MorningDepr+CDSS_EarlyWake+CDSS_Suicide+CDSS_ObservedDepr)

boris.clin.dat<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/BORIS%20relapse%20clinical.csv")

sessions.dat<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/sessions1.csv")
sci.dat<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/SCI.csv")

qa_ima.dat.1<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/qa_ima.csv")
qa_ima.dat<- qa_ima.dat.1 %>%
  group_by(MRI_session) %>%
  summarise(mean_pct_scrub = mean(PctScrub.5 ),mins_left = sum(MinsLeft.5),mean_DvarsPost = mean(MeanDvarsPost),mean_CorrelFDPost = mean(CorrelFDPost),mean_FD=mean(FDpostScrub))

data1<- bamm.clin.dat %>% 
  full_join(boris.clin.dat) %>% 
  left_join(sessions.dat,by="GridNumber") %>% 
  left_join(sci.dat,by="MRI_session") %>% 
  left_join(qa_ima.dat,by="MRI_session")  %>% 
  mutate(study_group=replace(study_group,study_group=="Non-Adherent","APF")) %>% 
  select(GridNumber,STUDYID,ID,APPrescribed,sex,AgeAtExam,Risperidone_Lvl,Paliperidone_Lvl,Olanzapine_Lvl,Aripiprazole_Lvl,Haloperidol_Lvl,Fluphenazine_Lvl,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,HRLSI_Total,BMI,UrineToxStatus,LQ,CD_RISC_Total,MRI_session,`Subject ID`,study_group,SCI_GSR,SCI_NoGSR,mean_pct_scrub,mins_left,mean_DvarsPost,mean_CorrelFDPost,mean_FD) %>%
  filter(!is.na(MRI_session))  

hc.clin.dat.1<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/HC_clin2.csv")
hc.qa.dat.1<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/QA_HC.csv")
qa_ima.dat<- hc.qa.dat.1 %>%
  group_by(MRI_session) %>%
  summarise(mean_pct_scrub = mean(PctScrub.5 ),mins_left = sum(MinsLeft.5),mean_DvarsPost = mean(MeanDvarsPost),mean_CorrelFDPost = mean(CorrelFDPost),mean_FD=mean(FDpostScrub))

hc.clin.dat<-hc.clin.dat.1 %>%
  left_join(qa_ima.dat,by="MRI_session")

diag.dat<-read_csv("https://raw.githubusercontent.com/lorente01/psychosisrelapseRSFC/master/diag_BAMM.csv")

data<-hc.clin.dat %>%
  full_join(data1)%>%
  left_join(diag.dat,by= "GridNumber")

#Boxplot of SCI distributions
bp1 <- ggplot(data, aes(x=factor(study_group,levels = c('HC','APF','BAMM'),ordered = TRUE),y=SCI_GSR)) +
  geom_violin(aes(col=factor(study_group,levels = c('HC','APF','BAMM'),ordered = TRUE)))+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme(plot.title =element_text(size = 14,face = "bold"),
        axis.title =element_text(size = 12),
        legend.title =element_text(size = 10))+
  scale_x_discrete("Participant group")+
  scale_color_discrete(name="Participant group",labels=c("HC=Healthy Controls","APF=Antipsychotic free upon relapse","BAMM=Breakthrough psychosis"))+
  theme(legend.position = "bottom")+
  ylab("Striatal Connectivity Index (SCI)")+
  ggtitle("Comparison in SCI by antipsychotic adherence status upon relapse")+
  guides(linetype=guide_legend(nrow=2))
bp1

#Tabulation of sample characteristics by group
library(gtsummary)
tab1<- data %>%
  select(AgeAtExam,sex,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,CD_RISC_Total,HRLSI_Total,UrineToxStatus,APPrescribed,study_group,Primary_diagnosis) %>%
  mutate(sex=replace(sex,sex==1,"Male"))%>%
  mutate(sex=replace(sex,sex==2,"Female"))

tabt<- tab1 %>% select (AgeAtExam,sex,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,CD_RISC_Total,HRLSI_Total,UrineToxStatus,APPrescribed)

label<-list(sex ~ "Gender", AgeAtExam ~ "Age",APPrescribed ~ "LAI Antipsychotic prescribed", Total_Psychotic ~"Brief Psychopathology Rating Scale Psychotic Subscore",BPRS_Total_1_18 ~ "Brief Psychopathology Rating Scale Score",NSA_GlobalNegativeSymptoms ~ "Negative Symptom Assessment Score",CDSS_Total ~ "Calgary Depressive Symptom Scale Score",YMRS_Total ~ "Young Mania Rating Scale Score",CD_RISC_Total ~ "Connor Davidson Resiliency Scale Score",HRLSI_Total ~ "Holmes Rahe Stressful Life Events Scale Score",UrineToxStatus ~"Positive urine toxicology screen")
type<-list(NSA_GlobalNegativeSymptoms ~"continuous")
statistic <- list(all_continuous() ~ "{mean} ({sd})")
digits <- list(all_continuous() ~ 2)
#value = list(APPrescribed ~ "Aripiprazole",APPrescribed ~"Fluphenazine",APPrescribed ~"Haloperidol",APPrescribed ~"Paliperidone")
table1g<-tbl_summary(tab1,by=study_group,statistic = statistic,digits = digits,missing="no" ,label=label,type=type) %>% 
  add_p() %>% 
  bold_labels()

table1t<-tbl_summary(tabt,digits = digits,statistic = statistic,missing="no" ,label=label,type=type)
table1<-tbl_merge(tbls=list(table1t,table1g),tab_spanner = c("Total Sample","By Group")) 

table1 

#linear model adjusting for sex and age
#Rename outcome variable so reference condition (BAMM) is picked first in the model

data2<- data %>%
mutate(study_group=replace(study_group,study_group=="BAMM","1BAMM")) %>%
mutate(study_group=replace(study_group,study_group=="APF","2APF")) %>%
mutate(study_group=replace(study_group,study_group=="HC","3HC")) 
  
#Covariates of interest
SCI_GSR.v<-data2$SCI_GSR
SCI_NoGSR.v<-data2$SCI_NoGSR
group.v<-as.factor(data2$study_group)
sex.v<-data2$sex
age.v<-data2$AgeAtExam

#lm GSR
model1<-lm(SCI_GSR.v~group.v+sex.v+age.v)
summary(model1)

#lm No GSR
model2<-lm(SCI_NoGSR.v~group.v+sex.v+age.v)
summary(model2)

#Now removing BP-I diagnosis
data3<- data2 %>% 
  filter(Primary_diagnosis=="Schizophrenia"|Primary_diagnosis=="Schizoaffective"|Primary_diagnosis=="Psychosis NOS" |Primary_diagnosis=="No axis I disorder")

SCI_GSR.v<-data3$SCI_GSR
SCI_NoGSR.v<-data3$SCI_NoGSR
group.v<-as.factor(data3$study_group)
sex.v<-data3$sex
age.v<-data3$AgeAtExam

#lm GSR
model3<-lm(SCI_GSR.v~group.v+sex.v+age.v)
summary(model3)

#lm No GSR
model4<-lm(SCI_NoGSR.v~group.v+sex.v+age.v)
summary(model4)

#Calculate effect size of BAMM vs HC
library(effsize)
group<-data$study_group[data$study_group=="BAMM" | data$study_group=="HC"]

#GSR
SCI_GSR.c<-data$SCI_GSR[data$study_group=="BAMM" | data$study_group=="HC"]
cohen.d(SCI_GSR.c,group)

#No GSR
SCI_NoGSR.c<-data$SCI_NoGSR[data$study_group=="BAMM" | data$study_group=="HC"]
cohen.d(SCI_NoGSR.c,group)

#Calculate effect size of BAMM vs APF
group<-data$study_group[data$study_group=="BAMM" | data$study_group=="APF"]

#GSR
SCI_GSR.c<-data$SCI_GSR[data$study_group=="BAMM" | data$study_group=="APF"]
cohen.d(SCI_GSR.c,group)

#No GSR
SCI_NoGSR.c<-data$SCI_NoGSR[data$study_group=="BAMM" | data$study_group=="APF"]
cohen.d(SCI_NoGSR.c,group)

#Correlation between GSR and No GSR derived SCIs
cor(data$SCI_GSR,data$SCI_NoGSR)
