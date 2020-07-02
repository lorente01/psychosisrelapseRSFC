# Code for project: "Striatal functional connectivity in psychosis relapse: A comparison between antipsychotic adherent and non-adherent patients at the time of relapse"
# Code developed 3/6/20 to 4/30/20
# Jose M Rubio M.D.

library(tidyverse)
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

data<- bamm.clin.dat %>% 
  full_join(boris.clin.dat) %>% 
  left_join(sessions.dat,by="GridNumber") %>% 
  left_join(sci.dat,by="MRI_session") %>% 
  left_join(qa_ima.dat,by="MRI_session")  %>% 
  mutate(study_group=replace(study_group,study_group=="Non-Adherent","APF")) %>% 
  select(GridNumber,STUDYID,ID,APPrescribed,sex,AgeAtExam,Risperidone_Lvl,Paliperidone_Lvl,Olanzapine_Lvl,Aripiprazole_Lvl,Haloperidol_Lvl,Fluphenazine_Lvl,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,HRLSI_Total,BMI,UrineToxStatus,LQ,CD_RISC_Total,MRI_session,`Subject ID`,study_group,SCI_GSR,SCI_NoGSR,mean_pct_scrub,mins_left,mean_DvarsPost,mean_CorrelFDPost,mean_FD) %>%
  filter(!is.na(MRI_session))  
  
# DESCRIPTIVE STATISTICS
library(gtsummary)
tab1<- data %>%
  select(AgeAtExam,sex,LQ,BMI,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,CD_RISC_Total,HRLSI_Total,UrineToxStatus,APPrescribed,mins_left,study_group) %>%
  mutate(sex=replace(sex,sex==1,"Male"))%>%
  mutate(sex=replace(sex,sex==2,"Female"))
  

tabt<- tab1 %>% select (AgeAtExam,sex,LQ,BMI,Total_Psychotic,BPRS_Total_1_18,NSA_GlobalNegativeSymptoms,CDSS_Total,YMRS_Total,CD_RISC_Total,HRLSI_Total,UrineToxStatus,APPrescribed,mins_left)
  
label<-list(sex ~ "Gender", AgeAtExam ~ "Age",APPrescribed ~ "LAI Antipsychotic prescribed", LQ ~ "Laterality Quotient",BMI~"Body Mass Index", Total_Psychotic ~"Brief Psychopathology Rating Scale Psychotic Subscore",BPRS_Total_1_18 ~ "Brief Psychopathology Rating Scale Score",NSA_GlobalNegativeSymptoms ~ "Negative Symptom Assessment Score",CDSS_Total ~ "Calgary Depressive Symptom Scale Score",YMRS_Total ~ "Young Mania Rating Scale Score",CD_RISC_Total ~ "Connor Davidson Resiliency Scale Score",HRLSI_Total ~ "Holmes Rahe Stressful Life Events Scale Score",UrineToxStatus ~"Positive urine toxicology screen",mins_left~"Mean resting state minutes acquired")
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

#GROUP DIFFERENCES IN RSFC PER ROI

roi.dat<-data %>%
  select(MRI_session,study_group,sex,AgeAtExam,mean_pct_scrub ,mins_left ,mean_FD,mean_DvarsPost ,mean_CorrelFDPost) 

# SCI PER GROUP
SCI_GSR.v<-data$SCI_GSR
SCI_NoGSR.v<-data$SCI_NoGSR
group.v<-data$study_group
sex.v<-data$sex
age.v<-data$AgeAtExam

l.GSR<-lm(SCI_GSR.v~group.v+sex.v+age.v)
l.NoGSR<-lm(SCI_NoGSR.v~group.v+sex.v+age.v)

summary(l.GSR)
summary(l.NoGSR)

cor(SCI_GSR.v,SCI_NoGSR.v)

bp1 <- ggplot(data, aes(x=study_group,y=SCI_GSR)) +
  geom_boxplot(aes(col=factor(study_group)),outlier.size = 4, outlier.shape = 16, outlier.stroke = 2)+
  geom_point(aes(col=factor(study_group)))+
  scale_x_discrete()+
  scale_color_discrete(name="Participant group")+
  ylab("Striatal Connectivity Index (SCI)")+
  ggtitle("Global Signal Regression")

bp2 <- ggplot(data, aes(x=study_group,y=SCI_NoGSR)) +
    geom_boxplot(aes(col=factor(study_group)),outlier.size = 4, outlier.shape = 16, outlier.stroke = 2)+
    geom_point(aes(col=factor(study_group)))+
    scale_x_discrete()+
    scale_color_discrete(name="Participant group")+
  ylab("Striatal Connectivity Index (SCI)")+
  ggtitle("No Global Signal Regression")

boxplot1<-gridExtra::grid.arrange(bp1, bp2, ncol=2)
boxplot1



