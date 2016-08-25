# 40101_Neutropenia_Association

setwd("~/Desktop/R/Rpractice")

library(GenABEL)
library(gdata)
library(survival)
library(ggplot2)

gwasfile="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/R_data_files/GWA40101CACAUC759reducedauto.RData"
tools::md5sum(gwasfile)

phenofile="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/ac_arm_anc_40101jan32014_fmEdits.csv"
tools::md5sum(phenofile)

attach(gwasfile)
GDAT=GDATCAUCreducedauto

###Removing the two PatIDs later determined not to have given consent.
exc=c("batch2l_plate10_C11", "batch1l_plate02_A04") ##genotyped but did not give consent
GDAT=GDAT[!is.element(GDAT@phdata$id,exc),] ##Removing the two patients that did give consent

GDAT@gtdata@nids
GDAT@gtdata@nsnps

save(GDAT,file=paste("~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/R_data_files/GWA40101CACAUC757reducedauto.RData",sep=""))

pheno=read.csv(phenofile)
table(pheno$anc_event,exclude=NULL)

pdat=read.table("~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/C40101_dbGaP_SubjectPhenotypesDD_CJ.txt",header=TRUE)

GDAT=GDAT[is.element(GDAT@phdata$patid,as.character(pdat$patid[pdat$consent=="Yes"])),]

### ANC from on study which should be multiplied by 10^3
### wbc which is an indicator 1=WBC*10^3<5 
### bili is an indicator 1=bilirubin>1; 0=bilirubin<=1
### cc is  Cockroft and Gault formula=Creatinine Clearance = (140 - age) x weight (kg) x 0.85 (for females)/72 X Serum Creatinine (mg/dL);
### GFR is an indicator for cc<60,            1=cc<60;0=cc>=60
### Prior_chemo is an indicator for whether the patient had prior adjuvant chemotherapy for this malignancy.

##value ancevent 
###   1= 'Neutropenic Event' 
###   0= 'Control' 
###  21='Inadequate AE Documentation' 
###  22='Less than 4 cycles of treatment' 
###  23='Only Drug given on 1st cycle' 
###  24='Wrong Drug Given' 
###  25='No treatment drug information' 
###  26='Concurrent non-protocol chemotherapy';

GDAT@phdata=merge(GDAT@phdata,pheno[,c("patid","anc_event","ANC","wbc","bili","cc","GFR","G_CSF_cyc1","anc_event","prdi","age_reg")],all.x=TRUE)

GDAT@phdata=GDAT@phdata[match(GDAT@gtdata@idnames,GDAT@phdata$id),]
rownames(GDAT@phdata)=as.character(GDAT@phdata$id)

table(GDAT@phdata$anc_event,exclude=NULL)
GDAT@phdata$ANCevent=ifelse(GDAT@phdata$anc_event==1,1,ifelse(GDAT@phdata$anc_event==0,0,NA))
table(GDAT@phdata$ANCevent,GDAT@phdata$anc_event,exclude=NULL)

mod=glm(formula = ANCevent ~ G_CSF_cyc1 + GFR + bili + prdi + wbc + age_reg, family = "binomial", data = GDAT@phdata)
mod ## shows results

## The summary() function gives the 'Analysis of Maximum Likelihood Estimates' in the '40101 ANC 1000 analysis' doc that Flora generated
summary(mod)

## The following are commands that are useful to know when extracting data from the model output
coefficients(mod) # model coefficients
confint(mod, level=0.95) # CIs for model parameters 
fitted(mod) # predicted values
residuals(mod) # residuals
anova(mod) # anova table 
vcov(mod) # covariance matrix for model parameters 
influence(mod) # regression diagnostics

## odds ratios only
exp(coef(mod))


## AC induced Neutropenia
result=mlreg(ANCevent~1,trait.type="binomial",GDAT[!is.na(GDAT@phdata$ANCevent),])
RES=results(result)

result=mlreg(ANCevent~G_CSF_cyc1 + GFR + bili + prdi + wbc + age_reg,trait.type="binomial",GDAT[!is.na(GDAT@phdata$ANCevent),])
adjRES=results(result)

aa=summary(GDAT[!is.na(GDAT@phdata$ANCevent),])
RES4cys=cbind(RES[,c(6,7,10)],data.frame(adjpval=adjRES$P1df,logadjHR=adjRES$effB,logSE_HR=adjRES$se_effB))

RES4cys=merge(RES4cys,aa[,c(1:2,6:12)],by=0)
RES4cys=RES4cys[order(RES4cys$P1df),]

save(RES4cys, file="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/New/C40101ACarm_ANC4CYS_042616.RData")
write.csv(RES4cys[1:500,], file="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/New/C40101ACarm_ANC4CYS_Top500_042616.csv",row.names=FALSE)

RES4cyssub=RES4cys[RES4cys$Q.2>0.01,]
write.csv(RES4cyssub[1:500,],file="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/New/C40101ACarm_ANC4CYS_MAF01_Top500_042616.csv",row.names=FALSE)


## ANC baseline:
RES=mlreg(log(ANC)~1,dat=GDAT[!is.na(GDAT@phdata$ANC),],trait.type = "gaussian")

adjRES=mlreg(log(ANC)~GFR + bili + wbc,dat=GDAT[!is.na(GDAT@phdata$ANC),],trait.type = "gaussian")

ANCRES=cbind(RES@results[,c(1,2,5)],data.frame(adjpval=(adjRES@results)$P1df,logadjHR=(adjRES@results)$effB),logSE_HR=(adjRES@results)$se_effB)

aa=summary(GDAT[!is.na(GDAT@phdata$ANC),])
ANCRES=merge(ANCRES,aa[1,2,6:12],by=0)

ANCRES=ANCRES[order(ANCRES$P1df),]

save(ANCRES,file="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/New/C40101ACarm_LogANCbaseline_042616.RData")
write.csv(ANCRES[1:500,],file="~/documents/postdoc/Neutropenia_Data/Gen_Phen_files+Script/New/C40101ACarm_LogANCbaselineTop500_042616.csv",row.names=FALSE)

