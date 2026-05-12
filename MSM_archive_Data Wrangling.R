  # Load data and package, adjust format and order  
  pacman::p_load(tidyverse,data.table,readxl,ggplot2,writexl,RColorBrewer,lubridate,ggalluvial,ggsurvfit,survminer,lmtest,sandwich,survival,rms,mice)
  dx <- readRDS("E:/OneDrive - The University of Hong Kong - Connect/RA explore/RA/Diagnosis/Diagnosis.RDS")
  rx <- readRDS("E:/OneDrive - The University of Hong Kong - Connect/RA explore/RA/Prescription/rx_ra_20002023_20022023.RDS")
  death <- readRDS("E:/OneDrive - The University of Hong Kong - Connect/RA explore/RA/Death Cause/death_20002023.RDS")
  Demo <- readRDS("E:/OneDrive - The University of Hong Kong - Connect/RA explore/RA/Identification/Identification.RDS")
  lab <- readRDS("E:/OneDrive - The University of Hong Kong - Connect/RA explore/RA/Lab/Lab_RA_diagnosemarker.RDS")

  colnames(Demo) <- make.names(colnames(Demo))
  setorder(rx,PrescriptionStartDate)
  rx <- rx[!is.na(PrescriptionStartDate)]
  rx$PrescriptionEndDate <- as.Date(rx$PrescriptionEndDate)
  rx$PrescriptionStartDate <- as.Date(rx$PrescriptionStartDate)
  dx$date <- as_date(dx$date)
  setorder(dx,date)
  rx <- rx[!PrescriptionStartDate>as.Date("2024-12-31")]#Extreme records/missing in start date/missing in end date
  #Population:
  #1: Patients diagnosed with rheumatoid arthritis since 2002
  #2: Without other common autoimmune disease
  co <- dx[str_detect(code,"^714.0")][order(date),head(.SD,1),ReferenceKey]
  #Identify patients with other autoimmune disease 
  coad <- unique(dx[ReferenceKey%in%co$ReferenceKey][str_detect(code,"^55[56]|^710.0|^720.0|^340|^696")]$ReferenceKey)
  #Other drugs indicated for RITUXIMAB
  coritu <- unique(dx[ReferenceKey%in%co$ReferenceKey][str_detect(code,"^446|^20[02]|^204.1")]$ReferenceKey)
  co <- co[!ReferenceKey%in%coad][!ReferenceKey%in%coritu]

  druglist <- c("INFLIXIMAB","ADALIMUMAB","CERTOLIZUMAB","ETANERCEPT",
                "GOLIMUMAB","ABATACEPT","RITUXIMAB","SARILUMAB","TOCILIZUMAB",
                "BARICITINIB","TOFACITINIB","UPADACITINIB")
  alldrug<- rx[ReferenceKey%in%co$ReferenceKey][str_detect(DrugName,regex("Infliximab|Adalimumab|Certolizumab|Etanercept|Golimumab|Abatacept|Rituximab|Sarilumab|Tocilizumab|Baricitinib|Tofacitinib|Upadacitinib",ignore_case = T))]
  for (i in 1:length(druglist)) {
    alldrug[str_detect(DrugName,regex(druglist[i],ignore_case = T)),agent:=toupper(druglist[i])]
    print(i)
  }
   
  #rituximab (a B lymphocyte depleting agent), abatacept (targets T cell co-stimulation)
  #With at least two consecutive records
  singlerecord <- alldrug[,.N,.(ReferenceKey,agent)][N==1]
  alldrug <- anti_join(alldrug,singlerecord)
  
  ccc <- alldrug %>%
    arrange(ReferenceKey, PrescriptionStartDate) %>%
    mutate(PrescriptionStartDate = ymd(PrescriptionStartDate), 
           PrescriptionEndDate = ymd(PrescriptionEndDate)) %>%
    group_by(ReferenceKey,agent) %>%
    mutate(indx = c(0, cumsum(as.numeric(lead(PrescriptionStartDate)) > 
                                cummax(as.numeric(PrescriptionEndDate)))[-n()])) %>% select(ReferenceKey, start = PrescriptionStartDate,end = PrescriptionEndDate,agent) %>% distinct()

    calc_time <- function(refDate, endDate) {
    period <- as.period(interval(refDate, endDate), unit = "day")
    period$day
  }
  
  ccc.gap <- ccc %>% 
    arrange(ReferenceKey, start, agent) %>%
    group_by(ReferenceKey,agent) %>%
        mutate(gap = as.numeric(start - lag(end, default = min(start))),
           change = if_else(gap > 180, 1, 0),
           duration = calc_time(start, end)) %>% 
    group_by(ReferenceKey) %>%
    mutate(change = if_else(agent != lag(agent)|gap>180&agent!="RITUXIMAB"|gap>365&agent=="RITUXIMAB", 1, 0),
           regimen = cumsum(change))  
  
  setDT(ccc.gap) 
  
  ccc.gap[agent%in%c("TOFACITINIB","BARICITINIB","UPADACITINIB"),moa:="JAKi"]
  ccc.gap[agent%in%c("ETANERCEPT","INFLIXIMAB","GOLIMUMAB","CERTOLIZUMAB","ADALIMUMAB"),moa:="TNFi"]
  ccc.gap[agent%in%c("TOCILIZUMAB","SARILUMAB"),moa:="IL6"]
  ccc.gap[agent%in%c("RITUXIMAB","ABATACEPT"),moa:="Lymphocyte"]  

  ccc.gap[is.na(change),change:=0]
  ccc.gap[,switch:=cumsum(change),.(ReferenceKey)]
  ccc.new <- ccc.gap[,.(start=min(start),end=max(end),agent),.(ReferenceKey,moa,switch)]
  ccc.new <- distinct(ccc.new)
  
  #Clean the prescription records
  #Exclude patients who switched back to their initial b/ts DMARDs
  switchback <- unique(ccc.new[,.N,.(ReferenceKey,agent)][N>1]$ReferenceKey)
  ccc.back <- ccc.new[ReferenceKey%in%switchback]
  ccc.back[order(start),gap := as.numeric(start - lag(end, default = min(start))),.(agent,ReferenceKey)]
  ccc.back[,switch:=ifelse(gap>180,1,0)]
  ccc.back[,switch:=cumsum(switch),.(agent,ReferenceKey)]
  ccc.back <- distinct(ccc.back[switch==0,.(start=min(start),end=max(end),moa),.(ReferenceKey,agent)])
  ccc.new <- rbind(ccc.back, ccc.new[!ReferenceKey%in%switchback][,switch:=NULL])
  #Exclude records with overlap greater than 180 days
  ccc.dup <- copy(ccc.new)
  ccc.dup[order(start),start_s:=lead(start),ReferenceKey]
  ccc.dup <- ccc.dup[ReferenceKey%in%ccc.dup[start_s<end&start_s>start]$ReferenceKey,]
  ccc.dup[,gap:=as.numeric(end-start_s)]
  ccc.dup[gap>180]#180 days consistent with previous definition of treatment gap
  ccc.new <- ccc.new[!ReferenceKey%in%ccc.dup[gap>180]$ReferenceKey]
  #If previous prescription end date is later than subsequent prescription start date,change the end date
  ccc.new[,T_enddate_bts:=lead(start),.(ReferenceKey)]
  ccc.new[end>=T_enddate_bts,end:=T_enddate_bts-1] 
  ccc.new[,T_enddate_bts:=NULL]
  ccc.new <- merge(ccc.new,Demo[,.(ReferenceKey,sex=Sex,birthdate=as.Date(DateofBirth.yyyy.mm.dd.))],by = c("ReferenceKey"))
  ccc.new <- ccc.new[start>=as.Date("2009-01-01")&start<=as.Date("2022-12-31")]
  #add the lag period
  ccc.new[,end:=ifelse(agent!="RITUXIMAB",end+90,end+183)]
  ccc.new[end>as.Date("2023-12-31"),end:=as.Date("2023-12-31")]
  ccc.new <- merge(ccc.new,death[,.(ReferenceKey,DateofRegisteredDeath)],by = c("ReferenceKey"),all.x = T)
  ccc.new[!is.na(DateofRegisteredDeath)&DateofRegisteredDeath>=as.Date("2024-1-1"),DateofRegisteredDeath:=NA]
  ccc.new[,age:=as.numeric(round((start-birthdate)/365.25,1))]
  ccc.new <- ccc.new[age>=18]#only adult patients  
  ccc.new[,switch:=1][,switch:=cumsum(switch)-1,.(ReferenceKey)]
  ccc.new[,temname:=paste0(ReferenceKey,"_",switch)]
  #Baseline time-fixed comorbidities
  dxcode <- read_xlsx("Variable_safety.xlsx",sheet="CCI")
  dxcc <- merge(dx[ReferenceKey%in%ccc.new$ReferenceKey],ccc.new[,.(ReferenceKey,start,switch,temname)],by=c("ReferenceKey"),allow.cartesian = T)
  dxcc <- dxcc[date<start]
  
  for (i in 1:length(dxcode$Regex)) {
    l <- dxcc[str_detect(code,dxcode$Regex[i])]$temname
    ccc.new[,dxcode$Name[i]:=ifelse(temname%in%l,1,0)]
    print(i)
  } 
  #Split observation period using 3 months interval
  for (i in 0:19) {
    ccc.new[,paste0("t",i):=start%m+% months(i*3)]
    print(i)
  }
  
  ccc.long <- melt(ccc.new, id.vars = colnames(ccc.new[,!(t0:t19)]),
                measure.vars = colnames(ccc.new[,t0:t19]),
       variable.name = "visit",value.name = "timepoint",
       variable.factor = T)
  ccc.long$visit <- as.numeric(str_remove_all(ccc.long$visit,"t"))
  ccc.long[,trt:=ifelse(timepoint<end,as.numeric(1),as.numeric(0))]
  
  ccc.long[,count:=cumsum(trt),.(ReferenceKey,start,moa,switch)]
  ccc.long[,keep:=visit<=max(count),.(ReferenceKey,start,moa,switch)]
  ccc.long <-ccc.long[keep==TRUE]  
  ccc.long[,timepoint_bk:=timepoint%m-%months(3)]
  ccc.long[,outcome_end:=timepoint%m+% months(3)]
  ccc.long[DateofRegisteredDeath<outcome_end&DateofRegisteredDeath>=timepoint,death:=1]
  ccc.long[is.na(death),death:=0]
  ccc.long <- rbind(ccc.long[!timepoint>DateofRegisteredDeath],ccc.long[is.na(DateofRegisteredDeath)])#Some people die before the drug end date
#Time-varing variables including exposure within 3 month (binary) and prescription dosage of  
  
  calc_time_2 <- function(A1,A2,B1,B2) {
  period <- max(0, min(A2, B2) - max(A1, B1) + 1) 
  return(period)
  }
  
  #Create matching table
  
  rx_nsaid <-rx[ReferenceKey%in%ccc.long$ReferenceKey][str_detect(DrugName,"CELECOXIB|ETORICOXIB|DICLOFENAC|IBUPROFEN|MEFENAMIC|NAPROXEN|PIROXICAM|SULINDAC")]
  #Remove missing value and remove unwanted indications
  rx_nsaid <- rx_nsaid[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(PrescriptionStartDate))][str_detect(Route,"ORAL")&!str_detect(DrugStrength,"ML")]    
  nsaid_table <- data.table(Frequency=unique(rx_nsaid$DrugFrequency))
  
  rx_pred <-rx[ReferenceKey%in%ccc.long$ReferenceKey][str_detect(DrugName,"^PREDNISOLONE")]
  rx_pred <-rx_pred[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(PrescriptionStartDate))][Route=="ORAL"&!str_detect(DrugStrength,"ML")]
  pred_table <- data.table(Frequency=unique(rx_pred$DrugFrequency))
  
  rx_opioid <- rx[ReferenceKey%in%ccc.long$ReferenceKey][str_detect(DrugName,"CODEINE|DEXTROPROPOXYPHENE|DIHYDROCODEINE|METHADONE|MORPHINE|OXYCODONE|TRAMADOL")]
  rx_opioid <-rx_opioid[!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(PrescriptionStartDate))][Route=="ORAL"&!str_detect(DrugStrength,"ML")]
  opioid_table <- data.table(Frequency=unique(rx_opioid$DrugFrequency))
  
  #frequency_table <- distinct(rbind(nsaid_table,pred_table,opioid_table))
  # write_xlsx(frequency_table,"Frequency_table.xlsx")  

  frequency_table <- data.table(readxl::read_excel("frequency_table.xlsx"))
  for (i in 1:length(frequency_table$Frequency_text)) {
    rx_nsaid[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    rx_opioid[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    rx_pred[DrugFrequency==frequency_table$Frequency_text[i],Frequency:=frequency_table$Frequency_numeric[i]]
    print(i)
  }
  rx_nsaid <- rx_nsaid[!Frequency==0]
  rx_opioid <- rx_opioid[!Frequency==0]
  rx_pred <- rx_pred[!Frequency==0]
  
  nsaidname <- c("CELECOXIB","ETORICOXIB","DICLOFENAC","IBUPROFEN","MEFENAMIC","NAPROXEN","PIROXICAM","SULINDAC")
  for (i in 1:length(nsaidname)) {
    rx_nsaid[str_detect(DrugName,regex(nsaidname[i],ignore_case = T)),agent:=nsaidname[i]]
    print(i)    
  }    
  
  setorder(rx_nsaid,PrescriptionStartDate,PrescriptionEndDate)
  rx_nsaid <- distinct(rx_nsaid,ReferenceKey,PrescriptionStartDate,PrescriptionEndDate,agent,.keep_all = T)
  
  opioidname <- c("CODEINE","DEXTROPROPOXYPHENE","DIHYDROCODEINE","METHADONE","MORPHINE","OXYCODONE","TRAMADOL")
  for (i in 1:length(opioidname)) {
    rx_opioid[str_detect(DrugName,regex(opioidname[i],ignore_case = T)),agent:=opioidname[i]]
    print(i)    
  }   
  
  setorder(rx_opioid,PrescriptionStartDate,PrescriptionEndDate)
  rx_opioid <- distinct(rx_opioid,ReferenceKey,PrescriptionStartDate,PrescriptionEndDate,agent,.keep_all = T)
  
  ddd <- setDT(read_excel("Drug DDD.xlsx"))
  
  setorder(rx_pred,PrescriptionStartDate,PrescriptionEndDate)
  rx_pred <- distinct(rx_pred,ReferenceKey,PrescriptionStartDate,PrescriptionEndDate,.keep_all = T)

  rx_pred <- merge(ccc.long[,.(ReferenceKey,temname,timepoint,timepoint_bk)],rx_pred,all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
  rx_pred[,T_enddate:=lead(PrescriptionStartDate),.(timepoint,temname)]
  rx_pred[PrescriptionEndDate>=T_enddate,PrescriptionEndDate:=T_enddate-1]
  rx_pred[,overlap:=mapply(calc_time_2,timepoint_bk,timepoint,PrescriptionStartDate,PrescriptionEndDate)]
  rx_pred[,strength:=parse_number(DrugStrength)]
  rx_pred[,dosage:=parse_number(Dosage)]
  ddd_pred <- ddd[Type=="STERIOD"]$`WHO defined daily dose`
  rx_pred <- rx_pred[,.(prednisolone=sum(round((overlap*strength*Frequency*dosage)/(ddd_pred*(as.integer(timepoint-timepoint_bk)+1)),2),na.rm = T)),.(temname,timepoint,timepoint_bk)]

  # rx_pred <- rx_pred[,.(prednisolone=sum(overlap*truedose/ddd_pred,na.rm = T)),.(temname,timepoint,timepoint_bk)]
  # rx_pred[,prednisolone:=round(prednisolone/(as.numeric(timepoint-timepoint_bk)+1),2),.(temname,timepoint,timepoint_bk)]
  ccc.long <- merge(ccc.long,rx_pred,all.x=T,allow.cartesian = T,by=c("temname","timepoint","timepoint_bk"))
  ccc.long[is.na(prednisolone),prednisolone:=as.numeric(0)]
  
  
  
  
  rx_nsaid <- merge(ccc.long[,.(ReferenceKey,temname,timepoint,timepoint_bk)],rx_nsaid,all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
  rx_nsaid[,T_enddate:=lead(PrescriptionStartDate),.(timepoint,temname,agent)]
  rx_nsaid[PrescriptionEndDate>=T_enddate,PrescriptionEndDate:=T_enddate-1]
  rx_nsaid[,overlap:=mapply(calc_time_2,timepoint_bk,timepoint,PrescriptionStartDate,PrescriptionEndDate)]
  rx_nsaid[,strength:=parse_number(DrugStrength)]
  rx_nsaid[,dosage:=parse_number(Dosage)]
  rx_nsaid[dosage>5,dosage:=1]#some dosage was written in combination of drug strength and dosage
  rx_nsaid[,truedose:=strength*Frequency*dosage]  
  for (i in 1:length(ddd[Type=="NSAID"]$`WHO defined daily dose`)) {
    rx_nsaid[str_detect(DrugName,regex(ddd[Type=="NSAID"]$Description[i],ignore_case = T)),ddd_nsaid:=ddd[Type=="NSAID"]$`WHO defined daily dose`[i]]
    print(i)
  }
  rx_nsaid <- rx_nsaid[,.(nsaid=sum(round(overlap*truedose/ddd_nsaid/(as.numeric(timepoint-timepoint_bk)+1),2),na.rm = T)),.(temname,timepoint,timepoint_bk)]
  
  # rx_nsaid <- rx_nsaid[,.(nsaid=sum(overlap*truedose/ddd_nsaid,na.rm = T)),.(temname,timepoint,timepoint_bk)]
  # rx_nsaid[,nsaid:=round(nsaid/(as.numeric(timepoint-timepoint_bk)+1),2),.(temname,timepoint,timepoint_bk)]  
  ccc.long <- merge(ccc.long,rx_nsaid,all.x=T,allow.cartesian = T,by=c("temname","timepoint","timepoint_bk"))
  ccc.long[is.na(nsaid),nsaid:=as.numeric(0)]

   
  rx_opioid <- merge(ccc.long[,.(ReferenceKey,temname,timepoint,timepoint_bk)],rx_opioid,all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
  rx_opioid[,T_enddate:=lead(PrescriptionStartDate),.(timepoint,temname,agent)]
  rx_opioid[PrescriptionEndDate>=T_enddate,PrescriptionEndDate:=T_enddate-1]
  rx_opioid[,overlap:=mapply(calc_time_2,timepoint_bk,timepoint,PrescriptionStartDate,PrescriptionEndDate)]
  rx_opioid[,strength:=parse_number(DrugStrength)]
  rx_opioid[,dosage:=parse_number(Dosage)]
  rx_opioid[,truedose:=strength*Frequency*dosage]  
  for (i in 1:length(ddd[Type=="OPIOID"]$`WHO defined daily dose`)) {
    rx_opioid[str_detect(DrugName,regex(ddd[Type=="OPIOID"]$Description[i],ignore_case = T)),ddd_opioid:=ddd[Type=="OPIOID"]$`WHO defined daily dose`[i]]
    print(i)
  }
  rx_opioid <- rx_opioid[,.(opioid=sum(round(overlap*truedose/ddd_opioid/(as.numeric(timepoint-timepoint_bk)+1),2),na.rm = T)),.(temname,timepoint,timepoint_bk)]
  
  # rx_opioid <- rx_opioid[,.(opioid=sum(overlap*truedose/ddd_opioid,na.rm = T)),.(temname,timepoint,timepoint_bk)]
  # rx_opioid[,opioid:=round(opioid/(as.numeric(timepoint-timepoint_bk)+1),2),.(temname,timepoint,timepoint_bk)]  
  ccc.long <- merge(ccc.long,rx_opioid,all.x=T,allow.cartesian = T,by=c("temname","timepoint","timepoint_bk"))
  ccc.long[is.na(opioid),opioid:=as.numeric(0)]

  druglist_csdmard <- c("HYDROXYCHLOROQUINE","SULPHASALAZINE","METHOTREXATE","LEFLUNOMIDE","AZATHIOPRINE","CYCLOSPORIN")
  for (i in 1:length(druglist_csdmard)) {
  rx_csdmards <- merge(ccc.long[,.(ReferenceKey,temname,timepoint,timepoint_bk)], rx[ReferenceKey%in%ccc.long$ReferenceKey][str_detect(DrugName,druglist_csdmard[i])][!(is.na(Dosage)|is.na(DrugStrength)|is.na(DrugFrequency)|is.na(PrescriptionStartDate))],all.x=T,allow.cartesian = T,by=c("ReferenceKey"))
  rx_csdmards[,overlap:=mapply(calc_time_2,timepoint_bk,timepoint,PrescriptionStartDate,PrescriptionEndDate)]
  rx_csdmards <- rx_csdmards[,.(csdmards=sum(overlap,na.rm = T)),.(temname,timepoint,timepoint_bk)]
  rx_csdmards[is.na(csdmards),csdmards:=0]
  rx_csdmards[,druglist_csdmard[i]:=ifelse(csdmards>0,1,0)]
  rx_csdmards[,csdmards:=NULL]
  ccc.long <- merge(ccc.long,rx_csdmards,all.x=T,allow.cartesian = T,by=c("temname","timepoint","timepoint_bk"))
  }
  
  #Previous exposure to btsDMARDs
  pre_bts <- distinct(merge(ccc.long[,.(ReferenceKey,temname,start)],alldrug[,.(ReferenceKey,agent,PrescriptionStartDate)],all.x=T,allow.cartesian = T,by=c("ReferenceKey")))
  pre_bts <- pre_bts[PrescriptionStartDate<start,.(base_btsdmards=uniqueN(agent)),.(temname,start)]
  
  #Previous exposure to csDMARDs
  pre_csdmard <- rx[ReferenceKey%in%ccc.long$ReferenceKey][str_detect(DrugName,regex("HYDROXYCHLOROQUINE|SULPHASALAZINE|METHOTREXATE|LEFLUNOMIDE|AZATHIOPRINE|CYCLOSPORIN",ignore_case = T))]
  for (i in 1:length(druglist_csdmard)) {
    pre_csdmard[str_detect(DrugName,regex(druglist_csdmard[i],ignore_case = T)),agent:=toupper(druglist_csdmard[i])]
    print(i)
  }  
  pre_csdmard <- distinct(merge(ccc.long[,.(ReferenceKey,temname,start)],pre_csdmard[,.(ReferenceKey,agent,PrescriptionStartDate)],all.x=T,allow.cartesian = T,by=c("ReferenceKey")))
  pre_csdmard <- pre_csdmard[PrescriptionStartDate<start,.(base_csdmards=uniqueN(agent)),.(temname,start)]
  ccc.long <- merge(ccc.long,pre_bts,all.x=T,allow.cartesian = T,by=c("temname","start"))
  ccc.long <- merge(ccc.long,pre_csdmard,all.x=T,allow.cartesian = T,by=c("temname","start"))
  ccc.long[is.na(base_btsdmards),base_btsdmards:=as.numeric(0)]
  ccc.long[is.na(base_csdmards),base_csdmards:=as.numeric(0)]
  ccc.long[,sex:=ifelse(sex=="M",1,0)]
  
  #Missing value imputation for lab test results
  lab[,LISReferenceDatetime:=as_date(substr(LISReferenceDatetime,1,10))]
  lab[,LISResultNumericResult:=as.numeric(LISResultNumericResult)]
  ra_labs <- rbind(lab[grepl("crp|C-reactive Protein", LISTestDescription, ignore.case=T) & LISTestUnit=="mg/dL"][, test_type := "crp"], lab[grepl("esr|Erythrocyte", LISTestDescription, ignore.case=T)][, test_type := "esr"])
  ra_labs <- merge(ccc.long[, .(temname, ReferenceKey, visit, timepoint)], ra_labs, by="ReferenceKey", allow.cartesian=T)[as.Date(LISReferenceDatetime) >= (timepoint %m-% months(6)) & as.Date(LISReferenceDatetime) <= timepoint]
  ra_labs <- merge(ccc.long[, c(colnames(ccc.long[,cci.mi:cci.liver_modsev]),"temname","visit","age","sex","prednisolone","nsaid","opioid","HYDROXYCHLOROQUINE","SULPHASALAZINE","METHOTREXATE","LEFLUNOMIDE","AZATHIOPRINE","CYCLOSPORIN","base_btsdmards","base_csdmards"), with=F], dcast(ra_labs[, .(temname, visit, test_type, result=as.numeric(LISResultNumericResult))], temname + visit ~ test_type, fun.aggregate = mean, na.rm=T), by=c("temname","visit"), all.x=T)
  ra_labs$esr <- round(ra_labs$esr,2)
  ra_labs$crp <- round(ra_labs$crp,2)
  ra_labs_imp <- cbind(ra_labs[, .(temname, visit)], crp=complete(mice(copy(ra_labs)[, c("temname","visit","esr"):=NULL], m=1,method = "pmm",seed=123))$crp, esr=complete(mice(copy(ra_labs)[, c("temname","visit","crp"):=NULL], m=1,method = "pmm",seed=123))$esr)
  ccc.long <- merge(ccc.long, ra_labs_imp, by=c("temname","visit"))
  

  #Create baseline value for time-varing variables
  timevary_b <- ccc.long[visit==0,.(prednisolone_b=prednisolone,nsaid_b=nsaid,opioid_b=opioid,HYDROXYCHLOROQUINE_b=HYDROXYCHLOROQUINE,SULPHASALAZINE_b=SULPHASALAZINE,METHOTREXATE_b=METHOTREXATE,LEFLUNOMIDE_b=LEFLUNOMIDE,AZATHIOPRINE_b=AZATHIOPRINE,CYCLOSPORIN_b=CYCLOSPORIN,esr_b=esr,crp_b=crp),.(ReferenceKey,start,moa,switch)]
  ccc.long <- merge(ccc.long,timevary_b,all.x=T,allow.cartesian = T,by=c("ReferenceKey","start","moa","switch"))
  saveRDS(ccc.long,"msm_data_20250420_corrected.RDS")
  

