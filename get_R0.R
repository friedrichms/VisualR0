### Required Packages
library(mongolite)
require(jsonlite)
#library(ggplot2)
library(tidyverse)
#### Load Data from RKI

json <- "https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.geojson"
inp <- jsonlite::fromJSON(json)
df <- data.frame(inp$features$properties)

### Transform (Lankreise/Bundesland)
data.lk <-df %>% mutate(Datum=as.Date(Meldedatum)) %>% arrange(Landkreis,IdLandkreis,Datum) %>% group_by(Landkreis,IdLandkreis,Datum) %>% distinct() %>% ungroup() %>% arrange(Landkreis,IdLandkreis,Datum) %>% group_by(Landkreis,IdLandkreis,Datum) %>% summarise(ninfected=sum(AnzahlFall)) %>% mutate(infected=cumsum(ninfected))
#data.bl <-df %>% mutate(Datum=as.Date(Meldedatum)) %>% arrange(Bundesland,Datum) %>% group_by(Bundesland,Datum) %>% distinct() %>% ungroup() %>% arrange(Bundesland,Datum) %>% group_by(Bundesland,Datum) %>% summarise(ninfected=sum(AnzahlFall)) %>% mutate(infected=cumsum(ninfected))

### Ignore last Date because of lagged reporting
data.lk <- data.lk[!data.lk$Datum==max(data.lk$Datum),]

#ggplot(data.lk,aes(x=Datum,y=infected,fill=Landkreis))+geom_bar(stat="identity")+theme(legend.position = "none")

### Get a specific province
get_lk<-function(key,data){
  prov.data <- data[data$Landkreis==key,]
  #prov.data <- aggregate(prov.data$infected,FUN=sum,by=list(time=prov.data$time))
  tmin <- min(prov.data$Datum)
  prov.data$day <- as.numeric(prov.data$Datum+1-tmin)
  prov.data
}


### Function calculates the effective rt for a time series
calc_rt <- function(infected,a=3,b=0.3){
  tmax <- length(infected)
  ### Use Gamma dist as Heterogeneity with mean a/b
  rho <- dgamma(1:tmax,a,b)
  mn <- a/b
  ind <- sapply(c(1:tmax),function(x){1:x})
  infm <- infected[unlist(rev(ind))]
  mat <- matrix(0,tmax,tmax)
  mat[lower.tri(mat,diag=T)]<-infm
  mat <- t(mat)
  rt<-infected/colSums(mat*rho)
  rt[1:as.integer(mn)]<-NA
  rt
}

### Analyse RT for a given province
rt <- function(key,data){
  prov.data <- get_lk(key,data)
  if(dim(prov.data)[1]==0){
    return(data.frame(time=NA,day=NA,R0=NA))
  }
  days <- 1:max(prov.data$day)
  infected <- rep(0,max(days))
  infected[prov.data$day]<- prov.data$infected
  out<- calc_rt(infected)
  data.frame(time=prov.data$Datum,day=prov.data$day,R0=out[prov.data$day])
}

### Get present R0 (with a small lag to aviod outliers)
rt_last<-function(key,data,lag=2){
  r0 <- rt(key,data)$R0
  len <- length(r0)
  if(len < lag || all(is.na(r0))){return(data.frame(mr0=NA,tr0=NA))}
  df <- data.frame(r0=r0[(len-lag):len],t=(len-lag):len)
  data.frame(mr0=mean(df$r0),tr0=coef(lm(r0~t,df))[2])
}

### Get all Landkreise
get_r_all<-function(data){
  df <- data.frame()
  id <- unique(data$IdLandkreis)
  i = 1
  for(region in unique(data$Landkreis)){

    #print(id)
    #print(region)
    r<-rt_last(region,data)
    r <- data.frame(Landkreis=region,IdLandkreis=id[i],mr0=r$mr0,tr0=r$tr0)
    #print(r)
    df<-rbind(df,r)
    i=i+1
  }
  df
}
