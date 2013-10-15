library(plyr)
library(xtable)
library(ggplot2)
library(reshape2)
library(lubridate) # for quarter function. Q1 is jan-mar

memory.limit(size = 2500)
options(scipen = 8)

percent <- function(n, d){
  round((sum(n, na.rm = TRUE)/sum(d, na.rm = TRUE))*100,2)
}

percent2 <- function(n, d){
  round((n/d)*100,2)
}

simpleCap<-function(x){
  s<-as.character(x)
  s<-paste(toupper(substring(s,1,1)),tolower(substring(s,2)),sep="")
  #  s<-as.factor(s)
}

setwd("F:\antimicrobial resistance\Simon\thorax2\FIS\poster\R code")

df <- read.csv("F:\\antimicrobial resistance\\Simon\\thorax2SuperlrtallorgFINAL_fixed4.txt", 
               header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# should be 1 125 084 records long
names(df) <- tolower(names(df))
df$earliest.specimen.date <- as.Date(df$earliest.specimen.date, format = "%d/%m/%Y")

df$uniq <- 0 
df$uniq[!duplicated(df$opie.id)] <- 1

df$year <- as.character(format(df$earliest.specimen.date, "%Y"))
df$quarter <- quarter(df$earliest.specimen.date)
df$yrqtr <- paste(df$year, df$quarter, sep = "-")
df$age.grouping <- NULL
df$lab.name <- NULL
df$specimen.accession.number <- NULL
df$patient.hospital.number <- NULL
df$soundex <- NULL
df$earliest.specimen.week <- NULL

source("fig1.R", echo = FALSE)
