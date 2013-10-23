library(car) # for deltaMethod
library(plyr)
library(binom)
library(knitr)
library(xtable)
library(lmtest) # for coeftest, for robust std errors
library(survey) # for svycontrast, for interactions. 
library(ggplot2)
library(reshape2)
library(sandwich) # for sandwich function, for robust std errors
library(lubridate) # for quarter function. Q1 is jan-mar

memory.limit(size = 3500)
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

setwd("F:\\antimicrobial resistance\\Simon\\thorax2\\FIS\\poster\\R code")

df <- read.csv("F:\\antimicrobial resistance\\Simon\\thorax2\\SuperlrtallorgFINAL_fixed4_w_region.txt", 
               header = TRUE, sep = "\t", stringsAsFactors = FALSE)

names(df) <- tolower(names(df))
df$earliest.specimen.date <- as.Date(df$earliest.specimen.date, format = "%d/%m/%Y")

df$uniq <- 0 
df$uniq[!duplicated(df$opie.id)] <- 1

df$year <- as.character(format(df$earliest.specimen.date, "%Y"))
df$quarter <- quarter(df$earliest.specimen.date)
df$yrqtr <- paste(df$year, df$quarter, sep = "-")
df$age.grouping <- NULL
#df$lab.name <- NULL
df$specimen.accession.number <- NULL
df$patient.hospital.number <- NULL
df$soundex <- NULL
df$earliest.specimen.week <- NULL
df$earliest.specimen.quarter..based.on.date <- NULL
df$specimen.sequence.number <- NULL
df$earliest.specimen.month.and.year <- NULL
gc()

# Headline figures ####
# number of labs
length(levels(factor(df$lab.name)))
table(df$year, length(levels(factor(df$lab.name))) )
n.labs <- ddply(df, .(year, quarter), summarise, n = length(levels(factor(lab.name))) )
n.labs 
n.labs$yr.qtr <- paste(n.labs$year, n.labs$quarter, sep = "-")
qplot(x = yr.qtr, y = n, data = n.labs, geom = "point")

# specimen type
table(df$specimen.type.description, df$uniq)
sum(df$uniq)
sum(df$uniq[df$specimen.type.description == "SPUTUM"])
round( (sum(df$uniq[df$specimen.type.description == "SPUTUM"])/sum(df$uniq))*100 , 2)

# specimen source
table(df$specimen.source.type.description, df$uniq)
sum(df$uniq[df$specimen.source.type.description == "GENERAL PRACTITIONER"])
round( (sum(df$uniq[df$specimen.source.type.description == "GENERAL PRACTITIONER"])/sum(df$uniq))*100 , 2)

# figures and tables ####
source("fig1.R", echo = FALSE)

# Classify resistances
source("F:\\antimicrobial resistance\\Simon\\thorax2\\non_susceptibility.R")

# Resistance to any recommended ceph ####
abx <- subset(df, antimicrobial.name=="CEFTRIAXONE" | antimicrobial.name=="CEFOTAXIME" | 
                antimicrobial.name=="CEFUROXIME", select = c(opie.id, antimicrobial.name))
#head(abx)
abx <- subset(abx, !duplicated(opie.id))
abx$rec.cef.tested <- 1
sum(abx$rec.cef.tested)
abx <- abx[,c(1,3)]
df <- join(df, abx, type = "left", match = "all")

abx <- subset(df, antimicrobial.name=="CEFTRIAXONE" | antimicrobial.name=="CEFOTAXIME" | 
                antimicrobial.name=="CEFUROXIME", 
              select = c(opie.id, antimicrobial.name, susceptibility.result.description))

abx <- subset(abx, 
              susceptibility.result.description != "SUSCEPTIBLE" & susceptibility.result.description != "",
              select = c(opie.id, antimicrobial.name, susceptibility.result.description))
#head(abx)
abx$rec.cef.resistant <- 1
# abx$uniq <- 0
# abx$uniq[!duplicated(abx$opie.id)] <- 1
abx <- subset(abx, !duplicated(opie.id), select= c(opie.id, rec.cef.resistant))
df <- join(df, abx, type = "left", match = "all")
#head(df)
df$rec.cef.tested[df$uniq!=1] <- 0
df$rec.cef.resistant[df$uniq!=1] <- 0
rm(abx)
gc()

# fig 2 ####
source("fig2.R", echo = FALSE)

# table 1 #### incident rate ratios
source("table1.R", echo = FALSE)

