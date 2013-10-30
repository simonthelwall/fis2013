# Table 1: Incident rate ratios for resistance

# functions for extracting IRRs ####
orCalc <- function(outcome, riskf){
  # x should be outcome, y stratifying variable. 
  # cribbed from mhodds in epicalc package
  # with reference to p157 + p164 of Kirkewood and Sterne, Essential Medical Statistics. 2nd Ed
  tab <- table(riskf, outcome)
  #print(tab)
  or <- c(1:dim(tab)[1]) # create vector of same length as table rows. 
  se.log.or <- c(1:dim(tab)[1])
  lci <- c(1:dim(tab)[1])
  uci <- c(1:dim(tab)[1])
  z <- c(1:dim(tab)[1])
  p <- c(1:dim(tab)[1])
  for (i in 1:dim(tab)[1]){
    or[i] <- (tab[1,1]*tab[i,2])/(tab[1,2]*tab[i,1])
    se.log.or[i] <-sqrt(1/tab[1,1] + 1/tab[1,2] + 1/tab[i,1] + 1/tab[i,2])
    lci[i] <- or[i]/exp(1.96 * se.log.or[i])
    uci[i] <- or[i]*exp(1.96 * se.log.or[i])
    z[i] <- log(or[i])/se.log.or[i]
    p[i] <- 2*(1-pnorm(abs(z[i])))
  }
  m <- as.matrix(cbind(or, se.log.or, lci, uci, z, p))
  return(m)
}

lincom <- function(svycontrast_object){
  require(survey)
  if (class(svycontrast_object)=="svystat"){
    or <- exp(svycontrast_object[1])
    lci <- or / exp(1.96*sqrt(attributes(svycontrast_object)$var))
    uci <- or * exp(1.96*sqrt(attributes(svycontrast_object)$var))
    return(as.list(c(or, lci, uci)))
  } else {
    print("Requires object of class svystat")
  }
}

sehac<-function(fit,vcov=sandwich){ #Convenience function for robust standard errors
  coeftest(fit,vcov)
}

funinteff<-function(mod,var,vcov=sandwich){ #mod is an lm() object, var is the name of the main effect that was interacted, vcov is the type of variance covariance method you want to use 
  #Extract Coefficient names create 'beta names' to feed to deltaMethod()
  cnames<-coef(mod)
  pnams<-data.frame('b'=paste('b',0:(length(cnames)-1),sep=""),'est'=cnames) #assign parameter names so that deltaMethod does not throw an error
  
  #Extract the specific parameters of interest
  vars<-grep(var,names(cnames),value=T)
  var1<-vars[1]
  intvars<-vars[2:length(vars)]
  bi<-pnams[var1,'b']
  
  #--Create Data Frame to store Main Effect
  int<-sehac(mod,vcov=vcov)[var1,c('Estimate','Std. Error')]
  int<-as.data.frame(t(int))
  names(int)<-c('Estimate','SE')
  row.names(int)<-var1
  
  #Loop through and store the results in a data.frame
  for(i in 1:length(intvars)){
    bint<-pnams[intvars[i],'b']
    eq<-paste(bi,bint,sep="+")
    interac<-deltaMethod(mod,eq,parameterNames=pnams[,1],vcov=vcov)
    row.names(interac)<-intvars[i]
    
    int<-rbind(int,interac)
  }
  return(int)
}


#  regressions ####
ampamox <- subset(df, ampamox.tested == 1 & uniq == 1)
ampamox <- ampamox[complete.cases(ampamox[,c(2, 18, 14, 10, 43, 44)]),] 
ampamox$ampamox.resistant[is.na(ampamox$ampamox.resistant)==TRUE] <- 0
length(ampamox$ampamox.tested)
table(ampamox$ampamox.resistant, useNA = "ifany")

ampamox.m <- glm(ampamox.resistant ~ year * organism.name + age2 + factor(quarter) + postcode.derived.region.name +
                   specimen.source.type.description + male, data = ampamox, family = "poisson")

ampamox.m2 <- glm(ampamox.resistant ~ year * organism.name * age2 + factor(quarter) + postcode.derived.region.name +
                    specimen.source.type.description + male, data = ampamox, family = "poisson")
epicalc::lrtest(ampamox.m, ampamox.m2)

ampamox.robust <- funinteff(ampamox.m, "year") # worth noting that due to grep in funinteff estimates for 
# age in years gets pulled out too. These can be safely dropped. 
ampamox.robust$z <- ampamox.robust$Estimate/ampamox.robust$SE
ampamox.robust$p <- 2*pnorm(-abs(ampamox.robust$z))
rm(ampamox, ampamox.m, ampamox.m2)


# cla ####
cla <- subset(df, cla.tested == 1 & uniq == 1)
cla$cla.resistant[is.na(cla$cla.resistant)==TRUE] <- 0
length(cla$cla.tested)
table(cla$cla.resistant, useNA = "ifany")

# H. influenzae
cla.hi <- cla[cla$organism.name == "HAEMOPHILUS INFLUENZAE",] # not using subset due to weird R things happening with subset.
cla.hi <- cla.hi[complete.cases(cla.hi[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(cla.hi$cla.tested)
table(cla.hi$cla.resistant, useNA = "ifany")

cla.m <- glm(cla.resistant ~ year + age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = cla.hi, family = "poisson")

cla.m.simple <- glm(cla.resistant ~ year, data = cla.hi, family = "poisson")


cla.m2 <- glm(cla.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = cla.hi, family = "poisson")

cla.hi.robust <- funinteff(cla.m2, "year")
cla.hi.robust$z <- cla.hi.robust$Estimate/cla.hi.robust$SE
cla.hi.robust$p <- 2*pnorm(-abs(cla.hi.robust$z))
cla.hi.robust
rm(cla.m2, cla.hi)

# S. aureus
cla.sa <- cla[cla$organism.name == "STAPHYLOCOCCUS AUREUS",]
#cla.sa <- subset(cla, organism.name == "STREPTOCOCCUS PNEUMONIAE")
cla.sa <- cla.sa[complete.cases(cla.sa[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(cla.sa$cla.tested)
table(cla.sa$cla.resistant, useNA = "ifany")

cla.m2 <- glm(cla.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = cla.sa, family = "poisson")

cla.sa.robust <- funinteff(cla.m2, "year")
cla.sa.robust$z <- cla.sa.robust$Estimate/cla.sa.robust$SE
cla.sa.robust$p <- 2*pnorm(-abs(cla.sa.robust$z))
cla.sa.robust
rm(cla.m2, cla.sa)

# S. pneumoniae
cla.sp <- cla[cla$organism.name == "STREPTOCOCCUS PNEUMONIAE",]
#cla.sp <- subset(cla, organism.name == "STREPTOCOCCUS PNEUMONIAE")
cla.sp <- cla.sp[complete.cases(cla.sp[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(cla.sp$cla.tested)
table(cla.sp$cla.resistant, useNA = "ifany")

cla.m2 <- glm(cla.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = cla.sp, family = "poisson")

cla.sp.robust <- funinteff(cla.m2, "year")
cla.sp.robust$z <- cla.sp.robust$Estimate/cla.sp.robust$SE
cla.sp.robust$p <- 2*pnorm(-abs(cla.sp.robust$z))
cla.sp.robust
rm(cla.m2, cla.sp, cla)


# dox ####
dox <- subset(df, dox.tested == 1 & uniq == 1)
dox$dox.resistant[is.na(dox$dox.resistant)==TRUE] <- 0
length(dox$dox.tested)
table(dox$dox.resistant, useNA = "ifany")

# H. influenzae
dox.hi <- dox[dox$organism.name == "HAEMOPHILUS INFLUENZAE",]
#dox.hi <- subset(dox, organism.name == "HAEMOPHILUS INFLUENZAE")
dox.hi <- dox.hi[complete.cases(dox.hi[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(dox.hi$dox.tested)
table(dox.hi$dox.resistant, useNA = "ifany")

dox.m2 <- glm(dox.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = dox.hi, family = "poisson")

dox.hi.robust <- funinteff(dox.m2, "year")
dox.hi.robust$z <- dox.hi.robust$Estimate/dox.hi.robust$SE
dox.hi.robust$p <- 2*pnorm(-abs(dox.hi.robust$z))
dox.hi.robust
rm(dox.m2, dox.hi)

# S. aureus
dox.sa <- dox[dox$organism.name == "STAPHYLOCOCCUS AUREUS",]
#dox.sa <- subset(dox, organism.name == "STAPHYLOCOCCUS AUREUS")
dox.sa <- dox.sa[complete.cases(dox.sa[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(dox.sa$dox.tested)
table(dox.sa$dox.resistant, useNA = "ifany")

dox.m2 <- glm(dox.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = dox.sa, family = "poisson")

dox.sa.robust <- funinteff(dox.m2, "year")
dox.sa.robust$z <- dox.sa.robust$Estimate/dox.sa.robust$SE
dox.sa.robust$p <- 2*pnorm(-abs(dox.sa.robust$z))
dox.sa.robust
rm(dox.m2, dox.sa)

# S. pneumoniae
dox.sp <- dox[dox$organism.name == "STREPTOCOCCUS PNEUMONIAE",]
#dox.sp <- subset(dox, organism.name == "STREPTOCOCCUS PNEUMONIAE")
dox.sp <- dox.sp[complete.cases(dox.sp[,c(18, 14, 10, 43, 44)]),]  # complete for variables to be adjusted
length(dox.sp$dox.tested)
table(dox.sp$dox.resistant, useNA = "ifany")

dox.m2 <- glm(dox.resistant ~ year * age2 + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male, data = dox.sp, family = "poisson")

dox.sp.robust <- funinteff(dox.m2, "year")
dox.sp.robust$z <- dox.sp.robust$Estimate/dox.sp.robust$SE
dox.sp.robust$p <- 2*pnorm(-abs(dox.sp.robust$z))
dox.sp.robust
rm(dox.m2, dox.sp, dox)

# new combine ####
abx <- c("amp", rep("cla", 5), rep("dox",5))
age <- c("-", "$<45$", "45-64", "65-74", "$\\geq 75$", "Unknown", "$<45$", "45-64", "65-74", 
         "$\\geq 75$", "Unknown")
h.influenzae <- rep(NA, 11)
s.aureus <- rep(NA, 11)
s.pneumoniae <- rep(NA, 11)
col <- as.data.frame(cbind(abx, age, h.influenzae, s.aureus, s.pneumoniae), stringsAsFactors = FALSE)
col
rm(age, abx, h.influenze, s.aureus, s.pneumoniae)

# amp
col$h.influenzae[col$abx=="amp"] <- paste(sprintf("%.2f", round(exp(ampamox.robust$Estimate[1]), 2)), " (",
                sprintf("%.2f", round(exp(ampamox.robust$Estimate[1])/exp(1.96*ampamox.robust$SE[1]) ,2)), "-",
                sprintf("%.2f", round(exp(ampamox.robust$Estimate[1])*exp(1.96*ampamox.robust$SE[1]) ,2)), ")",
                                          sep = "")

col$s.aureus[col$abx=="amp"] <- paste(sprintf("%.2f", round(exp(ampamox.robust$Estimate[5]), 2)), " (",
                                          sprintf("%.2f", round(exp(ampamox.robust$Estimate[5])/exp(1.96*ampamox.robust$SE[5]) ,2)), "-",
                                          sprintf("%.2f", round(exp(ampamox.robust$Estimate[5])*exp(1.96*ampamox.robust$SE[5]) ,2)), ")",
                                          sep = "")

col$s.pneumoniae[col$abx=="amp"] <- paste(sprintf("%.2f", round(exp(ampamox.robust$Estimate[6]), 2)), " (",
                                      sprintf("%.2f", round(exp(ampamox.robust$Estimate[6])/exp(1.96*ampamox.robust$SE[6]) ,2)), "-",
                                      sprintf("%.2f", round(exp(ampamox.robust$Estimate[6])*exp(1.96*ampamox.robust$SE[6]) ,2)), ")",
                                      sep = "")

# cla
col$h.influenzae[col$abx=="cla" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(cla.hi.robust$Estimate[1]), 2)), " (",
                                          sprintf("%.2f", round(exp(cla.hi.robust$Estimate[1])/exp(1.96*cla.hi.robust$SE[1]) ,2)), "-",
                                          sprintf("%.2f", round(exp(cla.hi.robust$Estimate[1])*exp(1.96*cla.hi.robust$SE[1]) ,2)), ")",
                                          sep = "")

col$h.influenzae[col$abx=="cla" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(cla.hi.robust$Estimate[5]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[5])/exp(1.96*cla.hi.robust$SE[5]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[5])*exp(1.96*cla.hi.robust$SE[5]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="cla" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(cla.hi.robust$Estimate[6]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[6])/exp(1.96*cla.hi.robust$SE[6]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[6])*exp(1.96*cla.hi.robust$SE[6]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="cla" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(cla.hi.robust$Estimate[7]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[7])/exp(1.96*cla.hi.robust$SE[7]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.hi.robust$Estimate[7])*exp(1.96*cla.hi.robust$SE[7]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="cla" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(cla.hi.robust$Estimate[8]), 2)), " (",
                                                                    sprintf("%.2f", round(exp(cla.hi.robust$Estimate[8])/exp(1.96*cla.hi.robust$SE[8]) ,2)), "-",
                                                                    sprintf("%.2f", round(exp(cla.hi.robust$Estimate[8])*exp(1.96*cla.hi.robust$SE[8]) ,2)), ")",
                                                                    sep = "")
# cla s. aureus
col$s.aureus[col$abx=="cla" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(cla.sa.robust$Estimate[1]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[1])/exp(1.96*cla.sa.robust$SE[1]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[1])*exp(1.96*cla.sa.robust$SE[1]) ,2)), ")",
                                                               sep = "")

col$s.aureus[col$abx=="cla" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(cla.sa.robust$Estimate[5]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[5])/exp(1.96*cla.sa.robust$SE[5]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[5])*exp(1.96*cla.sa.robust$SE[5]) ,2)), ")",
                                                               sep = "")

col$s.aureus[col$abx=="cla" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(cla.sa.robust$Estimate[6]), 2)), " (",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[6])/exp(1.96*cla.sa.robust$SE[6]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(cla.sa.robust$Estimate[6])*exp(1.96*cla.sa.robust$SE[6]) ,2)), ")",
                                                               sep = "")

col$s.aureus[col$abx=="cla" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(cla.sa.robust$Estimate[7]), 2)), " (",
                                                                    sprintf("%.2f", round(exp(cla.sa.robust$Estimate[7])/exp(1.96*cla.sa.robust$SE[7]) ,2)), "-",
                                                                    sprintf("%.2f", round(exp(cla.sa.robust$Estimate[7])*exp(1.96*cla.sa.robust$SE[7]) ,2)), ")",
                                                                    sep = "")

col$s.aureus[col$abx=="cla" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(cla.sa.robust$Estimate[8]), 2)), " (",
                                                                 sprintf("%.2f", round(exp(cla.sa.robust$Estimate[8])/exp(1.96*cla.sa.robust$SE[8]) ,2)), "-",
                                                                 sprintf("%.2f", round(exp(cla.sa.robust$Estimate[8])*exp(1.96*cla.sa.robust$SE[8]) ,2)), ")",
                                                                 sep = "")

# cla s. pneumo
col$s.pneumoniae[col$abx=="cla" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(cla.sp.robust$Estimate[1]), 2)), " (",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[1])/exp(1.96*cla.sp.robust$SE[1]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[1])*exp(1.96*cla.sp.robust$SE[1]) ,2)), ")",
                                                           sep = "")

col$s.pneumoniae[col$abx=="cla" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(cla.sp.robust$Estimate[5]), 2)), " (",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[5])/exp(1.96*cla.sp.robust$SE[5]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[5])*exp(1.96*cla.sp.robust$SE[5]) ,2)), ")",
                                                           sep = "")

col$s.pneumoniae[col$abx=="cla" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(cla.sp.robust$Estimate[6]), 2)), " (",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[6])/exp(1.96*cla.sp.robust$SE[6]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(cla.sp.robust$Estimate[6])*exp(1.96*cla.sp.robust$SE[6]) ,2)), ")",
                                                           sep = "")

col$s.pneumoniae[col$abx=="cla" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(cla.sp.robust$Estimate[7]), 2)), " (",
                                                                sprintf("%.2f", round(exp(cla.sp.robust$Estimate[7])/exp(1.96*cla.sp.robust$SE[7]) ,2)), "-",
                                                                sprintf("%.2f", round(exp(cla.sp.robust$Estimate[7])*exp(1.96*cla.sp.robust$SE[7]) ,2)), ")",
                                                                sep = "")

col$s.pneumoniae[col$abx=="cla" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(cla.sp.robust$Estimate[8]), 2)), " (",
                                                             sprintf("%.2f", round(exp(cla.sp.robust$Estimate[8])/exp(1.96*cla.sp.robust$SE[8]) ,2)), "-",
                                                             sprintf("%.2f", round(exp(cla.sp.robust$Estimate[8])*exp(1.96*cla.sp.robust$SE[8]) ,2)), ")",
                                                             sep = "")


# dox
col$h.influenzae[col$abx=="dox" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(dox.hi.robust$Estimate[1]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[1])/exp(1.96*dox.hi.robust$SE[1]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[1])*exp(1.96*dox.hi.robust$SE[1]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="dox" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(dox.hi.robust$Estimate[5]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[5])/exp(1.96*dox.hi.robust$SE[5]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[5])*exp(1.96*dox.hi.robust$SE[5]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="dox" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(dox.hi.robust$Estimate[6]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[6])/exp(1.96*dox.hi.robust$SE[6]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.hi.robust$Estimate[6])*exp(1.96*dox.hi.robust$SE[6]) ,2)), ")",
                                                               sep = "")

col$h.influenzae[col$abx=="dox" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(dox.hi.robust$Estimate[7]), 2)), " (",
                                                                    sprintf("%.2f", round(exp(dox.hi.robust$Estimate[7])/exp(1.96*dox.hi.robust$SE[7]) ,2)), "-",
                                                                    sprintf("%.2f", round(exp(dox.hi.robust$Estimate[7])*exp(1.96*dox.hi.robust$SE[7]) ,2)), ")",
                                                                    sep = "")

col$h.influenzae[col$abx=="dox" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(dox.hi.robust$Estimate[8]), 2)), " (",
                                                                 sprintf("%.2f", round(exp(dox.hi.robust$Estimate[8])/exp(1.96*dox.hi.robust$SE[8]) ,2)), "-",
                                                                 sprintf("%.2f", round(exp(dox.hi.robust$Estimate[8])*exp(1.96*dox.hi.robust$SE[8]) ,2)), ")",
                                                                 sep = "")
# dox s. aureus
col$s.aureus[col$abx=="dox" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(dox.sa.robust$Estimate[1]), 2)), " (",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[1])/exp(1.96*dox.sa.robust$SE[1]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[1])*exp(1.96*dox.sa.robust$SE[1]) ,2)), ")",
                                                           sep = "")

col$s.aureus[col$abx=="dox" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(dox.sa.robust$Estimate[5]), 2)), " (",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[5])/exp(1.96*dox.sa.robust$SE[5]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[5])*exp(1.96*dox.sa.robust$SE[5]) ,2)), ")",
                                                           sep = "")

col$s.aureus[col$abx=="dox" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(dox.sa.robust$Estimate[6]), 2)), " (",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[6])/exp(1.96*dox.sa.robust$SE[6]) ,2)), "-",
                                                           sprintf("%.2f", round(exp(dox.sa.robust$Estimate[6])*exp(1.96*dox.sa.robust$SE[6]) ,2)), ")",
                                                           sep = "")

col$s.aureus[col$abx=="dox" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(dox.sa.robust$Estimate[7]), 2)), " (",
                                                                sprintf("%.2f", round(exp(dox.sa.robust$Estimate[7])/exp(1.96*dox.sa.robust$SE[7]) ,2)), "-",
                                                                sprintf("%.2f", round(exp(dox.sa.robust$Estimate[7])*exp(1.96*dox.sa.robust$SE[7]) ,2)), ")",
                                                                sep = "")

col$s.aureus[col$abx=="dox" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(dox.sa.robust$Estimate[8]), 2)), " (",
                                                             sprintf("%.2f", round(exp(dox.sa.robust$Estimate[8])/exp(1.96*dox.sa.robust$SE[8]) ,2)), "-",
                                                             sprintf("%.2f", round(exp(dox.sa.robust$Estimate[8])*exp(1.96*dox.sa.robust$SE[8]) ,2)), ")",
                                                             sep = "")

# dox s. pneumo
col$s.pneumoniae[col$abx=="dox" & col$age == "$<45$"] <- paste(sprintf("%.2f", round(exp(dox.sp.robust$Estimate[1]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[1])/exp(1.96*dox.sp.robust$SE[1]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[1])*exp(1.96*dox.sp.robust$SE[1]) ,2)), ")",
                                                               sep = "")

col$s.pneumoniae[col$abx=="dox" & col$age == "45-64"] <- paste(sprintf("%.2f", round(exp(dox.sp.robust$Estimate[5]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[5])/exp(1.96*dox.sp.robust$SE[5]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[5])*exp(1.96*dox.sp.robust$SE[5]) ,2)), ")",
                                                               sep = "")

col$s.pneumoniae[col$abx=="dox" & col$age == "65-74"] <- paste(sprintf("%.2f", round(exp(dox.sp.robust$Estimate[6]), 2)), " (",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[6])/exp(1.96*dox.sp.robust$SE[6]) ,2)), "-",
                                                               sprintf("%.2f", round(exp(dox.sp.robust$Estimate[6])*exp(1.96*dox.sp.robust$SE[6]) ,2)), ")",
                                                               sep = "")

col$s.pneumoniae[col$abx=="dox" & col$age == "$\\geq 75$"] <- paste(sprintf("%.2f", round(exp(dox.sp.robust$Estimate[7]), 2)), " (",
                                                                    sprintf("%.2f", round(exp(dox.sp.robust$Estimate[7])/exp(1.96*dox.sp.robust$SE[7]) ,2)), "-",
                                                                    sprintf("%.2f", round(exp(dox.sp.robust$Estimate[7])*exp(1.96*dox.sp.robust$SE[7]) ,2)), ")",
                                                                    sep = "")

col$s.pneumoniae[col$abx=="dox" & col$age == "Unknown"] <- paste(sprintf("%.2f", round(exp(dox.sp.robust$Estimate[8]), 2)), " (",
                                                                 sprintf("%.2f", round(exp(dox.sp.robust$Estimate[8])/exp(1.96*dox.sp.robust$SE[8]) ,2)), "-",
                                                                 sprintf("%.2f", round(exp(dox.sp.robust$Estimate[8])*exp(1.96*dox.sp.robust$SE[8]) ,2)), ")",
                                                                 sep = "")


col

# Tidy up
names(col) <- c("Antibiotic", "\\shortstack{Age \\\\group}", "\\textit{H. influenzae}", "\\textit{S. aureus}", 
                "\\textit{S. pneumoniae}")

col$Antibiotic <- c("\\shortstack{Ampicillin/\\\\amoxicillin}", "Clarithromycin", "", "", "", "", 
                    "Doxycycline", "", "", "", "" )
col[11, 5] <- "-" # remove NA for S. pneumo&dox as no observations.

sum(df$uniq[df$organism.name=="STREPTOCOCCUS PNEUMONIAE" & df$age2 == "Unknown"])
sum(df$dox.tested[df$organism.name=="STREPTOCOCCUS PNEUMONIAE" & df$age2 == "Unknown"], na.rm = TRUE)

col

col <- xtable(col, 
                caption = "Incident rate ratios for annual change in susceptibility to recommended antimicrobials, England, Wales and Northern Ireland 2008-2013.
   Adjusted for calendar quarter, laboratory region, specimen source, male sex and patient age.", 
                label = "table1")
sink("table1.txt")
print(col)
sink()

rm(ampamox.robust, cla.hi.robust, cla.m, cla.m.simple, cla.sa.robust, cla.sp.robust, col, dox.hi.robust, 
   dox.sa.robust, dox.sp.robust, funinteff, h.influenzae, lincom, orCalc, percent, percent2, sehac, simpleCap, 
   yr.trend)
