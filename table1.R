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
# complete for variables to be adjusted
# organism, quarter, region, specimen source, male sex + age2
ampamox <- ampamox[complete.cases(ampamox[,c(2, 16, 17, 14, 9, 46, 45)]),] 
ampamox$ampamox.resistant[is.na(ampamox$ampamox.resistant)==TRUE] <- 0
length(ampamox$ampamox.tested)
table(ampamox$ampamox.resistant, useNA = "ifany")
ampamox.m <- glm(ampamox.resistant ~ year + factor(quarter) + postcode.derived.region.name +
                   specimen.source.type.description + male + age2, data = ampamox, family = "poisson")
summary(ampamox.m)
#orWrapper(ampamox.m)

ampamox.m <- glm(ampamox.resistant ~ year + organism.name + factor(quarter) + postcode.derived.region.name +
                   specimen.source.type.description + male + age2, data = ampamox, family = "poisson")

ampamox.m2 <- glm(ampamox.resistant ~ year * organism.name + factor(quarter) + postcode.derived.region.name +
                    specimen.source.type.description + male + age2, data = ampamox, family = "poisson")
epicalc::lrtest(ampamox.m, ampamox.m2)

svycontrast(ampamox.m2, c("year" = 1, "year:organism.nameSTAPHYLOCOCCUS AUREUS" =1))

ampamox.robust <- funinteff(ampamox.m2, "year") # worth noting that due to grep in funinteff estimates for 
# age in years gets pulled out too. This can be safely dropped. 
ampamox.robust$z <- ampamox.robust$Estimate/ampamox.robust$SE
ampamox.robust$p <- 2*pnorm(-abs(ampamox.robust$z))
rm(ampamox)

# cla
cla <- subset(df, cla.tested == 1)
cla <- cla[complete.cases(cla[,c(2, 16, 17, 14, 9, 46, 45)]),]  # complete for variables to be adjusted
cla$cla.resistant[is.na(cla$cla.resistant)==TRUE] <- 0
length(cla$cla.tested)
table(cla$cla.resistant, useNA = "ifany")

cla.m <- glm(cla.resistant ~ year + organism.name + factor(quarter) + postcode.derived.region.name +
               specimen.source.type.description + male + age2, data = cla, family = "poisson")

cla.m2 <- glm(cla.resistant ~ year * organism.name + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male + age2, data = cla, family = "poisson")
epicalc::lrtest(cla.m, cla.m2)

cla.robust <- funinteff(cla.m2, "year")
cla.robust$z <- cla.robust$Estimate/cla.robust$SE
cla.robust$p <- 2*pnorm(-abs(cla.robust$z))
cla.robust
rm(cla.m, cla)

# dox
dox <- subset(df, dox.tested == 1)
dox <- dox[complete.cases(dox[,c(2, 16, 17, 14, 9, 46, 45)]),]  # complete for variables to be adjusted
dox$dox.resistant[is.na(dox$dox.resistant)==TRUE] <- 0
length(dox$dox.tested)
table(dox$dox.resistant, useNA = "ifany")

dox.m <- glm(dox.resistant ~ year + organism.name + factor(quarter) + postcode.derived.region.name +
               specimen.source.type.description + male + age2, data = dox, family = "poisson")

dox.m2 <- glm(dox.resistant ~ year * organism.name + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male + age2, data = dox, family = "poisson")
epicalc::lrtest(dox.m, dox.m2)

dox.robust <- funinteff(dox.m2, "year")
dox.robust$z <- dox.robust$Estimate/dox.robust$SE
dox.robust$p <- 2*pnorm(-abs(dox.robust$z))
dox.robust
rm(dox.m, dox)

# cef
cef <- subset(df, rec.cef.tested == 1)
cef <- cef[complete.cases(cef[,c(2, 16, 17, 14, 9, 46, 45)]),]  # complete for variables to be adjusted
cef$rec.cef.resistant[is.na(cef$rec.cef.resistant)==TRUE] <- 0
length(cef$rec.cef.tested)
table(cef$rec.cef.resistant, useNA = "ifany")

cef.m <- glm(rec.cef.resistant ~ year + organism.name + factor(quarter) + postcode.derived.region.name +
               specimen.source.type.description + male + age2, data = cef, family = "poisson")

cef.m2 <- glm(rec.cef.resistant ~ year * organism.name + factor(quarter) + postcode.derived.region.name +
                specimen.source.type.description + male + age2, data = cef, family = "poisson")
epicalc::lrtest(cef.m, cef.m2)

cef.robust <- funinteff(cef.m2, "year")
cef.robust$z <- cef.robust$Estimate/cef.robust$SE
cef.robust$p <- 2*pnorm(-abs(cef.robust$z))
cef.robust
rm(cef.m, cef)

# combine ####
ampamox.robust$abx <- "Ampicillin"
cef.robust$abx <- "Any recommended cephalosporin"
dox.robust$abx <- "Doxycycline"
cla.robust$abx <- "Clarithromycin"

# dropping mis-grepped rows
ampamox.robust <- ampamox.robust[c(1,5,6),]
cla.robust <- cla.robust[c(1,5,6),]
dox.robust <- dox.robust[c(1,5,6),]
cef.robust <- cef.robust[c(1,5,6),]

col <- rbind(ampamox.robust, cef.robust, dox.robust, cla.robust)
col$organism <- row.names(col)
row.names(col) <- seq(1,length(col$organism),1)
col$organism[col$organism == "year" | col$organism == "year1" | col$organism == "year2" | 
               col$organism == "year3"] <- "H. influenzae"
col$organism[col$organism == "year:organism.nameSTAPHYLOCOCCUS AUREUS" | 
               col$organism == "year:organism.nameSTAPHYLOCOCCUS AUREUS1" | 
               col$organism == "year:organism.nameSTAPHYLOCOCCUS AUREUS2" | 
               col$organism == "year:organism.nameSTAPHYLOCOCCUS AUREUS3"] <- "S. aureus"
col$organism[col$organism == "year:organism.nameSTREPTOCOCCUS PNEUMONIAE" | 
               col$organism == "year:organism.nameSTREPTOCOCCUS PNEUMONIAE1" | 
               col$organism == "year:organism.nameSTREPTOCOCCUS PNEUMONIAE2" | 
               col$organism == "year:organism.nameSTREPTOCOCCUS PNEUMONIAE3"] <- "S. pneumoniae"

col$irr <- exp(col$Estimate)
col$lci <- col$irr / exp(1.96*col$SE)
col$uci <- col$irr * exp(1.96*col$SE)
col <- col[,names(col)=="abx" | names(col)=="organism" | names(col)=="irr" | names(col)=="lci" | 
             names(col)=="uci" | names(col)=="p"]
col
col$irr95ci <- paste(sprintf("%.2f", round(col$irr, 2)), " (", sprintf("%.2f", round(col$lci, 2)),
                     " - ", sprintf("%.2f", round(col$uci, 2)), ")", sep = ""
)
col2 <- col[, names(col)=="abx" | names(col)=="organism" | names(col)=="irr95ci"]

c.col <- dcast(col2, organism ~ abx)
c.col
names(c.col)[1] <- "Organism"
c.col <- xtable(c.col, 
  caption = "Incident rate ratios for annual change in susceptibility to recommended antimicrobials, England, Wales and Northern Ireland 2008-2013.
   Adjusted for calendar quarter, laboratory region, specimen source, male sex and patient age.", 
                label = "table1")
sink("table1.txt")
print(c.col, include.rownames = FALSE, booktabs = TRUE)
sink()
rm(ampamox.robust, cef.robust, cla.robust, dox.robust, col, col2, c.col, cef.m2, cla.m2, ampamox.m2, dox.m2, 
   ampamox.m)
