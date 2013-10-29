library(foreign)
memory.limit(size = 3500)

#twelve <- subset(df, year == 2012)
region <- ddply(df, .(postcode.derived.region.name, organism.name), summarise, n = sum(uniq),
                amp.tested = sum(ampamox.tested, na.rm = TRUE), 
                amp.resistant = sum(ampamox.resistant, na.rm = TRUE), cla.tested = sum(cla.tested, na.rm = TRUE),
                cla.resistant = sum(cla.resistant, na.rm = TRUE), dox.tested = sum(dox.tested, na.rm = TRUE),
                dox.resistant = sum(dox.resistant, na.rm = TRUE))
head(region)

region$amp.pc <- binom.confint(region$amp.resistant, region$amp.tested, method = "exact")$mean
region$amp.lci <- binom.confint(region$amp.resistant, region$amp.tested, method = "exact")$lower
region$amp.uci <- binom.confint(region$amp.resistant, region$amp.tested, method = "exact")$upper

region$cla.pc <- binom.confint(region$cla.resistant, region$cla.tested, method = "exact")$mean
region$cla.lci <- binom.confint(region$cla.resistant, region$cla.tested, method = "exact")$lower
region$cla.uci <- binom.confint(region$cla.resistant, region$cla.tested, method = "exact")$upper

region$dox.pc <- binom.confint(region$dox.resistant, region$dox.tested, method = "exact")$mean
region$dox.lci <- binom.confint(region$dox.resistant, region$dox.tested, method = "exact")$lower
region$dox.uci <- binom.confint(region$dox.resistant, region$dox.tested, method = "exact")$upper

region <- region[, names(region) != "amp.tested" & names(region) != "amp.resistant" & 
               names(region) != "cla.tested" & names(region) != "cla.resistant" &
               names(region) != "dox.tested" & names(region) != "dox.resistant"]

region$postcode.derived.region.name[region$postcode.derived.region.name == "E MIDS"] <- "East Midlands"
region$postcode.derived.region.name[region$postcode.derived.region.name == "EAST"] <- "East of England"
region$postcode.derived.region.name[region$postcode.derived.region.name == "LONDON"] <- "London"
region$postcode.derived.region.name[region$postcode.derived.region.name == "N EAST"] <- "North East"
region$postcode.derived.region.name[region$postcode.derived.region.name == "N IRELAND"] <- "Northern Ireland"
region$postcode.derived.region.name[region$postcode.derived.region.name == "N WEST"] <- "North West"
region$postcode.derived.region.name[region$postcode.derived.region.name == "S EAST"] <- "South East"
region$postcode.derived.region.name[region$postcode.derived.region.name == "S WEST"] <- "South West"
region$postcode.derived.region.name[region$postcode.derived.region.name == "W MIDS"] <- "West Midlands"
region$postcode.derived.region.name[region$postcode.derived.region.name == "WALES"] <- "Wales"
region$postcode.derived.region.name[region$postcode.derived.region.name == "YORK&HUM"] <- "Yorkshire & Humber"


write.csv(region, file = "regional_resistance.csv", row.names = FALSE)
