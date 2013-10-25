library(foreign)
memory.limit(size = 3500)
pop <- read.dta("F:\\HPR\\UKpopulation1.dta")
head(pop)
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

write.csv(region, file = "regional_resistance.csv", row.names = FALSE)

# m.region <- melt(region, id.var = c("organism.name", "postcode.derived.region.name", "year"))
# 
# m.region$abx <- sub(".\\w+$", "", m.region$variable, perl = TRUE)
# 
# m.region$measure <- sub("^\\w\\w\\w.", "", m.region$variable, perl = TRUE)
# m.region$variable <- NULL
# head(m.region)
# 
# region2 <- dcast(m.region, organism.name + postcode.derived.region.name + year + abx ~ measure)
# head(region2)
# 
# p <- ggplot(region2, aes(year, pc, group = postcode.derived.region.name)) + 
#   geom_line(aes(colour = postcode.derived.region.name)) + facet_grid(abx ~ organism.name)
# p

library(spatial)
library(sp)
library(rgdal)
library(maptools)

#regions <- readShapeSpatial("F:\\HPR\\Enterococci\\2012\\merged_renamed.shp")
regions <- readShapeSpatial("F:\\HPR\\Enterococci\\2012\\merged_renamed_simplified_pr_100.shp")

fregions<-fortify(regions, region = "HPAREGION")
#fregions$country<-"England"
head(fregions)
fregions$piece <- as.numeric(fregions$piece)
#fregions2 <- ddply(fregions, .(piece), function(x)rbind(x, fregions[1, ]))
fregions <- fregions[order(fregions$order, fregions$piece),]
# rm(ewni_map)
# ewni_map <- ggplot(fregions, aes(long, lat, group = id)) + 
#   geom_polygon(colour = I("White"), aes(fill = id, group = id))
# ewni_map

rm(ewni_map)
ewni_map <- ggplot(fregions) + geom_map(map = fregions, aes(long, lat, map_id = id, fill = id))
ewni_map

# load in previously saved 
regions <- read.csv("regional_resistance.csv", header = TRUE, stringsAsFactors = FALSE)
table(regions$postcode.derived.region.name)
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "E MIDS"] <- "East Midlands"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "EAST"] <- "East of England"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "LONDON"] <- "London"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "N EAST"] <- "North East"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "N IRELAND"] <- "Northern Ireland"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "N WEST"] <- "North West"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "S EAST"] <- "South East"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "S WEST"] <- "South West"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "W MIDS"] <- "West Midlands"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "WALES"] <- "Wales"
regions$postcode.derived.region.name[regions$postcode.derived.region.name == "YORK&HUM"] <- "Yorkshire & Humber"

# reshape to long to allow multiple maps in ggplot
rm(regions)
m.regions <- melt(regions, id.vars = c("postcode.derived.region.name", "organism.name"))
head(m.regions)

m.regions <- subset(m.regions, variable != "n")
m.regions <- droplevels(m.regions)
m.regions$abx <- sub(".\\w+$", "", m.regions$variable, perl = TRUE)

m.regions$measure <- sub("^\\w\\w\\w.", "", m.regions$variable, perl = TRUE)
m.regions$variable <- NULL
head(m.regions)
table(m.regions$measure)

#merge geographic and resistance data
names(fregions)
names(m.regions)
names(m.regions)[names(m.regions)=="postcode.derived.region.name"] <- "id"
gc()
fregions2 <- join(fregions, m.regions, type = "left")

head(fregions2)

# testing small section of data
m.regions <- subset(m.regions, abx == "amp" & measure != "lci" & measure != "uci" & measure != "resistant")
regions2 <- gSimplify(regions, tol = 1)
fregions2 <- fortify(regions2, region = "HPAREGION")
ewni_map <- ggplot(fregions2) + geom_map(map = fregions2, aes(long, lat, map_id = id, fill = id))
ewni_map
fregions <- join(fregions, m.regions, type = "left")
head(fregions)
