# Figure 2:  trend in resistances by organism

yr.trend <- ddply(df, .(year, organism.name), summarise, 
                  ampamox.tested = sum(ampamox.tested, na.rm = TRUE), 
                  ampamox.resistant = sum(ampamox.resistant, na.rm = TRUE),
                  cla.tested = sum(cla.tested, na.rm = TRUE), cla.resistant = sum(cla.resistant, na.rm = TRUE),
                  dox.tested = sum(dox.tested, na.rm = TRUE), dox.resistant = sum(dox.resistant, na.rm = TRUE),
                  rec_cef.tested = sum(rec.cef.tested, na.rm = TRUE), 
                  rec_cef.resistant = sum(rec.cef.resistant, na.rm = TRUE))

m.yr.trend <- melt(yr.trend)
m.yr.trend$abx <- sub("\\.\\w+$","", m.yr.trend$variable)
m.yr.trend$tested <- grepl("tested", m.yr.trend$variable)
m.yr.trend$tested[m.yr.trend$tested == "TRUE"] <- "tested"
m.yr.trend$resistant <- grepl("resistant", m.yr.trend$variable)
m.yr.trend$resistant[m.yr.trend$resistant == "TRUE"] <- "resistant"
m.yr.trend$variable <- NULL

yr.trend.2 <- dcast(m.yr.trend, year + organism.name + abx ~ tested + resistant, value.var = "value")
names(yr.trend.2) <- c("Year", "Organism", "Antibiotic", "Resistant", "Tested")
yr.trend.2$Organism <- simpleCap(yr.trend.2$Organism)
yr.trend.2$pc.resistant <- binom.confint(yr.trend.2$Resistant, yr.trend.2$Tested, methods = "exact")$mean*100
yr.trend.2$lci <- binom.confint(yr.trend.2$Resistant, yr.trend.2$Tested, methods = "exact")$lower*100
yr.trend.2$uci <- binom.confint(yr.trend.2$Resistant, yr.trend.2$Tested, methods = "exact")$upper*100

# Very low numbers of S. aureus tested against cephalosporins anyway, no point including. 
yr.trend.2$drop <- 0
yr.trend.2$drop[yr.trend.2$Organism == "Staphylococcus aureus" & yr.trend.2$Antibiotic == "rec_cef"] <- 1
yr.trend.2 <- yr.trend.2[yr.trend.2$drop == 0,]
yr.trend.2$drop <- NULL
head(yr.trend.2)

yr.trend.2$Antibiotic[yr.trend.2$Antibiotic == "ampamox"] <- "Ampicillin"
yr.trend.2$Antibiotic[yr.trend.2$Antibiotic == "cla"] <- "Clarithromycin"
yr.trend.2$Antibiotic[yr.trend.2$Antibiotic == "dox"] <- "Doxycylcine"
yr.trend.2$Antibiotic[yr.trend.2$Antibiotic == "rec_cef"] <- "Any recommended \n cephalosporin"
#yr.trend.2$Antibiotic[yr.trend.2$Antibiotic == "Any recommended \\n cephalosporin"] <- "Any recommended \n cephalosporin" #oops

f2 <- ggplot(yr.trend.2, aes(x = Year, y = pc.resistant, group = Antibiotic)) + facet_wrap(~Organism) +
  geom_line(aes(colour = Antibiotic)) + 
  geom_errorbar(aes(ymax = uci, ymin = lci, colour = Antibiotic), width = 0.25)

png("figure2.png", width = 1530, height = 1912, res = 300)
f2 + theme_grey(base_size = 10) + 
  theme(strip.text.x = element_text(face = 'italic'), legend.position = c(0.84,0.86), 
        legend.background = element_rect(fill = "#ffffffaa", colour = NA), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  scale_y_continuous("Per cent resistant", limits = c(0,100))
dev.off()

rm(m.yr.trend, yr.trend, yr.trend.2, f2)
gc()
