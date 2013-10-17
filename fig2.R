# Figure 2:  trend in resistances by organism

yr.trend <- ddply(df, .(year, organism.name), summarise, n = sum(uniq), 
                   ampamox.tested = sum(ampamox.tested, na.rm = TRUE), 
                   ampamox.resistant = sum(ampamox.resistant, na.rm = TRUE),
                   cla.tested = sum(cla.tested, na.rm = TRUE), cla.resistant = sum(cla.resistant, na.rm = TRUE),
                   dox.tested = sum(dox.tested, na.rm = TRUE), dox.resistant = sum(dox.resistant, na.rm = TRUE),
                   rec.cef.tested = sum(rec.cef.tested, na.rm = TRUE), 
                   rec.cef.resistant = sum(rec.cef.resistant, na.rm = TRUE))

yr.trend$pc.ampamox.res <- percent2(yr.trend$ampamox.resistant, yr.trend$ampamox.tested)
yr.trend$pc.ampamox.res.lci <- binom.confint(yr.trend$ampamox.resistant, yr.trend$ampamox.tested, methods = "exact")$lower*100
yr.trend$pc.ampamox.res.uci <- binom.confint(yr.trend$ampamox.resistant, yr.trend$ampamox.tested, methods = "exact")$upper*100

yr.trend$pc.cla.res <- percent2(yr.trend$cla.resistant, yr.trend$cla.tested)
yr.trend$pc.cla.res.lci <- binom.confint(yr.trend$cla.resistant, yr.trend$cla.tested, methods = "exact")$lower*100
yr.trend$pc.cla.res.uci <- binom.confint(yr.trend$cla.resistant, yr.trend$cla.tested, methods = "exact")$upper*100

yr.trend$pc.dox.res <- percent2(yr.trend$dox.resistant, yr.trend$dox.tested)
yr.trend$pc.dox.res.lci <- binom.confint(yr.trend$dox.resistant, yr.trend$dox.tested, methods = "exact")$lower*100
yr.trend$pc.dox.res.uci <- binom.confint(yr.trend$dox.resistant, yr.trend$dox.tested, methods = "exact")$upper*100

yr.trend$pc.rec.cef.res <- percent2(yr.trend$rec.cef.resistant, yr.trend$rec.cef.tested)
yr.trend$pc.rec.cef.res.lci <- binom.confint(yr.trend$rec.cef.resistant, yr.trend$rec.cef.tested, methods = "exact")$lower*100
yr.trend$pc.rec.cef.res.uci <- binom.confint(yr.trend$rec.cef.resistant, yr.trend$rec.cef.tested, methods = "exact")$upper*100

yr.trend <- subset(yr.trend, select = c("year", "organism.name", "pc.ampamox.res", "pc.cla.res", "pc.dox.res", 
                                          "pc.rec.cef.res", "pc.ampamox.res.lci", "pc.ampamox.res.uci",
                                        "pc.cla.res.lci", "pc.cla.res.uci", "pc.dox.res.lci", "pc.dox.res.uci",
                                        "pc.rec.cef.res.lci", "pc.rec.cef.res.uci"))

yr.trend
m.yr.trend <- melt(yr.trend)
m.yr.trend$sa.cef <- 0
m.yr.trend$sa.cef[(m.yr.trend$organism.name == "STAPHYLOCOCCUS AUREUS" & m.yr.trend$variable == "pc.rec.cef.res") |
                    (m.yr.trend$organism.name == "STAPHYLOCOCCUS AUREUS" & m.yr.trend$variable == "pc.rec.cef.res.lci") |  
                    (m.yr.trend$organism.name == "STAPHYLOCOCCUS AUREUS" & m.yr.trend$variable == "pc.rec.cef.res.uci")] <- 1

m.yr.trend <- subset(m.yr.trend, sa.cef==0)
m.yr.trend$sa.cef <- NULL

f2 <- ggplot(m.yr.trend, aes(x = year, y = value, group = variable)) + facet_wrap(~organism.name) +
  geom_line(aes(colour = variable)) + geom_errorbar(ymax = )
f2
