qtr.trend <- ddply(df, .(organism.name, yrqtr), summarise, isolates = sum(uniq))
qtr.trend$organism.name <- simpleCap(qtr.trend$organism.name)
qtr.trend$period <- as.numeric(factor(qtr.trend$yrqtr))
trend.lm <- glm(isolates~factor(organism.name) + period, data = qtr.trend, family = "poisson")
qtr.trend$predicted <- predict(trend.lm, type = "response")
qtr.trend$se <- predict(trend.lm, type = "response", se.fit = TRUE)$se.fit
qtr.trend$lci <- qtr.trend$predicted - 1.96*qtr.trend$se
qtr.trend$uci <- qtr.trend$predicted + 1.96*qtr.trend$se
p <- ggplot(qtr.trend, aes(x = yrqtr, y = isolates, group = organism.name)) 
p <- p + facet_wrap(~organism.name, ncol = 1, scales = "free_y") + geom_line(aes(y = predicted)) + 
  geom_line(aes(y = lci), linetype = 2) + geom_line(aes(y = uci), linetype = 2)
png("figure1.png", width = 1530, height = 1912, res = 300)
p + geom_line(aes(colour = organism.name), size = 1) + 
  theme(legend.position = "none", #legend.position = "bottom", 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(face = 'italic')) +
  scale_x_discrete("Year-quarter") + scale_y_continuous("n. isolates")
dev.off()
rm(qtr.trend, trend.lm, p)
