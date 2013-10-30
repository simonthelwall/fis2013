qtr.trend <- ddply(df, .(organism.name, yrqtr), summarise, isolates = sum(uniq))
qtr.trend$organism.name <- simpleCap(qtr.trend$organism.name)
qtr.trend$period <- as.numeric(factor(qtr.trend$yrqtr))
trend.lm2 <- lm(isolates~factor(organism.name) * period, data = qtr.trend)
qtr.trend$predicted <- predict(trend.lm2, type = "response")
qtr.trend$se <- predict(trend.lm2, type = "response", se.fit = TRUE)$se.fit
qtr.trend$lci <- qtr.trend$predicted - 1.96*qtr.trend$se
qtr.trend$uci <- qtr.trend$predicted + 1.96*qtr.trend$se

cols <- c("Trend" = "#000000", "95 % CI" = "#000000")

p <- ggplot(qtr.trend, aes(x = yrqtr, y = isolates, group = organism.name)) 
p <- p + facet_wrap(~organism.name, ncol = 1, scales = "free_y") + geom_line(aes(y = predicted, colour = "Trend")) + 
  geom_line(aes(y = lci, colour = "95 % CI"), linetype = 2) + geom_line(aes(y = uci, colour = "95 % CI"), linetype = 2)

png("figure1.png", width = 1530, height = 1912, res = 300)
p + geom_line(aes(colour = organism.name), size = 1) + 
  scale_color_manual(values = c("#FF0000", "#00AE9E", "#98002E", "#EEB111", "#000000"),
                     name = "", 
                     breaks = c("Trend", "95 % CI"),
                     labels = c("Trend", "95 % CI")
  ) + 
  theme(#legend.position = "none", #legend.position = "bottom", 
    legend.position = c(0.95,0.7), legend.direction = "horizontal", legend.justification = c(1,0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.x = element_text(face = 'italic')) +
  scale_x_discrete("Year-quarter") + scale_y_continuous("n. isolates")
dev.off()
rm(qtr.trend, trend.lm, p)
