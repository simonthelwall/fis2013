# age structure figure

df2 <- subset(df, uniq==1)
df2$organism.name <- simpleCap(df2$organism.name)

levels(factor(df2$age.group))
df2$age.group[df2$age.group == "75+ years"] <- "> 75 years"
Encoding(df2$age.group) <- "UTF-8"
df2$age.group <- factor(df2$age.group, levels = c("<1 month", "1-11 months", "1-4 years", "5-9 years", 
                                                  "10-14 years", "15-44 years", "45-64 years", "65-74 years", 
                                                  "> 75 years", "Unknown"))


png("age_fig.png", width = 1530, height = 1530, res = 300)
p <- ggplot(df2, aes(x = age.group)) + geom_histogram() + facet_wrap(~organism.name, ncol = 1) + 
  theme_grey(base_size = 8) + theme(strip.text.x = element_text(face = 'italic'), 
                                    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_x_discrete("Age group", labels = c("<1 month", "1-11 months", "1-4 years", "5-9 years", 
                                           "10-14 years", "15-44 years", "45-64 years", "65-74 years", 
                                           expression(phantom(x) >=75), "Unknown")) + scale_y_continuous("n. isolates")
p
dev.off()
