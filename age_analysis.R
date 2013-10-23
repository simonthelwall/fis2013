table(df$age.group)
age <- ddply(df, .(age.group, yrqtr, organism.name), summarise, n = sum(uniq))
age
age$age.group <- factor(age$age.group, levels = 
                          c("<1 month", "1-11 months", "1-4 years", "5-9 years","10-14 years","15-44 years",
                            "45-64 years","65-74 years","75+ years","Unknown"))
p <- ggplot(age, aes(yrqtr, n, group = age.group)) + facet_wrap(~organism.name, ncol = 1, scales = "free_y") + 
  geom_line(aes(colour = age.group))
p

names(df)[names(df)=="age2"] <- "age.gp2"

age2 <- ddply(df, .(organism.name, age.gp2, yrqtr), summarise, amp.tested = sum(ampamox.tested, na.rm = TRUE),
              amp.resistant = sum(ampamox.resistant, na.rm = TRUE), cla.tested = sum(cla.tested, na.rm = TRUE),
              cla.resistant = sum(cla.resistant, na.rm = TRUE), dox.tested = sum(dox.tested, na.rm = TRUE),
              dox.resistant = sum(dox.resistant, na.rm = TRUE))
age2$amp.pc <- binom.confint(age2$amp.resistant, age2$amp.tested, method = "exact")$mean
age2$amp.lci <- binom.confint(age2$amp.resistant, age2$amp.tested, method = "exact")$lower
age2$amp.uci <- binom.confint(age2$amp.resistant, age2$amp.tested, method = "exact")$upper

age2$cla.pc <- binom.confint(age2$cla.resistant, age2$cla.tested, method = "exact")$mean
age2$cla.lci <- binom.confint(age2$cla.resistant, age2$cla.tested, method = "exact")$lower
age2$cla.uci <- binom.confint(age2$cla.resistant, age2$cla.tested, method = "exact")$upper

age2$dox.pc <- binom.confint(age2$dox.resistant, age2$dox.tested, method = "exact")$mean
age2$dox.lci <- binom.confint(age2$dox.resistant, age2$dox.tested, method = "exact")$lower
age2$dox.uci <- binom.confint(age2$dox.resistant, age2$dox.tested, method = "exact")$upper

age2 <- age2[, names(age2) != "amp.tested" & names(age2) != "amp.resistant" & 
               names(age2) != "cla.tested" & names(age2) != "cla.resistant" &
               names(age2) != "dox.tested" & names(age2) != "dox.resistant"]
m.age2 <- melt(age2, id.var = c("organism.name", "age.gp2", "yrqtr"))

m.age2$abx <- sub(".\\w+$", "", m.age2$variable, perl = TRUE)

m.age2$measure <- sub("^\\w\\w\\w.", "", m.age2$variable, perl = TRUE)
m.age2$variable <- NULL
head(m.age2)
age2 <- dcast(m.age2, organism.name + age.gp2 + yrqtr + abx ~ measure)
head(age2)
age2 <- age2[ age2$age.gp2 != "Unknown", ]
p <- ggplot(age2, aes(yrqtr, pc, group = age.gp2)) + geom_line(aes(colour = age.gp2)) + 
  facet_grid(abx ~ organism.name)
p
rm(p, m.age2, age2, age)
