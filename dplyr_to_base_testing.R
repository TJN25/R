df1 <- data.frame(a = c("1", "2", "3", "1", "2", "3"), b = c(1, 2, 3, 4, 5, 6))
df2 <- data.frame(a = c("1", "2", "3"), c = c(11, 12, 13))

##left_join
df3 <- merge(df1, df2,by.x = "a")
df3

##filter
df4 <- df3[df3$c > 11,]
df4
df5 <- df3

##mutate with combined columns
df5$score <- df3$b +df3$c
df5

##group_by and summarise
aggregate(df5$score, by=list(Category=df5$a), FUN=sum)
