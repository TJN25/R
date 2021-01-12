library(tidyverse)
dat <- read.table("~/Downloads/imdb_movies.tsv", sep = "\t", comment.char = "", quote = "")

colnames(dat) <- c("id", "type", "title", "title.2", "t1", "year", "t2", "runtime", "genera")

head(dat)
dat <- dat %>% filter(type == "movie")

recent <- dat %>% filter(as.numeric(as.character(year)) > 1990, grepl(pattern = "dult", genera) == F, grepl(pattern = "ocumentary", genera) == F)


englishtitles <- recent %>% filter(grepl(pattern = "El", title, ignore.case = T) == F,
                                   grepl(pattern = "L'", title, ignore.case = T) == F,
                                   grepl(pattern = "La", title, ignore.case = T) == F,
                                   grepl(pattern = " l'", title, ignore.case = T) == F,
                                   grepl(pattern = "Ich ", title, ignore.case = T) == F,
                                   grepl(pattern = "un  ", title, ignore.case = T) == F,
                                   grepl(pattern = "de ", title, ignore.case = T) == F,
                                   grepl(pattern = "Ku ", title, ignore.case = T) == F,
                                   grepl(pattern = "da ", title, ignore.case = T) == F,
                                   grepl(pattern = "Der ", title, ignore.case = T) == F,
                                   grepl(pattern = "le ", title, ignore.case = T) == F,
                                   grepl(pattern = "yuen", title, ignore.case = T) == F,
                                   grepl(pattern = "zh", title, ignore.case = T) == F,
                                   grepl(pattern = "zw", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¿", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F,
                                   grepl(pattern = "¡", title, ignore.case = T) == F)


englishtitles <- englishtitles %>% arrange(title)


write.table(x = englishtitles, file = "~/Downloads/imdb_movies_english.tab", sep = "\t", quote = F, row.names = F)














netflix <- read.csv("~/Downloads/netflix.csv")
imdb <- read.csv("~/Downloads/imdb.csv")
colnames(netflix)
colnames(imdb)

imdb <- imdb %>% select(Title, IMDb.Rating, Num.Votes, Year, Genres, Description, Directors)

imdb <- imdb %>% unique()
netflix <- netflix %>% unique()

imdb <- imdb %>% mutate(Year = as.numeric(as.character(Year)))


write.csv(x = imdb, file = "~/Downloads/imdb_unique.csv", row.names = F)
