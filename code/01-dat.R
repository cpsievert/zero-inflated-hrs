#get the top 5 HR leaders from 2012
library("rvest")
html("http://espn.go.com/mlb/history/leaders/_/breakdown/season/year/2012/sort/homeRuns") %>%
  html_table() -> top
top5 <- top[[1]][3:8, 2]

library("dplyr")
# If you don't have a PITCHf/x database, go get yourself one ;)
# https://baseballwithr.wordpress.com/2014/03/24/422/
db <- src_sqlite("~/pitchfx/pitchRx.sqlite3")
atbat12 <- tbl(db, "atbat") %>% filter(batter_name %in% top5) %>%
  filter(date >= '2012_01_01' & date <= '2013_01_01') %>%
  select(inning_side, event, batter_name, date, url) 
game <- tbl(db, "game") %>% filter(game_type == "R") %>% 
  select(game_type, url) 
atbats <- inner_join(atbat12, game, by = "url") %>% collect

#aggregate to the game level
dat <- atbats %>%
  group_by(batter_name, url) %>%
  summarise(home = sum(inning_side == "bottom"),
            y = sum(event == "Home Run"), n = n()) %>%
  mutate(home = home/n) %>%
  collect
#write.csv(dat, file="data/dat.csv", row.names=FALSE)

# The mcmc() function assumes these are in global environment
N <- length(dat$y)
n.players <- length(unique(dat$batter_name))
ni <- as.numeric(table(dat$batter_name))
zero <- dat$y == 0                        #indicates whether there was a 'success'
cum.ni <- cumsum(ni)                      #useful for computing likelihood at the player level
