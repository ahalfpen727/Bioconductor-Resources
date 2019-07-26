## ----loadLib, message=FALSE----------------------------------------------
library(tidyverse)

## ----ex1-----------------------------------------------------------------
x <- 3
log(x)

## ----ex1sol--------------------------------------------------------------
x <- 3
x %>% log()

## ----diamondDF-----------------------------------------------------------
data(diamonds)
set.seed(25091309)
sample1000 <- sample(1:nrow(diamonds), 1000, replace = FALSE)
diamonds <- diamonds[sample1000, ]

## ----select--------------------------------------------------------------
small_diam <- diamonds %>% select(cut, color, price)
small_diam

## ----arrange-------------------------------------------------------------
ordered_diams <- diamonds %>% arrange(desc(color))
ordered_diams

## ----filter--------------------------------------------------------------
vg_diam <- diamonds %>% filter(cut == "Ideal")
summary(vg_diam)
vg_diam

## ----filterb-------------------------------------------------------------
vg_diam <- diamonds %>% filter(cut == "Ideal")
summary(vg_diam)
vg_diam

## ----exTopExp------------------------------------------------------------
top_exp <- diamonds %>% filter(price >= quantile(price, probs = 0.9))
summary(top_exp$price)

## ----mutate--------------------------------------------------------------
large_diam <- diamonds %>% mutate(ratio = price / carat,
                                  nothing = NA,
                                  weird = paste(cut, color, sep = "-")) %>%
  select(carat, price, color, cut, ratio, nothing, weird)
large_diam

## ----mutateb-------------------------------------------------------------
large_diam <- diamonds %>% mutate(ratio = price / carat,
                                  nothing = NA,
                                  weird = paste(cut, color, sep = "-")) %>%
  select(carat, price, color, cut, ratio, nothing, weird)
large_diam

## ----exNewVars-----------------------------------------------------------
ld <- diamonds %>% filter(color == "D") %>%
  mutate(log10 = log10(price),
         combo = paste(cut, clarity, sep = "-")) %>%
  select(price, cut, clarity, color, log10, combo)
ld

## ----summarize-----------------------------------------------------------
new_diams <- diamonds %>% summarise(av_depth = mean(depth),
                                    sd_depth = sd(depth))
new_diams
new_diams <- diamonds %>% 
  group_by(color) %>%
  summarise(av_price = mean(price),
            sd_price = sd(price))
new_diams

## ----summEx--------------------------------------------------------------
new_diamb <- diamonds %>% filter(color %in% c("D", "E")) %>%
  group_by(color, cut) %>%
  summarise(av_price = mean(price),
            sd_price = sd(price),
            count = length(price))
new_diamb

## ----ggplotSumm----------------------------------------------------------
p <- ggplot(new_diamb, aes(x = cut, y = av_price, colour = color, group = color)) + geom_point() + 
  geom_line() + geom_errorbar(aes(ymin = av_price - sd_price / sqrt(count),
                                  ymax = av_price + sd_price / sqrt(count)))
p

## ----fakeData------------------------------------------------------------
grades <- tibble(
  Name = c("Tommy", "Mary", "Gary", "Cathy"),
  Sexage = c("m.15", "f.15", "m.16", "f.14"),
  Math = c(10, 15, 16, 14),
  Philo = c(11, 13, 10, 12),
  English = c(12, 13, 17, 10)
)
grades

## ----separate------------------------------------------------------------
grades <- grades %>% 
  separate(Sexage, into = c("Sex", "Age")) # default separator is any nonalphanumeric character

## ----gather--------------------------------------------------------------
modif_grades <- grades %>%
  gather(Math, Philo, English, key = Topic, value = Grade)
modif_grades

## ----ggplotGather--------------------------------------------------------
p <- ggplot(modif_grades, aes(x = Topic, y = Grade)) + geom_boxplot() + 
  theme_bw()
p

## ----gratherEx, echo=FALSE-----------------------------------------------
modif2 <- grades %>%
  gather(Math, Philo, English, key = Topic, value = Grade) %>%
  filter(Sex == "f") %>%
  group_by(Topic) %>%
  mutate(minval = min(Grade),
         maxval = max(Grade))
modif2

## ----gratherEx2, echo=FALSE----------------------------------------------
p <- ggplot(modif2, aes(y = (minval + maxval) / 2, x = Topic)) + geom_point() +
  geom_errorbar(aes(ymin = minval, ymax = maxval)) + theme_bw()
p

