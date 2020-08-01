# Title     : "COVID-19 mortality"
# Objective : TODO
# Created by: carloscaro
# Created on: 2020-08-01

library(tidyverse)

# Data preparation
d <- read_csv("covid19_mortality.csv") %>%
  mutate(non_deaths = cases - deaths)
d

counts <- matrix(c(23, 27, 103, 94), nrow = 2, ncol = 2)
fisher.test(counts)
# The fisher test is to check if there are non-random associations between 2 variables

get_M <- function(di) {
  deaths_f <- filter(di, gender == "f")$deaths
  non_deaths_f <- filter(di, gender == "f")$non_deaths
  deaths_m <- filter(di, gender == "m")$deaths
  non_deaths_m <- filter(di, gender == "m")$non_deaths
  matrix(c(non_deaths_f, non_deaths_m,
           deaths_f, deaths_m),
         byrow = TRUE,
         nrow = 2, ncol = 2,
         dimnames = list(dead = c('n', 'y'),
                         gender = c('f', 'm')))
}

get_M(filter(d, student == "Ling")) %>%
  prop.table(2)

d %>%
  filter(student == "Antoine") %>%
  get_M() %>%
  fisher.test() %>%
  broom::tidy()

broom::tidy(fisher.test(get_M(d)))
