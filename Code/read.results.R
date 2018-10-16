### Reading results file
library(tidyverse)

setwd("Code/Results")
df <- list.files() %>%
  map(readRDS)

#Log likelihood plot
plot(seq(-1,1,length=50),unlist(lapply(df, '[[', 2)))

