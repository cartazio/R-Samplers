
# Create document-term matrix of Associated Press articles
# The data set is an object of class "DocumentTermMatrix" provided by package tm. 
# It is a document-term matrix which contains the term frequency of 10473 terms in 2246 documents.
# Associated Press data from the First Text Retrieval Conference (TREC-1) 1992.

library(tidyverse)
#install.packages("tm")
library(tm)
#install.packages("topicmodels")
library(topicmodels)
#install.packages("tidytext")
library(tidytext)
data("AssociatedPress", package = "topicmodels")
AssociatedPress

press <- tidy(AssociatedPress) %>% spread(key = term, value = count) %>% 
  mutate(document = 1:2246) %>% 
  select(document, everything()) %>% 
  mutate_if(is.numeric , replace_na, replace = 0)

dim(press)
press[1:10,1:10]

# Want to remove words that are too frequent or too infrequent
# Top 0.1% words: first government i  last million   new people percent police president state states   two  year years
# Stop word: I, remove
16/10473
press %>%
  select_if(vars(sum(.) > 1000))

#only 1000 most common words
press <- press %>%
  select(-i) %>% 
  select_if(vars(sum(.) > 95)) #%>% 
  #as.matrix()

index <- sample(nrow(press), size = 1500, replace = FALSE)

training <- press[index,]
testing <- press[-index,]

K = 10
a <- rep(1, times = K)
g <- rep(1, times = 1003)
train <- lda_vi(n = training[,2:1004], K, a, g, max_iter = 100)

