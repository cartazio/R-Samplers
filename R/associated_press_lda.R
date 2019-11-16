
# Create document-term matrix of Associated Press articles
# The data set is an object of class "DocumentTermMatrix" provided by package tm. 
# It is a document-term matrix which contains the term frequency of 10473 terms in 2246 documents.
# Associated Press data from the First Text Retrieval Conference (TREC-1) 1992.

#install.packages("kableExtra")
library(kableExtra)
library(knitr)
library(stargazer)
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

#only 100 most common words
press <- press %>%
  select(-i, -document) %>% 
  select_if(vars(sum(.) > 425)) 
dim(press)

columns <- names(press)

# As matrix
press <- press %>% 
  as.matrix()

# Test / Training sets
index <- sample(nrow(press), size = 1000, replace = FALSE)
training <- press[index,]
testing <- press[-index,]
testing <- testing[1:1000,]

# Build model
K = 10
a <- rep(1/K, times = K)
g <- rep(1, times = 104)
dim(training)

train <- lda_vi(n = training, K, a, g, max_iter = 100)

# Normalize
train$pi_ik <- t(apply(train$pi_ik, 1, function(x)(x/(sum(x))))) %>% as_tibble()
train$b_vp <- apply(train$b_vp, 2, function(x)(x/(sum(x)))) %>% as_tibble()

## Plot ELBO

plot_elbo <- cbind(elbo = train$elbo, iter = train$iter) %>% as_tibble()

plot_elbo %>% 
  ggplot(aes(x = iter, y = elbo)) + geom_line() +
  theme_bw() +
  labs(x = "Evidence Lower Bound", y = "Iterations")

# Visualize

# Topic mixture of word loadings
word_mix <- train$b_vp %>% as_tibble() %>% 
  mutate_all(vars(./sum(.))) %>% as.matrix()

rownames(word_mix) <- columns
heatmap(word_mix)

# Words with high loadings to topics

topic_members <- word_mix %>% as_tibble() %>% 
  cbind(., max = colnames(word_mix)[max.col(word_mix,ties.method="first")]) %>% 
  as_tibble() %>% 
  cbind(., word = columns) %>% as_tibble()
  
topics <- topic_members %>% 
  select(max, word) %>% 
  arrange(max)

View(topics)

table <- topic_members %>% 
  spread(key = max, value = word)

apply(table, 2, sort) %>% knitr::kable("latex")

# Check model on held out data

norm((training), type = "F")
norm((testing), type = "F")
norm((as.matrix(train$pi_ik) %*% t(as.matrix(train$b_vp))), type = "F")

norm((training - as.matrix(train$pi_ik) %*% t(as.matrix(train$b_vp))), type = "F")
norm((testing - as.matrix(train$pi_ik) %*% t(as.matrix(train$b_vp))), type = "F")
