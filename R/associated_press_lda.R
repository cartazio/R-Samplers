
# Create document-term matrix of Associated Press articles
# The data set is an object of class "DocumentTermMatrix" provided by package tm. 
# It is a document-term matrix which contains the term frequency of 10473 terms in 2246 documents.
# Associated Press data from the First Text Retrieval Conference (TREC-1) 1992.

#install.packages("kableExtra")
library(kableExtra)
library(knitr)
#install.packages("stargazer")
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
  mutate_if(is.numeric , replace_na, replace = 0) %>% 
  select(-document)

dim(press)
press[1:10,1:10]

# 100 most common words
press <- press %>%
  select(-i) %>% 
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
#testing <- testing[1:10,1:50]

## STAN MODEL
#install.packages("rstan")
install.packages("brms")
library(rstan)
library(rstanarm)
library(bayesplot)
library(brms)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
#Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')
#pkgbuild::has_build_tools(debug = TRUE)

ap_dat <- list(
  K = 5,
  V = ncol(training),
  M = nrow(training),
  X = training,
  alpha = rep(1, times = 5),
  beta = rep(1, times = ncol(training))
)
ap_dat

fit <- stan(
  file = "lda.stan",  # Stan program
  data = ap_dat,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 2,              # number of cores (could use one per chain)
  verbose = TRUE
)

vi_model <- rstan::stan_model(
  file = "lda.stan")

fit_vi <- rstan::vb(
  object = vi_model,  # Stan program
  data = ap_dat,      # named list of data
  tol_rel_obj = 0.001 # convergence tolerance on the relative norm of the objective, default = 0.01.
)


#save(fit, file = "lda_out.RDA")
load(file = "lda_out.RDA")

print(fit)
pairs(fit)

la <- rstan::extract(fit, permuted = FALSE) # return a list of arrays 
### return an array of three dimensions: iterations, chains, parameters 

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
dim(m)
d <- as.data.frame(fit)

rstanarm::launch_shinystan(fit)

format(bayesplot::available_ppc()) # omit ppc_
format(bayesplot::available_mcmc()) # omit mcmc_

mcmc_rhat(rhat(fit))

###################
###################
###################

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
