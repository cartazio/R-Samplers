## Example
library(tidyverse)
library(janitor)
library(gridExtra)
library(stargazer)
source("./R/one_dim_ex.R")

## Read in data
undata <- read_csv("./data/country_profiles.csv") %>%
  clean_names() %>%
  mutate_at(vars(3:50), as.numeric) %>%
  dplyr::select(gdp = gdp_gross_domestic_product_million_current_us) %>% 
  na_if(-99) %>% 
  drop_na(.) %>% 
  mutate(gdp = scale(gdp)) %>% 
  as.matrix(.) %>% as.vector(.)

three_levels <- Gauss1dim(k = 3, x = undata)

df <- three_levels$TheProposalRecord %>% as_tibble(.name_repair = "unique")

all_mu_3 <- df %>% dplyr::select(grep("mu", colnames(.))) %>%
  mutate(mu_num = row.names(.)) %>% 
  gather(key = "t", value = "mu", -mu_num) %>% 
  spread(mu_num, mu) %>% 
  mutate(t = 1:nrow(.)) %>% 
  rename(mu_1 = 2, mu_2 = 3, mu_3 = 4)

all_pi_3 <- df %>% dplyr::select(grep("pi", colnames(.))) %>% 
  mutate(pi_num = row.names(.)) %>% 
  gather(key = "t", value = "pi", -pi_num) %>% 
  spread(pi_num, pi) %>% 
  mutate(t = 1:nrow(.)) %>% 
  rename(pi_1 = 2, pi_2 = 3, pi_3 = 4)

all_z_3 <- df %>% dplyr::select(grep("z", colnames(.))) %>% 
  mutate(z_num = row.names(.)) %>% 
  gather(key = "t", value = "z", -z_num) %>% 
  spread(z_num, z) %>% 
  mutate(t = 1:nrow(.)) %>% 
  rename(z_1 = 2, z_2 = 3, z_3 = 4)

all_post <- left_join(all_mu_3, all_pi_3) %>% 
  left_join(., all_z_3)

mu_1 <- all_mu_3 %>% 
  ggplot(aes(x = t, y = mu_1)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "MU 1")

mu_2 <- all_mu_3 %>% 
  ggplot(aes(x = t, y = mu_2)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "MU 2")

mu_3 <- all_mu_3 %>% 
  ggplot(aes(x = t, y = mu_3)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "MU 3")

pi_1 <- all_pi_3 %>% 
  ggplot(aes(x = t, y = pi_1)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "PI 1")

pi_2 <- all_pi_3 %>% 
  ggplot(aes(x = t, y = pi_2)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "PI 2")

pi_3 <- all_pi_3 %>% 
  ggplot(aes(x = t, y = pi_3)) + geom_line() +
  theme_bw() +
  labs(x = "Marcov Chain states", y = "PI 3")

pdf("hw1_plots.pdf")
grid.arrange(mu_1, mu_2, mu_3, pi_1, pi_2, pi_3, nrow = 2)
dev.off()

# deviance
# log p(x1:n, z1:n, µ1:K, pi1:K)
# prob( X, ClusterCenters, Z, PI )
# Prob (A , B , C , D) = P (A | B ,C ,D) * P(B | C ,D) * P (C | D) * P (D)

# post_pi <- all_post %>% 
#   dplyr::select(pi_1, pi_2, pi_3) %>% as.matrix()
#   
# sum(post_pi[1,])
# 
# MCMCpack::ddirichlet(post_pi, c(1/3, 1/3, 1/3))

uncountry <- read_csv("./data/country_profiles.csv") %>%
  clean_names() %>%
  dplyr::select(country, gdp_gross_domestic_product_million_current_us) %>% 
  na_if(-99) %>% 
  drop_na(.) %>% 
  dplyr::select(country)

labeled <- cbind(uncountry, t(three_levels$FinalState$z) %>% as_tibble())

labeled %>% filter(V1 == 1) %>% select(country) %>% stargazer(summary = F)
labeled %>% filter(V2 == 1) %>% select(country) %>% stargazer(summary = F)
labeled %>% filter(V3 == 1) %>% select(country) %>% stargazer(summary = F)

knitr::kable()
