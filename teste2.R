
library(tibble)
# data_frame
library(dplyr)


set.seed(12345)


##########################################
##########    GERANDO DADOS     ##########
##########################################


signal_to_noise <- .2
n               <- 750
prop_ones       <- .03

beta_1 <- sqrt(signal_to_noise)
#
dat <- data_frame(x1 = rnorm(n),
                  intercept = qnorm(prop_ones, sd = sqrt(beta_1 ^ 2 + 1)),
                  y_signal = beta_1 * x1,
                  y_star = intercept + y_signal + rnorm(n),
                  y = y_star >= 0,
                  true_prob = pnorm(y_signal + intercept))

dat %>%
  count(y) %>%
  mutate(prop = round(n/sum(n),2))


####################################
##########    FIT GLM     ##########
####################################


train <- dat %>%
  sample_frac(.66)

test <- dat %>%
  anti_join(train, by = c("x1", "y_star", "y", "true_prob"))

X_train <- train %>%
  select(x1)

X_test <- test %>%
  select(x1)

fit <- glm(as.factor(y) ~ x1, data = train, family = "binomial")


coefs <- fit %>% 
  coef() %>% 
  as.numeric()

V <- vcov(fit)


# MATRIZ DE INFORMACAO DE FISHER
V
summary(fit)$cov.scaled

########################################
##########      PREDITOS      ##########
########################################


logisticPred <- function(X, coef) {
  
  X %>%
    na.omit() %>%
    mutate(int = 1) %>%
    select(int, everything()) %>%
    as.matrix(.) %*% coef %>%
    as.vector() %>%
    (function(x) 1 / (1 + exp(-x)))
  
}

pred <- X_train %>%
  logisticPred(coefs)
hist(pred)

test_pred <- X_test %>%
  logisticPred(coefs)
hist(test_pred)

###########################################
##########      KING & ZENG      ##########
###########################################


# CONVERTENDO EM MATRIZ POR PRECISAR DE OPERACOES ALGEBRICAS

X_matrix <- X_train %>%
  mutate(bias = 1) %>%
  select(bias, everything()) %>%
  as.matrix()

X_test_matrix <- X_test %>%
  mutate(bias = 1) %>%
  select(bias, everything()) %>%
  as.matrix()


# CALCULANDO VIES

W <- diag(pred * (1 - pred))

Q <- X_matrix %*% solve(t(X_matrix) %*% W %*% X_matrix) %*% t(X_matrix)

e <- 0.5 * diag(Q) * (2 * pred - 1)

bias <- (solve(t(X_matrix) %*% W %*% X_matrix) %*% t(X_matrix) %*% W %*% e) %>%
  as.numeric()

# SUBTRAINDO O VIES

unbiased_coefs <- coefs - bias

updated_var <- (nrow(X_train) / (nrow(X_train) + ncol(X_train) + 1)) ^ 2 * V


##########################################################
##########    COMPARANDO COM O PACOTE ZELIG     ##########
##########################################################
install.packages("Zelig")
library(Zelig)

re_model <- zelig(y ~ x1, data = train, model = "relogit", cite = FALSE) 

unbiased_coefs
coef(re_model)
# estimativas batem

sqrt(diag(updated_var))
sqrt(diag(as.matrix(vcov(re_model)[[1]])))
# erros nao batem

zelig(y ~ x1, data = train, model = "relogit", cite = FALSE) %>% 
  summary

glm(y ~ x1, data = train, family = "binomial") %>%
  summary



#################   TESTE NA SW   #################





