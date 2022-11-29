source("functions.R")

# sample the data ----
set.seed(123)
n <- 1000
p <- 6
X <- matrix(runif(n*p, 0, 1), n, p)
beta <- c(0, 7.5, 4.5, 1.5, -1.5, -4.5, -7.5)
pY <- sigma(cbind(1,X) %*% beta)
Y <- rbinom(n, size=1,  prob=pY)

# set the label frequency and sample S ----
c <- 0.8
S <- Y*rbinom(n, size=1,  prob=c)

# set the hyperparameters and initial values ----
beta_b <- rep(0, p +1)
beta_B <- 10*diag(p + 1)
beta_start <- rep(0, p + 1)
c_alpha <- 1
c_beta <- 1
c_start <- 0.5

# simulation ----

results <- PUPG(S=S, X=cbind(1, X), 
                N=500, B=0,
                beta_b=beta_b, beta_B=beta_B, beta_start=beta_start,
                c_alpha=c_alpha, c_beta=c_beta, c_start=c_start)

# plots ----

library(dplyr)
library(tidyverse)
library(ggplot2)

### C ###

df_c <- data.frame(c=results$c_values, N=1:500)

plot_c <- ggplot(data=df_c, aes(x=N, y=c)) +
  geom_line() +
  theme_bw() +
  ylim(0, 1) +
  geom_hline(yintercept=c) +
  annotate("rect", xmin = 251, xmax = 500, ymin = 0, ymax = 1, 
           alpha = .25) +
  ggtitle("Chain of sampled values of C")
plot_c


df_c_hist <- data.frame(c=results$c_values[251:500])

plot_c_hist <- ggplot(df_c_hist, aes(x=c, y=..density..)) +
  geom_histogram(bins=nclass.FD(df_c_hist$c), fill="gray60" ,color="gray60") +
  geom_vline(xintercept = mean(results$c_values[251:500])) +
  theme_bw() + 
  ggtitle("Histogram of sampled values of C") +
  xlab(label="c")
plot_c_hist


### beta ###

df_beta <- data.frame(results$beta_values[,-1], N=1:500) %>% pivot_longer(1:p)

plot_beta <- ggplot(df_beta) +
  scale_colour_grey(labels=lapply(c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]"), function(x) parse(text = x)),
                    name="Coef. \n name") +
  theme_bw() +
  ggtitle(expression("The chain of sampled values of "~beta)) +
  geom_hline(data=data.frame(beta=beta[-1], name1=paste0("coef", 1:6)), aes(yintercept=beta), linetype="dashed") +
  geom_label(data=data.frame(beta0=beta[-1], 
                             name0=c("bar(beta)[1]", "bar(beta)[2]", "bar(beta)[3]", "bar(beta)[4]", "bar(beta)[5]", "bar(beta)[6]")),
             aes(x=520, y=beta0, label=name0), parse=T) +
  geom_line(aes(x=N, y=value, color=name)) +
  ylab(expression(beta)) +
  xlab("N")
plot_beta


### P(Y=1) ###

df_Y <- data.frame(Y=factor(Y), S=factor(S), pY=pY, pYpost=results$pY1_values[500,])

plot_Y <- ggplot(data=df_Y, aes(x=pY, y=pYpost, color=Y, shape=S)) +
  geom_point() +
  scale_color_manual(values=c("grey70", "black")) +
  scale_shape_manual(values=c(16, 3)) +
  theme_bw() +
  xlab(expression(sigma[0])) +
  ylab(expression(p[Y[i]==1])) +
  ggtitle("Comparison of probabilities of positive observation")
plot_Y
