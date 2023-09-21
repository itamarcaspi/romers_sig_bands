library(tidyverse)
library(broom)
library(haven)
library(dynlm)
library(sandwich)
library(lmtest)



h <- 17
p <- 0.05
lags <- 6
nwlag <- lags


df_raw <- read_dta("data/aggregatedata_final.dta")


df <- df_raw %>%
  mutate(
    date = observation_date,
    dlrgdp = 100 * dlrgdp,
    dlcpi  = 100* dlcpi,
    lcpi   = 100 * lcpi,
    lrgdp  = 100* lrgdp
  ) %>% 
  filter(
    date >= ymd("1985-01-01"),
    date <= ymd("2007-12-31")
  ) %>% 
  select(date, lcpi, rr_shock, dlcpi, dlrgdp, dstir) %>% 
  rename(
    y = lcpi,
    z = rr_shock
  )


lp_mat <- matrix(data = NA, ncol = 7, nrow = h+1)
lp_mat[,1] <- 0:h

for (i in 0:h) {
  
  lm_fit <- dynlm(L(ts(df$y), -i) - L(ts(df$y),1) ~ ts(df$z) + L(ts(df$dlrgdp), 1:lags) + L(ts(df$dlcpi), 1:lags) + L(ts(df$dstir), 1:lags))
  
  lm_coef <- coeftest(lm_fit, vcov. = NeweyWest(lm_fit, lag = i, prewhite = FALSE))
  lm_summary <- tidy(lm_coef)
  estimate   <- lm_summary %>% slice(2) %>% select(estimate) 
  se         <- lm_summary %>% slice(2) %>% select(std.error) 
  lp_mat[i, 2] <- estimate[[1]]
  lp_mat[i, 3] <- se[[1]]
  lp_mat[i, 4] <- estimate[[1]] + 1 * se[[1]]
  lp_mat[i, 5] <- estimate[[1]] - 1 * se[[1]]
  lp_mat[i, 6] <- estimate[[1]] + 2 * se[[1]]
  lp_mat[i, 7] <- estimate[[1]] - 2 * se[[1]]
}


model_y <- dynlm(ts(df$y) - L(ts(df$y),1) ~ L(ts(df$dlrgdp), 1:lags) + L(ts(df$dlcpi), 1:lags) + L(ts(df$dstir), 1:lags))
r_y <- model_y$residuals

model_z <- dynlm(ts(df$z) ~ L(ts(df$dlrgdp), 1:lags) + L(ts(df$dlcpi), 1:lags) + L(ts(df$dstir), 1:lags))
r_z <- model_z$residuals

w <- r_z^2
mw <- mean(w)
eta <- r_y*r_z

model_eta <- lm(eta ~ 1)
nw_se_eta <- sqrt(diag(vcovHAC(model_eta, lag = nwlag)))
seta <- nw_se_eta["(Intercept)"]

sbeta <- seta/mw

bju <- qnorm(1 - (p / (2 * (h + 1)))) * sbeta
bjd <- -bju



lp_tbl <- as_tibble(lp_mat)
colnames(lp_tbl) <- c("horizon", "mid", "se", "up1", "down1", "up2", "down2", "R2")

lp_tbl <- lp_tbl %>% 
  mutate(
    bju = bju,
    bjd = bjd
  ) %>% 
  drop_na()


lp_tbl %>% 
  ggplot(aes(x = horizon)) +
  geom_line(aes(y = mid), color = 'darkgreen', size = 1.5) +
  geom_line(aes(y = bju), color = "red", linetype = 2, size = 1.5) +
  geom_line(aes(y = bjd), color = "red", linetype = 2, size = 1.5) +
  geom_ribbon(aes(ymin = down1, ymax = up1),
              fill = "darkgreen", linetype = 0, alpha = 0.3) +
  geom_ribbon(aes(ymin = down2, ymax = up2),
              fill = "darkgreen", linetype = 0, alpha = 0.1) +
  scale_x_continuous(expand = c(0,0), breaks = seq(0, h, 5)) +
  scale_y_continuous(expand = c(0,0), limits = c(-6, 2)) +
  xlab("Quarter") +
  ylab("Percent") +
  geom_hline(aes(yintercept = 0)) +
  theme_light(12)
