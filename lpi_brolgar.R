###lpi + brolgar - visualising tsibbles
library(brolgar)
library(dplyr)
library(zoo)
library(tidyr)
library(tsibble)
library(ggplot2)

lpi <- read.csv("LPI_pops_20160523_edited.csv")

lpi_long <- lpi %>% 
  dplyr::select(ID, Class, Binomial, Latitude, Longitude, starts_with("X19"), starts_with("X20")) %>% 
  pivot_longer( -c(ID, Class, Binomial, Latitude, Longitude), names_to = "year", values_to = "abundance_est") %>% 
  mutate(abundance_est = as.numeric(as.character(abundance_est)), year = as.numeric(gsub("X", "", year)),log_abnd = log10(abundance_est + 1)) %>% 
  group_by(ID) %>% 
  arrange(year) %>% 
  filter(between(row_number(),min(which(!is.na(abundance_est))),max(which(!is.na(abundance_est))))) %>% 
  arrange(ID, year) %>% 
  ungroup()

           
lpi_na_fill <- function(lpi_pop, lpi_long = lpi_long){
  
  lpi_df <- lpi_long %>% 
    filter(ID == lpi_pop)
  
  if(sum(!is.na(lpi_df$log_abnd)) >= 6){
    
    SmoothParm <- round(sum(!is.na(lpi_df$log_abnd))/2) #K value for gam
    lpi_gam <-  mgcv::gam(log_abnd ~ s(year, k = SmoothParm), fx = TRUE, data = lpi_df)
    pred_lpi <- predict(lpi_gam, lpi_df, type = "response", se = TRUE)
    diff_pred <- diff(pred_lpi$fit)
    sum_lambda <- sum(diff_pred)
    mean_lambda <- mean(diff_pred)  
    Rsq <- summary(lpi_gam)$r.sq
    model <- "gam"
    
  } else {
      
    lpi_lm <- lm(log_abnd ~ year, data = lpi_df)
    pred_lpi <- na.approx(lpi_df$log_abnd)
    diff_pred <- diff(pred_lpi)
    sum_lambda <- sum(diff_pred)
    mean_lambda <- mean(diff_pred)
    Rsq <- summary(lpi_lm)$r.sq
    model <- "lm"
    
    }
  
  lpi_out <- data.frame(ID = lpi_pop,
                        class = lpi_df$Class,
                        binomial = lpi_df$Binomial,
                        year = lpi_df$year,
                        true_abnd = lpi_df$abundance_est,
                        pred_abnd = ifelse(model == rep("gam", nrow(lpi_df)), pred_lpi$fit, pred_lpi),
                        diff_out = c(NA, diff_pred), 
                        sum_lambda = sum_lambda,
                        mean_lambda = mean_lambda,
                        Rsq = Rsq,
                        model = model)
  
  print(lpi_out)
  return(lpi_out)
  }

lpi_fill <- lapply(X = unique(lpi_long$ID), lpi_na_fill, lpi_long = lpi_long)

lpi_filled <- do.call(rbind, lpi_fill)


lpi_filt <- lpi_filled %>% 
  group_by(ID) %>% 
  filter(Rsq > 0.5 | n() > 5)


years <- 1950:2015
id_vec <- rep(unique(lpi_filt$ID), each = length(years))

id_df <- data.frame(ID = id_vec, year = years)

lpi_merge <- merge(lpi_filt, id_df, by = c("ID", "year"), all = TRUE)

lpi_tsib <- as_tsibble(x = lpi_merge,
                       key = ID, 
                       index = year, 
                       regular = TRUE)


####make tsibble - plot dataaaaa

ggplot(lpi_tsib, aes(x = year, y = pred_abnd, group = ID))+
  geom_line()+
  facet_sample()



lpi_tsib %>%
  filter(sum(!is.na(pred_abnd))>20) %>% 
  sample_n_keys(size = 10) %>%
  ggplot(aes(x = year,
             y = pred_abnd,
             group = ID, 
             colour = as.factor(ID))) + 
  geom_line()

    

lpi_mono <- lpi_tsib %>% 
  filter(!is.na(pred_abnd)) %>% 
  features(pred_abnd, feat_monotonic) %>% 
  left_join(lpi_tsib, by  = "ID") 
           