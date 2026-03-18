load("data/byleaf2018.Rdata")
library(tidyverse)
library(survival)
library(contsurvplot)
library(survminer)


first_date = min(byleaf.2018$Survey.Date)

##### ------------- Crown Rust ------------------
#https://www.sthda.com/english/wiki/cox-proportional-hazards-model
byleaf.2018.rust = byleaf.2018 %>%arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1)),
         rust.transition = ifelse(rust_event == 1 & lag(rust_event, default = 0) == 0, 1, 0)) %>% drop_na() %>%
  filter(lag(rust_event, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date)) %>% data.frame()

byleaf.2018.rust$prev.anth = relevel(byleaf.2018.rust$prev.anth, ref = "0") #make sure factor reference levels are zero

#https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
cox_model.rust <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~  PC1 + PC2 + PC3 + PC4 + prev.anth + prev.bp,
  data = byleaf.2018.rust, x = TRUE, model = TRUE, cluster = ID
)

surv_fit <- survfit(cox_model.rust, data = byleaf.2018.rust)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.rust, data = byleaf.2018.rust)
ggsave("results/CoxPH.Rust.PCs.pdf", width = 6.5, height = 4)


##### ------ Anthracnose ---------------
byleaf.2018.anth = byleaf.2018 %>%arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         #prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1)),
         prev.rust = factor(lag(crown_rust, 1))) %>%
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date)) %>% drop_na() %>%data.frame()

cox_model.anth <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.bp,
  data = byleaf.2018.anth, x = TRUE, model = TRUE, cluster = ID
)


surv_fit <- survfit(cox_model.anth, data = byleaf.2018.anth)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.anth, data = byleaf.2018.anth)
ggsave("results/CoxPH.Anth.PCs.pdf", width = 6.5, height = 4)


####------------- brown patch ------------------
byleaf.2018.bp = byleaf.2018 %>%arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         #prev.bp = factor(lag(brown_patch, 1)),
         prev.rust = factor(lag(crown_rust, 1))) %>%
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date)) %>% drop_na() %>%data.frame()

cox_model.bp <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.anth,
  data = byleaf.2018.bp, x = TRUE, model = TRUE, cluster = ID
)


surv_fit <- survfit(cox_model.bp, data = byleaf.2018.bp)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.bp, data = byleaf.2018.bp)

ggsave("results/CoxPH.BP.PCs.pdf", width = 6.5, height = 4)

