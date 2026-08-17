load("data/weather.pcs.Rdata") #pc.roll
library(tidyverse)
library(survival)
library(contsurvplot)
library(survminer)
#byleaf.2018
dis.pal
### ERROR: Originally this data had the incorrect treatments included; let's not blanket exclude the treatments
#but instead we can look at it at both the leaf and the plant level for any treatments where *fungicide was not sprayed*
forest.plot <- function(cox_model, focal.dis){
  forest.pal = c("#1F78B4", "#33A02C", "#FEC601", "gray40", "gray40", "gray40", "gray40")

  names(forest.pal) = c("Prior ANTH", "Prior BP", "Prior CR", "PC1", "PC2", "PC3", "PC4")
  
  y = summary(cox_model)
  coef = y$coefficients
  cis = data.frame(y$conf.int)
  cis$var = rownames(cis)
  cis = cis %>% mutate(var = case_when(var == "prev.bp1"~"Prior BP",
                                       var == "prev.anth1" ~ "Prior ANTH", 
                                       var == "prev.rust1" ~ "Prior CR",
                                       .default = var))
  cis$var = factor(cis$var, levels = rev(sort(unique(cis$var))))
  
  cis = cis %>% mutate(bio = ifelse(grepl("PC[0-9]", var) == T, "abio", "bio"))
  ci.pal = forest.pal[levels(cis$var)]
  p <-  ggplot(cis, aes(x=log(exp.coef.), y=var, col = var)) + 
    geom_vline(xintercept =0, col = 'gray', lty = 'dashed') +
    geom_point() +
    geom_errorbarh(aes(xmin = log(lower..95), xmax = log(upper..95), width = .2)) + 
    theme_classic() + 
    xlab(paste("log(HR:",focal.dis, ")", sep="")) + ylab("") + 
    scale_color_manual(values=ci.pal) + 
    theme(text = element_text(size = 8), legend.position = "none")
  return(p)
} 

calc_ll_robust <- function(full, red){
  
  ll_full <- logLik(full)
  ll_reduced <- logLik(red)
  
  lrt_stat <- as.numeric(2 * (ll_full - ll_reduced))
  lrt_df <- attr(ll_full, "df") - attr(ll_reduced, "df")
  lrt_p <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
  return(data_frame(stat=lrt_stat, df = lrt_df, p = lrt_p))
}
#we need the 2017 and 2018 data
trtments = read.csv("data/Rita_longitudinal/plot_treatments.csv")
trtments = trtments %>% mutate(
  start_2017=case_when(
  fungicide == "7-months" ~ "05/15/2017",
  fungicide == "9-months" ~ "05/15/2017",
  fungicide == "Year-round" ~ "05/15/2017"),
  end_2017 = case_when(
    fungicide == "7-months" ~ "07/28/2017",
    fungicide == "9-months" ~ "09/18/2017",
    fungicide == "Year-round" ~ "02/18/2020"),
  start_2018 = case_when(
    fungicide == "7-months" ~ "01/15/2018",
    fungicide == "9-months" ~ "01/15/2018"),
  end_2018 = case_when(
    fungicide == "7-months" ~ "07/19/2018",
    fungicide == "9-months" ~ "09/19/2018")
  ) %>% mutate(across(start_2017:end_2018, function(x) as.Date(x, format = "%m/%d/%Y")))

grun2017 = read.csv("data/Rita_longitudinal/longitudinal_disease_survey_2017.csv") %>% left_join(trtments) %>%
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y")) %>%
  mutate(PlantID = outplant, 
         LeafID = paste("y17", plot, outplant, tiller_number, leaf, sep = "_"))

grun2018 = read.csv("data/Rita_longitudinal/longitudinal_disease_survey_2018.csv") %>% left_join(trtments) %>%
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y")) %>%
  mutate(PlantID = tiller, 
         LeafID = paste("y18", plot, tiller, leaf, sep = "_")) #these leaves are numbered correctly, numbers keep increasing, but leaves are not recorded after senesence
keepcols = intersect(colnames(grun2017), colnames(grun2018)) 

grun1718 = rbind(grun2017[,keepcols], grun2018[,keepcols])
#start with 2018

#1---- exclude any observations during a period of active fungicide spray
filt.obs = grun1718 %>% mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y")) %>%
  mutate(keep_obs = case_when(
    fungicide == "Never" ~ "y",
    fungicide == "Year-round" ~"n",
    year(survey_date) == 2017 & (survey_date > end_2017 | survey_date < start_2017) ~ "y",
    year(survey_date) == 2018 & (survey_date > end_2018 | survey_date < start_2018) ~ "y"
  ) ) %>% filter(keep_obs == "y")

filt.obs = filt.obs %>% left_join(pc.roll, by = c("survey_date" = "Survey.Date")) %>%
  mutate(exp_year = as.factor(year(survey_date))) %>% filter(year(survey_date) == 2018) #exclude the outplants, they look funky

##### ------ LEAF LEVEL CROWN RUST -----------
#for leaf level observations, it might also be good to use the cumsum of the prior infections
rust.leaf = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  arrange(LeafID, survey_date) %>%
  group_by(LeafID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1))) %>%
  filter(lag(crown_rust, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(LeafID =as.factor(LeafID), PlantID = as.factor(PlantID))

rust.leaf$prev.anth = relevel(rust.leaf$prev.anth, ref = "0") #make sure factor reference levels are zero


cox_model.rust.leaf <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~  PC1 + PC2 + PC3 + PC4 + prev.anth + prev.bp,
  data = rust.leaf, x = TRUE, model = TRUE, id = LeafID, cluster = PlantID
)

cox_model.rust.leaf

#what if we just grab one random leaf per plant?
set.seed(222)
keep.leaves = filt.obs %>% select(LeafID, PlantID) %>% distinct() %>%
  group_by(PlantID) %>%
  slice_sample(n = 1) %>%
  ungroup() 


cox_model.rust.leaf.red <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~ PC1 + PC2 + PC3 + PC4 + prev.anth + prev.bp,
  data = rust.leaf.red, x = TRUE, model = TRUE, id = LeafID, cluster = plot
)

cox_model.rust.leaf.red
cox_model.rust.leaf #basically the same
cox_model.rust.plant #the direction of observation is different 

##### ------ TILLER-LEVEL CROWN RUST -----------
rust.plant = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  group_by(PlantID, survey_date) %>%
  mutate(anthracnose = ifelse(sum(anthracnose) > 0, 1, 0),
         brown_patch = ifelse(sum(brown_patch) > 0, 1, 0), 
         crown_rust = ifelse(sum(crown_rust) > 0, 1, 0)) %>% 
  ungroup() %>%
  arrange(PlantID, survey_date) %>% 
  group_by(PlantID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1))) %>%
  filter(lag(crown_rust, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(PlantID =as.factor(PlantID), plot = as.factor(plot))

rust.plant$prev.anth = relevel(rust.plant$prev.anth, ref = "0") #make sure factor reference levels are zero


cox_model.rust.plant <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~   (PC1 + PC2 + PC3 + PC4 + prev.anth + prev.bp), 
  data = rust.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
) #different effects in different years 

cox_model.rust.plant

cox_model.rust.plant_noinf <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~   (PC1 + PC2 + PC3 + PC4 ), 
  data = rust.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
) #different effects in different years 

calc_ll_robust(cox_model.rust.plant, cox_model.rust.plant_noinf)
concordance(cox_model.rust.plant)
#the story so far: at the plant level, prior infection by anthracnose facilitates subsequent infection by crown rust
#however, at the leaf level a different story emerges; anthracnose inhibits leaf-level infections, which makes sense because
#it kills leaf tissue that an obligate biotroph would infect. 

##### --------- LEAF LEVEL ANTHRACNOSE ---------------

anth.leaf = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  arrange(LeafID, survey_date) %>%
  group_by(LeafID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.rust = factor(lag(crown_rust, 1)),
         prev.bp = factor(lag(brown_patch, 1))) %>%
  filter(lag(anthracnose, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(LeafID =as.factor(LeafID), PlantID = as.factor(PlantID))

cox_model.anth.leaf <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.bp,
  data = anth.leaf, x = TRUE, model = TRUE, id = LeafID, cluster = PlantID
) #but rust facilitates anth infection 
cox_model.anth.leaf


anth.leaf.red = anth.leaf %>% filter(LeafID %in% keep.leaves$LeafID)

cox_model.anth.leaf.red <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.bp,
  data = anth.leaf.red, x = TRUE, model = TRUE, id = LeafID, cluster = plot
)
cox_model.anth.leaf
cox_model.anth.leaf.red
#comparable to the model with all leaves that is clustered by plant


##### ------ TILLER-LEVEL ANTHRACNOSE -----------
anth.plant = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  group_by(PlantID, survey_date) %>%
  mutate(anthracnose = ifelse(sum(anthracnose) > 0, 1, 0),
         brown_patch = ifelse(sum(brown_patch) > 0, 1, 0), 
         crown_rust = ifelse(sum(crown_rust) > 0, 1, 0)) %>% 
  ungroup() %>%
  arrange(PlantID, survey_date) %>%
  group_by(PlantID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.rust = factor(lag(crown_rust, 1)),
         prev.bp = factor(lag(brown_patch, 1))) %>%
  filter(lag(anthracnose, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(PlantID =as.factor(PlantID), plot = as.factor(plot))

cox_model.anth.plant <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.bp,
  data = anth.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
)

cox_model.anth.plant_noinf <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4,
  data = anth.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
)

cox_model.anth.plant #we find similar patterns at the leaf and plant level for the effects of rust on anthracnose
calc_ll_robust(cox_model.anth.plant, cox_model.anth.plant_noinf)
concordance(cox_model.anth.plant)
##### --------- LEAF LEVEL BROWN PATCH ---------------

bp.leaf = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  arrange(LeafID, survey_date) %>%
  group_by(LeafID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.rust = factor(lag(crown_rust, 1)),
         prev.anth= factor(lag(anthracnose, 1))) %>%
  filter(lag(brown_patch, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(LeafID =as.factor(LeafID), PlantID = as.factor(PlantID))
bp.leaf$prev.anth = relevel(bp.leaf$prev.anth, ref = "0") #make sure factor reference levels are zero

cox_model.bp.leaf <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.anth,
  data = bp.leaf, x = TRUE, model = TRUE, id = LeafID, cluster = PlantID
) #but rust facilitates anth infection 
cox_model.bp.leaf


bp.leaf.red = bp.leaf %>% filter(LeafID %in% keep.leaves$LeafID)

cox_model.bp.leaf.red <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.anth,
  data = bp.leaf.red, x = TRUE, model = TRUE, id = LeafID, cluster = plot
)
cox_model.bp.leaf
cox_model.bp.leaf.red
#comparable to the model with all leaves that is clustered by plant
ph_test <- cox.zph(cox_model.bp.leaf)
print(ph_test)

##### ------ TILLER-LEVEL BROWN PATCH -----------
bp.plant = filt.obs %>%
  select(plot, exp_year, PlantID, survey_date, anthracnose, crown_rust, brown_patch, 
         PC1:PC4, LeafID) %>%
  group_by(PlantID, survey_date) %>%
  mutate(anthracnose = ifelse(sum(anthracnose) > 0, 1, 0),
         brown_patch = ifelse(sum(brown_patch) > 0, 1, 0), 
         crown_rust = ifelse(sum(crown_rust) > 0, 1, 0)) %>% 
  ungroup() %>%
  arrange(PlantID, survey_date) %>%
  group_by(PlantID) %>%
  mutate(first_date = min(survey_date)) %>%
  mutate(prev.time = lag(survey_date, 1),
         prev.rust = factor(lag(crown_rust, 1)),
         prev.anth = factor(lag(anthracnose, 1))) %>%
  filter(lag(brown_patch, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(survey_date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  drop_na() %>%
  mutate(PlantID =as.factor(PlantID), plot = as.factor(plot))
bp.plant$prev.anth = relevel(bp.plant$prev.anth, ref = "0") #make sure factor reference levels are zero

cox_model.bp.plant <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.anth,
  data = bp.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
)

cox_model.bp.plant #at the plant level it's not significant for rust or for anth. level is wrong for anth (nearly significant for rust)
concordance(cox_model.bp.plant)

cox_model.bp.plant_noinf <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4,
  data = bp.plant, x = TRUE, model = TRUE, id = PlantID, cluster = plot
)

calc_ll_robust(cox_model.bp.plant, cox_model.bp.plant_noinf)


########------ CREATE NICE FOREST PLOTS --------
forest.plot(cox_model.bp.plant, "brown patch")
ggsave(file="results/SurvBP.pdf", width = 6.5/3, height = 2, units = "in", dpi = 500)
forest.plot(cox_model.anth.plant, "anthracnose")
ggsave(file="results/SurvANTH.pdf", width = 6.5/3, height = 2, units = "in", dpi = 500)
forest.plot(cox_model.rust.plant, "crown rust")
ggsave(file="results/SurvRUST.pdf", width = 6.5/3, height = 2, units = "in", dpi = 500)


#first_date = min(long.withweather$Survey.Date)
#collapse to plant level
long.withweather = long.withweather %>%
  group_by(tiller, plot, Survey.Date, PC1, PC2, PC3, PC4) %>%
  summarise(anthracnose = ifelse(sum(anthracnose) > 0, 1, 0),
         crown_rust = ifelse(sum(crown_rust) > 0, 1, 0), 
         brown_patch = ifelse(sum(brown_patch) > 0, 1, 0))

long.withweather$ID = paste(long.withweather$plot, long.withweather$tiller)
##### ------------- Crown Rust ------------------
#https://www.sthda.com/english/wiki/cox-proportional-hazards-model
byleaf.2018.rust = long.withweather %>% arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(first_date = min(Survey.Date)) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1)),
         rust.transition = ifelse(crown_rust == 1 & lag(crown_rust, default = 0) == 0, 1, 0)) %>% drop_na() %>%
  filter(lag(crown_rust, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% data.frame() %>%
  mutate(ID =as.factor(ID))

byleaf.2018.rust$prev.anth = relevel(byleaf.2018.rust$prev.anth, ref = "0") #make sure factor reference levels are zero

#https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
cox_model.rust <- coxph(
  Surv(start.days, end.days, event = crown_rust) ~  PC1 + PC2 + PC3 + PC4 + prev.anth + prev.bp,
  data = byleaf.2018.rust, x = TRUE, model = TRUE, id = ID, cluster = plot
)

surv_fit <- survfit(cox_model.rust, data = byleaf.2018.rust)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.rust, data = byleaf.2018.rust)
ggsave("results/CoxPH.Rust.PCs.pdf", width = 6.5, height = 4)


##### ------ Anthracnose ---------------
byleaf.2018.anth = long.withweather %>%arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(first_date = min(Survey.Date)) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         #prev.anth = factor(lag(anthracnose, 1)),
         prev.bp = factor(lag(brown_patch, 1)),
         prev.rust = factor(lag(crown_rust, 1))) %>%
  filter(lag(anthracnose, default = 0) == 0) %>% #only keep the transition event, but not the following events
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date)+1) %>% drop_na() %>%data.frame()

cox_model.anth <- coxph(
  Surv(start.days, end.days, event = anthracnose) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.bp,
  data = byleaf.2018.anth, x = TRUE, model = TRUE, id = ID, cluster = plot
)


surv_fit <- survfit(cox_model.anth, data = byleaf.2018.anth)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.anth, data = byleaf.2018.anth)
ggsave("results/CoxPH.Anth.PCs.pdf", width = 6.5, height = 4)


####------------- brown patch ------------------
byleaf.2018.bp =  long.withweather %>%arrange(ID, Survey.Date) %>%
  group_by(ID) %>%
  mutate(first_date = min(Survey.Date)) %>%
  mutate(prev.time = lag(Survey.Date, 1),
         prev.anth = factor(lag(anthracnose, 1)),
         #prev.bp = factor(lag(brown_patch, 1)),
         prev.rust = factor(lag(crown_rust, 1))) %>%
  filter(lag(brown_patch, default =0) == 0 ) %>%
  mutate(end.days = as.numeric(Survey.Date - first_date), 
         start.days = as.numeric(prev.time - first_date) + 1) %>% drop_na() %>%data.frame()

byleaf.2018.bp$prev.rust = relevel(byleaf.2018.bp$prev.rust, ref = "0") #make sure factor reference levels are zero
byleaf.2018.bp$prev.anth = relevel(byleaf.2018.bp$prev.anth, ref = "0") #make sure factor reference levels are zero


cox_model.bp <- coxph(
  Surv(start.days, end.days, event = brown_patch) ~  PC1 + PC2 + PC3 + PC4 + prev.rust + prev.anth,
  data = byleaf.2018.bp, x = TRUE, model = TRUE, cluster = ID
)


surv_fit <- survfit(cox_model.bp, data = byleaf.2018.bp)
ggsurvplot(surv_fit, conf.int = TRUE)
ggforest(cox_model.bp, data = byleaf.2018.bp)

ggsave("results/CoxPH.BP.PCs.pdf", width = 6.5, height = 4)

