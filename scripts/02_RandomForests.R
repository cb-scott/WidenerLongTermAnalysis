#################################################################
####### BUILD AND RUN RF MODELS FOR EACH DISEASE ################
#################################################################
#Notes:
#. - Running these models takes a while (about 30 minutes with 11 cores) I've provided the final model outputs as Rdata files and you can load from those points in the script
#. - I've set up the parallelization on MacOS with M1 chip. I know that depending on your version of R and operating system (Windows, notiriously has conflicts with parallelization in R)
#    this parallelization scheme may not work for you. If this breaks the script you can run it in serial, but expect the run time to increase dramatically
#. - In the temporal splits, you will get some splits that don't have any positive observations. These are generally during winter, where disease is not present in the field
#.  You don't need to worry about the warning, unless it is most of your splits or those you know to have ample data.

library(readxl)
library(tidyverse)
library(themis)
library(randomForest)
library(caret)
library(tidymodels)
library(iml)
library(vip)
source("scripts/helper_functions.R")

load("data/Env+DiseaseData.Rdata")
############---------------- MACHINE LEARNING: CROWN RUST  -------------- ########################
set.seed(222)

metric <- metric_set(rmse)
var_only <- rf.data %>%
  mutate(crown_rust = as.numeric(as.character(crown_rust))) %>% select(!c(npoints)) %>%
  mutate(anthracnose = as.numeric(as.character(anthracnose)), brown_patch = as.numeric(as.character(brown_patch))) %>% arrange(survey_date) %>%
  select(!pc1:pc4)

#set up temporally explicit train and test sets to avoid temporal leakage.
split_data = initial_time_split(var_only, prop = .8)

library(themis)

#Set up temporally explicit cross-validation windows to avoid temporal leakage.
temporal_split <- sliding_period(
  data = training(split_data),
  index = "survey_date",  # The temporal index column, must be a date format
  period = "month",
  lookback = 12,            # Training window = previous 12 months
  assess_stop = 1)  #predict the next one month

#basic recipe
basic_rec <- recipe(crown_rust ~ ., data = training(split_data)) %>%
  update_role(survey_date, new_role = "SamplingDate") %>%
  update_role(plot_id, new_role = "plot_id") %>%
  step_zv(all_predictors())

recipe_mod <- basic_rec |>
  step_novel(plot_id) %>%
  step_dummy(all_nominal_predictors()) |>
  step_mutate(crown_rust = factor(crown_rust)
  ) |>
  step_upsample(crown_rust)

#set up a model
rf_spec <- rand_forest(
  mtry = tune(), #tune the mtry hyperparameter here
  trees = 3000, #up the number of trees to better stabilize variable importances
  mode = "classification") %>%
  set_engine("ranger", importance = "permutation")

#set up a workflow
rf_workflow <- workflow() %>%
  add_recipe(recipe_mod) %>%
  add_model(rf_spec) 

#Try mtry values from 4 to 32, this determines the number of variables each RF consideres
grid <- grid_regular(mtry(range = c(4, 32)), levels = 8)

#Run it with temporal resampling.... doing some hefty lifting here. 
#will get some errors for some splits that lack data/temporal resolution, I think this is okay
library(future)
plan(multisession, workers = 11) 
rf_results <- tune_grid(
  rf_workflow,
  resamples = temporal_split,
  grid = grid
) 

plan(sequential)

#  COLLECT THE BEST RESULTS
best_params <- select_best(rf_results, metric = "roc_auc")
best_params #best mtry is 4, this sounds more right!


#FINALIZE THE MODEL WITH THE BEST PARM SPACE
final_rf <- rf_workflow %>%
  finalize_workflow(select_best(rf_results, metric = "roc_auc")) %>%
  last_fit(split_data)

save(final_rf, file="results/crown_rust_rf_mtry4_upsample.Rdata")

load("results/crown_rust_rf_mtry4_upsample.Rdata")


#Extract best fit and metrics from final model
all_metrics = collect_metrics(final_rf)
predictions_tbl <- collect_predictions(final_rf)

test_metrics <- predictions_tbl %>%
  metrics(truth = crown_rust, estimate = .pred_class)
test_metrics #kappa greater than 0 = better than random chance. Good values start at about .21
test_metrics <- predictions_tbl %>%
  f_meas(truth = crown_rust, estimate = .pred_class) #High values = better!
test_metrics
predictions_tbl %>% 
  yardstick::precision(crown_rust, .pred_class) 
predictions_tbl %>% 
  yardstick::recall(crown_rust, .pred_class) 

predictions_tbl %>% accuracy(crown_rust, .pred_class)


library(caret)
out <- confusionMatrix(predictions_tbl$.pred_class, predictions_tbl$crown_rust)

plt <- as.data.frame(out$table)
plt$Prediction <- factor(plt$Prediction)

ggplot(plt, aes(Reference ,Prediction, fill= Freq)) +
  geom_tile(col = 'black') + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  cowplot::theme_cowplot()+
  theme(legend.position = "none") +
  labs(x = "Reference",y = "Prediction")# +

final_fit <- extract_workflow(final_rf)
final_rec <- extract_recipe(final_rf)
final_model <- extract_fit_parsnip(final_rf)
test_data_processed <- bake(final_rec, new_data = testing(split_data)) %>% select(!survey_date)

#Calculate Variable Importances
vi.model = vi(final_model) 
vi.model = vi.model %>% mutate(Type = case_when(
  grepl("anthracnose$|brown_patch|crown_rust$", Variable) ~ "Biological",
  grepl("env_", Variable) ~ "Environmental",
  TRUE ~ "Other")) %>%
  head(150) %>% arrange(Importance)
vi.model$Variable <- gsub("_", "\n", vi.model$Variable)
vi.model$Variable <- gsub("env\n", "", vi.model$Variable)
vi.model$Variable <- gsub("\nin|\nmicromolm2s|\nwm2|\nm3m3", "", vi.model$Variable)
vi.model$Variable = factor(vi.model$Variable, levels = vi.model$Variable)

#Fifteen most important variables
vi.model %>% tail(15) %>%
  ggplot(aes(x=Importance, y=Variable, fill = Type)) + 
  geom_bar(stat="identity", col = "grey20") + theme_bw() + ggtitle("Rust Observation") +
  theme(text=element_text(size = 8)) + 
  scale_fill_grey(start = .1, end = .9) + 
  theme(legend.position = c(.7, .2),
        text=element_text(size = 8), axis.text = element_text(size=6), legend.title = element_text(size = 8))+
  theme(legend.key.size = unit(0.35, "cm"))+
  scale_y_discrete(labels = scales::label_wrap(width = 35))+ ylab("") + ggtitle("") # Adjust the number of features
ggsave("results/RustRF_Final.pdf", width = 3.25, height =  4)

#check correlation between variables. take top 3 most important and non-correlated variables (at the .8 level)
# do this by hand

numeric_data <- test_data_processed %>% dplyr::select(where(is.numeric))
voi = vi(final_model) %>% head(15)
cor_rf = cor(numeric_data)[voi$Variable, voi$Variable] %>% data.frame()

colnames(cor_rf) <- gsub("_", " ", colnames(cor_rf))
colnames(cor_rf) <- gsub("env ", "", colnames(cor_rf))
colnames(cor_rf) <- gsub(" in| micromolm2s| wm2| m3m3", "", colnames(cor_rf))
colnames(cor_rf) <- substr(colnames(cor_rf), 1, 35)
#colnames(cor_rf) = str_wrap(colnames(cor_rf), width = 35)
rownames(cor_rf) = colnames(cor_rf)

pdf("results/RF_Rust_VarCorr.pdf",width = 6.5, height = 6.5)
pheatmap(abs(cor_rf), cluster_rows = T, cluster_cols = T)
dev.off()


#maximum photosynthetically active radiation
#lag anthracnose
#any interaction term
#average wind speed
#maximum relative humidity
#This is coming from the iml library
predictor <- Predictor$new(model=final_model$fit,
                           data = dplyr::select(test_data_processed,-c(crown_rust)),
                           y=test_data_processed$crown_rust)

get_ale = names(test_data_processed)[grepl("^env_maximum_photosynthetically|^env_maximum_relative_humid|^env_average_wind_speed|lag_anthracnose", names(test_data_processed))]
get_ale = names(test_data_processed)[names(test_data_processed) %in% voi$Variable]
#create a plot for each one
plist <- list()
for(var in 1:length(get_ale)){
  var1 <- gsub("_", " ", get_ale[var])
  var1 <- gsub("env ", "",var1)
  var1 = paste(strwrap(var1, 30), collapse = '\n')
  plist[[var]] = viz_ales_1d(predictor, get_ale[var]) + xlab(var1) + theme(text = element_text(size = 6), axis.text = element_text(size = 6))
  ggsave(filename = paste("results/CR_RF_ALE_", get_ale[var], ".pdf", sep = ""), plot = plist[[var]], width = (6.5-3.25)/2, height = 4/2)
}


library(cowplot)
plot_grid(plotlist = plist, nrow = 5, ncol = 3)
ggsave("results/allALEs_RUST.pdf", width = 6.5, height = 8)

#modify based on var of interest
get_ale_pathos <- names(test_data_processed)[grepl("^anth|^brown", names(test_data_processed))]
out.ale <- data.frame()
for(pred in 1:length(get_ale_pathos)){
  path.ale <- FeatureEffect$new(predictor,
                                feature = get_ale_pathos[pred],
                                method = "ale",
                                grid.size = 25)$results %>%
    filter(.class == "X1") %>% filter(!!sym(get_ale_pathos[pred]) == 1) %>%
    mutate(patho = get_ale_pathos[pred]) %>% select(!matches(get_ale_pathos[pred]))
  out.ale <- rbind(out.ale, path.ale)
}

#special barplot
dis.pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(1:2, 3:4)], "#FEE99A", "#FEC601")

out.ale %>% 
  ggplot(aes(x=patho, y = .value, fill = patho)) + geom_bar(stat = "identity") + xlab("Disease Presence") + ylab("ALE") + 
  theme_bw() + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c("#1F78B4", "#33A02C")) + 
  theme(legend.position = "none", text = element_text(size = 6))

ggsave(filename = paste("results/CR_RF_ALE_", "coinfection", ".pdf", sep = ""), width =(6.5-3.25)/2, height = 4/2)

###### ------------- MACHINE LEARNING FOR BROWN PATCH ---------------- ##############
set.seed(222)

metric <- metric_set(rmse)
var_only <- rf.data %>%
  mutate(crown_rust = as.numeric(as.character(crown_rust))) %>% select(!c(npoints)) %>%
  mutate(anthracnose = as.numeric(as.character(anthracnose)), brown_patch = as.numeric(as.character(brown_patch))) %>% arrange(survey_date) %>%
  select(!pc1:pc4)

split_data = initial_time_split(var_only, prop = .8)

rf_spec <- rand_forest(
  mtry = tune(), #tune the mtry hyperparameter here
  trees = 3000, #up the number of trees to better stabilize variable importances
  mode = "classification") %>%
  set_engine("ranger", importance = "permutation")

#set up a workflow
rf_workflow <- workflow() %>%
  add_recipe(recipe_mod) %>%
  add_model(rf_spec) 

#Try mtry values from 4 to 32, this determines the number of variables each RF consideres
grid <- grid_regular(mtry(range = c(4, 32)), levels = 8)

#Run it with temporal resampling.... doing some hefty lifting here. 
#will get some errors for some splits that lack data/temporal resolution, I think this is okay
library(future)
plan(multisession, workers = 11) 
rf_results <- tune_grid(
  rf_workflow,
  resamples = temporal_split,
  grid = grid
) 

plan(sequential)
#basic recipe
basic_rec <- recipe(brown_patch ~ ., data = training(split_data)) %>%
  update_role(survey_date, new_role = "SamplingDate") %>%
  update_role(plot_id, new_role = "PlotID") %>%
  step_zv(all_predictors())

recipe_mod <- basic_rec |>
  step_novel(plot_id) %>%
  step_dummy(all_nominal_predictors()) |>
  step_mutate(brown_patch = factor(brown_patch)
  ) |>
  step_upsample(brown_patch)

rf_spec <- rand_forest(
  mtry = 4,
  trees = 3000,
  mode = "classification") %>%
  set_engine("ranger", importance = "permutation")

#set up a workflow
rf_workflow <- workflow() %>%
  add_recipe(recipe_mod) %>%
  add_model(rf_spec) 

library(future)
plan(multisession, workers = 11) #not quire sure how to close this connection :( 

rf_results <- tune_grid(
  rf_workflow,
  resamples = temporal_split
  #control = ctrl_grid
) 

plan(sequential)
#  COLLECT THE BEST RESULTS
best_params <- select_best(rf_results, metric = "roc_auc")
best_params

#best mtry = 4
#FINALIZE THE MODEL WITH THE BEST PARM SPACE
final_rf <- rf_workflow %>%
  finalize_workflow(select_best(rf_results, metric = "roc_auc")) %>%
  last_fit(split_data)
save(final_rf, file="results/bp_rf_mtry4_roll7_upsample.Rdata")

load("results/bp_rf_mtry4_roll7_upsample.Rdata")

all_metrics = collect_metrics(final_rf)
all_metrics
predictions_tbl <- collect_predictions(final_rf)
predictions_tbl

test_metrics <- predictions_tbl %>%
  metrics(truth = brown_patch, estimate = .pred_class)
test_metrics #kappa greater than 0 = better than random chance. Good values start at about .21
test_metrics <- predictions_tbl %>%
  f_meas(truth = brown_patch, estimate = .pred_class) #High values = better!
test_metrics
predictions_tbl %>% 
  yardstick::precision(brown_patch, .pred_class)

predictions_tbl %>%
  yardstick::recall(brown_patch, .pred_class)

final_fit <- extract_workflow(final_rf)
final_rec <- extract_recipe(final_rf)
final_model <- extract_fit_parsnip(final_rf)
test_data_processed <- bake(final_rec, new_data = testing(split_data)) %>% select(!survey_date)


library(caret)
out <- confusionMatrix(predictions_tbl$.pred_class, predictions_tbl$brown_patch)

plt <- as.data.frame(out$table)
plt$Prediction <- factor(plt$Prediction)

ggplot(plt, aes(Reference ,Prediction, fill= Freq)) +
  geom_tile(col = 'black') + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  cowplot::theme_cowplot()+
  theme(legend.position = "none") +
  labs(x = "Reference",y = "Prediction")




# Plot variable importance case_when(grepl("_x_", Variable) ~"Interaction",
library(vip)
vi.model = vi(final_model) 
vi.model = vi.model %>% mutate(Type = case_when(
  grepl("crown_rust$|anthracnose|brown_patch$", Variable) ~ "Biological",
  grepl("env_", Variable) ~ "Environmental",
  TRUE ~ "Other")) %>%
  head(150) %>% arrange(Importance)
vi.model$Variable <- gsub("_", "\n", vi.model$Variable)
vi.model$Variable <- gsub("env\n", "", vi.model$Variable)
vi.model$Variable <- gsub("\nin|\nmicromolm2s|\nwm2|\nm3m3", "", vi.model$Variable)
vi.model$Variable = factor(vi.model$Variable, levels = vi.model$Variable)


vi.model %>% tail(15) %>%
  ggplot(aes(x=Importance, y=Variable, fill = Type)) + 
  geom_bar(stat="identity", col = "grey20") + theme_bw() + ggtitle("Brown Patch Observation") +
  theme(text=element_text(size = 8)) + 
  scale_fill_grey(start = .1, end = .9) + 
  theme(legend.position = c(.7, .2),
        text=element_text(size = 8), axis.text = element_text(size=6), legend.title = element_text(size = 8))+
  theme(legend.key.size = unit(0.35, "cm"))+
  scale_y_discrete(labels = scales::label_wrap(width = 35))+ ylab("") + ggtitle("") # Adjust the number of features
ggsave("results/BP_RF_Final.pdf", width = 3.25, height =  4)

#check correlation between variables. take top 3 non-correlated, #am I missing variables on my inrf?
numeric_data <- test_data_processed %>% dplyr::select(where(is.numeric))
voi = vi(final_model) %>% head(15)
cor_rf = cor(numeric_data)[voi$Variable, voi$Variable] %>% data.frame()

colnames(cor_rf) <- gsub("_", " ", colnames(cor_rf))
colnames(cor_rf) <- gsub("env ", "", colnames(cor_rf))
colnames(cor_rf) <- gsub(" in| micromolm2s| wm2| m3m3", "", colnames(cor_rf))
colnames(cor_rf) <- substr(colnames(cor_rf), 1, 35)
#colnames(cor_rf) = str_wrap(colnames(cor_rf), width = 35)
rownames(cor_rf) = colnames(cor_rf)

#par(lheight = 0.3)
#pdf("results/RF_Anth_VarCorr.pdf",width = 6.5, height = 6.5)
pheatmap(cor_rf, cluster_rows = T, cluster_cols = T)
#dev.off()

library(iml)
predictor <- Predictor$new(model=final_model$fit,
                           data = dplyr::select(test_data_processed,-c(brown_patch)),
                           y=test_data_processed$brown_patch)

get_ale = names(test_data_processed)[grepl("^env_maximum_soil_temp|^env_average_wind_speed|^env_maximum_relative_humidity|lag_anth|lag_crown", names(test_data_processed))]
get_ale = names(test_data_processed)[names(test_data_processed) %in% voi$Variable]

#create a plot for each one
plist <- list()
for(var in 1:length(get_ale)){
  var1 <- gsub("_", " ", get_ale[var])
  var1 <- gsub("env ", "",var1)
  var1 = paste(strwrap(var1, 30), collapse = '\n')
  plist[[var]] = viz_ales_1d(predictor, get_ale[var]) + xlab(var1) + theme(text = element_text(size = 6), axis.text = element_text(size = 6))
  ggsave(filename = paste("results/BP_RF_ALE_", get_ale[var], ".pdf", sep = ""), plot = plist[[var]], width = (6.5-3.25)/2, height = 4/2)
}



library(cowplot)
plot_grid(plotlist = plist, nrow = 5, ncol = 3)
ggsave("results/allALEs_BP.pdf", width = 6.5, height = 8)

#modify based on var of interest
get_ale_pathos <- names(test_data_processed)[grepl("^crown|^anth", names(test_data_processed))]
out.ale <- data.frame()
for(pred in 1:length(get_ale_pathos)){
  path.ale <- FeatureEffect$new(predictor,
                                feature = get_ale_pathos[pred],
                                method = "ale",
                                grid.size = 25)$results %>%
    filter(.class == "X1") %>% filter(!!sym(get_ale_pathos[pred]) == 1) %>%
    mutate(patho = get_ale_pathos[pred]) %>% select(!matches(get_ale_pathos[pred]))
  out.ale <- rbind(out.ale, path.ale)
}

#special barplot
dis.pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(1:2, 3:4)], "#FEE99A", "#FEC601")

out.ale %>% 
  ggplot(aes(x=patho, y = .value, fill = patho)) + geom_bar(stat = "identity") + xlab("Disease Presence") + ylab("ALE") + 
  theme_bw() + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c( "#1F78B4", "#FEC601")) + 
  theme(legend.position = "none", text = element_text(size = 6))

ggsave(filename = paste("results/BP_RF_ALE_", "coinfection", ".pdf", sep = ""), width =(6.5-3.25)/2, height = 4/2)


###### ------------- MACHINE LEARNING FOR ANTHRACNOSE ---------------- ##############
set.seed(222)

var_only <- rf.data %>%
  mutate(crown_rust = as.numeric(as.character(crown_rust))) %>% select(!c(npoints)) %>%
  mutate(anthracnose = as.numeric(as.character(anthracnose)), brown_patch = as.numeric(as.character(brown_patch))) %>% arrange(survey_date) %>%
  select(!pc1:pc4)

split_data = initial_time_split(var_only, prop = .8)

library(themis)

temporal_split <- sliding_period(
  data = training(split_data),
  index = "survey_date",  # The temporal index column
  period = "month",
  lookback = 12,            # Training window
  assess_stop = 1)  

#basic recipe
basic_rec <- recipe(anthracnose ~ ., data = training(split_data)) %>%
  update_role(survey_date, new_role = "SamplingDate") %>%
  update_role(plot_id, new_role = "PlotID") %>%
  step_zv(all_predictors())

recipe_mod <- basic_rec |>
  step_novel(plot_id) %>%
  step_dummy(all_nominal_predictors()) |>
  step_mutate(anthracnose = factor(anthracnose)
  ) |>
  step_upsample(anthracnose)

#set up a model

rf_spec <- rand_forest(
  mtry = tune(),
  trees = 3000, #up the number of trees to better stabilize variable importances
  mode = "classification") %>%
  set_engine("ranger", importance = "permutation")

#set up a workflow
rf_workflow <- workflow() %>%
  add_recipe(recipe_mod) %>%
  add_model(rf_spec) 


grid <- grid_regular(mtry(range = c(4, 32)), levels = 8)

library(future)
plan(multisession, workers = 11) #not quire sure how to close this connection :( 

rf_results <- tune_grid(
  rf_workflow,
  resamples = temporal_split,
  grid = grid
) 

plan(sequential)
#  COLLECT THE BEST RESULTS
best_params <- select_best(rf_results, metric = "roc_auc")
best_params #best mtry is 4

#FINALIZE THE MODEL WITH THE BEST PARM SPACE
final_rf <- rf_workflow %>%
  finalize_workflow(select_best(rf_results, metric = "roc_auc")) %>%
  last_fit(split_data)
save(final_rf, file="results/w7_anth_rf_mtry4_upsample.Rdata")

load("results/w7_anth_rf_mtry4_upsample.Rdata")

all_metrics = collect_metrics(final_rf)
all_metrics
predictions_tbl <- collect_predictions(final_rf)

test_metrics <- predictions_tbl %>%
  metrics(truth = anthracnose, estimate = .pred_class)
test_metrics #kappa greater than 0 = better than random chance. Good values start at about .21
test_metrics <- predictions_tbl %>%
  f_meas(truth = anthracnose, estimate = .pred_class) #High values = better!
test_metrics
predictions_tbl %>% 
  yardstick::precision(anthracnose, .pred_class)
predictions_tbl %>% 
  yardstick::recall(anthracnose, .pred_class)


final_fit <- extract_workflow(final_rf)
final_rec <- extract_recipe(final_rf)
final_model <- extract_fit_parsnip(final_rf)
test_data_processed <- bake(final_rec, new_data = testing(split_data)) %>% select(!survey_date)

library(caret)
out <- confusionMatrix(predictions_tbl$.pred_class, predictions_tbl$anthracnose)

plt <- as.data.frame(out$table)
plt$Prediction <- factor(plt$Prediction)

ggplot(plt, aes(Reference ,Prediction, fill= Freq)) +
  geom_tile(col = 'black') + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  cowplot::theme_cowplot()+
  theme(legend.position = "none") +
  labs(x = "Reference",y = "Prediction")


# Plot variable importance 
vi.model = vi(final_model) 
vi.model = vi.model %>% mutate(Type = case_when(
  grepl("crown_rust$|brown_patch|anthracnose$", Variable) ~ "Biological",
  grepl("env_", Variable) ~ "Environmental",
  TRUE ~ "Other")) %>%
  head(150) %>% arrange(Importance)
vi.model$Variable <- gsub("_", "\n", vi.model$Variable)
vi.model$Variable <- gsub("env\n", "", vi.model$Variable)
vi.model$Variable <- gsub("\nin|\nmicromolm2s|\nwm2|\nm3m3", "", vi.model$Variable)
vi.model$Variable = factor(vi.model$Variable, levels = vi.model$Variable)


vi.model %>% tail(15) %>%
  ggplot(aes(x=Importance, y=Variable, fill = Type)) + 
  geom_bar(stat="identity", col = "grey20") + theme_bw() + ggtitle("Rust Observation") +
  theme(text=element_text(size = 8)) + 
  scale_fill_grey(start = .5, end = .9) + 
  theme(legend.position = c(.72, .2),
        text=element_text(size = 8), axis.text = element_text(size=6), legend.title = element_text(size = 8))+
  theme(legend.key.size = unit(0.35, "cm"))+
  scale_y_discrete(labels = scales::label_wrap(width = 35))+ ylab("") +
  xlim(c(0, .05))+
  ggtitle("") # Adjust the number of features
ggsave("results/AnthRF_Final.pdf", width = 3.25, height =  4)

#check correlation between variables. take top 3 most important & non-correlated
numeric_data <- test_data_processed %>% dplyr::select(where(is.numeric))
voi = vi(final_model) %>% head(15)
cor_rf = cor(numeric_data)[voi$Variable, voi$Variable] %>% data.frame()

colnames(cor_rf) <- gsub("_", " ", colnames(cor_rf))
colnames(cor_rf) <- gsub("env ", "", colnames(cor_rf))
colnames(cor_rf) <- gsub(" in| micromolm2s| wm2| m3m3", "", colnames(cor_rf))
colnames(cor_rf) <- substr(colnames(cor_rf), 1, 35)
rownames(cor_rf) = colnames(cor_rf)

par(lheight = 0.3)
pdf("results/RF_Anth_VarCorr.pdf",width = 6.5, height = 6.5)
pheatmap(cor_rf, cluster_rows = T, cluster_cols = T)
dev.off()

predictor <- Predictor$new(model=final_model$fit,
                           data = dplyr::select(test_data_processed,-c(anthracnose)),
                           y=test_data_processed$anthracnose)

#these come from the correlation matrix, manually defined
get_ale = names(test_data_processed)[grepl("^env_average_soil_temp|^env_maximum_relative_humidity|^env_minimum_soil_moisture|^env_maximum_soil_moisture", names(test_data_processed))]
get_ale = names(test_data_processed)[names(test_data_processed) %in% voi$Variable]
#create a plot for each one
plist <- list()
for(var in 1:length(get_ale)){
  var1 <- gsub("_", " ", get_ale[var])
  var1 <- gsub("env ", "",var1)
  var1 = paste(strwrap(var1, 30), collapse = '\n')
  plist[[var]] = viz_ales_1d(predictor, get_ale[var]) + xlab(var1) + theme(text = element_text(size = 6), axis.text = element_text(size = 6))
  ggsave(filename = paste("results/ANTH_RF_ALE_", get_ale[var], ".pdf", sep = ""), plot = plist[[var]], width = (6.5-3.25)/2, height = 4/2)
}

library(cowplot)
plot_grid(plotlist = plist, nrow = 5, ncol = 3)
ggsave("results/allALEs_ANTH.pdf", width = 6.5, height = 8)

#modify based on var of interst
get_ale_pathos <- names(test_data_processed)[grepl("^crown|^brown", names(test_data_processed))]
out.ale <- data.frame()
for(pred in 1:length(get_ale_pathos)){
  path.ale <- FeatureEffect$new(predictor,
                                feature = get_ale_pathos[pred],
                                method = "ale",
                                grid.size = 25)$results %>%
    filter(.class == "X1") %>% filter(!!sym(get_ale_pathos[pred]) == 1) %>%
    mutate(patho = get_ale_pathos[pred]) %>% select(!matches(get_ale_pathos[pred]))
  out.ale <- rbind(out.ale, path.ale)
}

#special barplot
dis.pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(1:2, 3:4)], "#FEE99A", "#FEC601")

out.ale %>% 
  ggplot(aes(x=patho, y = .value, fill = patho)) + geom_bar(stat = "identity") + xlab("Disease Presence") + ylab("ALE") + 
  theme_bw() + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c( "#33A02C", "#FEC601")) + 
  theme(legend.position = "none", text = element_text(size = 6))

ggsave(filename = paste("results/ANTH_RF_ALE_", "coinfection", ".pdf", sep = ""), width =(6.5-3.25)/2, height = 4/2)


