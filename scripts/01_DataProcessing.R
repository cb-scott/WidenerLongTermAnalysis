##########################################
####### CLEAN INPUT DATA #################
##########################################
#NOTE: TO REPLICATE OUR ANALYSES YOU SHOULD NOT HAVE TO RUN THIS SCRIPT (UNLESS YOU WOULD LIKE TO REDO THE CLEANING)
# I'VE PROVIDED THE OUTPUTS AS INTERMEDIATE PRODUCTS IN THE GITHUB REPO THAT ARE LOADED IN SUBSEQUENT SCRIPTS
# YOU SHOULD ONLY NEED: "data/Env+DiseaseData.Rdata" and "

#Load in data from each of the past experiments. This is the long term plot surveys from widener 2021-2024,
#longitudinal data from 2018, the fungicide data from 2017-2019, and the treatment plot surveys from 2024

#For experiments where there was a treatment, keep only the control (untreated) plots. 
#For measurements at the leaf level, combine up to the plant level for consistency.

#Combine weather data from the NC Cardinal Station at the Durham Wastewater Treatment Plant
#The online NC Cardinal system only lets you download so many observations at a time, so they need to be manually combined here

#Calculate lagged disease values (disease prevalence in a given plot at the previous timepoint) and the 7-day rolling average of environmental varaiables.
#Combine all data together.



#this whole script uses the 7 day rolling average of environmental variables. Replace w <- 7 in the script to get a different rolling average.
library(readxl)
library(tidyverse)
library(pheatmap)
library(zoo)

#library(contsurvplot)
#source("scripts/helper_functions.R")

set.seed(222)
dis.pal <- c(RColorBrewer::brewer.pal(12, "Paired")[c(1:2, 3:4)], "#FEE99A", "#FEC601")
names(dis.pal) <- c("anth_n", "anth_y", "bp_n", "bp_y", "rust_n", "rust_y")

Widener21 <- list.files("data/Widener_2021_22", pattern = ".xlsx", full.names = T)
outWidener21 <- data.frame(read_xlsx(Widener21[2]))
for (file in 3:length(Widener21)){
  x <- read_xlsx(Widener21[file])
  x <- x[,intersect(colnames(x), colnames(outWidener21))]
  outWidener21 <- outWidener21[,intersect(colnames(x), colnames(outWidener21))]
  outWidener21 <- rbind(outWidener21, x)
}
Widener22.24 <- list.files("data/Widener (2022-24)", pattern = ".xlsx", full.names = T)
outWidener22.24 <- data.frame(read_xlsx(Widener22.24[1]))
for (file in 2:length(Widener22.24)){
  x <- read_xlsx(Widener22.24[file])
  x <- x[,intersect(colnames(x), colnames(outWidener22.24))]
  outWidener22.24 <- outWidener22.24[,intersect(colnames(x), colnames(outWidener22.24))]
  outWidener22.24 <- rbind(outWidener22.24, x)
}
keepcols = intersect(colnames(outWidener21), colnames(outWidener22.24))
outWidener22.24 <- rbind(outWidener21[,keepcols], outWidener22.24[,keepcols])
outWidener22.24 <- outWidener22.24 %>% mutate(Survey.Date=as.Date(Survey.Date, format = c("%Y-%m-%d")))
outWidener22.24 <- outWidener22.24 %>% drop_na() #check this later and see what's up


#widener long term plot data is at the tiller level, combine here
fungicide = read.csv("data/Rita_longitudinal/monthly_disease_survey_171819.csv", header = T) %>%
  filter(fungicide == "Never") %>% select(month, year, plot, survey_date, plant, leaf, anthracnose, crown_rust, brown_patch) %>%
  pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(plant, survey_date, plot, month, year, disease) %>%
  summarise(disease.sum = sum(presence)) %>% mutate(disease.sum = ifelse(disease.sum >= 1, 1, 0)) %>%
  pivot_wider(names_from = disease, values_from = disease.sum)

midmonth.control = read.csv("data/Rita_longitudinal/midmonthly_disease_survey_2019.csv", header = T) %>% filter(fungicide == "Never") %>%
  select(month, year, plot, survey_date, plant, leaf, anthracnose, crown_rust, brown_patch) %>%
  pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(plant, survey_date, plot, month, year, disease) %>%
  summarise(disease.sum = sum(presence)) %>% mutate(disease.sum = ifelse(disease.sum >= 1, 1, 0)) %>%
  pivot_wider(names_from = disease, values_from = disease.sum)

long2018 = read.csv("data/Rita_longitudinal/longitudinal_disease_survey_2018.csv", header = T)
merge2018 = long2018 %>% select(PlotLoc, ID, Leaf, anthracnose, crown_rust, brown_patch, Survey.Date) #only keep 2018
#get rid of individual plant ID's, we're not going to track them through time in the overall model
longitudinal =  merge2018 %>% pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(PlotLoc, ID, disease, Survey.Date) %>% summarise(disease.sum = sum(presence)) %>%
  mutate(disease.sum = ifelse(disease.sum >=1, 1, 0)) %>% pivot_wider(names_from = disease, values_from = disease.sum) %>%
  separate(PlotLoc, into = c("PlotID", "SubID")) %>% mutate(PlotID = as.numeric(PlotID)) %>% drop_na()

#2024 plot survey data
#want only unsprayed uninoculated plots
#cleaned in other script -- the dates were a mess!
plots2024 = read.csv("data/2024PlotSurveys/PlotSummary2024_CNTL_ONLY.csv", header = T) %>% select(!X)

#what are the datasets we have?
colnames(outWidener21)
colnames(outWidener22.24)
colnames(fungicide)
colnames(midmonth.control)
colnames(longitudinal)
colnames(plots2024)

widener2021 = outWidener21 %>% select(sub.array.ID, Rhiz.prev, Rust.Prev, Anth.prev, Survey.Date) %>%
  mutate(sub.array.ID = paste("WidenerLT", sub.array.ID, sep = "_"))
colnames(widener2021) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date")

widener2022.24 = outWidener22.24 %>% select(sub.array.ID, Rhiz.prev, Rust.Prev, Anth.prev, Survey.Date) %>%
  mutate(sub.array.ID = paste("WidenerLT", sub.array.ID, sep = "_"))
colnames(widener2022.24) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date")

fung171819 = fungicide %>% mutate(survey_date = as.Date(survey_date, format = "%m/%d/%y")) %>% ungroup() %>%
  mutate(plotID = paste(gsub("plotID", "fencedID", plot))) %>% select(plotID, brown_patch, crown_rust, anthracnose, survey_date)
colnames(fung171819) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date")  

midmonth19 = midmonth.control %>% mutate(survey_date = as.Date(survey_date, format = "%m/%d/%y")) %>% ungroup() %>%
  mutate(plotID = paste(gsub("plotID", "fencedID", plot))) %>% select(plotID, brown_patch, crown_rust, anthracnose, survey_date)
colnames(midmonth19) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date") 

long18 = longitudinal %>% mutate(Survey.Date = as.Date(Survey.Date, format = "%m/%d/%y")) %>% ungroup() %>%
  mutate(plotID = paste("fencedID", PlotID, sep = "_")) %>% select(plotID, brown_patch, crown_rust, anthracnose, Survey.Date)
colnames(long18) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date") 

plots24 = plots2024 %>% mutate(Survey.Date = as.Date(Survey.Date, format = "%Y-%m-%d")) %>% ungroup() %>%
  mutate(plotID = paste("WidenerPlot", Plot.ID, sep = "_")) %>% select(plotID, Rhiz.Prev, Rust.prev, Anth.Prev, Survey.Date)
colnames(plots24) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date")  

master.widener = rbind(widener2021, widener2022.24, fung171819, midmonth19, long18, plots24)
write.csv(master.widener, "data/MasterDataset.csv", quote = F, row.names = F)

#now we need the environmental data. 
#see what dates we need to pull from
range(master.widener$Survey.Date)

##############################
##### GET ALL WEATHER DATA ####
##############################
## Had to use the DURH cardinal station to get data that went back to 2017. THis is slightly further than CHAP
weather = list.files("data/EnviroData_DURH/", full.names = T)
readin = function(x) {
  tmp = read_excel(x, skip = 12, col_names = T)
  colnames(tmp) = make.names(colnames(tmp))
  tmp |> 
    mutate(across(everything(), as.character))
}

all_weather = bind_rows(lapply(weather, readin)) %>%
  mutate(across(!Date, as.numeric)) %>%
  mutate(across(contains(".F."), function(x) (x-32)/1.8)) %>%
  mutate(across(contains(".in"), function(x) x*2.54)) %>%
  mutate(across(contains("..mph."), function(x) x*1.609 ))
colnames(all_weather)
colnames(all_weather) = gsub("\\.F", ".C", colnames(all_weather))
colnames(all_weather) = gsub("\\.in", ".cm", colnames(all_weather))
colnames(all_weather) = gsub("\\.mph", ".kph", colnames(all_weather))

all_weather.noDate = all_weather %>% select(!Date) %>% data.frame()
pa.tmp = all_weather.noDate
pa.tmp[!is.na(pa.tmp)] = 1
pa.tmp[is.na(pa.tmp)] = 0
#get rid of bads
all_weather.noDate = all_weather.noDate[,colnames(pa.tmp[colSums(pa.tmp) > 2000])]
all_weather = drop_na(all_weather[,c("Date", colnames(pa.tmp[colSums(pa.tmp) > 2000]))])

nona.weather = drop_na(all_weather.noDate)

#let's princomp the read in and then take weekly averages over the MDS scores
library(vegan)
weather.stand = data.frame(decostand(nona.weather, method = "standardize"))
weather.pr = prcomp(weather.stand)
screeplot(weather.pr, bstick = T) #keep first 4 pcs
plot(scores(weather.pr)[,1:2])
plot(scores(weather.pr)[,c(1,3)])
weather.scores = data.frame(scores(weather.pr)[,1:4])
weather.scores$Survey.Date = as.Date(all_weather$Date)
weather.scores %>%
  ggplot(aes(x=PC1, y = PC2, col = month(Survey.Date))) + 
  geom_point() + theme_bw() + scale_color_viridis_c()

### Calculate Rolling Averages for Environmental Data ######
w <- 7  # window size in days
pc.roll = weather.scores %>%
  arrange(Survey.Date) %>%
  mutate(Npoints = 1:n() - findInterval(Survey.Date - w, Survey.Date))  %>%
  mutate(across(PC1:PC4, function(x) rollapplyr(x, Npoints, mean, partial = TRUE, fill = NA)))
env.roll = all_weather  %>%
  mutate(Survey.Date = as.Date(Date)) %>% select(!Date) %>% arrange(Survey.Date) %>% 
  mutate(Npoints = 1:n() - findInterval(Survey.Date - w, Survey.Date))  %>%
  mutate(across(!Survey.Date, function(x) rollapplyr(x, Npoints, mean, partial = TRUE, fill = NA))) %>% select(!Npoints)

#### Calculate PC Loadings ########
weather.loadings <- data.frame(weather.pr$rotation)
rownames(weather.loadings) = janitor::make_clean_names(rownames(weather.loadings))

library(pheatmap)
library(RColorBrewer)
library(pracma)
pdf("results/WeatherLoadings.pdf",width = 6.5, height = 6.5)
pheatmap(weather.loadings[,1:4], cluster_cols = F, color = colorRampPalette(brewer.pal(n = 7, name =
                                                                                         "BrBG"))(30))
dev.off()

#### Calculate lagged disease values #######

mean.prev = master.widener %>% group_by(PlotID, Survey.Date) %>% 
  arrange(Survey.Date, .by_group = T) %>%
  summarise(across(brown_patch:anthracnose, mean)) %>%
  ungroup()

lag.prev = mean.prev %>% group_by(PlotID) %>%
  mutate(across(brown_patch:anthracnose, lag)) %>%
  rename(lag.brown_patch= brown_patch,
         lag.crown_rust = crown_rust,
         lag.anthracnose = anthracnose)


####### Merge master dataset with environmental values and lagged prevalences #######
data.meta = left_join(left_join(master.widener, pc.roll), lag.prev) %>% select(!Npoints) %>%
  drop_na() %>% mutate(month = month(Survey.Date), year = year(Survey.Date))

rf.data = left_join(left_join(left_join(master.widener, pc.roll), lag.prev), env.roll) %>% 
  drop_na() %>% mutate(month = month(Survey.Date), year = year(Survey.Date))

colnames(rf.data) = 
  janitor::make_clean_names(colnames(rf.data))

colnames(rf.data)[14:45] = paste("env_", colnames(rf.data)[14:45], sep = "")

save(rf.data, file = "data/Env+DiseaseData.Rdata")
write.csv(rf.data, file = "data/Env+DiseaseData.csv")
