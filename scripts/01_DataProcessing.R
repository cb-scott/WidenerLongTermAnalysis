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
  filter(fungicide == "Never") %>%
  select(month, year, plot, survey_date, plant, leaf, anthracnose, crown_rust, brown_patch) %>%
  pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(plant, survey_date, plot, month, year, disease) %>%
  summarise(disease.sum = sum(presence)) %>% mutate(disease.sum = ifelse(disease.sum >= 1, 1, 0)) %>%
  pivot_wider(names_from = disease, values_from = disease.sum) %>%
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y"))

midmonth.control = read.csv("data/Rita_longitudinal/midmonthly_disease_survey_2019.csv", header = T) %>%
  filter(fungicide == "Never") %>%
  select(month, year, plot, survey_date, plant, leaf, anthracnose, crown_rust, brown_patch) %>%
  pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(plant, survey_date, plot, month, year, disease) %>%
  summarise(disease.sum = sum(presence)) %>% mutate(disease.sum = ifelse(disease.sum >= 1, 1, 0)) %>%
  pivot_wider(names_from = disease, values_from = disease.sum)%>%
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y"))

#I think that I accidentally kept the fungicide-treated plants here? Check
long2018 = read.csv("data/Rita_longitudinal/longitudinal_disease_survey_2018.csv", header = T)
#need the plot treatments
trtments = read.csv('data/Rita_longitudinal/plot_treatments.csv')
long2018 = long2018 %>% left_join(trtments) %>% filter(fungicide == "Never")
merge2018 = long2018 %>% select(plot, tiller, leaf, anthracnose, crown_rust, brown_patch, survey_date) #only keep 2018
#this isn't structured the way I thought it was
#get rid of individual plant ID's, we're not going to track them through time in the overall model
longitudinal =  merge2018 %>% 
  pivot_longer(anthracnose:brown_patch, names_to = "disease", values_to = "presence") %>%
  group_by(plot, tiller, disease, survey_date) %>% summarise(disease.sum = sum(presence)) %>%
  mutate(disease.sum = ifelse(disease.sum >=1, 1, 0)) %>% 
  pivot_wider(names_from = disease, values_from = disease.sum) %>%
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%Y"))


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

long18 = longitudinal %>% 
  mutate(survey_date = as.Date(survey_date, format = "%m/%d/%y")) %>% ungroup() %>%
  mutate(PlotID = paste("fencedID", gsub("plotID_", "", plot), sep = "_")) %>%
  select(PlotID, brown_patch, crown_rust, anthracnose, survey_date)
colnames(long18) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date") 

plots24 = plots2024 %>% mutate(Survey.Date = as.Date(Survey.Date, format = "%Y-%m-%d")) %>% ungroup() %>%
  mutate(plotID = paste("WidenerPlot", Plot.ID, sep = "_")) %>% 
  select(plotID, Rhiz.Prev, Rust.prev, Anth.Prev, Survey.Date)
colnames(plots24) = c("PlotID", "brown_patch", "crown_rust", "anthracnose", "Survey.Date")  

#something is wrong, why do my dates start at 2020
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

#just PC1
pc1 = data.frame(weather.loadings[c("average_soil_temperature_c", "average_air_temperature_c", "average_soil_moisture_m3_m3", "average_relative_humidity", "average_wind_speed_kph", "total_precipitation_cm"),1:4])
rownames(pc1)=c("Soil Temp", "Air Temp", "Soil Mositure", "Humidity", "Wind Speed", "Precipitation")

max_val <- max(abs(pc1), na.rm = TRUE)

# create breaks symmetric around zero
breaks <- seq(-max_val, max_val, length.out = 7)  # 31 breaks = 30 colors

dev.off()
pdf("results/simple_loadings.pdf", width = 1.6, height = 2.4)
pheatmap(pc1, cluster_cols = FALSE, cluster_rows = T,
         color = colorRampPalette(brewer.pal(n = 7, name = "BrBG"))(6),
         legend = T,
         fontsize = 8,
         breaks = breaks, 
         legend_breaks = c(-max_val,  0, max_val),
         legend_labels = c(round(-max_val,1),  "0", round(max_val,1)))
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


save(pc.roll, file = "data/weather.pcs.Rdata")
####### Merge master dataset with environmental values and lagged prevalences #######
data.meta = left_join(left_join(master.widener, pc.roll), lag.prev) %>% select(!Npoints) %>%
  drop_na() %>% mutate(month = month(Survey.Date), year = year(Survey.Date))

rf.data = left_join(left_join(left_join(master.widener, pc.roll), lag.prev), env.roll) %>% 
  drop_na() %>% mutate(month = month(Survey.Date), year = year(Survey.Date)) #why is this starting at 2020

colnames(rf.data) = 
  janitor::make_clean_names(colnames(rf.data))

colnames(rf.data)[14:45] = paste("env_", colnames(rf.data)[14:45], sep = "")

save(rf.data, file = "data/Env+DiseaseData.Rdata")
write.csv(rf.data, file = "data/Env+DiseaseData.csv")

long2018.cntl = read.csv("data/Rita_longitudinal/longitudinal_disease_survey_2018.csv", header = T)
#need the plot treatments
trtments = read.csv('data/Rita_longitudinal/plot_treatments.csv')
long2018.cntl = long2018.cntl %>% left_join(trtments) %>% filter(fungicide == "Never")

long.withweather = long2018.cntl %>% mutate(Survey.Date = as.Date(survey_date, format = "%m/%d/%Y"))%>% dplyr::select(!survey_date) %>%
  left_join(pc.roll)

save(long.withweather, file = "data/long2018weather.Rdata")

#######################################################
##### CREATE SUMMARY PLOTS OF DISEASE PREVALENCE ######
#######################################################

rf.data %>% mutate(year = year(survey_date), doy = as.Date(yday(survey_date), origin = "2020-01-01")) %>% group_by(plot_id, survey_date, doy, year) %>%
  summarise(across(brown_patch:anthracnose, mean, na.rm = T)) %>%
  pivot_longer(brown_patch:anthracnose, names_to = "disease", values_to = "prev") %>%
  ggplot(aes(x=doy, y=prev, group = year)) + geom_smooth(aes(col = year, fill = year), alpha = .2) + 
  theme_bw() + 
  scale_color_viridis_c(option = "viridis") +
  scale_fill_viridis_c(option = "viridis") +
  ylim(c(-.1,1.1)) +
  scale_x_continuous(breaks = c( seq(from=as.Date("01/15/2020", format = "%m/%d/%Y"),
                                     to= as.Date("12/15/2020", format = "%m/%d/%Y"),
                                     by = "month")),
                      labels = c("Jan", "Feb", "Mar", "April", "May", "Jun", "July", "Aug", "Sept", "Oct", "Nov", "Dec"))+
  theme(axis.text.x = element_text(angle = 90), text = element_text(size =8 )) + 
  xlab("") + 
  ylab("Prevalence") +
  facet_wrap(~disease, nrow = 3)
ggsave("Figures/ByYear_Prevalence.pdf", width = 6.5, height = 6.5)
