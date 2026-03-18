library(tidyverse)
library(nicheROVER)

load("data/Env+DiseaseData.Rdata")

nichein = rf.data %>% select(pc1:pc4, anthracnose, crown_rust, brown_patch) 
nichein = nichein %>% pivot_longer(brown_patch:anthracnose, names_to = "disease", values_to = "pres") %>% filter(pres == 1) %>%
  select(!pres)
nsamples <- 1e4
system.time({
  fish.par <- tapply(1:nrow(nichein), nichein$disease,
                     function(ii) niw.post(nsamples = nsamples, X = nichein[ii,1:4]))
})


# various parameter plots
clrs <- c("black", "red", "blue") # colors for each species

# mu1 (del15N), mu2 (del13C), and Sigma12
par(mar = c(4, 4, .5, .1)+.1, mfrow = c(1,3))
niche.par.plot(fish.par, col = clrs, plot.index = 1)
niche.par.plot(fish.par, col = clrs, plot.index = 2)
niche.par.plot(fish.par, col = clrs, plot.index = 1:2)
legend("topleft", legend = names(fish.par), fill = clrs)

niche.par.plot(fish.par, col = clrs, plot.mu = TRUE, plot.Sigma = FALSE)

#The plots represent the posterior probability that an individual from the species indicated by the 
#row will be found within the niche of the species indicated by the column header. 

clrs <- c("black", "red", "blue") # colors for each species
over.stat <- overlap(fish.par, nreps = nsamples, nprob = 1e4, alpha = .95)
pdf("results/NicheOverlap.pdf", width = 6, height = 6)
overlap.plot(over.stat, col = clrs, mean.cred.col = "turquoise", equal.axis = TRUE,
             xlab = "Overlap Probability (%) -- Niche Region Size: 95%")

dev.off()




# posterior distribution of niche size by species
fish.size <- sapply(fish.par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})


rbind(est = colMeans(fish.size),
      se = apply(fish.size, 2, sd))

clrs <- c("black", "red", "blue") # colors for each species
boxplot(fish.size, col = clrs, pch = 16, cex = .5,
        ylab = "Niche Size", xlab = "Species")


## try 2d plots
nsamples <- 10

# format data for plotting function
nichein = data.frame(nichein)
fish.par <- tapply(1:nrow(nichein), nichein$disease,
                   function(ii) niw.post(nsamples = nsamples, X = nichein[ii,1:4]))

fish.data <- tapply(1:nrow(nichein), nichein$disease, function(ii) X = nichein[ii,1:4])

niche.plot(niche.par = fish.par, niche.data = fish.data, pfrac = .05, ndens = 100, col = clrs)
legend("topleft", legend = names(fish.par), fill = clrs)



