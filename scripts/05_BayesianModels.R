library(ape)
library(brms)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(cowplot)

load("data/Env+DiseaseData.Rdata")
#####------------- CROWN RUST ---------------------

ci_colors = c("lightgray", RColorBrewer::brewer.pal(5, "YlOrRd")[1:5])
names(ci_colors) = c("<80%", "80-90%", "90-95%", "95-99%", "99-100%")



data.rust = rf.data %>% mutate(month = factor(month), anthracnose=factor(anthracnose), brown_patch = factor(brown_patch))
prior = c(set_prior("normal(0, 2)", class = "b")) #sin function and MDS3 are basically the same (70% correlated)
options(mc.cores = parallel::detectCores()-1)
month.brm.rust <- brm(
  crown_rust ~  lag_crown_rust +  month + year + 
    pc1 + pc2 + pc3 + pc4 + 
    anthracnose + brown_patch + anthracnose:brown_patch +
    anthracnose*(pc1 + pc2 + pc3 + pc4) + 
    brown_patch*(pc1 + pc2 + pc3 + pc4) + 
    (1 | plot_id),
  data = data.rust,
  family = "bernoulli",
  prior = prior,
  cores = 4,
  thin = 4, 
  iter = 6000,# parallel computation, adjust as needed
  threads = threading(2)                       # within-chain threads (cores × threads ≤ total CPUs)
) 

save(month.brm.rust, file = "results/bigdata.month.brm.rust.Rdata")
pp_check(month.brm.rust, ndraws = 50)
ggsave("results/ppcheck_rust.pdf", width = 3, height = 3)

plot(month.brm.rust)


sumrust = summary(month.brm.rust)
write.csv(signif(sumrust$fixed, digits = 3), file ="Rust_Brms_Summmary.csv")


## load from save: rust
load("results/bigdata.month.brm.rust.Rdata")
draws <- as_draws_array(month.brm.rust)
signif = summarise_draws(draws, default_summary_measures()) %>% filter(grepl("b_", variable)) %>%
  mutate(contain_zero = case_when(0 > q5 & 0 < q95 ~ "y",
                                  0 > q5 & 0 > q95 ~"n",
                                  0<q5 & 0<q95 ~ "n",
                                  0 > q5& 0>q95 ~ "n")) %>%
  mutate(pr_contain0 = 1 - ifelse(mean < 0,
                                  pnorm(0, mean, sd, lower.tail = F),
                                  pnorm(0, mean, sd, lower.tail = T))) %>%
  mutate(`Pr(Non-Zero)` = case_when(pr_contain0 >= .99 ~ "99-100%",
                                    pr_contain0 >= .95 & pr_contain0 < .99 ~ "95-99%",
                                    pr_contain0 >= .9 & pr_contain0 < .95 ~ "90-95%",
                                    pr_contain0 >= .8 & pr_contain0 < .9 ~ "80-90%",
                                    pr_contain0 < .8 ~ "<80%"))

posterior_beta <- month.brm.rust %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(intercept = str_detect(.variable, "Intercept"),
         month = str_detect(.variable, "month"), 
         pc = str_detect(.variable, "PC"), 
         ixn = str_detect(.variable, ":"),
         base_patho = str_detect(.variable, "b_anth|b_brown|b_crown_rust")) %>% left_join(signif, by = c(".variable"= "variable"))


#compare interaction and base effects
allinter = posterior_beta %>% filter(grepl("Intercept|month|year|lag|anthracnose1:brown", .variable) == F) %>%
  mutate(nonzero_prop = case_when(
    # CI entirely above zero
    q5 > 0 ~ 1,
    # CI entirely below zero
    q95 < 0 ~ 1,
    # CI crosses zero — proportion NOT overlapping zero
    mean > 0 ~ q95 / (q95 - q5),
    mean < 0 ~ abs(q5) / (q95 - q5)
  ),
  nonzero = case_when(
    nonzero_prop < 0.80 ~ "<80%",
    nonzero_prop < 0.90 ~ "80-90%",
    nonzero_prop < 0.95 ~ "90-95%",
    nonzero_prop < 0.99 ~ "95-99%",
    TRUE                ~ "99-100%"
  )) 
allinter %>% ggplot(aes(y=fct_rev(.variable), x=.value, fill = `Pr(Non-Zero)`)) + 
  geom_vline(xintercept = 0, col = "blue") +
  stat_slab(slab_alpha = .6,
            height = 2, color = "gray20", lwd = .5,
            expand = F, trim = F, density = "unbounded") +  #changed from stat halfeye #FFFFB2
  scale_fill_manual(values = ci_colors[unique(allinter$`Pr(Non-Zero)`)])+
  theme_bw() + theme(text = element_text(size = 8)) + ylab("") + xlab("Log-Odds of Rust Observation")
ggsave("results/logOdds_Rust.pdf", width = 6.5, height = 4)


###Conditional effects rust
rust.pos = data.rust %>% filter(crown_rust == 1) 
rust.neg = data.rust %>% filter(crown_rust == 0)



plist <- list()

for(patho in c("anthracnose", "brown_patch")){
  for(pc in c("pc1", "pc2", "pc3", "pc4")){
    ce = conditional_effects(month.brm.rust, effects = c(paste(pc, patho, sep = ":")), prob = 0)
    short.path = ifelse(patho == "anthracnose", "anth", 
                        ifelse(patho == "brown_patch", "bp", "rust"))
    plotting.colors = dis.pal[grepl(short.path, names(dis.pal))]
    names(plotting.colors) = NULL
    ce_range = max(ce[[1]][,"estimate__"]) - min(ce[[1]][,"estimate__"])
    p <- ggplot() +
      geom_line(data = ce[[1]], aes(x=!!sym(pc), y=estimate__, col = !!sym(patho)), lwd = 1)  +
      cowplot::theme_cowplot() +
      scale_color_manual(values=plotting.colors) + 
      geom_bin_2d(data= rust.pos, aes(x=!!sym(pc), y= (max(ce[[1]][,"estimate__"])+(.15*ce_range))), bins = list(y=10, x=150)) +
      scale_fill_gradient(low = "gray90", high = "black", transform = "identity")+
      geom_bin_2d(data= rust.neg, aes(x=!!sym(pc), y = (min(ce[[1]][,"estimate__"])-(.15*ce_range))), bins = list(y=10, x=150)) +
      theme(legend.position = "none",
            text=element_text(size = 6), 
            axis.text = element_text(size = 6)) + ylab("Pr(Rust)") #+
    plist[[length(plist) + 1]] = p
    ggsave(paste("results/CRbrms_ce_", patho, pc, ".pdf", sep = ""), plot = p, width = 1.5, height = 1.2)
  }
}

#####------------- Anthracnose ---------------------


data.anth = rf.data %>% mutate(month = factor(month), crown_rust=factor(crown_rust), brown_patch = factor(brown_patch))
prior = c(set_prior("normal(0, 2)", class = "b"))
options(mc.cores = parallel::detectCores()-1)
month.brm.anth <- brm(
  anthracnose ~  lag_anthracnose +  month + year + 
    pc1 + pc2 + pc3 + pc4 + 
    crown_rust + brown_patch + crown_rust:brown_patch +
    crown_rust*(pc1 + pc2 + pc3 + pc4) + 
    brown_patch*(pc1 + pc2 + pc3 + pc4) + 
    (1 | plot_id),
  data = data.anth,
  family = "bernoulli",
  prior = prior,
  cores = 4,
  thin = 4, 
  iter = 6000,
  threads = threading(2)                      
) 

save(month.brm.anth, file = "results/bigdata.month.brm.anth.Rdata")
pp_check(month.brm.anth, ndraws = 50)
ggsave("results/ppcheck_anth.pdf", width = 3, height = 3)

plot(month.brm.anth)

sumanth = summary(month.brm.anth)
write.csv(signif(sumanth$fixed, digits = 3), file ="Anth_Brms_Summmary.csv")


## load from save: anthracnose
load("results/bigdata.month.brm.anth.Rdata")
draws <- as_draws_array(month.brm.anth)
signif = summarise_draws(draws, default_summary_measures()) %>% filter(grepl("b_", variable)) %>%
  mutate(contain_zero = case_when(0 > q5 & 0 < q95 ~ "y",
                                  0 > q5 & 0 > q95 ~"n",
                                  0<q5 & 0<q95 ~ "n",
                                  0 > q5& 0>q95 ~ "n")) %>%
  mutate(pr_contain0 = 1 - ifelse(mean < 0,
                                  pnorm(0, mean, sd, lower.tail = F),
                                  pnorm(0, mean, sd, lower.tail = T))) %>%
  mutate(`Pr(Non-Zero)` = case_when(pr_contain0 >= .99 ~ "99-100%",
                                    pr_contain0 >= .95 & pr_contain0 < .99 ~ "95-99%",
                                    pr_contain0 >= .9 & pr_contain0 < .95 ~ "90-95%",
                                    pr_contain0 >= .8 & pr_contain0 < .9 ~ "80-90%",
                                    pr_contain0 < .8 ~ "<80%"))

posterior_beta <- month.brm.anth %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(intercept = str_detect(.variable, "Intercept"),
         month = str_detect(.variable, "month"), 
         pc = str_detect(.variable, "PC"), 
         ixn = str_detect(.variable, ":"),
         base_patho = str_detect(.variable, "b_crown|b_brown|b_anthracnose")) %>% left_join(signif, by = c(".variable"= "variable"))


#compare interaction and base effects
allinter = posterior_beta %>% filter(grepl("Intercept|month|year|lag|crown_rust1:brown", .variable) == F) %>%
  mutate(nonzero_prop = case_when(
    q5 > 0 ~ 1,
    q95 < 0 ~ 1,
    mean > 0 ~ q95 / (q95 - q5),
    mean < 0 ~ abs(q5) / (q95 - q5)
  ),
  nonzero = case_when(
    nonzero_prop < 0.80 ~ "<80%",
    nonzero_prop < 0.90 ~ "80-90%",
    nonzero_prop < 0.95 ~ "90-95%",
    nonzero_prop < 0.99 ~ "95-99%",
    TRUE                ~ "99-100%"
  )) 
allinter %>% ggplot(aes(y=fct_rev(.variable), x=.value, fill = `Pr(Non-Zero)`)) + 
  geom_vline(xintercept = 0, col = "blue") +
  stat_slab(slab_alpha = .6,
            height = 2, color = "gray20", lwd = .5,
            expand = F, trim = F, density = "unbounded") +
  scale_fill_manual(values = ci_colors[unique(allinter$`Pr(Non-Zero)`)])+
  theme_bw() + theme(text = element_text(size = 8)) + ylab("") + xlab("Log-Odds of Anthracnose Observation")
ggsave("results/logOdds_Anth.pdf", width = 6.5, height = 4)


###Conditional effects anthracnose
anth.pos = data.anth %>% filter(anthracnose == 1) 
anth.neg = data.anth %>% filter(anthracnose == 0)



plist <- list()

for(patho in c("crown_rust", "brown_patch")){
  for(pc in c("pc1", "pc2", "pc3", "pc4")){
    ce = conditional_effects(month.brm.anth, effects = c(paste(pc, patho, sep = ":")), prob = 0)
    short.path = ifelse(patho == "crown_rust", "rust", 
                        ifelse(patho == "brown_patch", "bp", "anth"))
    plotting.colors = dis.pal[grepl(short.path, names(dis.pal))]
    names(plotting.colors) = NULL
    ce_range = max(ce[[1]][,"estimate__"]) - min(ce[[1]][,"estimate__"])
    p <- ggplot() +
      geom_line(data = ce[[1]], aes(x=!!sym(pc), y=estimate__, col = !!sym(patho)), lwd = 1)  +
      cowplot::theme_cowplot() +
      scale_color_manual(values=plotting.colors) + 
      geom_bin_2d(data= anth.pos, aes(x=!!sym(pc), y= (max(ce[[1]][,"estimate__"])+(.15*ce_range))), bins = list(y=10, x=150)) +
      scale_fill_gradient(low = "gray90", high = "black", transform = "identity")+
      geom_bin_2d(data= anth.neg, aes(x=!!sym(pc), y = (min(ce[[1]][,"estimate__"])-(.15*ce_range))), bins = list(y=10, x=150)) +
      theme(legend.position = "none",
            text=element_text(size = 6), 
            axis.text = element_text(size = 6)) + ylab("Pr(Anthracnose)") 
    plist[[length(plist) + 1]] = p
    ggsave(paste("results/ANTHbrms_ce_", patho, pc, ".pdf", sep = ""), plot = p, width = 1.5, height = 1.2)
  }
}

#####------------- BROWN PATCH ---------------------

data.bp = rf.data %>% mutate(month = factor(month), crown_rust=factor(crown_rust), anthracnose = factor(anthracnose))
prior = c(set_prior("normal(0, 2)", class = "b"))
options(mc.cores = parallel::detectCores()-1)
month.brm.bp <- brm(
  brown_patch ~  lag_brown_patch +  month + year + 
    pc1 + pc2 + pc3 + pc4 + 
    crown_rust + anthracnose + crown_rust:anthracnose +
    crown_rust*(pc1 + pc2 + pc3 + pc4) + 
    anthracnose*(pc1 + pc2 + pc3 + pc4) + 
    (1 | plot_id),
  data = data.bp,
  family = "bernoulli",
  prior = prior,
  cores = 4,
  thin = 4, 
  iter = 6000,
  threads = threading(2)                      
) 

save(month.brm.bp, file = "results/bigdata.month.brm.bp.Rdata")
pp_check(month.brm.bp, ndraws = 50)
ggsave("results/ppcheck_bp.pdf", width = 3, height = 3)

plot(month.brm.bp)

sumbp = summary(month.brm.bp)
write.csv(signif(sumbp$fixed, digits = 3), file ="BP_Brms_Summmary.csv")


## load from save: brown patch
load("results/bigdata.month.brm.bp.Rdata")
draws <- as_draws_array(month.brm.bp)
signif = summarise_draws(draws, default_summary_measures()) %>% filter(grepl("b_", variable)) %>%
  mutate(contain_zero = case_when(0 > q5 & 0 < q95 ~ "y",
                                  0 > q5 & 0 > q95 ~"n",
                                  0<q5 & 0<q95 ~ "n",
                                  0 > q5& 0>q95 ~ "n")) %>%
  mutate(pr_contain0 = 1 - ifelse(mean < 0,
                                  pnorm(0, mean, sd, lower.tail = F),
                                  pnorm(0, mean, sd, lower.tail = T))) %>%
  mutate(`Pr(Non-Zero)` = case_when(pr_contain0 >= .99 ~ "99-100%",
                                    pr_contain0 >= .95 & pr_contain0 < .99 ~ "95-99%",
                                    pr_contain0 >= .9 & pr_contain0 < .95 ~ "90-95%",
                                    pr_contain0 >= .8 & pr_contain0 < .9 ~ "80-90%",
                                    pr_contain0 < .8 ~ "<80%"))

posterior_beta <- month.brm.bp %>% 
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(intercept = str_detect(.variable, "Intercept"),
         month = str_detect(.variable, "month"), 
         pc = str_detect(.variable, "PC"), 
         ixn = str_detect(.variable, ":"),
         base_patho = str_detect(.variable, "b_crown|b_anth|b_brown_patch")) %>% left_join(signif, by = c(".variable"= "variable"))


#compare interaction and base effects
allinter = posterior_beta %>% filter(grepl("Intercept|month|year|lag|crown_rust1:anth", .variable) == F) %>%
  mutate(nonzero_prop = case_when(
    q5 > 0 ~ 1,
    q95 < 0 ~ 1,
    mean > 0 ~ q95 / (q95 - q5),
    mean < 0 ~ abs(q5) / (q95 - q5)
  ),
  nonzero = case_when(
    nonzero_prop < 0.80 ~ "<80%",
    nonzero_prop < 0.90 ~ "80-90%",
    nonzero_prop < 0.95 ~ "90-95%",
    nonzero_prop < 0.99 ~ "95-99%",
    TRUE                ~ "99-100%"
  )) 
allinter %>% ggplot(aes(y=fct_rev(.variable), x=.value, fill = `Pr(Non-Zero)`)) + 
  geom_vline(xintercept = 0, col = "blue") +
  stat_slab(slab_alpha = .6,
            height = 2, color = "gray20", lwd = .5,
            expand = F, trim = F, density = "unbounded") +
  scale_fill_manual(values = ci_colors[unique(allinter$`Pr(Non-Zero)`)])+
  theme_bw() + theme(text = element_text(size = 8)) + ylab("") + xlab("Log-Odds of Brown Patch Observation")
ggsave("results/logOdds_BP.pdf", width = 6.5, height = 4)


###Conditional effects brown patch
bp.pos = data.bp %>% filter(brown_patch == 1) 
bp.neg = data.bp %>% filter(brown_patch == 0)



plist <- list()

for(patho in c("crown_rust", "anthracnose")){
  for(pc in c("pc1", "pc2", "pc3", "pc4")){
    ce = conditional_effects(month.brm.bp, effects = c(paste(pc, patho, sep = ":")), prob = 0)
    short.path = ifelse(patho == "crown_rust", "rust", 
                        ifelse(patho == "anthracnose", "anth", "bp"))
    plotting.colors = dis.pal[grepl(short.path, names(dis.pal))]
    names(plotting.colors) = NULL
    ce_range = max(ce[[1]][,"estimate__"]) - min(ce[[1]][,"estimate__"])
    p <- ggplot() +
      geom_line(data = ce[[1]], aes(x=!!sym(pc), y=estimate__, col = !!sym(patho)), lwd = 1)  +
      cowplot::theme_cowplot() +
      scale_color_manual(values=plotting.colors) + 
      geom_bin_2d(data= bp.pos, aes(x=!!sym(pc), y= (max(ce[[1]][,"estimate__"])+(.15*ce_range))), bins = list(y=10, x=150)) +
      scale_fill_gradient(low = "gray90", high = "black", transform = "identity")+
      geom_bin_2d(data= bp.neg, aes(x=!!sym(pc), y = (min(ce[[1]][,"estimate__"])-(.15*ce_range))), bins = list(y=10, x=150)) +
      theme(legend.position = "none",
            text=element_text(size = 6), 
            axis.text = element_text(size = 6)) + ylab("Pr(Brown Patch)") 
    plist[[length(plist) + 1]] = p
    ggsave(paste("results/BPbrms_ce_", patho, pc, ".pdf", sep = ""), plot = p, width = 1.5, height = 1.2)
  }
}
