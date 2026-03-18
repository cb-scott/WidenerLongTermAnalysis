lag.prev <- function(indata, Var.OI){
  invar <- rlang::sym(Var.OI)
  lag.name = rlang::sym(paste(Var.OI, "lag", sep = "_"))
  var.name = rlang::sym(paste(Var.OI, "mean0", sep = "_"))
  lagged_prev = indata %>% mutate(invar = as.numeric(as.character(get(Var.OI)))) %>% 
    select(Xmeters, SurveyDate, invar) %>%
    group_by(Xmeters, SurveyDate) %>%
    summarise(!!var.name := mean(invar)) %>% mutate(!!lag.name := lag(!!var.name))
  return(lagged_prev)
}

restandardize <- function(dis.lag){
  restand = data.frame(decostand(dis.lag[,grepl("_lag|_mean0", colnames(dis.lag))], 
                                             method = "standardize"))
  restand = cbind(dis.lag[,c("coAR","Anth.prev", "Rust.Prev", "Rhiz.prev", "Xmeters", "PC1", "PC2", "PC3", "Year", "DoY", "Month")], restand) %>%
    mutate(Month = factor(Month))
  return(restand)
}

gather_posterior <- function(brmsobject){
  posterior_beta <- brmsobject %>% 
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(intercept = str_detect(.variable, "Intercept"),
           month = str_detect(.variable, "Month"), 
           mds = str_detect(.variable, "PC"), 
           coinf = str_detect(.variable, "prev"))
  draws <- as_draws_array(brmsobject)
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
  
  posterior_beta.time = posterior_beta %>% filter(month == T) %>%
    mutate(enviro = ifelse(str_detect(.variable, "PC1"), "PC1",
                           ifelse(str_detect(.variable,"PC2"), "PC2",
                                  ifelse(str_detect(.variable, "PC3"), "PC3",
                                         ifelse(str_detect(.variable, "PC4"), "PC4",
                                                ifelse(str_detect(.variable, "PC5"), "PC5", "base")))))) %>%
    mutate(basemonth = str_split_i(.variable, "_|:", 2)) %>% left_join(signif, by = c(".variable" = "variable"))
  posterior_beta.time$basemonth = factor(posterior_beta.time$basemonth, levels=str_sort(unique(posterior_beta.time$basemonth), numeric = T))
  
  levels(posterior_beta.time$basemonth) = c("Mar.", "April", "May", "June", "July", "Aug.", "Sept.", "Oct.", "Nov.", "Dec.")
  
  outlist = list(posterior_beta, signif, posterior_beta.time, draws)
  names(outlist) = c("beta.MEP", "all_signif", "beta.ME.Inter", "alldraws")
  return(outlist)
}

stand_classes <- function(indata, response){
  get_fakes.pres <- indata[indata[,response]==1,]
  get_fakes.abs <- indata[indata[,response] ==0,]
  fake_data <- rbind(get_fakes.abs,
                     get_fakes.pres[sample(nrow(get_fakes.pres),
                                           nrow(get_fakes.abs), replace = T),]) %>%
    arrange(Survey.Date) #mutated
  return(fake_data)
}

tr.ts.by.year <- function(fake_data, holdout){
  trainyears = unique(year(fake_data$survey_date))[!(unique(year(fake_data$survey_date)) %in% holdout)]
  testyears = unique(year(fake_data$survey_date))[(unique(year(fake_data$survey_date)) %in% holdout)]
  train <- fake_data[year(fake_data$survey_date) %in% trainyears,]
  test <- fake_data[year(fake_data$survey_date) %in% testyears,]
  out.list = list(train, test)
  names(out.list) = c("train", "test")
  return(out.list)
}

viz_ales_1d <- function(predictor, pred){
  getale <- FeatureEffect$new(predictor,
                              feature = pred,
                              method = "ale",
                              grid.size = 25)$results %>%
    filter(.class == "X1") 
  voi = pred
  p = getale %>%
    ggplot(aes(x=get(voi), y=.value)) +
    geom_line(lwd = .9) + theme_bw() + geom_hline(yintercept = 0) + 
    xlab(pred) + ylab("ALE") 
  return(p)
}

# To predict values we did not observe, use splines with NO SMOOTHING to connect points
yearly.ale.splines <- function(yearlymeans, ale.obj, voi){
  ale.obj = data.frame(ale.obj)
  HumiditySpline = smooth.spline(yearlymeans$DoY, yearlymeans[,voi])
  Humid.df = data.frame(DoY=HumiditySpline$x, var = HumiditySpline$y)
  
  spline_result <- smooth.spline(ale.obj[,voi], ale.obj$.value) #add a spline with no smoothing to predict values not observed by model
  newvals <- data.frame(predict(spline_result, x = Humid.df$var)) #
  colnames(newvals) = c("var", "ale")
  out.pred = left_join(newvals, Humid.df)
  out.pred.smooth = smooth.spline(out.pred$DoY, out.pred$ale, spar = .7)
  return(out.pred.smooth)
}


