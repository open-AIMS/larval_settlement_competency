## ---- loadPackages
library(tidyverse)
#library(glmmTMB)
library(DHARMa)
library(emmeans)
library(brms)
library(ggeffects)
library(tidybayes)
library(rstan)
library(patchwork)
## ----end

## ---- paths
DATA_PATH <<- "../data/"
OUTPUT_PATH <<- "../output/"
if (!dir.exists(DATA_PATH)) dir.create(DATA_PATH)
if (!dir.exists(paste0(DATA_PATH, "processed"))) dir.create(paste0(DATA_PATH, "processed"))
if (!dir.exists(paste0(DATA_PATH, "modelled"))) dir.create(paste0(DATA_PATH, "modelled"))

if (!dir.exists(OUTPUT_PATH)) dir.create(OUTPUT_PATH)
if (!dir.exists(paste0(OUTPUT_PATH, "figures"))) dir.create(paste0(OUTPUT_PATH, "figures"))
if (!dir.exists(paste0(OUTPUT_PATH, "tables"))) dir.create(paste0(OUTPUT_PATH, "tables"))
## ----end


## ----functionEDA
## The following function generates two EDA plots:
## 1. the raw data against larval age
## 2. the settlement (binary) data avainst larval age
EDA <- function(dat) {
    dat <- dat %>% filter(!is.na(Settlement)) %>% droplevels()
    g1<-dat %>% ggplot(aes(y = NoSet/(NoSet + NoNotSet), x = LarvalAge, colour = SpecificTreatment)) +
        geom_point(aes(shape = factor(Settlement)), alpha = 0.5) +
        geom_smooth(method = "gam",
                formula = y ~ s(x, bs = "cs"),
                method.args = list(family = 'binomial'),
                se = FALSE) +
        geom_line(data = dat,
                  stat = 'summary',
                  fun.y = 'mean',
                  aes(y = NoSet/(NoSet + NoNotSet),
                      x = LarvalAge,
                      colour = SpecificTreatment,
                      group = interaction(SpecificTreatment, factor(Plate))),
                  inherit.aes = FALSE,
                  alpha = 0.5) + 
        scale_y_continuous('Settlement (%)', label = scales::label_percent()) +
        scale_x_continuous('Larval age (d)') +
        scale_colour_discrete('Treatment') +
        scale_shape_discrete('Above threshold\nsettlement (cumm.)', breaks = c(0,1), labels = c('FALSE','TRUE')) +
        theme_classic()

    g2<-dat %>% ggplot(aes(y = Settlement, x = LarvalAge, colour = SpecificTreatment)) +
        geom_point(aes(shape = factor(Settlement))) +
        geom_smooth(method = 'glm',
                    method.args = list(family = 'binomial'),
                    se = FALSE) +
        geom_line(data = dat,
                  stat = 'summary',
                  fun.y = 'mean',
                  aes(y = Settlement,
                      x = LarvalAge,
                      colour = SpecificTreatment,
                      group = interaction(SpecificTreatment, factor(Plate))),
                  inherit.aes = FALSE,
                  alpha = 0.5) + 
        theme_classic()
    g1/g2
}
## ----end

## ----functionFitModel
## The following function fits the logistic regression model
fitModel <- function(dat) {
    dat <- dat %>% dplyr::select(Settlement, SpecificTreatment, LarvalAge, Plate) %>%
        na.omit() %>%
        droplevels()
    ## p <- round(max(2.5/model.matrix(~SpecificTreatment*LarvalAge, dat)[,-1] %>% apply(2, sd)),2)
    p <- round(2.5/model.matrix(~SpecificTreatment*LarvalAge, dat)[,-1] %>% apply(2, sd),2)
    priors <- prior(normal(0, 2.5), class = 'Intercept')
    for (i in 1:length(p)) {
        priors <- priors +
            set_prior(paste0("normal(0, ", p[i],")"), class = "b", coef = names(p)[i])
            ## prior(normal(0, p[i]), class = 'b', coef = !sym(names(p)[i]))
        }
    ## prior(normal(0,p), class = 'b') +
    priors <- priors + prior(student_t(3,0,2.5), class = "sd")
    form <- bf(Settlement|trials(1) ~ SpecificTreatment * LarvalAge + (1|Plate),
               family = 'binomial')
    mod <- brm(form,
               prior = priors,
               data = dat,
               iter = 6000,
               warmup = 2000,
               chains = 3, cores = 3,
               thin = 10,
               control = list(adapt_delta = 0.99, max_treedepth = 20),
               backend = 'cmdstanr')
    mod
}
## ----end

## ----functionMCMCdiagnostics
## The following functions produce MCMC diagnostic plots
MCMCdiagnostics_trace <- function(mod) {
    vars <- mod %>% get_variables()
    pars <- vars %>% stringr::str_subset("^b_.*|^sd_.*")
    mod$fit %>% rstan::stan_trace(pars = pars) 
}
MCMCdiagnostics_ac <- function(mod) {
    vars <- mod %>% get_variables()
    pars <- vars %>% stringr::str_subset("^b_.*|^sd_.*")
    mod$fit %>% rstan::stan_ac(pars = pars) 
}
MCMCdiagnostics_rhat <- function(mod) {
    mod$fit %>% rstan::stan_rhat() 
}
MCMCdiagnostics_ess <- function(mod) {
    mod$fit %>% rstan::stan_ess() 
}
## ----end

## ----functionDHARMa
## The following function produces DHARMa based diagnostics
DHARMa <- function(mod) {
    preds <- mod %>% posterior_predict(ndraws = 250, summary = FALSE)
    mod.resids <- createDHARMa(
        simulatedResponse = t(preds),
        observedResponse = mod$data$Settlement,
        fittedPredictedResponse = apply(preds, 2, median),
        integerResponse = FALSE
    )
    wrap_elements(~plot(mod.resids))
}
## ----end

## ----functionSumTable
## The following function generates a summary table of the
## regression parameters
sumTable <- function(mod) {
    mod %>% summarise_draws() %>%
        filter(str_detect(variable, "^b_.*|^s_.*"))
}
## ----end

## ----functionPartials
## The following function calculates the cell means from the model
partials <- function(mod) {
    dat <- mod$data
    newdata <-  with(dat,
                     list(LarvalAge = modelr::seq_range(LarvalAge, n = 100),
                          SpecificTreatment = levels(SpecificTreatment)))
    
    mod %>% emmeans(~LarvalAge|SpecificTreatment, at = newdata, type = "response") %>%
        data.frame() 
}
## ----end

## ----functionPartialPlot
## The following function produces a simple partial plot
partialPlot <- function(pred) {
    pred %>%
        ggplot(aes(y = prob,
                   x = LarvalAge,
                   colour = SpecificTreatment,
                   fill = SpecificTreatment)) +
        geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.3, colour = NA) +
        geom_line() +
        theme_classic()
}
## ----end

## ----functionLD50
## The following function calculates the LD50 for each SpecificTreatment
LD50 <- function(mod) {
    ## mod %>% spread_draws(`b_.*`, regex = TRUE) %>% 
    ##     dplyr::select(-.chain, -.iteration, -.draw) %>%
    ##     as.matrix() -> coefs
    ## ## start by creating intercept and slope columns
    ## Xmat <- cbind("CCA_intercept" =     c(1,0,0,0,0,0,0,0),
    ##               "CCA_slope" =         c(0,0,0,0,1,0,0,0),
    ##               "control_intercept" = c(1,1,0,0,0,0,0,0),
    ##               "control_slope" =     c(0,0,0,0,1,1,0,0),
    ##               "disc_intercept" =    c(1,0,1,0,0,0,0,0),
    ##               "disc_slope" =        c(0,0,0,0,1,0,1,0),
    ##               ## "extract_intercept=c(1,0,0,0,0,0,0,0),
    ##               ## "extract_slope" =  c(0,0,0,0,1,0,0,0),
    ##              ## "peptide_intercept"=c(1,0,0,0,0,0,0,0),
    ##               ## "peptide_slope" =  c(0,0,0,0,1,0,0,0),
    ##               "rubble_intercept" =  c(1,0,0,1,0,0,0,0),
    ##               "rubble_slope" =      c(0,0,0,0,1,0,0,1)
    ##               )
    ## fit <- coefs %*% (Xmat)
    e1 <- mod %>% emmeans(~SpecificTreatment|LarvalAge, at = list(LarvalAge =0)) %>%
        tidy_draws() %>%
        rename_with(~str_replace(.,"(SpecificTreatment) (.*),.*","\\1 \\2"), starts_with("Specific"))
    e2 <- mod %>% emtrends(specs = 'SpecificTreatment', var = 'LarvalAge') %>%
        tidy_draws() %>% 
        rename_with(~str_replace(.,"SpecificTreatment","Slope"), starts_with("Specific"))
    ee <- e1 %>% full_join(e2)
    ee %>% as.data.frame() %>%
        mutate(.draw = 1:n()) %>%
        pivot_longer(cols = c(-.chain,-.iteration,-.draw)) %>%
        separate(col = name, into = c("Parameter", "Variable"), sep=" ") %>%
        group_by(.draw, Variable) %>%
        summarise(LD50 =-1*value[1]/value[2]) #%>%
        ## ungroup() %>%
        ## group_by(Variable) %>%
        ## median_hdci()
} 
## ----end

## ---- functionLD50dif
LD50dif <- function(d) {
    d <- d %>% ungroup 
    coefs <- d %>%
        pivot_wider(names_from = "Variable", values_from = "LD50") %>%
        dplyr::select(-.draw) %>%
        as.matrix()
    varnames <- d %>%
        dplyr::select(Variable) %>%
        distinct() %>%
        pull(Variable)
    k <- varnames %>% length()
    pairwiseTests(coefs, exponentiate = FALSE, k=k, varnames = varnames, comp.value = 0)
}
## ----end

## ---- function pairwise tests
pairwiseTests <- function(coefs, exponentiate = FALSE, k, varnames, comp.value = 0) {
    Xmat <- NULL
    rnames <- NULL
    
    kindx <- 1:k
    for (i in 1:(k-1)) {
        for (j in (i+1):k) {
            Xmat <- rbind(Xmat,
                          as.numeric(kindx == j) - as.numeric(kindx == i))
            rnames <- c(rnames, paste(varnames[j], "-", varnames[i])) 
        }
    }
    x <- coefs %*% t(Xmat)
    if (exponentiate) x <- x %>% exp()
    colnames(x) <- rnames
    x %>% as.data.frame() %>%
        summarise_draws(median,
                        HDInterval::hdi,
                        `P>0` = ~ sum(.x > comp.value)/length(.x),
                        `P<0` = ~ sum(.x < comp.value)/length(.x))
}
## ----end

## ---- functionoddsQ1
oddsRatio <- function(mod) {
   mod %>% emtrends(specs = 'SpecificTreatment', var = 'LarvalAge') %>%
        tidy_draws() %>% 
        rename_with(~str_replace(.,"SpecificTreatment","Slope"), starts_with("Specific")) %>%
        mutate(across(c(-.chain,-.iteration,-.draw), exp)) 
}
oddsRatioTab <- function(odds) {
    odds %>%
        summarise_draws(median, HDInterval::hdi)
}

oddsRatiodif <- function(odds) {
    odds <- odds %>% ungroup 
    coefs <- odds %>%
      dplyr::select(-.chain,-.iteration,-.draw) %>%
      as.matrix() %>%
      log()
    varnames <- colnames(coefs) 
    k <- varnames %>% length()
    pairwiseTests(coefs, exponentiate = TRUE, k=k, varnames=varnames, comp.value = 1)
}
## ----end
