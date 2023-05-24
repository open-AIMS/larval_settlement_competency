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
        ## geom_line(data = dat,
        ##           stat = 'summary',
        ##           fun.y = 'mean',
        ##           aes(y = NoSet/(NoSet + NoNotSet),
        ##               x = LarvalAge,
        ##               colour = SpecificTreatment,
        ##               group = interaction(SpecificTreatment, factor(Plate))),
        ##           inherit.aes = FALSE,
        ##           alpha = 0.5) + 
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
        ## geom_line(data = dat,
        ##           stat = 'summary',
        ##           fun.y = 'mean',
        ##           aes(y = Settlement,
        ##               x = LarvalAge,
        ##               colour = SpecificTreatment,
        ##               group = interaction(SpecificTreatment, factor(Plate))),
        ##           inherit.aes = FALSE,
        ##           alpha = 0.5) + 
        scale_x_continuous('Larval age (d)') +
        scale_colour_discrete('Treatment') +
        scale_shape_discrete('Above threshold\nsettlement (cumm.)', breaks = c(0,1), labels = c('FALSE','TRUE')) +
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
    pars <- vars %>% stringr::str_subset("^b_.*|^bs_.*|^sd_.*")
    mod$fit %>% rstan::stan_trace(pars = pars) 
}
MCMCdiagnostics_ac <- function(mod) {
    vars <- mod %>% get_variables()
    pars <- vars %>% stringr::str_subset("^b_.*|^bs_.*|^sd_.*")
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
    wrap_elements(~testUniformity(mod.resids)) +
        wrap_elements(~plotResiduals(mod.resids)) +
        wrap_elements(~testDispersion(mod.resids))
}
## ----end

## ----functionDHARMaQ2
## The following function produces DHARMa based diagnostics
DHARMaQ2 <- function(mod) {
    preds <- mod %>% posterior_predict(ndraws = 250, summary = FALSE)
    mod.resids <- createDHARMa(
        simulatedResponse = t(preds),
        observedResponse = mod$data$NoSet,
        fittedPredictedResponse = apply(preds, 2, median),
        integerResponse = FALSE
    )
    wrap_elements(~testUniformity(mod.resids)) +
        wrap_elements(~plotResiduals(mod.resids)) +
        wrap_elements(~testDispersion(mod.resids))
}
## ----end
## ----functionSumTable
## The following function generates a summary table of the
## regression parameters
sumTable <- function(mod) {
    mod %>% summarise_draws() %>%
        filter(str_detect(variable, "^b_.*|^bs_.*|^s_.*"))
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


## ----functionPartialsQ2
## The following function calculates the cell means from the model
## and retains all draws
partialsQ2 <- function(mod) {
    dat <- mod$data
    newdata <-  with(dat,
                     list(LarvalAge = modelr::seq_range(LarvalAge, n = 100),
                          SpecificTreatment = levels(SpecificTreatment)))
    
    mod %>% emmeans(~LarvalAge|SpecificTreatment, at = newdata, type = "link") %>%
        ## contrast(by = "LarvalAge",
        ##          method = 'pairwise') %>%
        gather_emmeans_draws()
}
## ----end

## ----functionAUCQ2
## this function calculates the areas under curves
AUC <- function(mod) {
    partials <- mod %>% partialsQ2()
    partials %>%
        mutate(.value = plogis(.value),
               T = 1) %>%
        ungroup() %>%
        group_by(.draw, SpecificTreatment) %>%
        summarise(Area = sum(.value),
                  Total = sum(T))
}
## ----end

## ---- functionAUCPlot
AUCplot <- function(x, spec_treat_levels) {
    x %>%
        ggplot(aes(colour = SpecificTreatment,
                   fill = SpecificTreatment,
                   x = Area)) +
        stat_halfeye(alpha = 0.3, normalize = "groups") +
        scale_x_continuous('Area under curve') +
        scale_fill_discrete('Inducer', limits = spec_treat_levels) +
        scale_colour_discrete('Inducer', limits = spec_treat_levels) +
        theme_classic() +
        ## scale_y_continuous('',expand = c(0,0)) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank())
}
## ----end

## ---- functionAUC_diff
AUC_diff <- function(x, v = c('CCA','rubble','disc','control')) {
    ## fix this
    ## want to work out a way to control the comparisons,
    ## but dont always have all the of the SpecificTreatments
    
    xx <- x %>% pull(SpecificTreatment) %>% unique 
    vv <- v[str_which(string = v, paste(xx, collapse='|'))]
    vs <- combn(vv, 2, simp = FALSE)
    x %>%
        compare_levels(variable = Area,
                       ## comparison = 'pairwise',
                       ## comparison = list(
                       ##     c('CCA','control'),
                       ##     c('CCA','disc'),
                       ##     c('CCA','rubble'),
                       ##     c('rubble','disc'),
                       ##     c('rubble','control'),
                       ##     c('disc','control')),
                       comparison = vs,
                       by = SpecificTreatment) 
}
## ----end

## ---- functionAUC_ratio
AUC_ratio <- function(x, v = c('CCA','rubble','disc','control')) {
    ## fix this
    ## want to work out a way to control the comparisons,
    ## but dont always have all the of the SpecificTreatments
    
    xx <- x %>% pull(SpecificTreatment) %>% unique 
    vv <- v[str_which(string = v, paste(xx, collapse='|'))]
    vs <- combn(vv, 2, simp = FALSE)
    x %>%
        mutate(lArea = log(Area)) %>%
        compare_levels(variable = lArea,
                       ## comparison = 'pairwise',
                       ## comparison = list(
                       ##     c('CCA','control'),
                       ##     c('CCA','disc'),
                       ##     c('CCA','rubble'),
                       ##     c('rubble','disc'),
                       ##     c('rubble','control'),
                       ##     c('disc','control')),
                       comparison = vs,
                       by = SpecificTreatment) 
}
## ----end
## ---- functionAUC_diff_sums
AUC_diff_sums <- function(x) {
    x %>% ungroup() %>%
        ## mutate(SpecificTreatment = factor(SpecificTreatment,
        ##                                   levels = c('control', 'disc', 'rubble', 'CCA'))) %>%
        group_by(SpecificTreatment) %>%
        summarise_draws(median,
                        HDInterval::hdi,
                        `P<0`=function(x) sum(x<0)/length(x),
                        `P>0`=function(x) sum(x>0)/length(x))
}
## ----end

## ---- functionAUC_ratio_sums
AUC_ratio_sums <- function(x) {
    x %>% ungroup() %>%
        mutate(Area = exp(lArea)) %>%
        group_by(SpecificTreatment) %>%
        summarise_draws(median,
                        HDInterval::hdi,
                        `P<1`=function(x) sum(x<1)/length(x),
                        `P>1`=function(x) sum(x>1)/length(x))
}
## ----end

## ---- functionAUC_diff_plot
AUC_diff_plot <- function(x) {
    x %>%
        ## mutate(SpecificTreatment = factor(SpecificTreatment,
        ##                                   levels = c('disc - control',
        ##                                              'rubble - control',
        ##                                              'CCA - control',
        ##                                              'CCA - disc',
        ##                                              'rubble - disc',
        ##                                              'CCA - rubble'))) %>%
        ggplot() +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        stat_halfeye(aes(y = SpecificTreatment, x = Area), alpha = 0.4,
                     fill_type = 'gradient') +
        ## geom_density() +
        ## facet_grid(SpecificTreatment~., scales = 'free')
        ## geom_density_ridges(aes(y = SpecificTreatment, x = Area),
        ##                     alpha = 0.4) +
        scale_x_continuous('Difference in area under curve') +
        theme_classic() +
        theme(axis.title.y = element_blank())
}
## ----end

## ---- functionAUC_ratio_plot
AUC_ratio_plot <- function(x) {
    x %>%
        ggplot() +
        geom_vline(xintercept = 1, linetype = 'dashed') +
        stat_halfeye(aes(y = SpecificTreatment, x = exp(lArea)), alpha = 0.4,
                     fill_type = 'gradient') +
        scale_x_continuous('Fractional difference in area under curve',
                           trans = scales::log2_trans(),
                           breaks = scales::breaks_log(base = 2)) +
        theme_classic() +
        theme(axis.title.y = element_blank())
}
## ----end

## ----functionPartialPlot
## The following function produces a simple partial plot
partialPlot <- function(pred, species, T=NULL, limits = c('CCA','control','disc', 'rubble')) {
    load(paste0(DATA_PATH, "processed/Species abbreviations.RData"))
    SP <- speciesAbbrev %>% filter(Abbreviation == species) %>% pull(Species) %>% unique()

    pred %>%
        ggplot(aes(y = prob,
                   x = LarvalAge,
                   colour = SpecificTreatment,
                   fill = SpecificTreatment)) +
        geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.3, colour = NA) +
        scale_fill_discrete('Inducer',limits = limits) +
        scale_colour_discrete('Inducer',limits = limits) +
        scale_x_continuous('Larval age (days)') +
        {if(!is.null(T)) scale_y_continuous(str_wrap(paste0('Cohort settlement prob. (P>',T,')'),25)) } + 
        {if(is.null(T)) scale_y_continuous(str_wrap(paste0('Cohort settlement prob'),25)) } + 
        geom_line() +
        ggtitle(SP) +
        theme_classic()
}
## ----end

## ----functionCompilation
partial_plot_compilations <- function(path, g, ncol = 3, dpi = 100) {
     gw <- g %>% 
         suppressMessages() %>%
         suppressWarnings() %>%
         wrap_plots() + guide_area() + plot_layout(guides = 'collect', ncol = ncol) &
         guides(
             fill = "none",
             colour = guide_legend(override.aes = list(shape = NA, size = 0.7))) 
     n_patches <- length(gw$patches$plots) + 1
     dims <- wrap_dims(n_patches, ncol = ncol, nrow = NULL)
     ggsave(path, gw, width = 4*dims[2], height = 3*dims[1], dpi = dpi)
}
## ----end

## ---- functionCompilationNew
partial_plot_compilations_new <- function(path, dat.mod, ncol = 3, dpi = 100,
                                          legend.position = "bottom",
                                          legend.direction = "horizontal",
                                          legend.justification = c(0.5, 0.5)) {
    load(paste0(DATA_PATH, "processed/Species abbreviations.RData"))
    
    data.compilation <- dat.mod %>%
        dplyr::select(Species, Partials) %>%
        unnest(Partials) %>%
        ungroup() %>%
        left_join(data.q1.mod %>%
                  dplyr::select(LD50tab) %>%
                  unnest(LD50tab) %>%
                  dplyr::select(SpecificTreatment = Variable,
                                LD50)) %>%
        left_join(speciesAbbrev %>%
                  dplyr::select(Species1 = Species,
                                Species = Abbreviation)) %>%
        dplyr::select(-Species) %>%
        dplyr::rename(Species = Species1)
    ld50 <- data.compilation %>%
        group_by(Species, SpecificTreatment) %>%
        summarise(LD50 = mean(LD50)) %>%
        ungroup() %>%
        group_by(Species) %>%
        summarise(SpecificTreatment = SpecificTreatment[which.min(LD50)],
                  LD50 = min(LD50))
    
        
    
    limits = c('CCA','control','disc', 'rubble')
    g <- data.compilation %>%
        ggplot(aes(y = prob,
                   x = LarvalAge,
                   colour = SpecificTreatment,
                   fill = SpecificTreatment)) +
        geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.3, colour = NA) +
        scale_fill_discrete('Inducer',limits = limits) +
        scale_colour_discrete('Inducer',limits = limits) +
        scale_x_continuous('Larval age (days)') +
        scale_y_continuous(paste0('Cohort settlement probability (P>',thresholdProp,')')) + 
        geom_line() +
        geom_vline(data = ld50, aes(xintercept = LD50, colour = SpecificTreatment),
                   linetype = 'dashed', show.legend = FALSE) + 
        facet_wrap(~Species, labeller = label_bquote(col = italic(.(Species))),
                   scales = 'free_x', ncol = ncol) +
        theme_classic() +
        theme(
            legend.position = legend.position,
            legend.direction = legend.direction,
            legend.justification = legend.justification,
            axis.title.y = element_text(margin = margin(0, r = 10, unit = "pt")),
            axis.title.x = element_text(margin = margin(0, t = 10, unit = "pt")),
            strip.background = element_blank()) +
        guides(
            fill = "none",
            colour = guide_legend(override.aes = list(shape = NA, size = 0.7))) 
    ## n_patches <- length(g$patches$plots) + 1
    n_panels <- length(unique(ggplot_build(g)$data[[1]]$PANEL))
    dims <- wrap_dims(n_panels, ncol = ncol, nrow = NULL)
    ggsave(path, g, width = 2*dims[2], height = 1.5*dims[1], dpi = dpi)
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
        summarise(LD50 =-1*value[1]/value[2]) %>%
        mutate(LD50 = ifelse(LD50<0, NA, LD50))                                #%>%
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
        summarise_draws(median = ~median(.x, na.rm = TRUE),
                        HDInterval::hdi,
                        ## `P>0` = ~ sum(.x > comp.value)/length(.x),
                        ## `P<0` = ~ sum(.x < comp.value)/length(.x))
                        `P>0` = ~ mean(.x > comp.value, na.rm = TRUE),
                        `P<0` = ~ mean(.x < comp.value, na.rm = TRUE))
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

## ---- functionfitModelQ2
fitModelQ2 <- function(dat) {
    dat <- dat %>%
        mutate(AgePlate = factor(paste(LarvalAge, Plate))) %>%
        dplyr::select(NoSet, NoNotSet, SpecificTreatment, LarvalAge, AgePlate) %>%
        na.omit() %>%
        droplevels()

    form <- bf(NoSet|trials(NoSet + NoNotSet) ~ s(LarvalAge, by = SpecificTreatment, k = 5) + (1|AgePlate),
                   family = "binomial")
    mod <- brm(form,
                data = dat,
                iter = 6000,
                warmup = 2000,
                chains = 3, cores = 3,
                thin = 2,
                control = list(adapt_delta = 0.99, max_treedepth = 20),
                backend = 'cmdstanr')
    mod
}
## ----end

## ---- functionCompilationQ2Area
partial_plot_compilations_Q2_Area <- function(path, dat.mod, ncol = 3, dpi = 100,
                                          legend.position = "bottom",
                                          legend.direction = "horizontal",
                                          legend.justification = c(0.5, 0.5)) {
    load(paste0(DATA_PATH, "processed/Species abbreviations.RData"))
    
    data.compilation <- dat.mod %>%
        dplyr::select(Species, Partials) %>%
        unnest(Partials) %>%
        ungroup() %>%
        left_join(speciesAbbrev %>%
                  dplyr::select(Species1 = Species,
                                Species = Abbreviation)) %>%
        dplyr::select(-Species) %>%
        dplyr::rename(Species = Species1)
    
    ## limits = c('CCA','control','disc', 'rubble')
    limits = c('control', 'CCA', 'disc', 'extract', 'peptide', 'rubble')

    g <- data.compilation %>%
        ggplot(aes(y = prob,
                   x = LarvalAge,
                   colour = SpecificTreatment,
                   fill = SpecificTreatment)) +
        geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.3, colour = NA) +
        geom_line() +
        scale_x_continuous('Larval age (days)') +
        scale_y_continuous('Cohort settlement probability') +
        scale_fill_discrete('Inducer', limits = limits) +
        scale_colour_discrete('Inducer', limits = limits) +
        facet_wrap(~Species, labeller = label_bquote(col = italic(.(Species))),
                   scales = 'free_x', ncol = ncol) +
        theme_classic() +
        ## scale_y_continuous('',expand = c(0,0)) +
        theme(legend.position = legend.position,
              legend.direction = legend.direction,
              legend.justification = legend.justification,
              axis.title.y = element_text(margin = margin(0, r = 10, unit = "pt")),
              axis.title.x = element_text(margin = margin(0, t = 10, unit = "pt")),
              strip.background = element_blank()) +
        guides(
            fill = "none",
            colour = guide_legend(override.aes = list(shape = NA, size = 0.7))) 
    ## n_patches <- length(g$patches$plots) + 1
    n_panels <- length(unique(ggplot_build(g)$data[[1]]$PANEL))
    dims <- wrap_dims(n_panels, ncol = ncol, nrow = NULL)
    if (ncol == 5) {
        ggsave(path, g, width = 2*dims[2], height = (1.5*dims[1]) + 2, dpi = dpi)
    } else {
        ggsave(path, g, width = 2*dims[2], height = 1.5*dims[1], dpi = dpi)
    }
}
## ----end

## ---- functionCompilationQ2AreaPosteriors
partial_plot_compilations_Q2_Area_posteriors <- function(path, dat.mod, ncol = 3, dpi = 100,
                                          legend.position = "bottom",
                                          legend.direction = "horizontal",
                                          legend.justification = c(0.5, 0.5)) {
    load(paste0(DATA_PATH, "processed/Species abbreviations.RData"))
    
    data.compilation <- dat.mod %>%
        dplyr::select(Species, AreaUnderCurve) %>%
        unnest(AreaUnderCurve) %>%
        ungroup() %>%
        left_join(speciesAbbrev %>%
                  dplyr::select(Species1 = Species,
                                Species = Abbreviation)) %>%
        dplyr::select(-Species) %>%
        dplyr::rename(Species = Species1)
    
    ## limits = c('CCA','control','disc', 'rubble')
    limits = c('control', 'CCA', 'disc', 'extract', 'peptide', 'rubble')

    g <- data.compilation %>%
        ggplot(aes(x = Area,
                   colour = SpecificTreatment,
                   fill = SpecificTreatment)) +
        stat_halfeye(alpha = 0.3, normalize = "groups") +
        scale_x_continuous('Area under curve') +
        scale_fill_discrete('Inducer', limits = limits) +
        scale_colour_discrete('Inducer', limits = limits) +
        facet_wrap(~Species, labeller = label_bquote(col = italic(.(Species))),
                   scales = 'free_x', ncol = ncol) +
        theme_classic() +
        ## scale_y_continuous('',expand = c(0,0)) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = legend.position,
            legend.direction = legend.direction,
            legend.justification = legend.justification,
            axis.title.x = element_text(margin = margin(0, t = 10, unit = "pt")),
            strip.background = element_blank()) +
        guides(
            ## fill = "none",
            colour = guide_legend(override.aes = list(shape = NA, size = 0.7))) 
    ## n_patches <- length(g$patches$plots) + 1
    n_panels <- length(unique(ggplot_build(g)$data[[1]]$PANEL))
    dims <- wrap_dims(n_panels, ncol = ncol, nrow = NULL)
    if (ncol == 5) {
        ggsave(path, g, width = 2*dims[2], height = (1.5*dims[1]) + 2, dpi = dpi)
    } else {
        ggsave(path, g, width = 2*dims[2], height = 1.5*dims[1], dpi = dpi)
    }
}
## ----end
