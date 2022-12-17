source('functions.R')

## Question 2: cue comparison within species
{
    ## ---- Analysis loop
    { 
        ## ---- loadData Q2
        load(file = paste0(DATA_PATH, "processed/data.q1.RData"))
        ## ----end
        
        ## ---- fitModelQ2
        dat <- data.q1[1,'data.sub'][[1]][[1]]
        data.q2.mod <- data.q1 %>%
            mutate(Mod = purrr::map(.x = data.sub, .f = fitModelQ2))
        ## ----end

    }
    ## ----end
}

fitModelQ2 <- function(dat) {
    dat <- dat %>% dplyr::select(NoSet, NoNotSet, SpecificTreatment, LarvalAge, Plate) %>%
        na.omit() %>%
        droplevels()

    form <- bf(NoSet|trials(NoSet + NoNotSet) ~ s(LarvalAge, by = SpecificTreatment, k = 4) + (1|Plate),
                   family = "binomial")
    mod <- brm(form,
                data = dat,
                ## iter = 6000,
                ## warmup = 2000,
                iter = 1000,
                warmup = 200,
                chains = 3, cores = 3,
                thin = 2,
                control = list(adapt_delta = 0.99, max_treedepth = 20),
                backend = 'cmdstanr')

    summary(mod)
    plot(conditional_smooths(mod, effects = 'LarvalAge:SpecificTreatment'), ask = FALSE) %>% wrap_plots()
    dev.off()
    plot(conditional_effects(mod, effects = 'LarvalAge:SpecificTreatment'), ask = FALSE) %>% wrap_plots()
    dev.off()
    plot(conditional_effects(mod, effects = 'LarvalAge:SpecificTreatment', method = 'posterior_predict', transform = NULL), ask = FALSE) %>% wrap_plots()
    dev.off()

    plot(conditional_smooths(mod, transform = plogit), ask = FALSE, plot = FALSE) -> a
    conditional_smooths(mod) -> a
    plot(conditional_effects(mod, transform = plogit), ask = FALSE, plot = FALSE) -> a
    conditional_smooths(mod) -> a

    a <- emmeans(mod, ~ LarvalAge|SpecificTreatment,
            at = list(LarvalAge = seq(0, 0.9, by = 0.01)),
            type = 'response') %>% as.data.frame()
    a %>% ggplot(aes(y = prob, x = LarvalAge)) +
        geom_line() +
        facet_wrap(~SpecificTreatment)
    dev.off()
    
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
