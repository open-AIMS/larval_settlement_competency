source('functions.R')

## Question 1: age to settlement
{
    ## ---- Analysis loop
    for (thresholdProp in c(0.1, 0.3, 0.5)) {
        ## ---- load Models Summaries
        load(file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod.RData"))
        ## ----end

        ## ----sumTableQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(sumTable = purrr::map(.x = Mod, .f = sumTable))
        ## ----end
        ## ----partialPlotQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(Partials = purrr::map(.x = Mod, .f = partials))
        ## ----end
        ## ----partialPlotQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(PartialPlot = purrr::map(.x = Partials, .f = partialPlot))
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 300)
        ## ----end

        
        ## ----LD50Q1
        mod = data.q1.mod[1,'Mod'][[1]][[1]]
        data.q1.mod <- data.q1.mod %>%
            mutate(LD50 = purrr::map(.x = Mod, .f = LD50),
                   LD50tab = purrr::map(.x = LD50,
                                        .f = function(x) x %>% ungroup %>%
                                                         group_by(Variable) %>%
                                                         median_hdci()
                                        ),
                   LD50dif = purrr::map(.x = LD50,
                                        .f = LD50dif) 
                   )
        ## ----end

        ## ---- oddsQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(Odds = purrr::map(.x = Mod, .f = oddsRatio),
                   Oddstab = purrr::map(.x = Odds,
                                        .f = function(x) x %>% ungroup %>%
                                                         summarise_draws(median,
                                                                         HDInterval::hdi)
                                        ),
                   Oddsdif = purrr::map(.x = Odds,
                                        .f = oddsRatiodif)
                   )
        ## ----end


        ## ----save Model Summaries
        data.q1.mod <- data.q1.mod %>%
            dplyr::select(-Mod)
        save(data.q1.mod, file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod_sum.RData"))
        ## ----end

    }
    ## ----end
}

