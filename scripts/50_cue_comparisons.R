source('functions.R')

## Question 2: cue comparison within species
{
    ## ---- Analysis loop
    { 
        ## ---- loadData Q2
        load(file = paste0(DATA_PATH, "processed/data.q1.RData"))
        ## ----end
        
        ## ---- fitModel Q2
        ## dat <- data.q1[1,'data.sub'][[1]][[1]]
        data.q2.mod <- data.q1 %>%
            ungroup() %>%
            mutate(Mod = purrr::map(.x = data.sub, .f = fitModelQ2))
        ## ----end
        ## ---- save Models Q2
        save(data.q2.mod, file = paste0("../data/modelled/data.q2.mod.RData"))
        ## ----end

        ## ---- MCMCdiagnostics Q2
        data.q2.mod <- data.q2.mod %>%
            mutate(Trace = purrr::map(.x = Mod, .f = MCMCdiagnostics_trace),
                   AC  = purrr::map(.x = Mod, .f = MCMCdiagnostics_ac),
                   Rhat  = purrr::map(.x = Mod, .f = MCMCdiagnostics_rhat),
                   ess  = purrr::map(.x = Mod, .f = MCMCdiagnostics_ess)
                   )
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q2.mod$Species,"_.pdf"),
             data.q2.mod$Trace, ggsave)
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q2.mod$Species,"_.png"),
             data.q2.mod$Trace, ggsave, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q2.mod$Species,"_large_.png"),
             data.q2.mod$Trace, ggsave, dpi = 300)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q2.mod$Species,"_.pdf"),
             data.q2.mod$AC, ggsave)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q2.mod$Species,"_.png"),
             data.q2.mod$AC, ggsave, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q2.mod$Species,"_large_.png"),
             data.q2.mod$AC, ggsave, dpi = 300)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q2.mod$Species,"_.pdf"),
             data.q2.mod$Rhat, ggsave, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q2.mod$Species,"_.png"),
             data.q2.mod$Rhat, ggsave, dpi = 72, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q2.mod$Species,"_large_.png"),
             data.q2.mod$Rhat, ggsave, dpi = 300, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q2.mod$Species,"_.pdf"),
             data.q2.mod$ess, ggsave, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q2.mod$Species,"_.png"),
             data.q2.mod$ess, ggsave, dpi = 72, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q2.mod$Species,"_large_.png"),
             data.q2.mod$ess, ggsave, dpi = 300, width = 6, height = 3)
        ## ----end
        ## ---- DHARMa Q2
        data.q2.mod <- data.q2.mod %>%
            mutate(DHARMa = purrr::map(.x = Mod, .f = DHARMaQ2))
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q2.mod$Species,"_.pdf"),
             data.q2.mod$DHARMa, ggsave, width = 12, height = 6)
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q2.mod$Species,"_.png"),
             data.q2.mod$DHARMa, ggsave, width = 12, height = 6, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q2.mod$Species,"_large_.png"),
             data.q2.mod$DHARMa, ggsave, width = 12, height = 6, dpi = 300)
        ## ----end
        ## ---- save Model Validations Q2
        data.q2.mod <- data.q2.mod %>%
            dplyr::select(-Mod)
        save(data.q2.mod,
             file = paste0("../data/modelled/data.q2",".mod_valid.RData"))
        ## ----end
        ## ---- partialsQ2
        load(file = paste0("../data/modelled/data.q2.mod.RData"))
        data.q2.mod <- data.q2.mod %>%
            mutate(Partials = purrr::map(.x = Mod, .f = partials))
        ## ----end
        ## ---- partialPlotQ2
        data.q2.mod <- data.q2.mod %>%
            mutate(
                PartialPlot = purrr::map(.x = Partials,
                                            .f = ~ partialPlot(.x, T = NULL)))
        save(data.q2.mod, file = paste0('../data/modelled/data.q2__.partials.RData'))
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$PartialPlot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q2.mod$Species,"__.png"),
             data.q2.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 300)
        ## ----end
        ## ---- areaQ2
        load(file = paste0('../data/modelled/data.q2__.partials.RData'))
        data.q2.mod <- data.q2.mod %>%
            mutate(AreaUnderCurve = purrr::map(.x = Mod,
                                               .f = ~ AUC(.x)),
                   Area_sums = purrr::map(.x = AreaUnderCurve,
                                          .f = ~ .x %>%
                                              ungroup() %>%
                                              group_by(SpecificTreatment) %>%
                                              summarise_draws(median,
                                                              HDInterval::hdi)
                                          ),
                   Area_plot = purrr::map(.x = AreaUnderCurve,
                                          .f = ~ AUCplot(.x))
                   )
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__.png"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4, dpi = 300)
        save(data.q2.mod, file = paste0('../data/modelled/data.q2__.area.RData'))
        ## ----end
        ## x <- data.q2.mod[2,'AreaUnderCurve'][[1]][[1]]
        ## ---- compare_Areas Q2
        load(file = paste0('../data/modelled/data.q2__.area.RData'))
        data.q2.mod <- data.q2.mod %>%
            ## ungroup() %>% slice(2) %>%
            mutate(Areas_diff = purrr::map(.x = AreaUnderCurve,
                                           .f = ~ AUC_diff(.x) 
                                           ),
                   Areas_ratio = purrr::map(.x = AreaUnderCurve,
                                            .f = ~ AUC_ratio(.x))
                                           )
        save(data.q2.mod, file = paste0("../data/modelled/data.q2",".mod_areas_diff.RData"))
        ## ----end
        ## ---- compare_Areas_sums Q2
        load(file = paste0("../data/modelled/data.q2",".mod_areas_diff.RData"))
        data.q2.mod <- data.q2.mod %>%
            mutate(Areas_diff_sums = purrr::map(.x = Areas_diff,
                                                .f = ~ AUC_diff_sums(.x)),
                   Areas_ratio_sums = purrr::map(.x = Areas_ratio,
                                                 .f = ~ AUC_ratio_sums(.x))
                   )
        save(data.q2.mod, file = paste0("../data/modelled/data.q2",".mod_areas_sums.RData"))
        ## ----end
        ## ---- compare_Areas_plots Q2
        load(file = paste0("../data/modelled/data.q2",".mod_areas_sums.RData"))
        data.q2.mod <- data.q2.mod %>%
            mutate(Areas_diff_plot = purrr::map(.x = Areas_diff,
                                                .f = ~ AUC_diff_plot(.x)),
                   Areas_ratio_plot = purrr::map(.x = Areas_ratio,
                                                 .f = ~ AUC_ratio_plot(.x))
                   )
        
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__.png"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4, dpi = 300)
        
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaRatios_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$Areas_ratio_plot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaRatios_",data.q2.mod$Species,"__.png"),
             data.q2.mod$Areas_ratio_plot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaRatios_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$Areas_ratio_plot, ggsave, width = 6, height = 4, dpi = 300)
        save(data.q2.mod, file = paste0("../data/modelled/data.q2",".mod_areas_plots.RData"))
        ## ----end
        
    }
    ## ----end
}
