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
            slice(1) %>% 
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
                                          .f = function(x) {
                                              x %>%
                                                  ## areas %>%
                                                  ggplot(aes(colour = SpecificTreatment,
                                                             fill = SpecificTreatment,
                                                             x = Area)) +
                                                  stat_halfeye(alpha = 0.3) +
                                                  scale_x_continuous('Area under curve') +
                                                  scale_fill_discrete('Inducer') +
                                                  scale_colour_discrete('Inducer') +
                                                  theme_classic() +
                                                  ## scale_y_continuous('',expand = c(0,0)) +
                                                  theme(axis.text.y = element_blank(),
                                                        axis.ticks.y = element_blank(),
                                                        axis.title.y = element_blank())
                                          }
                                          )
                   )
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__.png"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialArea_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$Area_plot, ggsave, width = 6, height = 4, dpi = 300)
        ## ----end
        ## ---- compare_Areas Q2
        data.q2.mod <- data.q2.mod %>%
            mutate(Areas_diff = purrr::map(.x = AreaUnderCurve,
                                           .f = function(x) {
                                               x %>%
                                                   compare_levels(variable = Area,
                                                                  ## comparison = 'pairwise',
                                                                  comparison = list(
                                                                      c('CCA','control'),
                                                                      c('CCA','disc'),
                                                                      c('CCA','rubble'),
                                                                      c('rubble','disc'),
                                                                      c('rubble','control'),
                                                                      c('disc','control')),
                                                                  by = SpecificTreatment) 
                                           }),
                   Areas_diff_sums = purrr::map(.x = Areas_diff,
                                                .f = function(x) {
                                                   x %>% ungroup() %>%
                                                       mutate(SpecificTreatment = factor(SpecificTreatment,
                                                                                         levels = c('control', 'disc', 'rubble', 'CCA'))) %>%
                                                       group_by(SpecificTreatment) %>%
                                                       summarise_draws(median,
                                                                       HDInterval::hdi,
                                                                       `P<0`=function(x) sum(x<0)/length(x),
                                                                       `P>0`=function(x) sum(x>0)/length(x))

                                                }))
        
        data.q2.mod <- data.q2.mod %>%
                   mutate(Areas_diff_plot = purrr::map(.x = Areas_diff,
                                                .f = function(x) {
                                                    x %>%
                                                        mutate(SpecificTreatment = factor(SpecificTreatment,
                                                                                          levels = c('disc - control',
                                                                                                     'rubble - control',
                                                                                                     'CCA - control',
                                                                                                     'CCA - disc',
                                                                                                     'rubble - disc',
                                                                                                     'CCA - rubble'))) %>%
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
                                                    })
                   )
                   
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__.pdf"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__.png"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialAreaDiffs_",data.q2.mod$Species,"__large_.png"),
             data.q2.mod$Areas_diff_plot, ggsave, width = 6, height = 4, dpi = 300)
        ## ----end
        ## ---- save Model Validations Q2
        data.q2.mod <- data.q2.mod %>%
            dplyr::select(-Mod)
        save(data.q2.mod,
             file = paste0("../data/modelled/data.q2",".mod_areas.RData"))
        ## ----end
        
    }
    ## ----end
}
