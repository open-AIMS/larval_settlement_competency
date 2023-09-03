source('functions.R')

## Question 1: age to settlement
{
    ## ---- Analysis loop
    for (thresholdProp in seq(0.1, 0.9, by = 0.1)) {
        ## ---- load Models Summaries
        load(file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod.RData"))
        ## ----end

        ## ----sumTableQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(sumTable = purrr::map(.x = Mod, .f = sumTable))
        ## ----end
        ## ----partialsQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(Partials = purrr::map(.x = Mod, .f = partials))
        ## ----end
        ## ----partialPlotQ1
        ## data.q1.mod <- data.q1.mod %>%
        ##     mutate(PartialPlot = purrr::map(.x = Partials, .f = partialPlot))
        data.q1.mod <- data.q1.mod %>%
            mutate(PartialPlot = purrr::map2(.x = Partials, .y = Species,
                                             .f = ~ partialPlot(.x, .y, T = thresholdProp,
                                                                limits = c('CCA','disc', 'rubble'))))
        save(data.q1.mod, file = paste0('../data/modelled/data.q1__', thresholdProp, '.partials.RData'))
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/PartialPlot_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$PartialPlot, ggsave, width = 6, height = 4, dpi = 300)
        ## ----end

        
        ## ----LD50Q1
        load(file = paste0('../data/modelled/data.q1__', thresholdProp, '.partials.RData'))
        ## mod = data.q1.mod[1,'Mod'][[1]][[1]]
        data.q1.mod <- data.q1.mod %>%
            mutate(LD50 = purrr::map(.x = Mod, .f = LD50),
                   LD50tab = purrr::map(.x = LD50,
                                        .f = function(x) x %>% ungroup %>%
                                                         group_by(Variable) %>%
                                                         median_hdci(na.rm = TRUE)
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
        ## ---- sumTabsQ1
        sum.tabs.q1 <- data.q1.mod %>%
            dplyr::select(Species, LD50tab, LD50dif, Oddstab, Oddsdif)
        save(sum.tabs.q1, file = paste0(DATA_PATH, "modelled/sum.tabs.q1__", thresholdProp, ".RData"))
        ## ----end

        ## ---- compiliationNewQ1
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_",
                                              thresholdProp, "_.png"),
                                  dat.mod = data.q1.mod,
                                  ncol = 4, dpi = 72,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_",
                                              thresholdProp, "_.pdf"),
                                  dat.mod = data.q1.mod,
                                  ncol = 4, dpi = 72,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_",
                                              thresholdProp, "_large_.png"),
                                  dat.mod = data.q1.mod,
                                  ncol = 4, dpi = 300,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_3cols",
                                              thresholdProp, "_.png"),
                                  dat.mod = data.q1.mod,
                                  ncol = 3, dpi = 72,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_3cols",
                                              thresholdProp, "_.pdf"),
                                  dat.mod = data.q1.mod,
                                  ncol = 3, dpi = 72,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))
        partial_plot_compilations_new(path=paste0(OUTPUT_PATH,
                                              "figures/partialCompilation_3cols",
                                              thresholdProp, "_large_.png"),
                                  dat.mod = data.q1.mod,
                                  ncol = 3, dpi = 300,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('CCA','disc', 'rubble'))

        ## ----end 
        ## ---- compiliationQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(PartialPlot = purrr::map2(.x = PartialPlot,
                                             .y = Species,
                                             .f = ~ .x + ggtitle(.y)))
        partial_plot_compilations(path=paste0(OUTPUT_PATH,
                                              "figures/compilation_",
                                              thresholdProp, "_.png"),
                                  g=data.q1.mod$PartialPlot,
                                  ncol = 5, dpi = 72)
        partial_plot_compilations(path=paste0(OUTPUT_PATH,
                                              "figures/compilation_",
                                              thresholdProp, "_.pdf"),
                                  g=data.q1.mod$PartialPlot,
                                  ncol = 5)
        partial_plot_compilations(path=paste0(OUTPUT_PATH,
                                              "figures/compilation_",
                                              thresholdProp, "_large_.png"),
                                  g=data.q1.mod$PartialPlot,
                                  ncol = 5, dpi = 300)
        ## ----end

        ## ----save Model Summaries
        data.q1.mod <- data.q1.mod %>%
            dplyr::select(-Mod)
        save(data.q1.mod, file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod_sum.RData"))
        ## ----end

    }
    ## ----end
    
    ## ---- days_vs_thresholdQ1
    files <- list.files(path = paste0(DATA_PATH, "modelled"),
                        pattern = "data.q1.*mod_sum.RData",
                        full.names = TRUE)
    load(files[1])
    SPECIES <- data.q1.mod %>% pull(Species) %>% unique()
    species_list <- setNames(vector('list', length = length(SPECIES)),
                             SPECIES)
    thress <- seq(0.1, 0.9, by = 0.1)
    for (i in 1:length(SPECIES)) {
        print(paste0("Species = ", SPECIES[i]))
        thresh_list <- setNames(vector('list', length = length(thress)),
                                thress)
        for (thresholdProp in seq(0.1, 0.9, by = 0.1)) {
            load(paste0(DATA_PATH, "modelled/data.q1__",thresholdProp,".mod_sum.RData"))
            thresh_list[[as.character(thresholdProp)]] <- data.q1.mod[i,'LD50tab'][[1]][[1]] %>%
                mutate(Threshold = thresholdProp,
                       Species = SPECIES[i])
        }
        species_list[[SPECIES[i]]] <- do.call('rbind', thresh_list)
        maxDays <- data.q1.mod[i, 'data'][[1]][[1]] %>% pull(LarvalAge) %>% max()
        species_list[[SPECIES[i]]] <- species_list[[SPECIES[i]]] %>%
            mutate(MaxDays = maxDays)
    }
    all_LD50 <- do.call('rbind', species_list) %>%
        group_by(Species) %>%
        nest() %>%
        mutate(data = map(.x = data,
                          .f = ~ .x %>%
                              ## mutate(across(c(LD50, .lower, .upper),
                              ##              ~ pmin(.x, MaxDays))) %>%
                              filter(LD50<MaxDays) %>%
                              mutate(Flag = LD50<MaxDays) %>%
                              group_by(Threshold) %>%
                              bind_rows(group_by(., Threshold) %>%
                                        summarise(
                                            BestVar = Variable[which.min(LD50)],
                                            Variable = 'Best',
                                            .lower = .lower[which.min(LD50)],
                                            .upper = .upper[which.min(LD50)],
                                            LD50 = LD50[which.min(LD50)],
                                           )) %>%
                              arrange(Threshold, Variable) %>%
                              ungroup() %>%
                              mutate(MaxY = max(.upper[which(Flag)]))
                          )
               )
    
    load(file = paste0(DATA_PATH, "processed/TreatmentPalette.RData"))
    load(file = paste0(DATA_PATH, "processed/TreatmentOrder.RData"))
    treatment_palette <- c("Best" = "#000000", TreatmentPalette)
    treatment_palette <- treatment_palette[c("Best","CCA", "disc", "rubble")]

    all_LD50 <- all_LD50 %>%
        mutate(g = map(.x = data,
                       .f = ~ .x %>%
                           ## filter(!is.na(BestVar)) %>%
                           ggplot(aes(y = LD50, x = Threshold)) +
                           geom_blank(data = NULL, aes(y = 0, x = 0.1)) +
                           ## geom_point(data = .x %>% filter(Variable == 'control') %>% droplevels(),
                           ##            aes(colour = Variable), position = position_dodge()) +
                           geom_line(data = . %>% filter(Variable != 'control') %>% droplevels(),
                                     aes(colour = Variable), position = position_dodge(width = 0.05)) +
                           geom_pointrange(data = . %>% filter(Variable != 'control') %>% droplevels(),
                                           aes(ymin = .lower, ymax = .upper, colour = Variable),
                                           position = position_dodge(width = 0.05)) +
                           scale_colour_manual("Inducer", breaks = names(treatment_palette),
                                               values = treatment_palette) +
                           theme_classic() +
                           scale_x_continuous('Cohort settlement threshold',
                                              breaks = thress) +
                           ## scale_y_continuous('Days to >0.5 settlement probability') +
                           scale_y_continuous('Days to >0.5 settlement probability',
                                              expand = c(0,0)) 
                           ## coord_cartesian(ylim = c(min(.$.lower), unique(.$MaxY)), clip = "off")
                       )
               )
    save(all_LD50, file = paste0("../data/modelled/all_LD50.RData"))
    ## a[1,'g'][[1]][[1]]
    ## dev.off()
    pwalk(list(Species = all_LD50$Species,
               g = all_LD50$g),
               .f = function(Species, g) {
                   ggsave(filename = paste0(OUTPUT_PATH, "figures/all_LD50__", Species,"_.png"),
                          g,
                          width = 6, height = 6/1.6,
                          dpi = 72
                          )
               }
               )
    ## ----end

    ## ---- compiliationLD50Q1
    load(file = paste0("../data/modelled/all_LD50.RData"))
    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_.png"),
                                  LD50 = all_LD50,
                                  ncol = 4, dpi = 72,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))
    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_.pdf"),
                                  LD50 = all_LD50,
                                  ncol = 4,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))
    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_large_.png"),
                                  LD50 = all_LD50,
                                  ncol = 4, dpi = 300,
                                  legend.position = c(0.625, 0.03),
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))

    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_3cols_.png"),
                                  LD50 = all_LD50,
                                  ncol = 3, dpi = 72,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))
    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_3cols_.pdf"),
                                  LD50 = all_LD50,
                                  ncol = 3,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))
    LD50_compilations(path=paste0(OUTPUT_PATH,
                                  "figures/LD50Compilation_3cols_large_.png"),
                                  LD50 = all_LD50,
                                  ncol = 3, dpi = 300,
                                  legend.position = 'bottom',
                                  legend.direction = 'horizontal',
                                  legend.justification = c(0.5,0),
                                  limits = c('Best','CCA','disc', 'rubble'))
}

