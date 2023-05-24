source('functions.R')

## Question 1: age to settlement
{
    ## ---- Analysis loop
    for (thresholdProp in seq(0.1, 0.9, by = 0.1)) {
        ## ---- load Model Validations
        load(file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod.RData"))
        ## ----end

        ## ---- MCMCdiagnosticsQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(Trace = purrr::map(.x = Mod, .f = MCMCdiagnostics_trace),
                   AC  = purrr::map(.x = Mod, .f = MCMCdiagnostics_ac),
                   Rhat  = purrr::map(.x = Mod, .f = MCMCdiagnostics_rhat),
                   ess  = purrr::map(.x = Mod, .f = MCMCdiagnostics_ess)
                   )
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$Trace, ggsave)
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$Trace, ggsave, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/Trace_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$Trace, ggsave, dpi = 300)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$AC, ggsave)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$AC, ggsave, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/AC_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$AC, ggsave, dpi = 300)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$Rhat, ggsave, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$Rhat, ggsave, dpi = 72, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/Rhat_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$Rhat, ggsave, dpi = 300, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$ess, ggsave, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$ess, ggsave, dpi = 72, width = 6, height = 3)
        map2(paste0(OUTPUT_PATH, "figures/ess_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$ess, ggsave, dpi = 300, width = 6, height = 3)
        ## ----end

        ## ---- DHARMaQ1
        data.q1.mod <- data.q1.mod %>%
            mutate(DHARMa = purrr::map(.x = Mod, .f = DHARMa))
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q1.mod$Species,"__",thresholdProp,"_.pdf"),
             data.q1.mod$DHARMa, ggsave, width = 12, height = 4)
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q1.mod$Species,"__",thresholdProp,"_.png"),
             data.q1.mod$DHARMa, ggsave, width = 12, height = 4, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/DHARMa_",data.q1.mod$Species,"__",thresholdProp,"_large_.png"),
             data.q1.mod$DHARMa, ggsave, width = 12, height = 4, dpi = 300)
        ## ----end

        ## ---- save Model Validations
        data.q1.mod <- data.q1.mod %>%
            dplyr::select(-Mod)
        save(data.q1.mod,
             file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod_valid.RData"))
        ## ----end
    }
    ## ---end
}
