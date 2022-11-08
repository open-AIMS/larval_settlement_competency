source('functions.R')

## ---- loadData
load(file = '../data/processed/data.RData')
## ----end

## Question 1: age to settlement
{
    ## ---- Analysis loop
    for (thresholdProp in c(0.1, 0.3, 0.5)) {

        ## Prepare the data
        ## 1. Score an observation as either 1 (proportion
        ##    settled greater than threshold) or 0 (not greater)
        ## 2. For each species, treatment and plate
        ##     - sort by Larval age
        ##     - score settlement as 1 (has currently or previously
        ##       obtained a settlement proportion greater than the
        ##       threshold) or 0 (has not)
        ## 3. Establish analysis nesting (on Species)

        ## ---- prepareDataQ1
        data.q1 <- data %>%
            mutate(Settle = as.numeric(NoSet/(NoSet + NoNotSet) > thresholdProp)) %>% 
            group_by(Species, SpecificTreatment, Plate) %>%
            arrange(LarvalAge) %>%
            mutate(cumsumSettle = cumsum(Settle),
                   Settlement = ifelse(cumsumSettle>0, 1, 0)) %>%
            ungroup() %>%
            group_by(Species) %>%
            nest()

        ## remove peptile and extract - A.N. suggestion
        data.q1 <- data.q1 %>%
            mutate(data.sub = map(.x = data,
                                  .f = function(x) {
                                      x %>%
                                          filter(!SpecificTreatment %in% c('peptide','extract')) %>%
                                          droplevels()
                                  }))
        ## ----end

        ## ----EDAQ1
        data.q1 <- data.q1 %>%
            ## filter(Species == 'Agla') %>%  droplevels() %>% 
            mutate(EDA = purrr::map(.x = data.sub, .f = EDA))
        ## data.q1[1,'EDA'][[1]][[1]]
        ## dev.off()
        save(data.q1, file = paste0(DATA_PATH, "processed/data.q1.RData"))
        map2(paste0(OUTPUT_PATH, "figures/eda_",data.q1$Species,"__",thresholdProp,"_.pdf"),
             data.q1$EDA, ggsave)
        map2(paste0(OUTPUT_PATH, "figures/eda_",data.q1$Species,"__",thresholdProp,"_.png"),
             data.q1$EDA, ggsave, dpi = 72)
        map2(paste0(OUTPUT_PATH, "figures/eda_",data.q1$Species,"__",thresholdProp,"_large_.png"),
             data.q1$EDA, ggsave, dpi = 300)
        ## ----end

        ## ---- fitModelQ1
        data.q1.mod <- data.q1 %>%
            ## filter(Species == 'Agla') %>%  droplevels() %>% 
            mutate(Mod = purrr::map(.x = data.sub, .f = fitModel))
        ## ----end

        ## ---- save Models
        save(data.q1.mod, file = paste0("../data/modelled/data.q1","__",thresholdProp,".mod.RData"))
        ## ----end
    }
    ## ----end
}