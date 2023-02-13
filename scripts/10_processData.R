source('functions.R')

## ---- readData
data <- read_csv(file = paste0(DATA_PATH, "primary/AllData_Murray.csv"), trim_ws = TRUE)
data <- data %>% mutate(Plate = factor(paste0(LarvalAge, Plate)))
glimpse(data)
## ----end

## ---- processData

## We collectively decided that when Larval age is 0, it is physically
## impossible for the larvae to settle, and therefore even though it
## is not actually measured, we will set the number of settled larvae
## to be 0 at this age.  Although there has been an attempt to do this
## is the dataset, it needs to be done at the Plate level for each
## Species and SpecificTreatment.  Therefore I will remove the current
## attempt and add a more complete version.

data.lookup <- data %>%
    filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
    droplevels() %>%
    group_by(Species, SpecificTreatment) %>%
    tidyr::expand(Plate) %>% 
    mutate(LarvalAge = 0, NoSet = 0, NoNotSet = 10)

data <- data %>%
    filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
    droplevels() %>%
    full_join(data.lookup) %>% 
    mutate(Species = factor(Species),
           Family = factor(Family),
           SpecificTreatment = factor(SpecificTreatment))
save(data, file = '../data/processed/data.RData')
## ----end

## Test
data %>% filter(Species == 'Agla',
                SpecificTreatment == 'CCA') %>%
    mutate(Settlement = NoSet/(NoSet + NoNotSet)) %>%
    dplyr::select(-`...1`, -Species, -Family, - SpecificTreatment) %>%
    filter(Settlement ==1)



## ---- EDA
g <- data %>%
    ## filter(Species == "Aaus",
    ##        !SpecificTreatment %in% c("control")) %>%
    ## droplevels() %>%
    ## ggplot(aes(y = NoSet/NoAlive, x = LarvalAge, colour = SpecificTreatment)) +
    ggplot(aes(y = NoSet/(NoSet + NoNotSet), x = LarvalAge, colour = SpecificTreatment)) +
    geom_point() +
    ## geom_smooth(method = lm, formula = y ~ splines::ns(x, 5)) +
    ## geom_smooth(method = "lm", formula = y ~ poly(x, 4), se = FALSE) +
    ## geom_smooth(se = FALSE) +
    geom_smooth(method = "gam",
                formula = y ~ s(x, bs = "cs"),
                method.args = list(family = 'binomial'),
                se = FALSE) +
    facet_wrap(~Species, scales = "free")
ggsave(filename = paste0(OUTPUT_PATH, "figures/EDA.pdf"), g, width = 12, height = 12)
ggsave(filename = paste0(OUTPUT_PATH, "figures/EDA.png"), g, width = 12, height = 12)
## ----end


