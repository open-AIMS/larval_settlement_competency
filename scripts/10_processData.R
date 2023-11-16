source('functions.R')

## ----readSpeciesAbbrev
speciesAbbrev <- read_csv(file = paste0(DATA_PATH, "primary/Species abbreviations.csv"), trim_ws = TRUE)
speciesAbbrev <- speciesAbbrev %>%
    mutate(Species = ifelse(Species == "Acropora tenuis",
                            "Acropora cf. kenti",
                            Species))
save(speciesAbbrev, file = paste0(DATA_PATH, "processed/Species abbreviations.RData"))
## ----end

## ---- readData
## data1 <- read_csv(file = paste0(DATA_PATH, "primary/AllData_Murray.csv"), trim_ws = TRUE)
## data1 <- data1 %>% mutate(Plate = factor(paste0(LarvalAge, Plate)))
data <- read_csv(file = paste0(DATA_PATH, "primary/AllData_Murray 2.csv"), trim_ws = TRUE)
glimpse(data)
data <- data %>% mutate(PLATE = factor(DateAgePlate))
## ----end

## ---- processData

## We collectively decided that when Larval age is 0, it is physically
## impossible for the larvae to settle, and therefore even though it
## is not actually measured, we will set the number of settled larvae
## to be 0 at this age.  Although there has been an attempt to do this
## is the dataset, it needs to be done at the Plate level for each
## Species and SpecificTreatment.  Therefore I will remove the current
## attempt and add a more complete version.

data.1 <- data %>%
    mutate(Plate = factor(paste0(SpawnDate, LarvalAge, Plate))) %>%
    filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
    droplevels()
data.2 <- data.1 %>%
    mutate(
        ReadDate = SpawnDate,
        LarvalAge = 0,
        NoSet = 0, NoNotSet = 10, NoAlive = 10
    )
data <- data.1 %>%
    bind_rows(data.2) %>%
    mutate(Species = factor(Species),
           Family = factor(Family),
           SpecificTreatment = factor(SpecificTreatment))

TreatmentOrder <- c("control", "CCA", "rubble", "disc", "extract", "peptide")
save(TreatmentOrder, file = paste0(DATA_PATH, "processed/TreatmentOrder.RData"))
TreatmentPalette <- c(control = "gray", CCA = "pink", rubble = "blue", disc = "green",
                      extract = "purple", peptide = "yellow")
TreatmentPalette <- scales::hue_pal()(6)[c(2,1,4,3,6,5)]
names(TreatmentPalette) <- TreatmentOrder
save(TreatmentPalette, file = paste0(DATA_PATH, "processed/TreatmentPalette.RData"))
## data <- data %>%
##     mutate(Plate = paste0(SpawnDate, LarvalAge, Plate)) %>%
##     filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
##     droplevels() %>%
    
##     group_by(Species, SpecificTreatment, SpawnDate) %>%
##     summarise(data = list(cur_data_all()), .groups = "drop") %>%
##     mutate(data = map(.x = data,
##                       .f = ~ {
##                           d <- .x %>%
##                               group_by(LarvalAge) %>%
##                               mutate(Count = n()) %>%
##                               dplyr::select(-ReadDate, -Plate, -DateAgePlate,
##                                             -AgePlate, -Well, -Rep, -NoSet, -NoNotSet,
##                                             -NoAlive, -TargetTotal, Q1, Q2) %>%
##                               distinct() %>%
##                               ungroup() %>% 
##                               ##filter(LarvalAge == min(LarvalAge)) %>%
##                               mutate(ReadDate = SpawnDate,
##                                      ## use distinct LarvalAge as a proxy for plate
##                                      Plate = paste0(ReadDate, " _ ", 0, " _ ", LarvalAge),
##                                      LarvalAge = 0,
##                                      NoSet = 0, NoNotSet = 10
##                                      )
##                           .x %>% bind_rows(d)
##                           }
##                       )) %>%
##     dplyr::select(data) %>%
##     unnest(c(data)) %>% 
##     ungroup() %>%
##     mutate(Species = factor(Species),
##            Family = factor(Family),
##            SpecificTreatment = factor(SpecificTreatment))
    
## Old attemp
## data.lookup <- data %>%
##     filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
##     droplevels() %>%
##     group_by(Species, SpecificTreatment) %>%
##     tidyr::expand(Plate) %>% 
##     mutate(LarvalAge = 0, NoSet = 0, NoNotSet = 10)

## data <- data %>%
##     filter(LarvalAge != 0) %>% # exclude previous attempts to set settlement at 0 when age is 0
##     droplevels() %>%
##     full_join(data.lookup) %>% 
##     mutate(Species = factor(Species),
##            Family = factor(Family),
##            SpecificTreatment = factor(SpecificTreatment))
save(data, file = '../data/processed/data.RData')
## ----end

## Test
data %>% filter(Species == 'Agla',
                SpecificTreatment == 'CCA') %>%
    mutate(Settlement = NoSet/(NoSet + NoNotSet)) %>%
    dplyr::select(-`...1`, -Species, -Family, - SpecificTreatment) %>%
    filter(Settlement ==1)



## ---- EDA
load(file = '../data/processed/data.RData')
load(file = paste0(DATA_PATH, "processed/TreatmentOrder.RData"))
load(file = paste0(DATA_PATH, "processed/TreatmentPalette.RData"))
g <- data %>%
    mutate(SpecificTreatment = factor(SpecificTreatment,
                                      levels = TreatmentOrder)) %>%
    left_join(speciesAbbrev %>%
              dplyr::select(Species1 = Species,
                            Species = Abbreviation)) %>%
    dplyr::select(-Species) %>%
    dplyr::rename(Species = Species1) %>%
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
    ## scale_colour_viridis_d() +
    scale_colour_manual("Treatment", breaks = TreatmentOrder, values = TreatmentPalette) +
    facet_wrap(~Species, scales = "free") +
    theme_classic()
ggsave(filename = paste0(OUTPUT_PATH, "figures/EDA.pdf"), g, width = 12, height = 12)
ggsave(filename = paste0(OUTPUT_PATH, "figures/EDA.png"), g, width = 12, height = 12)
## ----end


