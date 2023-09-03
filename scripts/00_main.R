source('functions.R')

source('10_processData.R')
source('20_fitModels.R')
source('30_modelValidations.R')
source('40_summariseModels.R')
source('50_cue_comparisons.R')

zip("../output/analysis.zip",
    flags = "-rj9X",
    files = c("../docs/analysis.html",
              "../output/figures/partialCompilation_3cols0.3_.pdf",
              "../output/figures/LD50Compilation_3cols_.pdf",
              "../output/figures/partialAreaCompilation_.pdf",
              "../output/figures/partialAreaPosteriorCompilation_.pdf",
              "../output/figures/partialAreaCompilation_Carly.pdf"
              )
    )
