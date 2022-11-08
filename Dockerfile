FROM rocker/r-ver:4.2.1

## Install the os packages
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    lbzip2 \
    libssl-dev \
    libudunits2-dev \
    imagemagick \
    pandoc \
    pandoc-citeproc \
    make \
    ghostscript \
    poppler-utils \
    zip \
    wget \
    fonts-dejavu-extra \
    curl \
    tk \
    libglpk-dev \ 
  && rm -rf /var/lib/apt/lists/*

## Install R package versions from MRAN (based on a date - YYYY-MM-DD)
RUN R -e "install.packages('stringi');	\
"

# RUN R -e "\ #options(repos = \
#     # list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-10-04/'));\
#   install.packages('stringi');	\
# "

RUN R -e "install.packages('tidyverse');	\
"

# RUN R -e "options(repos = \
#     list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-10-04/'));\
#   install.packages('tidyverse');	\
# "

RUN R -e " \
  install.packages('DHARMa');   \
  install.packages('patchwork');   \
"

RUN R -e " \
  install.packages('emmeans');   \
"

RUN R -e "  \
  install.packages('remotes');  \
  install.packages('rstan');  \ 
  install.packages('brms');   \
  install.packages('tidybayes'); 	\
  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos'))); \
  remotes::install_github('stan-dev/cmdstanr'); \
  library(cmdstanr); \
  check_cmdstan_toolchain(); \
  install_cmdstan(cores = 2); \
"

# RUN R -e "options(repos = \
#     list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2022-10-04/'));\
#   install.packages('tidybayes');	\
#   install.packages('emmeans');  \
# "

# RUN R -e "install.packages('checkpoint'); \
# library(checkpoint); \
# checkpoint('2022-10-01'); \
# install.packages('tidybayes'); \
# "

  #install.packages('bookdown'); \
  #install.packages('emmeans'); \
  #install.packages('DHARMa'); \
  #install.packages('patchwork'); \
#############################################################################
## NOTE: we could opt for installing specific versions of packages
## RUN R -e "install.packages('remotes');"
## RUN R -e "remotes::install_version('dplyr', '1.0.5');"
#############################################################################
