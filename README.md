# Larval settlement competency

This repository represents the code and environment required to
reproduce coral larval settlement competency analyses for Carly
Randall (AIMS).  **Note, this repo does not include the data.**

Whilst it is possible to run the R code (starting with `00_main.R`),
to ensure full compatibility, it is recommended that the code be run
through either a [docker](https://www.docker.com) image created from
the included `Dockerfile` or else a
[apptainer/singularity](https://apptainer.org/) image on a HPC,
similarly built from the aforementioned docker image.  These
approaches will safeguard the codebase from inevitable changes to the
underlying runtime environment(s).

# To run the code without `docker` or `singularity`

1. clone this repo (into the current working directory)

   `git clone https://github.com/open-AIMS/larval_settlement_competency.git`
   
   or if you use ssh keys
   
   `git clone git@github.com:open-AIMS/larval_settlement_competency.git`

2. start a new R session and run the `00_main.R` script in the
   `scripts/` folder

# To run within docker

## To build the Docker image

1. ensure [docker](https://www.docker.com) is installed on your local
   machine

2. navigate to this repositories root

3. build the docker image (with a tag to help distinguish this from
   other docker images you might have)
   
   `docker build . --tag R_brms`

4. check that he image has been created
   
   `docker images --all`

5. test whether the docker image works in a superficial way (does it
   drop you into an interactive R session - `q()` to quit)

   `docker run --rm -it R_brms R`

## To run the code in the docker container

   `docker run --rm -v <path to repo>:<path within the container> -w <working director> R_brms Rscript 00_main.R`
   
   where:
   
   - `<path to repo/scripts>` the path (relative or full) to the
   location of the local repository.  This can be simplified if you
   issue the `docker run` command from the root of the repo and use
   the bash variable `$PWD`

   - `<path within container>` is the path (full) of a location to
   mount the local folder to within the container.  This location will
   be created on the fly.

   - `<working director>` is the full path within the container to
     become the working directory.
   
   For example:
   `docker run --rm -v "$PWD":~/home/larval_settlement -w /home/larval_settlement/scripts R_brms Rscript 00_main.R`

# To run within apptainer/singularity on a HPC

## To build the apptainter/singularity image

1. ensure [apptainer/singularity](https://apptainer.org/) is installed on your local
   machine

2. save a local copy of the docker image as an archive (to the current
   working directory)
   
   `docker save R_brms -o R_brms.tar`
   
   This will create a zip (tar) of the docker image.
   
3. build the apptainer/singularity image from this archive

   `singularity build R_brms.sif docker-archive://R_brms.tar`

4. test whether the apptainer/singularity image works in a superficial
   way (does it drop you into an interactive R session - `q()` to
   quit)
   
   `singularity exec -B <path to repo>:<path within the container> R_brms.sif R` 

## Run the apptainer/singularity image on a HPC

1. clone the repository to the HPC (must be `ssh`'d into the HPC first)

   `git clone git@github.com:open-AIMS/larval_settlement_competency.git`


2. scp singularity image to the HPC

   `scp R_brms.sif <user>@hpc-l001.aims.gov.au:<path on HPC>`
   
   where:
   
   - `<user>` is your username of the HPC
   - `<path on HPC>` is the path of the repository on the HPC
   
3. run the apptainer/singularity container

   `singularity exect -B <path on HPC>:<path within container> R_brms.sif --pwd <working directory> Rscript 00_main.R`
   
   where:
   
   - `<path on HPC>` is the path of the repository on the HPC

   - `<path within container>` is the path (full) of a location to
      mount the host folder to within the container.  This location
      will be created on the fly.
	  
   - `<working directory>` is the full path within the container to
     become the working directory.

4. run the documentation (first navigate to the repo folder)

   `singularity exect -B <path on HPC>:<path within container> R_brms.sif --pwd <working directory> make -i`
   
   this should create a file called `analysis.html` in the `docs` folder
