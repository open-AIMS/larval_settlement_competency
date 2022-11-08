# Larval settlement competency

This repository represents the code and environment required to
reproduce coral larval settlement competency analyses for Carly
Randall (AIMS).  **Note, this repo does not include the data.**

Whilst it is possible to run the R code (starting with `00_main.R`),
to ensure full compatibility, it is recommended that the code be run
through either a `docker` image created from the included `Dockerfile`
or else a `singularity` image, similarly built from the aforementioned
docker image.  These approaches will safeguard the codebase from
inevitable changes to the underlying runtime environment(s).

# To build the Docker image

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

# To run the code in the docker container

   `docker run --rm -v <path to repo>:<path within the container> -w <working director> Rscript 00_main.R`
   
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
   `docker run --rm -v "$PWD":~/home/larval_settlement -w /home/larval_settlement/scripts Rscript 00_main.R`
   

