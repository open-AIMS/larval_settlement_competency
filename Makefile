# help:
# 	@echo "Usage: make -i SRC=<path/file> -> to make a specific file"
# 	@echo "       make -i                 -> to make all altered files"

SHELL:=/bin/bash

.PHONY: build build_singularity run code_container docs_container R_container code_local docs_local code_singularity docs_singularity code_slurm docs_slurm

build:
	docker build . --tag larval_settlement

build_singularity:
	docker save larval_settlement -o larval_settlement.tar 
	singularity build larval_settlement.sif docker-archive://larval_settlement.tar

# Run interactive R session in docker container
R_container:
	docker run --rm -it -v "$(shell pwd)":/home/Project larval_settlement R

code_container:
	docker run --rm -v "$(shell pwd)":/home/Project larval_settlement $(MAKE) -f scripts/Makefile

docs_container:
	docker run --rm -v "$(shell pwd)":/home/Project larval_settlement $(MAKE) -f docs/Makefile

code_local:
	$(MAKE) -f scripts/Makefile

docs_local:
	$(MAKE) -f docs/Makefile

code_singularity:
	@echo "Transfer to scripts/Makefile"
	module load singularity; \
	singularity exec -B .:/home/Project larval_settlement.sif $(MAKE) -f scripts/Makefile

docs_singularity:
	@echo "Transfer to docs/Makefile"
	module load singularity
	singularity exec -B .:/home/Project larval_settlement.sif $(MAKE) -f docs/Makefile

code_slurm:
	@echo "Submit slurm job to run code"
	sbatch code.slurm

docs_slurm:
	@echo "Submit slurm job to compile docs"
	sbatch docs.slurm

clean:
	rm -f *.log *.aux *.out texput.log *.stderr

# # Use this makefile to compile either the code or the documents
# # > make -i code
# # > make -i docs

# MODULES    := code docs
# DOCS_DIR   := $(addprefix ,docs)
# CODE_DIR   := $(addprefix ,scripts)

# SRC        := $(foreach sdir, $(DOCS_DIR), $(wildcard $(sdir)/*.Rmd))
# HTML_FILES := $(patsubst %.Rmd, %.html, $(SRC))


# $(info $(SRC))
# $(info $(HTML_FILES))

# docs: $(SRC) $(HTML_FILES)

# $(HTML_FILES): %.html: %.Rmd
# 	$(info Source = $<; Destination = $@)
# 	echo "library(rmarkdown); render(\"$<\", output_format = \"html_document\")" | R --no-save --no-restore;


# # echo "library(rmarkdown); render(\"$<\", output_format = \"html_document\", output_options = list(clean = FALSE)" | R --no-save --no-restore;
