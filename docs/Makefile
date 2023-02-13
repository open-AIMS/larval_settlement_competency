DOCS_SRC_DIR  = $(addprefix , docs)
ifdef FILE
DOCS_FILE = $(DOCS_SRC_DIR)/$(FILE)
else
# DOCS_FILE ?= $(foreach sdir, $(DOCS_SRC_DIR), $(wildcard $(sdir)/*.qmd))
DOCS_FILE ?= $(foreach sdir, $(DOCS_SRC_DIR), $(wildcard $(sdir)/*.Rmd))
endif

# HTML_FILE := $(patsubst %.qmd, %.html, $(DOCS_FILE))
HTML_FILE := $(patsubst %.Rmd, %.html, $(DOCS_FILE))

$(info ************************************)
$(info DOCS Source directory:     $(DOCS_SRC_DIR))
$(info DOCS Source files:         $(DOCS_FILE))
$(info ************************************)

# $(HTML_FILE) : %.html : %.qmd
$(HTML_FILE) : %.html : %.Rmd
	@echo "Compiling documntation for analyses"
	$(info Source = $<; Destination = $@)
	# echo "library(quarto); quarto_render(\"$<\")" | R --no-save --no-restore;
	echo "library(rmarkdown); render(\"$<\", output_format = \"html_document\")" | R --no-save --no-restore;
