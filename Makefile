# Use this makefile to compile either the code or the documents
# > make -i code
# > make -i docs

MODULES    := code docs
DOCS_DIR   := $(addprefix ,docs)
CODE_DIR   := $(addprefix ,scripts)

SRC        := $(foreach sdir, $(DOCS_DIR), $(wildcard $(sdir)/*.Rmd))
HTML_FILES := $(patsubst %.Rmd, %.html, $(SRC))


$(info $(SRC))
$(info $(HTML_FILES))

docs: $(SRC) $(HTML_FILES)

$(HTML_FILES): %.html: %.Rmd
	$(info Source = $<; Destination = $@)
	echo "library(rmarkdown); render(\"$<\", output_format = \"html_document\")" | R --no-save --no-restore;


# echo "library(rmarkdown); render(\"$<\", output_format = \"html_document\", output_options = list(clean = FALSE)" | R --no-save --no-restore;
