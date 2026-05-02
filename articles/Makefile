# Makefile — Petrification / articles
#
# Usage:
#   make <paper>        compile a specific paper (e.g. make z2_symmetry)
#   make all            compile every paper listed in PAPERS
#   make template       compile the template as a sanity check
#   make clean          remove all LaTeX auxiliary files
#   make clean-all      remove aux files AND generated PDFs
#
# Each paper lives in its own subdirectory with a .tex file
# of the same name and a symlink (or copy) of refs.bib.
# The Makefile automatically creates that symlink on first build.

LATEX  = pdflatex
BIBTEX = bibtex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error

# ---- Papers --------------------------------------------------
# Add a new entry here when starting a new article.
# The name must match both the subdirectory and the .tex filename.
PAPERS = \
    z2_symmetry \
    dyson_universality \
    perturbation_probing \
    chaos_control

# ---- Default target -----------------------------------------
.PHONY: all
all: $(PAPERS)

# ---- Template (sanity-check build) --------------------------
.PHONY: template
template: template.pdf

template.pdf: template.tex refs.bib
	$(LATEX) $(LATEXFLAGS) template
	$(BIBTEX) template
	$(LATEX) $(LATEXFLAGS) template
	$(LATEX) $(LATEXFLAGS) template

# ---- Per-paper rule -----------------------------------------
# Pattern: build <paper>/<paper>.tex using a local refs.bib symlink.
.PHONY: $(PAPERS)
$(PAPERS): %:
	@mkdir -p $*
	@if [ ! -f $*/$*.tex ]; then \
	    echo "WARNING: $*/$*.tex not found — copy template.tex to get started."; \
	    exit 0; \
	fi
	@if [ ! -e $*/refs.bib ]; then \
	    ln -s ../refs.bib $*/refs.bib; \
	fi
	cd $* && $(LATEX) $(LATEXFLAGS) $*
	cd $* && $(BIBTEX) $* || true
	cd $* && $(LATEX) $(LATEXFLAGS) $*
	cd $* && $(LATEX) $(LATEXFLAGS) $*

# ---- Clean --------------------------------------------------
AUX_EXTS = aux bbl blg log out toc lof lot synctex.gz fls fdb_latexmk

.PHONY: clean
clean:
	@for ext in $(AUX_EXTS); do \
	    find . -name "*.$$ext" -delete; \
	done
	@echo "Auxiliary files removed."

.PHONY: clean-all
clean-all: clean
	find . -name "*.pdf" -delete
	@echo "PDFs removed."
