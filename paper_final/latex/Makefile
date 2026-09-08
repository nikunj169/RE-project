# Makefile for LaTeX compilation

LATEX = pdflatex
BIBTEX = bibtex

.PHONY: all clean main supp

all: main supp

main: main.pdf

main.pdf: main.tex references.bib
	$(LATEX) main
	$(BIBTEX) main
	$(LATEX) main
	$(LATEX) main

supp: supplementary.pdf

supplementary.pdf: supplementary.tex references.bib
	$(LATEX) supplementary
	$(BIBTEX) supplementary
	$(LATEX) supplementary
	$(LATEX) supplementary

clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz
