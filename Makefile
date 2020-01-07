REFS = data/references
FIGS = results/figures
TABLES = results/tables
PROC = data/process
FINAL = submission/

# utility function to print various variables. For example, running the
# following at the command line:
#
#	make print-BAM
#
# will generate:
#	BAM=data/raw_june/V1V3_0001.bam data/raw_june/V1V3_0002.bam ...
print-%:
	@echo '$*=$($*)'


################################################################################

results/figures/Figure1.pdf : code/finish_trees.R\
 		code/edit_astro_tree.R\
		code/edit_mito_tree.R
	Rscript code/finish_trees.R

results/figures/Figure1.eps : results/figures/Figure1.pdf
	pdf2ps results/figures/Figure1.pdf
	mv Figure1.ps results/figures/Figure1.eps

submission/manuscript.md submission/manuscript.tex submission/manuscript.pdf : \
						submission/american-society-for-microbiology.csl\
						submission/references.bib\
						submission/manuscript.Rmd\
						results/figures/Figure1.pdf
	R -e 'render("submission/manuscript.Rmd", clean=FALSE)'
	mv submission/manuscript.knit.md submission/manuscript.md
	rm submission/manuscript.utf8.md


write.paper : submission/manuscript.md\
				submission/manuscript.tex submission/manuscript.pdf
