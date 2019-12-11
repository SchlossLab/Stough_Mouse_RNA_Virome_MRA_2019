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

#make a pdf to embed figure in manuscript pdf for Pat to review
results/figures/Figure1.pdf : results/figures/Figure1.tiff
	tiff2pdf -o results/figures/Figure1.pdf results/figures/Figure1.tiff

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
