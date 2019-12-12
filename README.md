## Mouse gut RNA virus genomes assembled from metatranscriptomes

Efforts to catalogue viral diversity in the gut microbiome have largely ignored RNA viruses. To address this, we screened assemblies of previously published mouse gut metatranscriptomes for the presence of RNA viruses. We identified the complete genomes of a Astrovirus and 5 Mitovirus-like viruses.



### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- manuscript.Rmd    # executable Rmarkdown for this study, if applicable
	| |- manuscript.md     # Markdown (GitHub) version of the *.Rmd file
	| |- manuscript.tex    # TeX version of *.Rmd file
	| |- manuscript.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- american-society-for-microbiology.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* R (v. 3.6.1) should be located in the user's PATH
* R packages:
  * `knitr`
  * `rmarkdown`
  * `ape`
  * `phangorn`
  * `ggtree`
  * `tidyverse`
  * `cowplot`
* SPAdes v3.13.1
* Trimmomatic v0.39
* Bowtie 2 v2.3.4.3
* ccontigs
* BBMap v36.68
* NCBI blast v2.9.0
* Prodigal v2.6.3
* Interproscan v5.39-77.0
* MAFFT v7.453
* Trimal v1.4
* IQ-TREE 1.6.12

#### Running analysis

```
git clone https://github.com/JMAStough/Stough_Mouse_RNA_Virome_MRA_2019.git
make write.paper
```
