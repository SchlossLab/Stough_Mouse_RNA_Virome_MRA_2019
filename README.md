## Mouse gut RNA virus genomes assembled from metatranscriptomes

Efforts to catalogue viral diversity in the gut microbiome have largely ignored RNA viruses. To address this, we screened assemblies of previously published mouse gut metatranscriptomes for the presence of RNA viruses. We identified the complete genomes of a Astrovirus and 5 Mitovirus-like viruses. The viral fraction of the mammalian gut microbiome forms a crucial component in the relationship between microbe and host. Bacterial viruses serve as an important source of genetic diversity and population control for the microbiota, driving its ecology and evolution [@ogilvie_human_2015]. Mammalian viruses disrupt the gut environment through infection and the response of the host immune system [@legoff_eukaryotic_2017]. Bacterial and mammalian viruses make significant contributions to host health and disease. Current efforts to describe the diversity of viruses present in the gut have focused on using shotgun metagenomics to identify double-stranded DNA viruses, predominantly bacteriophage and host pathogens [**INSERT REFERENCE**]. However, this method ignores viruses with RNA genomes, which make up a considerable portion of the environmental viromes [@culley_new_2018]. We re-analyzed deeply-sequenced metatranscriptome data produced by our lab for the study of microbiome dynamics in a mouse model for *Clostridioides difficile* infection [@jenior_clostridium_2017; @jenior_clostridium_2018]. Briefly, C57Bl/6 mice from a breeding colony we maintain at the University of Michigan were treated with one of three different antibiotics (clindamycin, streptomycin, or cefoperazone). After a 24 hour recovery period, the mice were infected with *C. difficile* strain 630. Cecal contents were removed from each animal 18 hours post infection and frozen for RNA extraction and sequencing. RNA sequences from each sample were assembled individually using rnaSPAdes v3.13.1 [@bankevich_spades:_2012] and concatenated for dereplication, resulting in **R 70,779** contigs longer than 1 kb. Contigs were then screened against a custom RefSeq database of viral RNA-dependent RNA polymerase (RdRP) protein sequences with a maximum e-value of 10^-20^, resulting in **R 22** contig hits. RdRP is conserved amongst almost all RNA viruses without a DNA stage in genome replication. These contigs were then annotated with Interproscan v5.39-77.0 [@hoang_ufboot2:_2018; @kalyaanamoorthy_modelfinder:_2017]. We constructed phylogenetic trees from RdRP protein sequences using IQ-TREE v1.6.12 [@nguyen_iq-tree:_2015]. Two classes of RNA viruses were assembled with high coverage with sequences originating from most of the mouse treatment groups, including germ-free mice. First, a **R 6811** base-long astrovirus genome was obtained with **R 1683**-fold coverage (Figure 1A). The genome contained 3 predicted open reading frames encoding a capsid, RdRP, and a trypsin-like peptidase. Second, 5 distinct, but closely related RNA virus genomes ranging in length from **R 2309 to 2447** bases with **R XXXXX** to **R XXXXXX**-fold coverage belonged to a previously undescribed clade of Narnaviridae adjacent to the Mitoviruses (Figure 1B). These RNA virus genomes will facilitate future studies of RNA virus biology in the murine microbiome.



### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
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
