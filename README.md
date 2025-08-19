# EvoTrace
A collection of R functions to work on Evolutionary Trace (ET) related jobs.

## Install
Vignettes requires pymol to build.
```
remotes::install_github("LichtargeLab/EvoTrace", build_vignettes = TRUE)
```

To install without vignettes:
```
remotes::install_github("LichtargeLab/EvoTrace", build_vignettes = FALSE)
```

## Vignettes
```
browseVignettes("EvoTrace")
```

## References
ET:
* [Lichtarge O, Bourne HR, Cohen FE. An evolutionary trace method defines 
binding surfaces common to protein families. J Mol Biol. 1996 
Mar 29;257(2):342-58.](https://pubmed.ncbi.nlm.nih.gov/8609628/)
* [Lua RC, Wilson SJ, Konecki DM, Wilkins AD, Venner E, Morgan DH, Lichtarge O. 
UET: a database of evolutionarily-predicted functional determinants of 
protein sequences that cluster as functional sites in protein structures. 
Nucleic Acids Res. 2016 Jan 4;44(D1):D308-12.](https://pubmed.ncbi.nlm.nih.gov/26590254/)

EA:
* [Katsonis P, Lichtarge O. A formal perturbation equation between
genotype and phenotype determines the Evolutionary Action of
protein-coding variations on fitness. Genome Res. 2014
Dec;24(12):2050-8.](https://genome.cshlp.org/content/24/12/2050.long)
* [Wang C, Govindarajan H, Katsonis P, Lichtarge O. ShinyBioHEAT: an
interactive shiny app to identify phenotype driver genes in E.coli and
B.subtilis. Bioinformatics. 2023 Aug
1;39(8):btad467.](https://academic.oup.com/bioinformatics/article/39/8/btad467/7234070)

SCW Z scores:
* [Mihalek I, Res I, Yao H, Lichtarge O. Combining inference from 
evolution and geometric probability in protein structure evaluation. 
J Mol Biol. 2003 Aug 1;331(1):263-79.](https://www.sciencedirect.com/science/article/abs/pii/S0022283603006636)
