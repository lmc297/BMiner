# BMiner
## A companion application for analyzing and viewing multiple BTyper output files in aggregate

## Overview

BMiner is BTyper's companion application for analyzing and viewing data from multiple BTyper output files in aggregate. BMiner can be used to:

* Analyze data generated from hundreds of *Bacillus cereus* group genomes

* Create reproducible interactive plots using concise, easily-shareable text files

* Produce publication-quality PDFs that incorporate categorical metadata from *B. cereus* group isolates

BMiner CANNOT:

* Determine if a particular statistical analysis is appropriate for a particular data set; it merely provides you (the user) with options. If you're unsure of whether an analysis (or your interpretation of it) is appropriate or not, consult a statistician.

* Definitively prove that an isolate or set of isolates is pathogenic or not. As always, interpret your results with caution. 

BMiner is a shiny application that can be run locally on your computer through R Studio. The application can be downloaded from https://github.com/lmc297/BMiner.

Post issues at https://github.com/lmc297/BMiner/issues

As always, some statistical analyses may not be appropriate for certain data sets. If you're unsure of what you're doing, consult a statistician.

Download and learn more about BTyper here: https://github.com/lmc297/BTyper.

### Citation

If you found the BMiner app and/or its source code to be useful, please cite:

Carroll, Laura M., Jasna Kovac, Rachel A. Miller, Martin Wiedmann. 2017. Rapid, high-throughput identification of anthrax-causing and emetic *Bacillus cereus* group genome assemblies using BTyper, a computational tool for virulence-based classification of *Bacillus cereus* group isolates using nucleotide sequencing data. Submitted to *Applied and Environmental Microbiology*.

------------------------------------------------------------------------

## Launching BMiner

1. Download R, if necessary: https://www.r-project.org/

2. Dowload R Studio, if necessary: https://www.rstudio.com/products/rstudio/download/

3. Open R Studio, and install the following packages, if necessary, by typing the following commands into R Studio's console:

```
install.packages("shiny")
install.packages("ggplot2")
install.packages("readr")
install.packages("stringr")
install.packages("vegan")
install.packages("plyr")
install.packages("dplyr")
install.packages("cluster")
install.packages("magrittr")
install.packages("ggrepel")
```

4. Launch the app by typing the following command into R Studio's console:
```
runGitHub("BMiner","lmc297")
```

You're ready to go!

------------------------------------------------------------------------

# Uploading Files

1. Once you've launched BMiner, be sure to click the **Open in Browser** button once the interface loads so you can select multiple files at a time.

2. Upload files by clicking the "Browse" button under "Choose BTyper final results files to analyze" in the left panel.

3. Select at least 2 BTyper final results files; final results files are text files with the extension "_final_results.txt" that BTyper deposits in the "btyper_final_results" directory that it creates within your specified output directory.

4. Optional: if you have a categorical metadata file, you can upload it by selecting the "Browse" button under "Upload a metadata file to accompany your BTyper final results files" in the left panel (to learn more about the metadata files that BMiner accepts, see the "Upload Categorical Metadata" section below).

5. Click on the tabs in the center panel to toggle between various typing schemes, and select different analyses using the menus and buttons in the left panel. Happy mining!

------------------------------------------------------------------------

# Upload Categorical Metadata

Overlay some categorical metadata on your plots! A metadata file is a 2-column tab-separated (.txt) file with no header; the first column should contain the name of each final results file, while the second should contain the corresponding metadata. To create a metadata file, follow these steps:

1. Make sure you're dealing with categorical metadata; currently, BMiner doens't support continuous metadata. Examples of categorical metadata might include host species (human, bovine, sheep, etc.), isolation source (food, clinical, environemnt, etc.), country or state of isolation (USA, France, etc.), illness-causing vs. non-illness-causing, etc.

2. Open up your favorite text editor.

3. Make a list with the name of each final_results.txt file that you are going to upload. After each file name, type **Tab**, and then type the metadata category that corresponds with that results file.

4. Save this file with a .txt extension, and upload!

A quick example: I have 4 *B. cereus* group isolates; strain "foo" was isolated from a cow, strain "bar" was isolated from a pig, and strains "spam" and "eggs" were isolated from humans. They have assemblies, each named foo.fasta, bar.fasta, spam.fasta, and eggs.fasta, respectively. I run BTyper, which gives me 4 final results files named foo_final_results.txt, bar_final_results.txt, spam_final_results.txt, and eggs_final_results.txt, respectively, so my metadata file should resemble the following, with a **Tab** separating the file names from the isolate sources:

```
foo_final_results.txt cow
bar_final_results.txt pig
spam_final_results.txt  human
eggs_final_results.txt  human
```

------------------------------------------------------------------------

# Analyses Using BMiner

## Analyses using data obtained from virulence gene detection

### Virulence Gene Distribution (Bar Plot)

This analysis produces a bar plot with virulence genes detected in your selected files on the x-axis, and the (1) total number of isolates in your uploaded results files that possess each gene, or (2) the percent (%) of isolates in your uploaded files that possess each gene, depending on your selection in y-axis display. 

To download the plot as a PDF, click the "Download this plot" in the left panel.

### Non-Metric Multidimensional Scaling (NMDS)

This performs NMDS using the 'metaMDS' function in R's 'vegan' package and produces an interactive, 2-dimensional plot using 'ggplot2'.

NMDS is a method of ordination in which input data (in this case, a virulence gene presence/absence matrix) are fit to *k* dimensions (in BMiner, k is set to 2 so that the results can be plotted in 2 dimensions). 

NMDS is valuable because it makes few assumptions about the nature of your data, and it allows any distance/dissimilarity measure to be used. At this time, I've given users the option of using (1) Jaccard, or (2) Raup-Crick dissimilarity for their data, as both are suitable for presence-absence data.

Potential drawbacks of NMDS are that it is more compuationally intensive and, thus, slower than other ordination methods. For very large data sets (i.e. hundreds of assemblies), I would recommend using PCA (described below), as it is much faster. Also, depending on the nature of your data, `metaMDS` may never be able to find an optimal solution. 
Note: I have hard-coded `metaMDS` as run through BMiner to try a maximum of 10000 random starts (`trymax=10000`). `metaMDS` will stop once this limit is reached.

In BMiner, NMDS plots display NMDS axis 1 on the x-axis (NMDS1) and NMDS axis 2 on the y (NMDS2). Points represent your input assemblies, while the location of virulence genes on the plot are shown by lines drawn using `ggrepel`. Clicking a point will give you a list of the isolates associated with it under the plot, while dragging your mouse and double-clicking allows you to zoom in on the plot (double-click to zoom back out). If you have input categorical data, you can draw convex hulls around corresponding assemblies by checking the "Overlay Metadata" box. To export your plot, click "Download this plot" in the left panel.

### Principal Component Analysis (PCA)

This performs PCA using a virulence gene presence/absence matrix based on the input data and the 'prcomp' function in R's 'stats' package. An  interactive plot is produced using 'ggplot2'.

By default, BMiner plots the first, second, and third principal components (PC1, PC2, and PC3, respectively) one the x-axis, y-axis, and as the point size, respectively. If you want to change which PCs are viewed, you can select them from the lists in the left panel. You can determine which assemblies are associated with which points by clicking on them, and you can zoom in by dragging your mouse and double-clicking on an area of the plot (double-click to zoom back out). You can color your points using your categorical metadata by checking the "Overlay Metadata" box. To download the plot, click "Download this plot" in the left panel.

Note: I have hard-coded `prcomp` to scale and center the data.

### Virulence Gene Presence/Absence Table

This produces the virulence gene presence/absence matrix that is used for the aforementioned virulence analyses. You can download this matrix by clicking "Download this table" in the left panel and use it for any post-hoc analyses you want (including the ones outlined in the BTyper manuscript!).

## Analyses using multi-locus sequence typing (MLST) data

### Sequence Type (ST) Distribution (Bar Plot)

This analysis produces a bar plot with sequence types (STs) found in your selected files on the x-axis, and the (1) total number of isolates in your uploaded results files that possess each ST, or (2) the percent (%) of isolates in your uploaded files that possess each ST, depending on your selection in y-axis display. 

To download the plot as a PDF, click the "Download this plot" in the left panel.

### Minimum Spanning Tree

This analysis produces a minimum spanning tree using your MLST data, the `spantree` function in R's `vegan` package, and `ggplot2`. Points represent MLST sequence types, while point size corresponds to the number of isolates with a particular sequence type. To download the plot, click the "Download this plot" in the left panel. To zoom in on an area of the plot, drag your mouse and double-click (double-click to zoom out).

The miminum spanning tree created by BMiner is based on the 7 *B. cereus* allelic types in a file (not nucleotide sequences themselves). As a result, the allelic types are treated as factors. For example, ST 1 is represented by a vector `c("1","1","1","1","1","1","1")`, while ST 3 is represented by a vector `c("2","1","1","1","1","1","1").` These 2 STs are more similar than they are to ST 778, which is represented by vector `c("56","1","93","1","193","37","160")`. By calculating the tree in this fashion, you can rapidly compare the allelic profiles of unknown sequence types (represented by "?" with an arbitrary number after them, to distinguish them from other unknown sequence types with different allelic profiles) to known ones to visually assess their similarity. Currently, I have hard-coded BMiner to produce minimum spanning trees using "gower" distance (as calculated using the `daisy` function in R's `cluster` package). I chose Gower's diatance because it can be used with categorical data, and allelic types are treated as categorical in BMiner.

Note: the "Overlay Metadata" and interactive plot-click functions for minimum spanning trees are currently under development. For small data sets, you can check the "Overlay Metadata", and pie charts will appear in lieu of points. However, this has not been optimized for larger data sets yet, so you may get an error message. To get the ST information of a particular point, look at the minimum spanning tree plot in the "Plots" tab of R Studio, or download a PDF of your plot by clicking the "Download this plot" button in the left panel (this will output a simple plot with ST information overlaid on points, as well as the `ggplot` version displayed by BMiner).

## Analyses using *rpoB* allelic typing data

### *rpoB* Allelic Type (AT) Distribution (Bar Plot)

This analysis produces a bar plot with *rpoB* allelic types (ATs) found in your selected files on the x-axis, and the (1) total number of isolates in your uploaded results files that possess each AT, or (2) the percent (%) of isolates in your uploaded files that possess each AT, depending on your selection in y-axis display. 

To download the plot as a PDF, click the "Download this plot" in the left panel.

## Analyses using *panC* clade typing data

### *panC* Clade Distribution (Bar Plot)

This analysis produces a bar plot with *panC* clade assignments found in your selected files on the x-axis, and the (1) total number of isolates in your uploaded results files that belong to each *panC* clade, or (2) the percent (%) of isolates in your uploaded files that belong to each *panC* clade, depending on your selection in y-axis display. 

To download the plot as a PDF, click the "Download this plot" in the left panel.

## Analyses using 16S rDNA data

### Closest 16S Sequence Distribution (Bar Plot)

This analysis produces a bar plot with the 16S assignment found in your selected files on the x-axis, and the (1) total number of isolates in your uploaded results files that were closest to each of 9 type strain 16S sequences, or (2) the percent (%) of isolates in your uploaded files that are closest to that 16S sequence, depending on your selection in y-axis display. 

To download the plot as a PDF, click the "Download this plot" in the left panel.

Note: Interpret 16S rDNA results with extreme caution! The 16S gene is not optimal for assessing *B. cereus* group taxonomy or pathogenic potential.

------------------------------------------------------------------------


# References

### R Packages

H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2009.

Hadley Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. Journal of Statistical Software, 40(1), 1-29. URL http://www.jstatsoft.org/v40/i01/.

Hadley Wickham (2017). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.2.0. https://CRAN.R-project.org/package=stringr

Hadley Wickham, Jim Hester and Romain Francois (2017). readr: Read Rectangular Text Data. R package version 1.1.0. https://CRAN.R-project.org/package=readr

Hadley Wickham and Romain Francois (2016). dplyr: A Grammar of Data Manipulation. R package version 0.5.0. https://CRAN.R-project.org/package=dplyr

Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2017). vegan: Community Ecology Package. R package version 2.4-2. https://CRAN.R-project.org/package=vegan

Kamil Slowikowski (2016). ggrepel: Repulsive Text and Label Geoms for 'ggplot2'. R package version 0.6.5. https://CRAN.R-project.org/package=ggrepel

Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K.(2017).  cluster: Cluster Analysis Basics and Extensions. R package version 2.0.6.

Stefan Milton Bache and Hadley Wickham (2014). magrittr: A Forward-Pipe Operator for R. R package version 1.5. https://CRAN.R-project.org/package=magrittr

Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2017). shiny: Web Application Framework for R. R package version 1.0.1. https://CRAN.R-project.org/package=shiny

### Journals

J.B. Kruskal. Multidimensional scaling by optimizing goodness
of fit to a nonmetric hypothesis. 1964. *Psychometrika* vol. 29, no. 1.

J.B. Kruskal. Nonmetric multidimensional scaling: A numerical method. 1964. *Psychometrika* vol. 29, issue 2, 115-129.

Jonathan M. Chase, Nathan J. B. Kraft, Kevin G. Smith, Mark Vellend, Brian D Inouye. Using null models to disentangle variation in community dissimilarity from variation in α-diversity. 2011. *Ecosphere* vol. 2, issue 2, 1-11.

### Other Materials

Steven M. Holland. Non-Metric Multidimensional Scaling (MDS). May 2008. http://strata.uga.edu/software/pdf/mdsTutorial.pdf

------------------------------------------------------------------------

Disclaimer: BTyper and BMiner are pretty neat! However, no tool is perfect, and BTyper and BMiner cannot definitively prove whether a *B. cereus* group isolate is pathogenic or not. As always, interpret your results with caution. We are not responsible for taxonomic misclassifications, misclassifications of an isolate's pathogenic potential, and/or misinterpretations (biological, statistical, or otherwise) of BTyper and/or BMiner results.
