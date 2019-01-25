# olfactory-expression

This is a Shiny app for quickly examining expression data. 

To run the app with Binder (in your browser), click this button: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/diyadas/olfactory-expression/master?filepath=rstudio), open app.R in the RStudio tab that it launches, and then click the "Run this app" button.

To run the app locally, either clone this GitHub repo or download a zip file. Then open app.R in RStudio, and then click the "Run this app" button.


By default, the app plots data from our analysis of stem cell-mediated regeneration in the olfactory epithelium.

Gadye L*, Das D*, Sanchez MA*, Street KN, Baudhuin A, Wagner A, Cole MB, Choi YG, Yosef N, Purdom E, Dudoit S, Risso D, Ngai J, Fletcher RB. Injury Activates Transient Olfactory Stem Cell States with Diverse Lineage Capacities. (2017). Cell Stem Cell _21_, 775-790.e9.

Read it [here](https://doi.org/10.1016/j.stem.2017.10.014).

But you can upload your own datasets too. It requires cluster labels for the heatmap, but you also have the option of uploading batch and experimental condition information too.

A sample gene list `oe_markers.txt` is provided for the default dataset.