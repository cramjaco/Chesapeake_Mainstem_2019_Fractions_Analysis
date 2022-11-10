# Analyis of size fractionated amplicon data from the Chesapeake Bay Mainstem 2019

Author: Jacob Cram

Date: 2022 November 09

This project contains the analysis of the amplicon sequencing data from Jacob Cram's lab's
size fractionation data from a 2019 RV Carson cruise on the Chesapeake Bay. 
These data complement the manuscript: "Microbial diversity and abundance vary along latitudinal and particle size gradients in the Chesapeake Bay".

Authors: Jacob Cram, Ashley Hollins, Alexandra McCarty, Grace Martinez, Clara Fuchsman
Data were collected by Jacob Cram, Ashley Hollins, Emily Dougherty, Mekayla Reynolds

Directory Structure follows:

 *  `RScripts`: Some produce intermediate output used by other things and some produce figures.
    * `InitalProcessing_3.R`: Data wrangling and pre-processing of ASV and environmental data. This is slow, and so generally I used its saved environment output in most of the notebooks.
    * `Brigandine_2.R `: Abundances of bacteria on different sizes of particles. Phyla, Planktomycetes ASVs, abundant ASVs on > 5 micron particles
    * `CBMap_Transect.R`: Make a map of all stations.
    * `AmpliconsBiogeochem`: Identify amplicons with names or closest sequenced relatives involved in key biogeochemical cycles. Returns a plot of the organisms' abundances.
    
* `ActiveNotebooks`: Notebook files in which much of the analysis in the manuscript took place.
  * `BreakawayAlphaDiversity`: Calculations of alpha diversity using the breakaway package. Called by `ConventioanlAlpha.Rmd`
  * `ConventionalAlpha.Rmd` All of the alpha diversity calculations, including the ones for the manuscript.
     Note that we settled on sharing Breakaway richness, Shannon Diversity and Evenness (Breakawy/Shannon)
  * `Exemplars.Rmd` Gives an example of the abundance patterns of one microbe at one site. Used in proposals but not the manuscript.
  * `AllMicrobesFigures.qmd`: Plots of total microbial abundance per L and per mg of particles.
  * `MockExplore.Rmd`: Provides a sanity check of the mock community data. In the manuscript, I mostly go so far as to say "looks good" but anyone is welcome to explore these more.
  * `OstensiblyFree`: Exploration of the bacteria that are "free living" in that they comprise at least 1% of the microbial community. Here I show that *all* of these microbes are also found on some particle size fractions.
  * `SummaryDataOverview.Rmd`: Some overview data about amplicon sequence quality
* Archive
