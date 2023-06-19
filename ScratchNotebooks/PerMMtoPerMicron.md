# Statement of the problem
I discovered, after submission to a journal, that when I divide by bin-width, my units are generally in #/L/μm but I've been reporting them as #/L/mm. I need to correct this in the code and in the manuscript.

Fortunately, this problem affects only units that I report things in and not any of the conclusions. I will fix everything here and in a draft of the manuscript. I may or may not ask the journal to send an update to reviewers.

The following figures need to be updated

Figure 2, panels A&C -- MassAndAbundanceNoMap.png -- AllMicrobesFigures.qmd; (Resaved)
Figure S2, all, MCN_Molarity.png -- ParticulateMolarity.Rmd; (Resaved)
Figure S23, all, bottomNS.png, oxyclineNS.png, surfaceNS.png -- OstensiblyFree.Rmd; (Resaved)

I need to comment the code so that when I define something as /L/mm, I clarify in the comments that the units are wrong.
Actually, it turns out that we devide by bin width when plotting happens and not before. So things are correct as written.

I need to update the text in two places.

For reasons I don't understand, the order of the plots are now wrong. I had to go in and manually reorder factors.
Also, in the new verison a few more POC/Mass points somehow appeared. I'm not sure why.