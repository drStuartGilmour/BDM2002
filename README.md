# BDM2002
This repository contains R and Stata code that can be used in support of my post where I re-assess the 2002 paper by Bertrand, Duflo and Muralaithan, *HOW MUCH SHOULD WE TRUST DIFFERENCES-IN-DIFFERENCES ESTIMATES?* (referred to throughout as *BDM2002*). The blog post can be found [here]().

The R code in this repository contains the function for randomly simulating data and the three simulation tasks described in the blog post.
The Stata code downloads Current Population Survey (CPS) data for analysis and filters it to match the conditions in BDM2002, then produces graphs for the blog post, runs a basic analysis, and applies the (ridiculously complicated) interaction model at the end of the post.

Note that because I'm putting this on the internet and I don't know anything about cybersecurity, I have removed filepaths from the cd/setwd() commands, so you need to put text in those for your own path names (or comment them out).

Good luck!
