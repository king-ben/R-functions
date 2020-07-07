# R-functions
Miscellaneous R functions for phylogenetics that I use regularly

Archive Contents
----------------

* `convert.lognormal.R` : converts a mean and standard deviation in real space to log space
* `geoplot.epoch.rates.R` : plots outpus of get.epoch.rates against the geological timescale
* `get.epoch.rates.R` : extracts mean (weighted branch length) evolutionary rates in defined epochs from phylogenetic tree samples
* `hpd.R` : calculates highest posterior density interval (default is 95%) from a vector of values.
* `/monophyl.multi.R` : estimates the posterior probability of a group (or several groups) being monophyletic from a posterior tree sample. Includes the option to exclude rogue taxa.
* `/operators.sa.xml.R` : Converts table of upper and lower age bounds for taxa into beast tip date operators for the sampled-ancestros package.
* `tipheight.R` : returns height of a taxon in the sense used by beast2 (i.e. age minus age of youngest taxon).
