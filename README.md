# Virkler-Data
Virkler et al's metal fatigue crack growth data derived from Figure 4.5.3 in "Probabilistic Models of Cumulative Damage" Bogdanoff and Kozin (1985)

The file 'VirklerData-Raw.csv' was obtained using a photograph Figure 4.5.3 in "Probabilistic Models of Cumulative Damage" Bogdanoff and Kozin (1985)
The photo was slightly warped, so the number had to be adjusted to account for the slight warping.  The details are found in the 'DataCleanRCode.R' which takes the raw data and converts it into the file used for degradation analysis 'VirklerData.csv'
One issue with the dervied data is that paths do not cross, however, in the original data (Figure 4.5.1 in Bogdanoff and Kozin (1985)) it is clear some paths do cross.

The file 'VirklerHierarchicalAnalysis.R' analyzes the data using a Bayesian hierarchical model that allows for a specimen specific random effect.  The MCMC computation is accomplished in JAGS.
