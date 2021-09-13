#' mycobacteria dataset
#' This data set was generated with a proton-transfer-reaction
#' quadruple time-of-flight mass spectrometer (Ionicon Analytik GmbH,
#' Innsbruck, Austria). The headspace has been
#' analyzed after 1 week of cell culture.
#'
#' exhaled air
#' Two volunteers blew into the mass spectrometer (Ionicon Analytik GmbH,
#' Innsbruck, Austria) equipped with a BET-med (buffer end tidal).
#'
#' simulated PTR-TOF-MS data from exhaled breath
#' One simulated file, generated as follows: first, peak clusters were
#' generated around nominal masses 21 to 400, with an asymmetric sech2 peak
#' shape distribution. Second, temporal evolutions were exacted from a large
#' in-house database of patient acquisitions (> 10,000 expiration and ambient
#' air profiles), after normalization and Savitzky-Golay smoothing. Third,
#' simulation parameters were randomly selected for each nominal mass: i.e.,
#' the asymmetry coefficient, the number of overlapping peaks (1 to 3), the peak
#' proximity, the intensity of the highest peak, the ratio of neighbor peaks, and
#' the class of temporal profile (“expiration”, “ambient air” or “constant”). The
#' exact m/z value of the first peak was selected from the formula library CxHyOzNt
#' used by PTRwid (Holzinger, 2015). Fourth, the peak width was chosen with a
#' Gaussian random sampling of the mean as a function of the mass (m/resolution;
#' with resolution set to 5,000). Finally, background noise  was added
#' by using a Poisson stochastic process described in Gundlach-Graham et al.
#' (2018), with a Gaussian distribution to model the single ion Pulse-Height.
#'


