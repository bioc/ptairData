# ptairData

## Description
This package contains two volatolomics raw datasets acquired with a Proton-Transfert-Reaction Time-Of-Flight mass spectrometer (PTR-TOF-MS; Ionicon Analytik GmbH, Innsbruck, Austria):

- exhaled air from two healthy individual, with three acquisitions (i.e. three files) per individual on distinct days [6 files]; each acquisition records two consecutive expirations

- mycobateria culture headspace: two replicates measurements from two different species and one control (culture medium only) [6 files]

To reduce the file size, the mass axis was truncated to the following m/z ranges: [20.4, 21.6] U [50, 150] for the human data set, and [20.4, 21.6] U [56.4, 90.6] for the mycobacteria.

## Installation
`devtools::install_github("camilleroquencourt/ptairData")`

## Preprocessing

The ptairMS package (also available on Bioconductor) is devoted to the processing of PTR-TOF-MS data.
