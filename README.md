# RangeShiftR <img src="man/figures/RSRlogo.png" align="right" height = 150/>

The RangeShiftR package implements the RangeShifter simulation platform for R.

[RangeShifter](https://rangeshifter.github.io/)
is a state-of-the-art eco-evolutionary modelling platform that is becoming 
increasingly used worldwide for both theoretical and applied purposes
[(Bocedi et al. 2014)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12162).

RangeShifter is a spatially-explicit, individual-based simulation platform that 
allows modelling species’ range dynamics, such as expansion and shifting, and 
patch connectivity by linking complex local population dynamics and dispersal 
behaviour, while also taking into account inter-individual variability and 
evolutionary processes. RangeShifter is highly flexible in terms of the spatial 
resolution and extent, and regarding the complexity of the considered ecological 
processes. Due to its modular structure, the level of detail in demographic and 
dispersal processes can be easily adapted to different research questions and 
available data.


## Installation

RangeShiftR is only available from this github repository.
(It may move to CRAN in the future.)

RangeShiftR has to be built from source and requires the package `Rcpp` as
well as a functional C++ compiler toolchain.

```r
# Install RangeShiftR from GitHub:
devtools::install_github("RangeShifter/RangeShiftR-package", ref="main")
```

## Usage

Please refer to our [website](https://rangeshifter.github.io/) for more information about RangeShifter simulation 
platform. 

A range of tutorials with theoretical and applied examples introduce you to 
the package's functionality and syntax. They can be found here:
https://rangeshifter.github.io/RangeshiftR-tutorials/


## References

 - Bocedi G, Palmer SCF, Pe’er G, Heikkinen RK, Matsinos YG, Watts K, Travis JMJ (2014). 
 *RangeShifter: A Platform for Modelling Spatial Eco-Evolutionary Dynamics and 
 Species’ Responses to Environmental Changes.* Methods in Ecology and Evolution 5: 388–96. 
