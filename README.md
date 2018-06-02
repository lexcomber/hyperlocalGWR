# Hyper-local GWR
The code and data in this repository are for the GW-ELNR simulation study described in the _Journal of Journal of Spatial
Information Science_ submission. The code is described in the file `Comber_hyperlocal.R`. To run the simulation data you will need to download the zip file from this repository and unzip the files to a directory of your choice. This should then be set as the working directory in R using the `setwd()` function. 

Please contact Lex Comber [a.comber@leeds.ac.uk](a.comber@leeds.ac.uk) if you have any questions.

## Paper title: Hyper-local geographically weighted regression: extending GWR through local model selection and local bandwidth optimisation
Alexis Comber<sup>1</sup> and Yunqiang Wang<sup>2</sup> Yihe Lü<sup>3</sup> Xingchang Zhang<sup>4</sup> Paul Harris<sup>5</sup>

Alexis Comber1, Yunqiang Wang2, Yihe Lü3, Xingchang Zhang4, and Paul Harris5

<sup>1</sup>School of Geography, University of Leeds, Leeds, UK LS2 9JT
<sup>2</sup>State Key Laboratory of Loess and Quaternary Geology, Institute of Earth Environment, Chinese Academy of Sciences, Xi'an, China\
<sup>3</sup>State Key Laboratory of Urban and Regional Ecology, Research Center for Eco-Environmental Sciences, Chinese Academy of Sciences; Joint Center for Global Change Studies; University of Chinese Academy of Sciences, Beijing, China\
<sup>4</sup>Institute of Soil and Water Conservation, Chinese Academy of Sciences and Ministry of Water Resources, Yangling, China\
<sup>5</sup>Rothamsted Research, North Wyke, Okehampton, Devon, UK EX20 2SB\

## Abstract
Geographically weighted regression (GWR) is an inherently exploratory technique for examining process non-stationarity in data relationships. This paper develops and applies a _hyper-local_ GWR which extends such investigations further. The hyper-local GWR simultaneously optimizes both local model selection (which covariates to include in each local regression) and local kernel bandwidth specification (how much data should be included locally). These are evaluated using a measure of model fit. By allowing models and bandwidths to vary locally, it extends the 'whole map model' and 'constant bandwidth calibration' under standard GWR. It provides an alternative and complementary interpretation of localized regression. The method is illustrated using a case study modeling soil total nitrogen (STN) and soil total phosphorus (STP) from data collected at 689 locations in a watershed in Northern China. The analysis compares linear regression, standard GWR and hyper-local GWR models of STN and STP and highlights the different locations at which covariates are identified as significant predictors of STN and STP by the different GWR approaches and the spatial variation in optimal bandwidths. The hyper-local GWR results indicate that the STN relationship processes are more non-stationary and localised than found via a standard application of GWR. By contrast, the results for STP are more confirmatory (i.e. similar) between the two GWR approaches providing extra assurance to the nature of the moderate non-stationary relationships observed. The overall benefits of hyper-local GWR are discussed, particularly in the context of the original investigative aims of standard GWR. Some areas of further work are suggested.

# Acknowledgements
This work was supported by National Natural Science Foundation of China (NSFC) and the Natural Environment Research Council (NERC) Newton Fund through the China-UK collaborative research on critical zone science (No. 41571130083 and NE/N007433/1), the NSFC (No.41530854) and UK Biotechnology and Biological Sciences Research Council grant (BBSRC BB/J004308/1). All of the analyses and mapping were undertaken in R 3.3.2 the open source statistical software. The GWR analyses used the GWmodel package, v2.0-1 \cite{gollini2015gwmodel}. 
