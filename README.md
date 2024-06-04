# How to compute the Mosquito basic reproduction number, $R_M$
This code computes the Mosquito basic reproduction number, $R_M$, for $Aedes$ $albopictus$ and $Aedes$ $aegypti$ following the methodology in:
https://www.biorxiv.org/content/10.1101/2024.05.31.596775v1

The Mosquito basic reproduction number is given by:
$$R_{M} =  \sqrt[3]{f\frac{a}{\delta_A}\frac{d_Eh}{(d_Eh + \delta_E)}     \frac{d_I}{(d_I + \delta_I)}} = \sqrt[3]{f\frac{a}{\delta_A}p_{EL}p_{LA}}.$$

Each parameter in $R_M$ could be found in funcR0.R for both species: $Aedes$ $albopictus$ and $Aedes$ $aegypti$. All of them depend on temperature, rainfall or human density.
In order to compute the Mosquito basic reproduction number, $R_M$, for a particular location you need average daily rainfall (mm), average temperature (Cº), and human density (number of habitants per km2). 
This will produce a number between 0 to 6.64 and 5.46 for $Aedes$ $albopictus$ and $Aedes$ $aegypti$ respectively. If this number is greater than 1 it means suitability for this species. This number represents the average number of female mosquitoes than one female mosquito will produce in her entire lifespan.

Example how to compute the Mosquito basic reproduction number in R:
```
source("~/RM_mosquito/funcR0.R") # Load the file with the functions
Te = 20 # Example of temperature in Cº
rain = 4 # Example of rainfall in mm
hum = 1000 # Example of human density in inhabitants per km2
R0_func_alb(Te,rain,hum) # R_M for Aedes albopictus
R0_func_aeg(Te,rain,hum) # R_M for Aedes aegypti
```
This number could be compute it at different scale but always taken the average rainfall. The units must be daily rainfall.
So if we do it for a month there should be the average monthly rainfall NOT the accumulated rainfall.
