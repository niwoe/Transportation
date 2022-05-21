# Transportation Economics

## Purpose
Simulate the evolution of the COVID-19 epidemic with one variant : demandSupplyList.m
Simulate the evolution of the COVID-19 epidemic with one variant : demandSupplyList2Variants.m


## Ecosystem
MATLAB trial version : (https://fr.mathworks.com/campaigns/products/trials.html?prodcode=ml)

## Usage
Open demandSupplyList.m or demandSupplyList2Variants.m with MATLAB
Run it, the Input window will open
The default values are already entered, if no changes are desired, press the "ok" button

## Possible issues
For the code with the 2 variants, a too big difference between the contagiousness or the number of contacts between the two variants can lead to stability problems. This is due to the fact that one of the variants will relatively quickly, compared to the necessary study period of the second variant, "go" close to 0. After a while this "0" will reach the precision limit of the computer and NaN values appear.
For instance : "probabilty of covid-19 transmission between an Infected and a Susceptible during a contact at peak infectiousness (Variant 2)" = 0.35 with all others parameters unchanged work fine. 

## Results
See the report of group 4

## License
Permissive

## Authors
Group 4, EPFL

## Acknowledgments
We would like to thank the teams of the course "Transportation Economics" for their teaching and for the advice they gave us which allowed the realization of this Laboratory
