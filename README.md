[![DOI](https://zenodo.org/badge/307830945.svg)](https://zenodo.org/badge/latestdoi/307830945)

# Source Code for

Luke V Colosi, Sarah T Gille, and Ana B Villas Bôas. The seasonal cycle of significant waveheight in the ocean: Local vs remote forcing. Journal of Geophysical Research: Oceans, submitted.

# Abstract

Significant wave height (SWH) stems from a combination of locally generated "wind-sea" and remotely generated "swell" waves. In the Northern and Southern Hemispheres, wave heights typically undergo a sinusoidal annual cycle, with larger SWH in winter in response to seasonal changes in high-latitude storm patterns that generate equatorward propagating swell.  However, some locations deviate from this hemispheric-scale seasonal pattern in SWH. For example, in the California coastal region, local wind events occur in boreal spring and summer, leading to a wind speed (WSP) annual cycle with a distinct maximum in boreal spring and a corresponding local response in SWH. Here ocean regions with a WSP annual cycle reaching a maximum in late spring, summer, or early fall are designated as seasonal wind anomaly regions (SWARs). Intra-annual variability of surface gravity waves is analyzed globally using two decades of satellite-derived SWH and WSP data. The phasing of the WSP annual cycle is used as a metric to identify SWARs. Global maps of probability of swell based on wave age confirm that during the spring and summer months, locally forced waves are more statistically more likely in SWARs than in surrounding regions. The magnitude of the deviation in the SWH annual cycle is determined by the exposure to swell and characteristics of the wave field within the region. Local winds have a more identifiable impact on Northern Hemisphere SWARs than on Southern Hemisphere SWARs due to the larger seasonality of Northern Hemisphere winds. 

# Plain Language Summary

At the ocean surface, wave height can give insight into ocean-atmosphere interactions. Storms generate waves, which are known as swell when they propagate away from their point of origin. Swell waves account for most of the global ocean's surface waves. They vary annually, with large waves in the winter and small waves in the summer, due to seasonal changes in high-latitude storm systems. In some coastal areas, including the coast of California, local wind effects cause exceptionally high wind speeds in late spring. These strong local winds result in large waves in springtime, separate from the global-scale winter maximum in swell waves. Places with strong local winds during the late spring, summer, and early fall, here referred to as seasonal wind anomaly regions (SWARs), are identified using global satellite observations of wave height and wind speed. Details vary by location. SWAR wave fields depend on the exposure to swell as well as the strength of the local winds. Compared with Southern Hemisphere storms, Northern Hemisphere storms have a stronger winter peak, which means that local winds have a larger influence in Northern Hemisphere SWARs than in Southern Hemisphere SWARs.

# Authors 
* [Luke Colosi](https://lcolosi.github.io/)<<lcolosi@ucsd.edu>>
* [Bia Villas Bôas](https://biavillasboas.github.io/) <<avillasboas@ucsd.edu>>
* [Sarah T. Gille](http://www-pord.ucsd.edu/~sgille/) <<sgille@ucsd.edu>>

# Data
CCMP Version-2.0 vector wind analysis data produced by Remote Sensing Systems is available at http://www.remss.com/measurements/ccmp/. Registration for an ftp account is required. The satellite altimetry significant wave height data produced by the French Research Institute for Exploitation of the Sea (IFREMER) is publically available at ftp://ftp.ifremer.fr/ifremer/cersat/products/swath/altimeters/waves. WAVEWATCH 3 hindcast significant wave height, surface wind velocity, and peak wave frequency is also produced by IFREMER and is publically available at  ftp://ftp.ifremer.fr/ifremer/ww3/HINDCAST. Intermediate data products including processed IFREMER, CCMP2, and WW3 data along with decorrelation time scales, weighted least-squares fit parameters and uncertainties, and probability of swell produced from CCMP2, IFREMER, and WW3 products are available through the [University of California, San Diego library](insert collection url).

# Funding
This work was supported by the NASA SWOT (awards NNX16AH67G and 80NSSC20K1136) and Ocean Vector Winds Science Teams (award 80NSSC19K0059), by a NASA Earth and Space Science Fellowship awarded to Ana Villas Bôas, and by the Hiestand Scholars program.

# How to use this repository

All figures in Colosi et al. (2021) can be reproduced using the Python scripts from this repository and the processed SWH, WSP and fp [data](insert collection url). To do so, follow these steps:

1. Make a local copy of this repository by either cloning or downloading it.

2. Download the processed SWH, WSP and fp [data](insert collection url), untar the files, and move all seven directories to `data` in the project root. After doing so, your directory tree should look like this:

```
WaveClimatology/
├── data
│   ├── ifremer_swh
│   ├── ccmp2_wsp
|   ├── ww3_swh
│   ├── ww3_wsp
│   ├── lsf_parameters
│   ├── decor_scales
│   └── prob_swell
├── figs
├── src
└── tools
```

3. Make sure that you create an environment with the package versions specified in `environment.yml`. If you are using [Conda](https://docs.conda.io/en/latest/) you can run 

`conda env create -f environment.yml`

from the project root to create the environment from the .yml file and run `conda activate waveclimate` to activate the environment.

4. If you follow the steps above you should be able to reproduce all figures, by running `python figXX.py` from the `src` directory without having to adjust any paths.

# How to cite this code

If you wish to use the code from this repository, you may cite it as: 

Colosi, Luke V. (2021, January 20). Source code for: 'The Seasonal Cycle of Significant Wave Height in the Ocean:  Local vs Remote Forcing'. Zenodo. (insert zenodo url)
