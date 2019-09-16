## About ACOLITE_MR
ACOLITE_MR contains the Dark Spectrum Fitting (DSF) atmospheric correction algorithms for aquatic applications of metre scale satellites developed at RBINS.

The algorithm was presented in Vanhellemont and Ruddick 2018, Atmospheric correction of metre-scale optical satellite data for inland and coastal water applications (https://www.sciencedirect.com/science/article/pii/S0034425718303481) and Vanhellemont 2019, Daily metre-scale mapping of water turbidity using CubeSat imagery (https://doi.org/10.1364/OE.27.0A1372).

ACOLITE_MR development was funded by the Belgian Science Policy Office STEREO program under contract SR/00/325 (PONDER project).

**ACOLITE is provided by RBINS as an experimental tool, without explicit or implied warranty. Use of the program is at your own discretion and risk.**

## Distribution
ACOLITE_MR is supported on the [ACOLITE forum](https://odnature.naturalsciences.be/remsem/acolite-forum/viewforum.php?f=15). 

## Dependencies
ACOLITE is coded in Python 3, and requires the following Python 
packages to run with all functionality:`to do!`

## Installation
* cd into a suitable directory and clone the git repository: `git clone 
https://github.com/acolite/acolite_mr`
* cd into the new acolite directory `cd acolite_mr`
* run `python acolite_mr.py --input $in --output $out`

## General use
At the least an inputfile $in and output directory $out has to be provided. $out is any path with write access, new directories will be created if necessary. 

The processor will output a NetCDF file with top-of-atmosphere reflectance (rhot_*), Rayleigh corrected reflectance (rhorc_*) and surface reflectance (rhos_*) in four bands with wavelengths derived from the band averaged spectral response function. Latitude and longitude are only  approximated from the scene bounding polygon vertex coordinates, and should be used with care. The NetCDF metadata contains the scene geometry and used atmospheric correction parameters (e.g. aerosol model, optical thickness, and derived path reflectance and transmittances per band.) For water applications the rhos_* data can be used as being rhow_* (rhow=pi*Rrs) although no sun glint correction has been applied.

For Pléiades $in is the full path to the extracted data obtained from Airbus. This is a bundle directory that typically contains a separate directory with MultiSpectral (MS) and/or Panchromatic (P) data. Some scenes are Pansharpened MultiSpectral (PMS). 

For example, bundle FCGC600220033 contains three directories (MS, P and LIBRARY) and three files:

/storage/Pleiades/FCGC600220033/IMG_PHR1A_MS_002
/storage/Pleiades/FCGC600220033/IMG_PHR1A_P_001
/storage/Pleiades/FCGC600220033/LIBRARY
/storage/Pleiades/FCGC600220033/DELIVERY.PDF
/storage/Pleiades/FCGC600220033/INDEX.HTM
/storage/Pleiades/FCGC600220033/VOL_PHR.XML

To process the full path to the bundle directory has to be specified as input file: /storage/Pleiades/FCGC600220033

SPOT bundles seem to have an extra level, and the directory containing the SPOT_VOL.XML file has to be provided. E.g. I have a SPOT7 scene, /storage/SPOT/ExampleScene, which contains:

/storage/SPOT/ExampleScene/LIBRARY
/storage/SPOT/ExampleScene/PROD_SPOT7_001
/storage/SPOT/ExampleScene/DELIVERY.PDF
/storage/SPOT/ExampleScene/INDEX.HTM
/storage/SPOT/ExampleScene/LICENSE.PDF
/storage/SPOT/ExampleScene/SPOT_LIST.XML

The directory /storage/SPOT/ExampleScene/PROD_SPOT7_001 contains:

/storage/SPOT/ExampleScene/PROD_SPOT7_001/LIBRARY
/storage/SPOT/ExampleScene/PROD_SPOT7_001/VOL_SPOT7_001_A
/storage/SPOT/ExampleScene/PROD_SPOT7_001/INDEX.HTM
/storage/SPOT/ExampleScene/PROD_SPOT7_001/SPOT_PROD.XML

and the directory /storage/SPOT/ExampleScene/PROD_SPOT7_001/VOL_SPOT7_001_A contains:

/storage/SPOT/ExampleScene/20150709/PROD_SPOT7_001/VOL_SPOT7_001_A/IMG_SPOT7_PMS_001_A
/storage/SPOT/ExampleScene/20150709/PROD_SPOT7_001/VOL_SPOT7_001_A/LIBRARY
/storage/SPOT/ExampleScene/20150709/PROD_SPOT7_001/VOL_SPOT7_001_A/INDEX.HTM
/storage/SPOT/ExampleScene/20150709/PROD_SPOT7_001/VOL_SPOT7_001_A/SPOT_VOL.XML

This level contains the SPOT_VOL.XML file and hence has to be passed as input file: /storage/SPOT/ExampleScene/20150709/PROD_SPOT7_001/VOL_SPOT7_001_A
