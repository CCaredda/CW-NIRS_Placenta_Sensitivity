# A Digital Instrument Simulator Platform to Support the Development of Non-Invasive Optical NIR Device for placenta monitoring


Charly Caredda, Frédéric Lange, Niccole Ranaei-Zamani, Uzair Hakim, Olayinka Kowobari, Dimitrios Siassakos, Sara Hillman, Anna L. David, Subhabrata Mitra, Ilias Tachtsidis, "Digital instrument simulator platform to support the development of noninvasive optical NIR device for placenta monitoring," J. Biomed. Opt. 31(2) 027003 (20 February 2026) \url{https://doi.org/10.1117/1.JBO.31.2.027003}


## Light propagation with Monte Carlo simulations

The code is contained in the folder "Data_generation_MCX"
For data generation Install:
- the latest version of mcxlab (http://mcx.space/) for the computation of Monte Carlo simulations 
- Iso2mesh toolbox: https://iso2mesh.sourceforge.net/cgi-bin/index.cgi

Please change the path in matlab scripts

## Light propagation with Redbird analytical model

The code is contained in the folder "Data_generation_redbird"
For data generation Install:
- the latest version of redbird (https://github.com/fangq/redbird)
- Iso2mesh toolbox: https://iso2mesh.sourceforge.net/cgi-bin/index.cgi

Please change the path in matlab scripts

## Data processing
The code is contained in the folder "Processing"
This folder contained Python scripts for calculating placental sensitivity, detection probabilities, calibration procedures

## Comparison between Monte Carlo and Redbird simulations
The code is contained in the folder "Compare_MCX_Redbird"

## Functions
This folder contains matlab functions used by RedBird or MCX

## Spectra
This folder contains aborption spectra (water, fat) and molar extinction coefficients (hemoglobin, cytochromes).




