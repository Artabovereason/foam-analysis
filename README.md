
# Foam Analysis to Detect Fibres effect on Foam Topology

The present scripts were developped during my summer internship at Institut Charles Sadron in the Mechanics of Interfaces and Multiphase Systems team (MIM). The aim of those codes is to find a quantity to describe the change of the Foam topology under the addition of fibres to the foam.


## Usage
* The script `FoamAnalysisFromCSV.m` require a `.csv` file that comes from the *Foam Analysis* software, the columns delimiter should be `;`, the comma for non-integers `.`. Simply run the script, you are then asked to select a `.csv` file and in a newly created folder where the `.csv` file is, the output figures will be into.

* The script `read_scale.m` require a `.tif` tomography picture with a visible scale. Simply make a two-point selection of the scale bar, where the first selection is the upper-right and the second selection is the bottom-left. You are then asked to input the written size of the scale (in mm). The corresponding resolution of the picture (in mm per pixel) is displayed in the terminal.

* The script `binarisation_fibers.m` require a folder in which all the `.tif` tomography pictures are into. It then create a three-dimensional reconstruction of the fibers position. You are required to manually select the fibres to *individualise* them. The way to individualise them is simply to make a two-point selection over three plane-cut of the three-dimensional representation. Once again, the first point is the upper-right and the second one is the bottom-left. The output will be a `.mat` file where a Matlab cell named `save_selected_cell` will have all the fibres individualised. There is also two times four figures, one with the fibres before the user-selection, and one after with all the possible interesting cuts : three-dimensional then three planes : xz, xy and yx. 

* The script `clean_up_fibers.m` require the fibres isolated Matlab file that was generated by `binarisation_fibers.m`. The aim of this script is to make the user clean-up manually all the possible noise that still exist within the fibres data, but one has to know that unless the sample is very noisy, the changes in results for fibres position is around 5%. For every fiber of the sample, the user is asked to give how many selection of noise area he wants to make, this has to be a positive integer value. Then it is required to make two-points selections again, the first point being the upper-right and teh second one the bottom-left. After this **lengthy** process the script will output a Matlab data file `fibers_cleaned_isolated.mat` where the cell `cleaned_selected_cell` is saved onto to make calculations latter on.

* The script `FoamAnalysisFromNodesAndStruts.m` require a `.mat` file from the Nodes and Struts extraction. This script will display some data on the terminal and output lot's of figures, saved onto a newly created folder with the same name as the input `.mat` file.
## Authors

- Pierre Guichard
- Jean Farago


## License
[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
