# Vercini_spt_analysis

 INFORMATION
 
This repository contains code involved in single-particle tracking analysis of vertically-oriented bacteria using VerCINI (Vertical Cell Imaging by Nanostructured Immobilisation). Much of it is associated with the publication: KD Whitley et al. Nature Microbiology (2024). It contains several functions and scripts used for different purposes. This software is provided as is and without warranty.

OVERVIEW

1. fit_gauss_sum.m: A MATLAB function that takes an array of speeds and fits them to a sum of independent Gaussian functions. It then calculates the expected values of the processive population.
2. plot_gauss_sum.m: A MATLAB function plots a histogram of speeds along with the associated fits to a sum of independent Gaussian functions (outputted from fit_gauss_sum.m).
3. fit_plot_gauss_sum.m: A simple MATLAB script that runs fit_gauss_sum.m and plot_gauss_sum.m one after the other. Open the file itself to add in your array of speeds and fitting parameters.
4. getCellLengths.m: A MATLAB function to measure the lengths of cells whose membranes have been stained with fluorescent dye, using the intensity profiles running lengthwise along them.
5. spt_tirf_trackmate_off_axis.m: A MATLAB function that measures the off-axis motion of single-molecule tracks along a defined axis and returns a structure variable containing information on the field of view.
6. spt_tirf_trackmate_off_axis_batch2.m: A MATLAB function that measures the off-axis motion of single-molecule tracks along a defined axis, using all files in a directory.
7. plot_spt_tirf_trackmate_off_axis2.m: A MATLAB function that plots the off-axis motion (sum of residuals to a linear fit) of single-molecule tracks, using output from spt_tirf_trackmate_off_axis.m or spt_tirf_trackmate_off_axis_batch2.m.

8. verciniSPT_v1.m: A MATLAB function that plots various things for single-particle tracking data in vertical cells. It plots tracks in polar coordinates, time traces in polar coordinates, and mean-squared displacements iwth fits. The output is an interactive figure, with several options for futher analysis.

SYSTEM REQUIREMENTS

This package does not require any special hardware. It can be performed on a standard computer.

This package development system has been tested on Windows 10, but should be compatible with other operating systems.

The package requires the following software:

- MATLAB (>=R2018b)
- MATLAB Image processing toolbox
- MATLAB Optimization toolbox

INSTALLATION

Add the 'Vercini_spt_analysis' directory to your MATLAB path, and save the updated path so that it will remain installed next time you start MATLAB. Typical install time should be a couple of minutes.

USAGE INSTRUCTIONS

To run the software, please see detailed help documentation in each function or script, which can be accessed by typing 'help [function or script name]' on the MATLAB command line. Each function should take between 10 s and 10 min to run on a 'standard' computer, depending on the size of the dataset and which function or script is being used.

LICENSING INFORMATION

All files are distributed under the GPLv3 and (c) 2024 Kevin Whitley, Newcastle University unless otherwise stated. See LICENSE.txt for full terms.
