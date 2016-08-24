# Sparse VFA T<sub>1</sub> by OMP

Sparse MRI Acquisition Variable Flip Angle (VFA) T1 Mapping using an Orthogonal Matching Pursuit (OMP) algorithm. 

This code investigates the effect of Orthogonal Matching Pursuit compressed censing processing applied to noisy 2D VFA data 
and its effect on T<sub>1</sub> values. 

The algorithm uses complete dictionnaries; the dictionnary is NOT sparsified using a K-SVD algorithm on training data (but 
feel free to add this feature if you want/need it).

It currently uses a modified Shepp-Logan phantom, but the could should be easily modifyable to any T1 maps. The measurements 
magnitude are calculated from T1 maps, inverse fourier transformed to K-Space and retroactively undersampled prior to any
OMP reconstruction.

## Requirements

* MATLAB_R2012b or later

* [OMP Toolbox](http://www.cs.technion.ac.il/~ronrubin/software.html)

## Installation

The OMP Toolbox needs to be compiled before using the package. After downloading and  unzipping the file, open the directory
and follow the instructions in "readme.txt", and add the OMP Toolbox folders to the Matlab path.

Prior to compiling the OMP Toolbox, MEX needs to be setup in Matlab. This can be done using `mex -setup`. Visit the MATLAB
page on [setting up MEX](https://www.mathworks.com/help/matlab/matlab_external/changing-default-compiler.html) for more 
information.

## Usage

The demo script provides a detailed example of how to use the package. To run the demo, execute from the MATLAB command line:

`mainDemo`

The expected results of this demo are contained in the `html/` directory. If you produce similar results on your system, it
should be an indication that the software package was installed successfully.

## Attributions

Not all the code was written by the author of this repository.

1. utils/fft2c.m, utils/iff2c.m, utils/genPDF.m and utils/genSampling originate from [Michael Lustigs Sparse MRI software](http://www.eecs.berkeley.edu/~mlustig/Software.html)

2. utils/ricernd.m originates from [Ged Ridgway's Matlab Rician package](https://www.mathworks.com/matlabcentral/fileexchange/14237-ricerician-distribution)

The idea of reconstructing sparsely acquire T1 mapping MRI data using OMP was originally published in the following papers:

*Doneva, M., Börnert, P., Eggers, H., Stehning, C., Sénégas, J. and Mertins, A. (2010), Compressed sensing reconstruction for magnetic resonance parameter mapping. Magn Reson Med, 64: 1114–1120. doi: 10.1002/mrm.22483*

*Li, W., Griswold, M. and Yu, X. (2012), Fast cardiac T1 mapping in mice using a model-based compressed sensing method. Magn Reson Med, 68: 1127–1134. doi: 10.1002/mrm.23323*

## About me

**Mathieu Boudreau** is a PhD Candidate at McGill University in the Department of Biomedical Engineering.
He holds a BSc in Physics from the Universite de Moncton ('09), and a MSc in Physics from the University 
of Western Ontario ('11).