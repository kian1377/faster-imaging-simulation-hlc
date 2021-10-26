# faster-imaging-simulation-hlc

This repo is meant to demostrate how to use interpolation of PSFs to run image simulations for complex systems with field dependent PSFs. The example used is that of the Roman HLC imaging mode.

In order to use all the code within this repository, all required python modules should be installed. This includes the modules PROPER and wfirst_phaseb-proper which can be found at https://sourceforge.net/projects/proper-library/files/ and https://github.com/kian1393/proper-models/tree/master/wfirst_cgi/models_phaseb/python respectively. All the Phase-B data files required can be found at on the IPAC website https://roman.ipac.caltech.edu/. 

These files include the 1D off-axis PSFs from IPAC, however, the PROPER PSFs can be created using the given create_hlc_proper-psfs.ipynb file. 

The interpolated arrays of PSFs can then be created using the create-interpolated_psfs_array.ipynb file. 

Once those steps are done, the simulation can be run with the run_simulations.ipynb file. 

If you use this repository in research, we ask that you please cite the Milani et al 2019 SPIE conference proceedings which is freely available from the [arXiv](https://arxiv.org/abs/2106.09122) and the [UArizona institutional repository]( https://repository.arizona.edu/handle/10150/658138). Please also consider directly citing the repository via Zenodo, [![DOI](https://zenodo.org/badge/283349486.svg)](https://zenodo.org/badge/latestdoi/283349486).

