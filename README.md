# faster-imaging-simulation-hlc

This repo is meant to demostrate how to use interpolation of PSFs to run image simulations for comlex systems with field dependent PSFs. The example used is that of the Roman HLC imaging mode.

In order to use all the code within this repository, all required python modules should be installed. This includes the modules PROPER and wfirst_phaseb-proper which can be found at https://sourceforge.net/projects/proper-library/files/ and https://github.com/kian1393/proper-models/tree/master/wfirst_cgi/models_phaseb/python respectively. 

These files include the 1D offaxis PSFs from IPAC, however, the PROPER PSFs can be created using the given create_hlc_proper-psfs.ipynb file. 

The interpolated arrays of PSFs can then be created using the create-interpolated_psfs_array.ipynb file. 

Once those steps are done, the simulation can be run with the run_simulations.ipynb file. 
