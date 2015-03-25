The files listed in this directory are to support the creation and testing of the XRS10.py algorithm.

figure1_EW py-idl - plot comparing EW locations from python and IDL algorithms

figure2_NS py-idl - plot comparing NS locations from python and IDL algorithms

figure3_X_flares - plot comparing EW and NS locations from python algorith with NOAA flare location for X-class flares

figure4_M_flares - plot comparing EW and NS locations from python algorith with NOAA flare location for M-class flares

figure5_C_flares - plot comparing EW and NS locations from python algorith with NOAA flare location for C-class flares

flarefile.dat - NOAA flare data information. Columns are flare start, flare peak, flare location and flare class

idl_results.dat - flare locations as determined by IDL routine. Columns are EW location and NS location.

README.txt - this file

Space Weather DR - XRS.10.pptx Powerpoint of delivery review presentation.

sxi_auxfile.dat - proxy SXI data. Columns are timestep, exposure time, roll angle, xoffset, yoffset, yaw_flip_flag, fov_eclipse and fov_offpoint. 

sxi_quadvals.dat - proxy quad-diode data derived from SXI images. Used to test XRS10.py. 5 columns, one for timestep and one for each quad-diode.

XRS10_output.txt - output from XRS10.py algorithm. Columns are time and flare location (either heliographic coordinate or position angle and radius).

XRS10_test.py - python routine used to test XRS10.py algorithm, compare the results with IDL algorithm results and with NOAA flare locations and to create the figures in this directory.

XRS10.py - the XRS10 science quality algorithm

XRS10.pyc - compiled form of XRS10 algorithm creatd by python.  

