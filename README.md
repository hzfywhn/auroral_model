Spatial modeling of aurora using Lattice Kriging

Before running the code, you should have:
1. SSUSI auroral EDR data grouped by the satellite code, this can be found at https://ssusi.jhuapl.edu/data_availability?type=edr-aur
2. THEMIS data with derived energy flux and characteristic energy as Gabrielse et al. (2021), contact UCLA THEMIS group for the required data
3. Kp index required for the Zhang and Paxton (2008) model can be obtained from SPDF OMNIWeb (https://omniweb.gsfc.nasa.gov/)

Then run codes in the following sequence:
1. combine_ssusi, to combine all orbits from one single satellite (f17_20140220.nc)
2. prepare_input, to combine SSUSI and THEMIS auroral observations and produce formated input for modeling (north_in.nc)
3. downsampling, to adjust the weight of each data source useful for modeling (north_in.nc -> north_flux_in.nc)
4. auroral_model, to predict energy flux or mean energy from existing observations (north_flux_in.nc -> north_flux_out_5d.nc)

The final output north_flux_out_5d.nc is the one containing predicted auroral flux/energy at high latitudes