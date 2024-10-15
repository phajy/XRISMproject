using SpectralFitting
using Plots

# Load the XRISM early release data for Perseus cluster
DATADIR = "data"
spec_path = joinpath(DATADIR, "xa_merged_p0px1000_Hp.pi")
response_path = joinpath(DATADIR, "xa_merged_p0px1000_HpS.rmf")
ancillary_path = joinpath(DATADIR, "rsl_standard_GVclosed.arf")
data = OGIPDataset(spec_path, response=response_path, ancillary=ancillary_path)

# Plot the spectrum
# plot(data, xrange=[1,10], yrange=[0.1,300])
plot(data, xrange=[6.25,6.75])
