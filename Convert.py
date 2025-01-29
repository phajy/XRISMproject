from astropy.io import fits
import numpy as np
import os

# Define file paths
input_file = "xa_merged_p0px1000_Hp.pi"  # Replace with actual file path
output_file = "xa_merged_p0px1000_Hp.fits"

    # Check if the input file exists
if not os.path.exists(input_file):
    raise FileNotFoundError(f"Input file not found: {input_file}")

    # Load data from .pi file (assuming binary float32 format)
data = np.fromfile(input_file, dtype=np.float32)

    # Optional: Reshape the data if necessary (modify based on expected format)
    # Example: If data should be a 100x100 array, uncomment the next line and adjust
    # data = data.reshape((100, 100))  

    # Create a FITS file
hdu = fits.PrimaryHDU(data)
hdu.writeto(output_file, overwrite=True)

print(f"Successfully converted {input_file} to {output_file}")

