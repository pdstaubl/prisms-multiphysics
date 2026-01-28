import numpy as np

# Create a pseudo-2D bicrystal structure
# Domain size, in voxels
size_x = 64
size_y = 256
size_z = 1

# The grain ID values to use, as strings (must be strings! use quotes)
upper_grain_ID = "1"
lower_grain_ID = "2"

# Name for the file to create
# Warning! The file will be overwritten if it already exists!
output_filename = "grainID_bicrystal_2D.txt"


# Begin the calculation and file creation (do not edit below this point)
midpoint_y = size_y // 2

with open(output_filename, "w") as f:
  # Write header
  f.write("Grain ID file, pseudo-2D bicrystal, size {}x{}x{} voxels\n".format(size_x, size_y, size_z))
  for x in np.arange(size_x):
    for y in np.arange(size_y):
      row = []
      for z in np.arange(size_z):
        if y < midpoint_y:
          row.append(lower_grain_ID)
        else:
          row.append(upper_grain_ID)
      f.write(" ".join(row) + "\n")

