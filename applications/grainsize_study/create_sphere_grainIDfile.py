import numpy as np

# The grain ID values to use, as strings (must be strings! use quotes)
inner_grain_ID = "1"
outer_grain_ID = "2"

# Grid domain, in voxels
# Lui 2021 has domain of 120 um, so for 32 voxels, one voxel is 120/32 = 3.72 um
size_x = 32
size_y = 32
size_z = 32

# Spherical grain diameter, in voxels
#sphere_diameter = 24.0    # 90 um
sphere_diameter = 21.0       # ~78 um
#sphere_diameter = 21.3333333333333333333333  # 80 um
#sphere_diameter = 16.0    # 60 um

# Name for the file to create
# Warning! The file will be overwritten if it already exists!
output_filename = "grainID_sphere_{:.1f}diam.txt".format(sphere_diameter)


# Values for the calculation (do not edit below this point)
sphere_radius_squared = (sphere_diameter / 2)**2
center_x = size_x / 2
center_y = size_y / 2
center_z = size_z / 2

with open(output_filename, "w") as f:
  # Write header
  f.write("Grain ID file, sphere of diameter {:.6f} voxels\n".format(sphere_diameter))
  for x in np.arange(size_x):
    for y in np.arange(size_y):
      row = []
      for z in np.arange(size_z):
        dist_squared = (x - center_x)**2 + (y - center_y)**2 + (z - center_z)**2
        if dist_squared < sphere_radius_squared:
          row.append(inner_grain_ID)
        else:
          row.append(outer_grain_ID)
      f.write(" ".join(row) + "\n")

