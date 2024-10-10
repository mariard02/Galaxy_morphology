from galaxymorph import galaxy
import matplotlib.pyplot as plt

z_f = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/multiwavelength/3379_z.dat"

z = galaxy("NGC 3379 (Z band)", 9.8, z_f)
z.galaxy_info()
z.galaxy_fit2(10, 2, 4, 10, 0.1, 1)

print("****************************************************************************************************")


u_f = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/multiwavelength/3379_u.dat"

u = galaxy("NGC 3379 (U band)", 9.8, u_f)
u.galaxy_info()
u.galaxy_fit2(0.3, 3, 4, 10, 0.1, 1)

print("****************************************************************************************************")

r_f = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/multiwavelength/3379_r.dat"

r = galaxy("NGC 3379 (R band)", 9.8, r_f)
r.galaxy_info()
r.galaxy_fit2(2, 3, 4, 20, 0.1, 1)

plt.show()