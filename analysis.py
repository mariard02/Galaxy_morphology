from galaxymorph import galaxy
import matplotlib.pyplot as plt

print("****************************************************************************************************")
print("M87")

file = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/m87.dat"

m87 = galaxy("M87", 16.4, file)
m87.galaxy_info()
m87.galaxy_fit1()

print("****************************************************************************************************")
print("M101")

file = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/m101.dat"

m101 = galaxy("M101", 6.4, file)
m101.galaxy_info()
m101.galaxy_fit1()

print("****************************************************************************************************")
print("M81")

file = "/Users/maria/Desktop/Máster/M1/S1/Galaxias/Morphology/Data/m81.dat"

m81 = galaxy("M81", 3.6, file)
m81.galaxy_info()
m81.galaxy_fit2()

plt.show()
