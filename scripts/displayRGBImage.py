#get a list of complete objects and display an rgb image for a given index
#imports
import os
from astropy.visualization import make_lupton_rgb
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
dataAddress = r"data/images/"

objects = list(os.listdir(dataAddress))

#list of complete objects
completeObjects = []

for object in objects:
    #we expect 10 files, 5 psfs and 5 images
    filesInObjectDirectory = len(list(os.listdir(f"{dataAddress}/{object}/")))
    if filesInObjectDirectory == 10:
        completeObjects.append(object)
    else:
        continue

#now, create colour images of each source and visually inspect for interfering foreground sources
index = 145
R = pyfits.open(f"{dataAddress}/{completeObjects[index]}/{completeObjects[index]}_HSC-I.fits")[1].data
G = pyfits.open(f"{dataAddress}/{completeObjects[index]}/{completeObjects[index]}_HSC-R.fits")[1].data
B = pyfits.open(f"{dataAddress}/{completeObjects[index]}/{completeObjects[index]}_HSC-G.fits")[1].data
image = make_lupton_rgb(R,G,B,stretch=0.5,Q=10)


plt.title(f"{completeObjects[index]}")
plt.imshow(image,origin="lower")
plt.show()