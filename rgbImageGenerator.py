#now saving images of all objects
import numpy as np
import matplotlib
from astropy.visualization import make_lupton_rgb
import astropy.io.fits as pyfits
from matplotlib import pyplot as plt
matplotlib.use('Agg')
dataAddress = r"images/"
imageAddress = r"rgbimages/"
contamCSV = r"contaminationCheck.csv"
objects,contaminated = np.loadtxt(contamCSV,unpack=True,delimiter=",",skiprows=1,dtype=str)

for i, object in enumerate(objects):
    if contaminated[i] == " 1":
        continue
    else:
        #make image for each object
        R = pyfits.open(f"{dataAddress}/{object}/{object}_HSC-I.fits")[1].data
        G = pyfits.open(f"{dataAddress}/{object}/{object}_HSC-R.fits")[1].data
        B = pyfits.open(f"{dataAddress}/{object}/{object}_HSC-G.fits")[1].data
        image = make_lupton_rgb(R,G,B,stretch = 0.5,Q=10)
        plt.imshow(image)
        plt.savefig(f"{imageAddress}/{object}.png")