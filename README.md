# RedQSOYoRiS
Ellie Saga's YoRiS project, testing the limitations of image decomposition using ground-based telescopes in preparation for LSST.

Included is the environment .yml file for recreating the environment used to run this code. A full list of packages is available in dependencies.txt

Importing the environment:
```
conda env create -f YoRiS_Environment.yml
```
Using the environment:
```
conda activate galightMain
```
To run the scripts, activate the environment as shown above and include an "images/" folder containing .fits files of fov and psf images. To run rgbImageGenerator.py, include a .csv with object names in column 0 and a 1 or 0 in column 1 indicating presence of contamination. To run decomposition.py, a csv is needed in the same directory including a target ID, right ascension and declination.
