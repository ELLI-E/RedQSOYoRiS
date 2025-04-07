#Data Aggregator Script
#!/usr/bin/env python3
from attr import dataclass
from scripts.QSOFunctions import CollectSpectroscopicData
from scripts.QSOFunctions import ImageGalaxyFraction
from scripts.QSOFunctions import MedianSN
import pandas as pd
import pickle
import numpy as np

"""
need to combine:
key elements of pkl files
S/N - spectral
L5100
Key elements of spectral decomposition and image decomposition
Solution: function which defines an object containing all relevant information about a given target

Class should contain the pkl file
Contents of _data.csv, _cont.csv and _line.csv
Redshift
"""
class TargetData:
    
    
    def __init__(self,targetName,specDataDirectory,imgResultDirectory,specResultDirectory,catalogueAddress):
        #This will collect all relevant data about a particular object
        #Also collect spectral data using CollectSpectroscopicData
        self.pklFiles = {"G":None,"R":None}
        self.specFGal = {"G":None,"R":None}
        self.imgFGal = {"G":None,"R":None}
        self.specSN = None
        self.name = targetName
        #first collect redshift
        try:
            catalogue = pd.read_csv(catalogueAddress)
            self.objNo = list(catalogue["DESI_ID"]).index(int(self.name))
            self.z = list(catalogue["Z_DESI"])[self.objNo]
        except FileNotFoundError:
            print(f"Attempted to collect redshift data for {self.name}, but the catalogue file is missing")
        #then image data
        try:
            file = open(f"{imgResultDirectory}/{self.name}/fitting_results/{self.name}-result-band-g.pkl","rb")
            self.pklFiles["G"] = pickle.load(file)
            self.imgFGal["G"] = ImageGalaxyFraction(self.pklFiles["G"])
            file.close()
        except FileNotFoundError:
            print(f"Attempted to collect g band image decomposition data for {targetName}, but the file is missing")
        try:
            file  = open(f"{imgResultDirectory}/{self.name}/fitting_results/{self.name}-result-band-r.pkl","rb")
            self.pklFiles["R"] = pickle.load(file)
            self.imgFGal["R"] = ImageGalaxyFraction(self.pklFiles["R"])
            file.close()
        except FileNotFoundError:
            print(f"Attempted to collect r band image decomposition data for {self.name}, but the file is missing")
        #Now attempt to collect spectral data
        try:
            self.wavelength,self.qsoFlux,self.hostFlux = np.loadtxt(f"{specResultDirectory}/{self.name}_data.csv",delimiter=",",unpack=True)
            self.contValue = np.loadtxt(f"{specResultDirectory}/{self.name}_cont.csv",delimiter=",",usecols=(2),skiprows=1)
            self.specFGal["G"] = CollectSpectroscopicData(self.name,"g")
            self.specFGal["R"] = CollectSpectroscopicData(self.name,"r")
            self.specSN = MedianSN(f"{specDataDirectory}/{self.name}_spec.csv")
        except FileNotFoundError:
            print(f"Attempted to collect spectral decomposition data for {self.name}, but file/s are missing")