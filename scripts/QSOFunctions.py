#!/usr/bin/env python3
import numpy as np
import copy
def toFlux(magnitude,zeroflux):
    return (zeroflux*(10**(-(magnitude)/(2.5))))
def toMagnitude(flux,zeroflux):
    return (-2.5)*np.log10(flux/zeroflux)
def CollectSpectroscopicData(target,band):
    wavelength,qso,host = np.loadtxt(f"results/PyQSOFit_test/{target}_data.csv",delimiter=",",unpack=True)
    GBandWavelength,GBandTransmissivity = np.loadtxt(f"data/filters/{band}_band.csv",skiprows=1,unpack=True)

    #first trim the spectrum data to only that within the bounds of the g band
    abswavedifferenceMin,abswavedifferenceMax = np.abs(np.subtract(wavelength,min(GBandWavelength))),np.abs(np.subtract(wavelength,max(GBandWavelength)))
    minIndex = np.where(abswavedifferenceMin == min(abswavedifferenceMin))[0][0]
    maxIndex = np.where(abswavedifferenceMax == min(abswavedifferenceMax))[0][0]

    wavelength = wavelength[minIndex:maxIndex]
    qso = qso[minIndex:maxIndex]
    host = host[minIndex:maxIndex]
    #calculate the galaxy fraction from these values:
    integratedqso = np.trapz(qso,wavelength)
    integratedhost = np.trapz(host,wavelength)
    spectroscopicGalaxyFraction = integratedhost/(integratedhost+integratedqso)
    #now iterate over the trimmed wavelength to find the nearest transmissivity
    interpolatedTransmissivity = []
    for i,value in enumerate(wavelength):
        absdifference = abs(np.subtract(GBandWavelength,value))
        nearestIndex = np.where(absdifference == min(absdifference))[0][0]
        #spectroscopic data has a higher resolution than filter data - so find out if the target wavelength is larger or smaller
        #then, find the next or previous point which the target value should lie between. Interpolate linearly for the average value
        if value < GBandWavelength[nearestIndex]:
            if i != 0:
                secondIndex = nearestIndex-1
            else:
                secondIndex = copy.deepcopy(nearestIndex)
        else:
            if i != (len(wavelength)-1):
                secondIndex = nearestIndex+1
            else:
                secondIndex = copy.deepcopy(nearestIndex)
                
        interpolatedTransmissivity.append(np.average([GBandTransmissivity[nearestIndex],GBandTransmissivity[secondIndex]]))

    #now adjust for transmissivity for what HSC would see
    transAdjustedFluxQSO = np.multiply(qso,interpolatedTransmissivity)
    transAdjustedFluxHost = np.multiply(host,interpolatedTransmissivity)

    #now need to integrate over the transadjusted specific flux - using trapezoid rule
    integratedQSOFlux = np.trapz(transAdjustedFluxQSO,wavelength)
    integratedHostFlux = np.trapz(transAdjustedFluxHost,wavelength)
    galaxyFraction = integratedHostFlux/(integratedHostFlux+integratedQSOFlux)
    #to find the magnitude, need to convert fluxes to photon counts
    photonCountsByWavelengthQSO = []
    photonCountsByWavelengthHost = []
    for i,x in enumerate(transAdjustedFluxQSO):
        energyPerPhotonSI = ((6.626e-34)*(3e8))/(wavelength[i]/1e10) # energy per photon in si units at this wavelength, converted wavelength to m from angstroms
        energyPerPhotonErgs = energyPerPhotonSI*1e7
        photonCountsByWavelengthQSO.append((transAdjustedFluxQSO[i]*10e-17)/energyPerPhotonErgs)
        photonCountsByWavelengthHost.append((transAdjustedFluxHost[i]*10e-17)/energyPerPhotonErgs)
    photonCountsQSO = np.trapz(photonCountsByWavelengthQSO,wavelength)
    photonCountsHost = np.trapz(photonCountsByWavelengthHost,wavelength)
    #currently in photons/s/cm^2
    fluxMag0 = 63095734448.01944
    QSOSpecMagnitude = toMagnitude(photonCountsQSO*10000,fluxMag0) #multiplying by 10,000 to go from photons/s/cm^2 to photons/s/m^2
    HostSpecMagnitude = toMagnitude(photonCountsHost*10000,fluxMag0)
    return {"host magnitude":HostSpecMagnitude,"qso magnitude":QSOSpecMagnitude,"galaxy fraction":galaxyFraction}
def ImageGalaxyFraction(pklFile):
    #returns galaxy fraction computed from magnitudes of host and qso
    fluxMag0 = 63095734448.01944
    HSCGalaxyMagnitude = pklFile.final_result_galaxy[0]["magnitude"]
    HSCQSOMagnitude = pklFile.final_result_ps[0]["magnitude"]
    HSCHostFlux = toFlux(HSCGalaxyMagnitude,fluxMag0)
    HSCQSOFlux = toFlux(HSCQSOMagnitude,fluxMag0)
    HSCGalaxyFraction = HSCHostFlux/(HSCHostFlux+HSCQSOFlux)
    return HSCGalaxyFraction
def MedianSN(specFilePath):
    flux,flux_err = np.loadtxt(specFilePath,delimiter=",",skiprows=1,usecols=(2,3),unpack=True)
    SN=[]
    for i, f in enumerate(flux):
        SN.append(f/flux_err[i])
    return np.median(np.array(SN))