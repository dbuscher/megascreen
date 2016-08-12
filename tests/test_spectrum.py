import scipy
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.table import Table
import time
import MegaScreen as ms

def interf_spectrum_quad(baseline, frequencies, r0, L0, eps1=1.0):
    return np.array([8 * (scipy.integrate.quad(lambda u, f, b, r0, L0:
                                    ms.VonKarmanSpectrum(np.sqrt(u**2+f**2), r0, L0) *
                                    (1.0 - np.cos(2 * np.pi * u * b)),
                                    0, eps1/baseline, args=(f, baseline, r0, L0), epsrel=1e-3,
                                    limit=400)[0]
                        + scipy.integrate.quad(lambda u, f, b, r0, L0:
                                    ms.VonKarmanSpectrum(np.sqrt(u**2+f**2), r0, L0),
                                    eps1/baseline, np.inf, args=(f, baseline, r0, L0), epsrel=1e-4,
                                    limit=100)[0]
                          -scipy.integrate.quad(lambda u, f, b, r0, L0:
                                    ms.VonKarmanSpectrum(np.sqrt(u**2+f**2), r0, L0),
                                    eps1/baseline, np.inf, args=(f, baseline, r0, L0), epsrel=1e-4,
                                    weight='cos',wvar=2*np.pi*baseline,
                                    limit=100)[0])

            for f in frequencies])


def SaveTable(results,prefix,args):
    for key in sorted(args):
        results.meta[key]=args[key]
    results.write(time.strftime("tmp/"+prefix+"%y%m%d-%H%M.dat"),
                format='ascii.ecsv')

def average_spectrum(screen):
    spectra=[]
    for line in screen:
        spectra.append(scipy.signal.periodogram(line,window="hanning"))
    spectra=np.array(spectra)
    return np.mean(spectra,axis=0)

def screen_spectra(generator,decimate,numScreen):
    spectra=[]
    for i in range(numScreen):
        screen = next(generator)
        spectra.append([average_spectrum(screen[:,::decimate]),
                        average_spectrum(screen.transpose()[:,::decimate])])
    spectra=np.array(spectra)
    return np.mean(spectra,axis=0)

def interf_spectrum(r0=5,L0=3000,baseline=250,step=32,numScreen=160000,nperseg=65536,freqOverlap=4.0,
                    nfftOuter=256,nfftInner=256,
                    prefix="interf_spec",plot=True,plotPhases=False):
    args=locals()
    args["MegaScreenVersion"]=ms.__version__
    g=ms.MegaScreen(r0, L0, windowShape=[step, 2], windowOrigins=[(0.0,0.0), (0.0, baseline)],
                 dx=step, frequencyOverlap=freqOverlap, nfftWoofer=nfftInner, nfftTweeter=nfftOuter)
    phase=np.concatenate([next(g) for i in range(numScreen)],axis=1)
    phase=phase[:,:,0]
    diff=phase[1]-phase[0]
    if plotPhases:
        plt.plot(phase[0])
        plt.plot(phase[1])
        plt.plot(diff)
        plt.show()
    results=Table(scipy.signal.welch(diff,nperseg=nperseg),names=["frequency","power"])
    SaveTable(results,prefix,args)
    results=Table(scipy.signal.welch(phase[1]+phase[0],nperseg=nperseg),names=["frequency","power"])
    SaveTable(results,prefix+"_add",args)
    results = Table(scipy.signal.welch(phase[0], nperseg=nperseg), names=["frequency", "power"])
    SaveTable(results, prefix + "_0", args)
    results = Table(scipy.signal.welch(phase[1], nperseg=nperseg), names=["frequency", "power"])
    SaveTable(results, prefix + "_1", args)
    if plot:
        plotit(results)
    return results


def component_spectrum(r0=10,L0=2000,step=32,numScreen=40000,which=2,nperseg=16384,
                       freqOverlap=4.0,fractionalSupport=1.0,
                       nfftOuter=256, nfftInner=256,
                       prefix="component_spec",plot=True):
    args=locals()
    args["MegaScreenVersion"]=ms.__version__
    g=ms.MegaScreen(r0, L0, windowShape=[step, 2], dx=step,
                 debug=True,
                 frequencyOverlap=freqOverlap, nfftWoofer=nfftInner, nfftTweeter=nfftOuter,
                 fractionalSupport=fractionalSupport)
    phase=np.concatenate([next(g)[which][:,0] for i in range(numScreen)])
    results=Table(scipy.signal.welch(phase,nperseg=nperseg),names=["frequency","power"])
    SaveTable(results,prefix+str(which)+"-",args)
    if plot:
        plotit(results)
    return results


def plotit(t):
    f,p=t["frequency"],t["power"]
    plt.loglog(f[1:],p[1:])
    plt.show()

# Generic command-line interface to run test function or give full control
if __name__ == '__main__':
    if len(sys.argv) == 1:
        interf_spectrum()
    else:
        exec(sys.argv[1])