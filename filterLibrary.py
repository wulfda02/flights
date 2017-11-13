# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:55:47 2014
All the definitions needed/used in filterAndFit
@author: dwulf
"""

import numpy as np
import matplotlib.pyplot as mpl
from scipy import signal as sig
import astropy.io.fits as fits
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
import time
import os
import h5py


##############################################################################
# Object Classes
##############################################################################

class fitParameters(object): # Object for storing fit numbers
    def __init__(self,templateLength,sampleForPulseArrival,segmentLength,
                 overlap):
        self.tmpltLen = templateLength
        self.plsArvl = sampleForPulseArrival
        self.segLen = segmentLength
        self.ovrlp = overlap
        self._params = np.array([], dtype=float) # fit parameters
        self._xPoints = np.array([], dtype=int) # fit points      
        self._badPulses = np.array([], dtype=int) # bad fit pulse numbers
    def setTemplateLength(self, templateLength):
        self.tmpltLen = templateLength
    def setSampleForPulseArrival(self, sampleForPulseArrival):
        self.plsArvl = sampleForPulseArrival
    def setSegmentLength(self, segmentLength):
        self.segLen = segmentLength
    def setOverlap(self, overlap):
        self.ovrlp = overlap
    def params(self):
        return self._params
    def numPulses(self):
        return len(self._params)/2
    def amplitudes(self):
        return self._params[::2]
    def positions(self):
        return self._params[1::2]+self.plsArvl
    def setXPoints(self):
        xP = np.array([],dtype=int)
        for pos in self.positions():
            plsLoc = int(np.around(pos))
            if (plsLoc>=0)and(plsLoc<self.segLen):
                # two xPoints, one for ph, one for pt
                xP = np.append(xP, [plsLoc, plsLoc+self.segLen]) 
            elif plsLoc<0:
                xP = np.append(xP,[0,self.segLen])
            elif plsLoc>=self.segLen:
                xP = np.append(xP,[self.segLen-1,2*self.segLen-1])
        self._xPoints = xP
    def addPulse(self, amplitude, position):
        self._params = np.append(self.params(),
                                 [amplitude,position-self.plsArvl])
        self.setXPoints()
    def setParams(self, params):
        self._params = params
        self.setXPoints()
    def xPoints(self):
        return self._xPoints
    def setBadPulses(self, badPulses):
        self._badPulses = np.unique(badPulses)
    def badPulses(self):
        return self._badPulses 
    def addBadPulse(self, pulseNumber):
        self._badPulses = np.unique(np.append(self.badPulses(), pulseNumber))
    def clearBadPulses(self):
        self._badPulses = np.array([], dtype=int) 
    def removeBadPulses(self): 
        rmIdx = np.append(2*self.badPulses(),2*self.badPulses()+1)
        self.setParams(np.delete(self.params(),rmIdx))
        self.clearBadPulses() 
    def clearParams(self): # clear stored parameters
        self._xPoints = np.array([], dtype=int)
        self._params = np.array([], dtype=float)
        self.clearBadPulses()
    def singlePulse(self, pulseNumber): # return params of selected pulse
        retIdx = np.append(2*pulseNumber,2*pulseNumber+1)
        return self.params()[retIdx]
    def complement(self, pulseNumber): # returns params for all pulses except one
        rmIdx = np.append(2*pulseNumber,2*pulseNumber+1)
        otherParams = np.delete(self.params(),rmIdx)
        return otherParams
    def lastStoredPulsePosition(self):
        return min(self.segLen-self.ovrlp+self.plsArvl,self.segLen)
    def storeParams(self): # return parameters for storing that won't be passed 
        retArray = np.array([],dtype=float)     
        for pulse in range(self.numPulses()):
            amp,plsStart = self.singlePulse(pulse)
            LHLimit = plsStart+self.tmpltLen>0
            RHLimit = plsStart+self.plsArvl<self.lastStoredPulsePosition()
            if LHLimit and RHLimit:
                retArray = np.append(retArray,[amp,plsStart])
        return retArray
    def passParams(self): # return parameters for passing to next segment
        retArray = np.array([],dtype=float)     
        for pulse in range(self.numPulses()):
            amp,plsStart = self.singlePulse(pulse)
            LHLimit = plsStart+self.plsArvl>=self.lastStoredPulsePosition()
            RHLimit = plsStart<self.segLen
            if LHLimit and RHLimit:
                retArray = np.append(retArray,[amp,plsStart])
        return retArray
    def copy(self,newParams=None): # copy object with choise of new parameters
        newInstance = fitParameters(self.tmpltLen,self.plsArvl,self.segLen,
                                    self.ovrlp)
        if newParams is not None:
            newInstance.setParams(newParams)
        else:
            newInstance.setParams(self.params())
        return newInstance
    def shiftParams(self): # shifts positions relative to next segment
        self._params[1::2] -= self.segLen-self.ovrlp
        self.setXPoints()
    def sort(self,order='position'): # sorts pulses keeping track of bad pulses
        if order=='position':
            sortKey = np.argsort(self.positions())
        elif order=='amplitude':
            sortKey = np.argsort(self.amplitudes())[::-1]
        else:
            print "Unrecognized Order"
            sortKey = np.arange(self.numPulses())
        badPos = self.positions()[self.badPulses()]
        newAmp = self._params[::2][sortKey]
        newPos = self._params[1::2][sortKey]
        newBad = np.where(badPos==newPos)[0]
        newParams = np.ones(2*self.numPulses())
        newParams[::2] = newAmp
        newParams[1::2] = newPos
        self.setParams(newParams)
        self.setBadPulses(newBad)
        
class fitSegements(object): # Object for storing data streams and residuals
    def __init__(self,filteredData,unfilteredData,mask=None,firstSample=0):
        self.segLen = len(unfilteredData)
        self._filtDat = filteredData
        self._filtRes = filteredData
        self._ufiltDat = unfilteredData
        self._ufiltRes = unfilteredData
        self._mask = mask
        self._i1 = firstSample
    def firstSample(self):
        return self._i1
    def segmentLength(self):
        return self.segLen
    def filteredData(self,which=None):
        if which is not None:
            i = int(which)
            retArr = self._filtDat[i*self.segLen:(i+1)*self.segLen]
        else: 
            retArr = self._filtDat
        return retArr
    def filteredResidual(self,which=None):
        if which is not None:
            i = int(which)
            retArr = self._filtRes[i*self.segLen:(i+1)*self.segLen]
        else: 
            retArr = self._filtRes
        return retArr
    def unfilteredData(self):
        return self._ufiltDat
    def unfilteredResidual(self):
        return self._ufiltRes
    def mask(self):
        return self._mask
    def setFilteredData(self,filteredData):
        self._filtDat = filteredData
        self.resetFilteredResidual()
    def setFilteredResidual(self,filteredResidual):
        self._filtRes = filteredResidual
    def setUnfilteredData(self,unfilteredData):
        self._ufiltDat = unfilteredData
        self.resetUnfilteredResidual()
    def setUnfilteredResidual(self,unfilteredResidual):
        self._ufiltRes = unfilteredResidual
    def setMask(self,mask):
        self._mask = mask
    def resetFilteredResidual(self):
        self._filtRes = self.filteredData()
    def resetUnfilteredResidual(self):
        self._ufiltRes = self.unfilteredData()
    def resetResiduals(self):
        self.resetFilteredResidual()
        self.resetUnfilteredResidual()
    def calcFilteredResidual(self,model):
        try:
            assert len(model)==len(self._filtDat)
            self._filtRes = self.filteredData() - model
        except AssertionError:
            print "Len Model Does Not Match Len Data"
    def calcUnfilteredResidual(self,model):
        try:
            assert len(model)==len(self._ufiltDat)
            self._ufiltRes = self.unfilteredData() - model
        except AssertionError:
            print "Len Model Does Not Match Len Data"
            
class filterObject(object): # object for making, storing, and using filters and templates
    def __init__(self,avgSignal,avgNoise,sampleRate):
        self.avgSig = avgSignal
        self.avgNos = avgNoise
        self.fs = sampleRate
        self.calEnergy = 3313.8
        self._analogNorm = np.abs(np.min(self.avgSig))
        self._phFilter = None
        self._ptFilter = None
        self._sfpa = None
        self._phNorm = None
        self._ptNorm = None
        self._plsTemplates = None
        self._phTemplates = None
        self._ptTemplates = None
        self.make() # fills in all the 'None' above
    def setCalenergy(self,calEnergy):
        self.calEnergy = calEnergy
        self.make()
    def setSampleRate(self,sampleRate):
        self.fs = sampleRate
        self.make()
    def filterLength(self):
        return len(self.avgSig)
    def templateLength(self):
        return 2*self.filterLength() - 1
    def frequency(self):
        df = self.fs/self.filterLength()
        return np.arange(len(self.avgNos))*df
    def NEP(self):  
        N = self.filterLength()
        df = self.fs/N
        ce = self.calEnergy*1.602E-19 # calenergy in joules
        n_f = np.sqrt(self.avgNos) # Vrms/sqrt(Hz)
        hcf = 0.1
        wndw = halfCosWindow(np.ones(N),hcf)  
        waf = np.sqrt(1./(np.mean(np.square(wndw)))) # window adjustment factor
        wndwdSig = halfCosWindow(self.avgSig, hcf)
        # Divide by sampling rate times energy in K k-alpha pulse in joules
        # Vrms / Hz * J == Vrms / W
        responsivity = (np.sqrt(2)*waf*np.fft.rfft(wndwdSig)/(N*df*ce))
        nep = n_f / responsivity
        return nep
    def bslnRes(self):
        df = self.fs/self.filterLength()
        nep = self.NEP()
        # Multiply by eV/J
        deltaErms = (1./np.sqrt(np.sum(4*df/np.square(np.abs(nep)))))*6.24E18 
        return 2.35482*deltaErms 
    def makeFilters(self):
        N = self.filterLength()
        df = self.fs/N
        n2_f = self.avgNos # Vrms^2/Hz
        hcf = 0.1
        wndw = halfCosWindow(np.ones(N),hcf)  
        waf = np.sqrt(1./(np.mean(np.square(wndw)))) # window adjustment factor
        wndwdSig = halfCosWindow(self.avgSig, hcf)
        s_f = np.sqrt(2)*waf*np.fft.rfft(wndwdSig)/(N*df) # Vrms/Hz
        # Equationb 28 of McCammon/Enss Cryogenic Particle Detection [Vrms^-1]
        phFilter_f = s_f.conjugate()/n2_f
        # Taken from Gwen's code, maximizes slope relative to noise [Hz / Vrms]
        ptFilter_f = -2*np.pi*np.complex(0,1)*self.frequency()*phFilter_f
        phFilter_t = np.fft.irfft(phFilter_f/np.sqrt(2))*N # V^-1
        ptFilter_t = np.fft.irfft(ptFilter_f/np.sqrt(2))*N/df # V^-1
        phFilter_t = np.roll(phFilter_t, N/2) 
        ptFilter_t = np.roll(ptFilter_t, N/2)
        self._phFilter = phFilter_t
        self._ptFilter = ptFilter_t 
    def filterAvgPulse(self,filt):
        hcf = .1
        wndwdSig = halfCosWindow(self.avgSig, hcf)/self._analogNorm
        fltrdSig = sig.fftconvolve(wndwdSig,filt,mode='full')
        assert(self.templateLength()==len(fltrdSig))
        return fltrdSig
    def phFilterTemplate(self):
        fltrdPls = self.filterAvgPulse(self._phFilter)
        self._phNorm = np.abs(np.max(fltrdPls))
        self._sfpa = np.argmax(fltrdPls)
        return fltrdPls/self._phNorm
    def ptFilterTemplate(self):
        fltrdPls = self.filterAvgPulse(self._ptFilter)
        self._ptNorm = np.abs(np.max(fltrdPls))
        return fltrdPls/self._ptNorm
    def pulseTemplate(self):
        tl = self.templateLength()
        fl = self.filterLength()
        sfpa = self._sfpa
        hcf = .1
        wndwdPls = halfCosWindow(self.avgSig, hcf)/self._analogNorm
        paddedPulse = np.append(wndwdPls,np.zeros(tl-fl))
        paddedPulse = np.roll(paddedPulse,sfpa-fl/4)
        return paddedPulse
    def sampleForPulseArrival(self):
        return self._sfpa
    def tmpltDict(self,tmplt,precision=10):
        sfpa = self._sfpa
        templateDict = {}
        for i in range(0,precision+1):
            decimal = float(i)*1.0/precision
            shiftedTemplate = fftShift(tmplt,decimal,sfpa)
            templateDict[decimal] = shiftedTemplate
        return templateDict 
    def make(self):
        self.makeFilters()
        self._phTemplates = self.tmpltDict(self.phFilterTemplate())
        self._ptTemplates = self.tmpltDict(self.ptFilterTemplate())
        self._plsTemplates = self.tmpltDict(self.pulseTemplate())
        print "Baseline Resolution: %.3f eV FWHM" % self.bslnRes()
    def normalizeData(self,data):
        return self.calEnergy*data/self._analogNorm
    def convolveData(self,data,filt):
        fl = self.filterLength()
        return longConvolution(data,filt,fl/4)
    def convolvePhFilter(self,data):
        return self.convolveData(data,self._phFilter)/self._phNorm
    def convolvePtFilter(self,data):
        return self.convolveData(data,self._ptFilter)/self._ptNorm
    def filtList(self):
        return [self._phTemplates,self._ptTemplates]
    def plsList(self):
        return [self._plsTemplates]
        
class fitter(object): # Details of fitting stored here
    def __init__(self,filterObj):
        self.filtObj = filterObj
        self._sfpa = self.filtObj.sampleForPulseArrival()
        self._tmpLen = self.filtObj.templateLength()
        self._ovrlp = int(1.1*self._tmpLen)
        self._mxTry = 10 # max number of fitting attemps per segment
        self._initThrshld = 2500. # approx. eV of initial threshold
        self._redFctr = 0.1 # should be >.06 to avoid ringing 
        self._nIter = 3
        self._minSep = 100 # minimum separation for finding new pulses
        self._dx = 20 # maximum difference between guess and fit sample
        self._ds = 0.05 # maximum difference between pulse shape parameter and 1
        self._i1 = 0
        self.segObj = None
        self._segLen = None
        self._mxPls = None # maw pulses per segement
        self.paramsObj = fitParameters(self._tmpLen, self._sfpa, self._segLen, 
                                       self._ovrlp)  
        self._currSegGood = False
        self._prvSegGood = False
    def lastStoredPulsePosition(self):
        return self.paramsObj.lastStoredPulsePosition()
    def sampleForPulseArrival(self):
        return self._sfpa
    def overlap(self):
        return self._ovrlp
    def setOverlap(self,overlap):
        self._ovrlp = int(overlap)
        self.paramsObj.setOverlap(overlap)
    def resetOverlap(self):
        self.setOverlap(int(1.1*self._tmpLen))
    def templateLength(self):
        return self._tmpLen
    def segmentLength(self):
        return self._segLen
    def params(self):
        return self.paramsObj.params()
    def xPoints(self):
        return self.paramsObj.xPoints()
    def numPulses(self):
        return self.paramsObj.numPulses()
    def fitSuccess(self):
        return self._currSegGood
    def previousFitSuccess(self):
        return self._prvSegGood
    def setSegmentLength(self,segmentLength):
        self._segLen = segmentLength
        self._mxPls = int(8*self._segLen/float(self._tmpLen)) 
        self.paramsObj.setSegmentLength(self._segLen)
    def setSegment(self,segmentObj):
        self.segObj = segmentObj
        self.setSegmentLength(self.segObj.segmentLength())
        self._i1 = self.segObj.firstSample()
    def setParams(self,params):
        self.paramsObj.setParams(params)
    def mask(self):
        return self.segObj.mask()
    def filteredData(self):
        return self.segObj.filteredData()
    def unfilteredData(self):
        return self.segObj.unfilteredData()
    def filteredResidual(self):
        return self.segObj.filteredResidual()
    def unfilteredResidual(self):
        return self.segObj.unfilteredResidual()
    def updateResiduals(self):
        self.segObj.calcFilteredResidual(self.filteredModel())
        self.segObj.calcUnfilteredResidual(self.unfilteredModel())
    def setDataToResidual(self):
        self.segObj.setFilteredData(self.filteredResidual())
        self.segObj.setUnfilteredData(self.unfilteredResidual())
    def copy(self,filterObj=None):
        if filterObj is None:
            return fitter(self.filtObj)
        else:
            return fitter(filterObj)
    def copyParams(self,params=None):
        if params is None:
            return self.paramsObj.copy()
        else:
            return self.paramsObj.copy(params)    
    def nextSegment(self,segmentObj):
        storedParams = self.copyParams(self.paramsObj.storeParams())
        passedParams = self.copyParams(self.paramsObj.passParams())
        storedParams.shiftParams()
        passedParams.shiftParams()
        self.setSegment(segmentObj)
        self.setParams(storedParams.params())
        self.updateResiduals()
        self.setDataToResidual()
        self.setParams(passedParams.params())
        self.updateResiduals()
        self._prvSegGood = self._currSegGood
        self._currSegGood = False
    def firstSample(self):
        return self._i1
    def nextSegmentStart(self):
        return self._i1 + self._segLen - self._ovrlp
    def threshold(self,iteration):
        return self._initThrshld*(self._redFctr**iteration)
    def addPulses(self,iteration):
        trig = self.threshold(iteration)
        residual = self.filteredResidual()
        i1 = 0
        if iteration<self._nIter-1:
            i2 = self._segLen
        else:
            i2 = self._segLen - self._sfpa
        maxLoc = findMultiMax(trig,residual[i1:i2],self._minSep)
        for j in maxLoc:
            ampGuess = residual[j+i1]
            self.paramsObj.addPulse(ampGuess,j+i1)
    def shapes(self,p=None):
        if p is None:
            pObj = self.paramsObj
        else:
            pObj = self.paramsObj.copy(p)
        sps = np.array([],dtype=float)
        for i in range(pObj.numPulses()):
            isoPls = (self.unfilteredData() 
                      -dataFunction(self._segLen,self.filtObj.plsList(),
                                    *pObj.complement(i)))
            sglPlsShp = dataFunction(self._segLen,self.filtObj.plsList(),
                                     *pObj.singlePulse(i))
            p0 = int(max(0,pObj.singlePulse(i)[1]))
            p1 = min(self._segLen,p0+self._tmpLen)
            sp = innerProduct(isoPls[p0:p1],sglPlsShp[p0:p1])
            sps = np.append(sps,sp) 
        return sps
    def badShape(self,p):
        return np.where(np.abs(self.shapes(p)-1)>self._ds)[0]
    def badPos(self,p):
        return np.where((np.abs(p[1::2]-self.params()[1::2])>self._dx))[0]
    def badAmp(self,p,i):
        return np.where(p[::2]<self.threshold(i))[0]
    def badPulses(self,p,i):
        badPulses = np.block([self.badAmp(p,i),self.badPos(p),self.badShape(p)])
        self.paramsObj.setBadPulses(badPulses) 
    def filteredModel(self):
        return dataFunction(self._segLen,self.filtObj.filtList(),
                            *self.params())
    def unfilteredModel(self):
        return dataFunction(self._segLen,self.filtObj.plsList(),
                            *self.params())
    def clearFit(self):
        self.paramsObj.clearParams()
        self.segObj.resetResiduals()
        self._currSegGood= False
    def lowResFit(self):
        i2 = self.lastStoredPulsePosition()        
        lowResPos = trigger(self.unfilteredResidual()[:i2])
        lowResPos = lowResPos[np.logical_not(self.mask())[lowResPos]]
        lowResAmp = pulseHeights(self.unfilteredResidual()[:i2],lowResPos)
        params = np.ones(2*len(lowResPos))
        params[::2] = lowResAmp
        params[1::2] = lowResPos
        return params
    def refit(self):
        fitFailed = False
        f = lambda x, *p: fitFunction(x, self._segLen, 
                                      self.filtObj.filtList(), *p)
        if self.numPulses()>0:
                fitx = self.xPoints()
                fity = self.filteredData()[fitx]
                guess = self.params()
                try:
                    popt,pcov = curve_fit(f, fitx, fity, p0=guess, xtol=1e-5)
                    self.setParams(popt)
                    self.updateResiduals()
                    self._currSegGood = True
                except RuntimeError:
                    fitFailed = True
                    self.clearFit()
        return fitFailed
    def fit(self):
        f = lambda x, *p: fitFunction(x, self._segLen, 
                                      self.filtObj.filtList(), *p) 
        i = 0 
        n = 0
        fitFailed = False
        self.addPulses(i)
        while i<self._nIter:
            if (n>self._mxTry) or (self.numPulses()>self._mxPls) or fitFailed:
                self.clearFit()
                break
            if self.numPulses()>0:
                fitx = self.xPoints()
                fity = self.filteredData()[fitx]
                guess = self.params()
                try:
                    popt,pcov = curve_fit(f, fitx, fity, p0=guess, xtol=1e-5)
                    self.badPulses(popt,i)
                    self.setParams(popt)
                except RuntimeError:
                    fitFailed = True
            if fitFailed:
                print "Fit Failed"
                pass
            elif len(self.paramsObj.badPulses())>0:
                self.paramsObj.removeBadPulses()
                n+=1
            elif i<self._nIter-1:
                i+=1
                n+=1
                self.updateResiduals()
                self.addPulses(i)
            else:
                self.updateResiduals()
                self._currSegGood = True
                i+=1
    def activeFit(self):
        stayInteractive = True
        mag=100
        mpl.figure(num=None, figsize=(16, 8), dpi=80)
        mpl.xlim([0,self._segLen])
        mpl.ylim([-90,90])
        threshold = self.threshold(self._nIter-1)
        mpl.plot(threshold*self.mask(),ls=':')
        mpl.axhline(0)
        mpl.axvline(self._sfpa)
        mpl.axvline(self.lastStoredPulsePosition())
        mpl.ylabel('Approx. eV (x%d)' % mag)
        mpl.xlabel('Sample')
        mpl.plot(self.filteredData()/mag+40) 
        mpl.plot(self.filteredModel()/mag+40) 
        fitx = self.xPoints()[::2]
        mpl.scatter(fitx,self.filteredModel()[fitx]/mag+40)
        loResP = self.lowResFit()
        loResX = loResP[1::2].astype(int)
        loResE = 0.5*loResP[::2]
        loResY = self.unfilteredData()[loResX]-loResE
        mpl.errorbar(loResX,loResY/mag-40,loResE/mag,fmt='.')                                      
        mpl.plot(self.unfilteredData()/mag-40)  
        mpl.plot(self.unfilteredModel()/mag-40)
        mpl.plot(self.filteredResidual()) 
        mpl.scatter(fitx,self.filteredResidual()[fitx])
        a = 0
        for b,c in zip(fitx,self.filteredResidual()[fitx]):
            mpl.annotate('%s' % str(a), (b,c), fontsize=14)
            a += 1
        mpl.ion()
        mpl.show()
        action = raw_input("Next Action:").lower()
        if action == 'a':
            print "Click on Pulse Position." 
            newPulsePosition = int(mpl.ginput(1)[0][0]) 
            if self.mask()[newPulsePosition]:
                ampGuess = self.filteredResidual()[newPulsePosition]
                self.paramsObj.addPulse(ampGuess,newPulsePosition)
            else:
                print "Unable To Add Pulse During Deadtime"
            mpl.close()
            fitFailed = self.refit()
            if not fitFailed:
                self.activeFit()
        elif action == 'r':
            badPulse = -999
            while (badPulse<0)or(badPulse>=self.numPulses()):
                while True:
                    try:
                        badPulse = raw_input("Enter Pulse Number:")
                        badPulse = int(badPulse)
                        break
                    except ValueError:
                        print "ERROR: INVALID INPUT."
                        print "ENTER NUMERIC VALUE"
                if (badPulse<0)or(badPulse>=self.numPulses()):
                    msg = "ERROR: PULSE NUMBER %d" % badPulse 
                    msg += "OUT OF RANGE FOR %d PULSES" % self.numPulses()
                    print msg
            self.paramsObj.addBadPulse(badPulse)                     
            self.paramsObj.removeBadPulses() 
            mpl.close()
            fitFailed = self.refit()
            if not fitFailed:
                self.activeFit()  
        elif action == 'p':
            print "Pulse\tAmplitude\tPosition"
            for pulse in range(self.numPulses()): 
                print ("%d\t%.1f\t\t%.1f" 
                       % (pulse,self.paramsObj.amplitudes()[pulse],
                          self.paramsObj.positions()[pulse]))  
            mpl.close()
            self.activeFit()
        elif action == 'c':
            mpl.close()
        elif action == 'm':
            mpl.close()
            stayInteractive = False
        else:
            if action == 'h':
                print "Help Menu"
            else:
                print "ERROR: UNRECOGNIZED OPTION" 
            print "a: add pulse, r: remove pulse"
            print "c: continue"
            print "p:  print fit parameters"
            print "m: toggle out of interactive mode"
        return stayInteractive
    def storeParams(self):
        self.paramsObj.sort()
        return self.paramsObj.storeParams()
              
class singlePixel(object):
    def __init__(self,run,pixel,interactive=True):
        self.run = run
        self.pixel = pixel
        self.interactive = interactive
        self._sampleRate = None 
        self._T0Sample = None 
        self._avgPulse = None
        self._noise2 = None
        self._data = None
        self._mask = None
        self._normData = None
        self._phData = None
        self._ptData = None
        self._pulses = self.resetPulses()
    def setData(self,data):
        self._data = data
    def dataInVolts(self):
        scale = 13.33/4096.
        return scale*self._data
    def setSampleRate(self,sampleRate):
        self._sampleRate = sampleRate
    def setT0Sample(self,T0Sample):
        self._T0Sample = T0Sample
    def setAvgPulse(self,avgPulse):
        self._avgPulse = avgPulse
    def setNoise2(self,noise2):
        self._noise2 = noise2
    def setPulses(self,pulses):
        self._pulses = pulses
    def avgPulse(self):
        return self._avgPulse
    def noise2(self):
        return self._noise2
    def setMask(self,mask):
        self._mask = mask
    def mask(self):
        return self._mask
    def pulses(self):
        return self._pulses
    def fileName(self, suffix):
        fileName = '%s.%s' % (self.run, suffix)
        return fileName
    def sampleRate(self):
        return self._sampleRate
    def T0Sample(self):
        return self._T0Sample
    def TtimeForSample(self, sampleNumber):
        return (sampleNumber-self.T0Sample()) / self.sampleRate()
    def sampleForTTime(self, TTime):
        return int(round(TTime * self.sampleRate() + self.T0Sample()))
    def dt(self):
        return 1./self.sampleRate()
    def loadFromHdf(self, suffix):
        fn = self.fileName('%s.hdf5' % suffix)
        f = h5py.File(fn,'r') 
        group = f['p%02d' % self.pixel]
        if 'sampleRate' in group.attrs:
            self.setSampleRate(group.attrs['sampleRate'])
        else:
            self.setSampleRate(1./96E-6)
        if 'T0Sample' in group.attrs:
            self.setT0Sample(group.attrs['T0Sample'])
        if 'data' in group.keys():
            self.setData(np.copy(group['data']))
        if 'avgPulse_t' in group.keys():
            self.setAvgPulse(np.copy(group['avgPulse_t']))
        if 'noise2_f' in group.keys():
            self.setNoise2(np.copy(group['noise2_f']))
        if 'mask' in group.keys():
            self.setMask(np.copy(group['mask']))
        elif self._data is not None:
            self.setMask(np.ones(len(self._data),dtype=bool))
        if 'eventList' in group.keys():
            self.setPulses(np.copy(group['eventList']))
        f.close()
    def writeHdfTiming(self, suffix): # intednded to modify existing file
        fn = self.fileName('%s.hdf5' % suffix)
        f = h5py.File(fn,'r+')
        group = f['p%02d' % self.pixel]
        if 'sampleRate' in group.attrs:
            group.attrs['sampleRate'] = self.sampleRate()
        else:
            group.attrs.create('sampleRate', self.sampleRate())
        if 'T0Sample' in group.attrs:
            group.attrs['T0Sample'] = self.T0Sample()
        else:
            group.attrs.create('T0Sample', self.T0Sample())
        f.close()
        return fn
    def writeHdfFilter(self, suffix): # intednded to modify existing file
        fn = self.fileName('%s.hdf5' % suffix)
        f = h5py.File(fn,'r+')
        group =  f['p%02d' % self.pixel]
        if 'avgPulse_t' in group.keys():
            del group['avgPulse_t']
        group.create_dataset('avgPulse_t',data=self.avgPulse())
        if 'noise2_f' in group.keys():
            del group['noise2_f']
        group.create_dataset('noise2_f',data=self.noise2())
        f.close()
        return fn
    def writeHdfMask(self, suffix):
        fn = self.fileName('%s.hdf5' % suffix)
        f = h5py.File(fn,'r+')
        group =  f['p%02d' % self.pixel]
        if 'mask' in group.keys():
            del group['mask']
        group.create_dataset('mask',data=self.mask(),dtype=bool)
        f.close()
        return fn
    def writeHdfEvents(self, suffix):
        fn = self.fileName('%s.hdf5' % suffix)
        f = h5py.File(fn,'r+')
        group = f['p%02d' % self.pixel]
        if 'pulseHeader' in group.attrs:
            group.attrs['pulseHeader'] = self.dtypeLabels()
        else:
            group.attrs.create('pulseHeader', self.dtypeLabels())
        if 'eventList' in group.keys():
            del group['eventList']
        group.create_dataset('eventList',data=self.pulses(),dtype=self.dtype())
        f.close()
        return fn  
    def writeHdfFit(self, suffix):
        self.writeHdfMask(suffix)
        self.writeHdfFilter(suffix)
        self.writeHdfEvents(suffix)
    def findSegInData(self,segment,sampleGuess):
        di = 0    
        iStart = sampleGuess + di
        segLen = len(segment)
        data =  self.dataInVolts()
        dataLen = len(data)
        segFound = False
        while (np.abs(iStart-sampleGuess)<1E6)and(not segFound):
            firstValue = data[iStart]
            offset = segment[0] - firstValue
            differences = np.sort(np.square(segment - offset 
                                            - data[iStart:iStart+segLen]))
            difference = np.sum(differences[:-5])
            if difference>1E-6:
                if di<=0:
                    di = -1*di + 1
                else:
                    di = -1*di
                if sampleGuess+di<0:
                    di = -1*di + 1
                if sampleGuess+di>dataLen-segLen:
                    di = -1*di
                iStart = sampleGuess + di
            else:
                segFound = True
        if segFound:
            return iStart
        else:
            return -1
    def timingFromSwp(self,T0=0):
        fn = "%s_swp.fits" % self.run
        if os.path.isfile(fn):
            swpFile = fits.open(fn)
            pixelDataHDU = swpFile['Pixel data']
            swpIdx = (pixelDataHDU.data['Pixel_Number']==self.pixel)
            fsGuess = 10416.7
            lenTrace = 256
            j = 0
            attempts = 0
            pulseFound = False
            while (not pulseFound) and (attempts<3):
                while pixelDataHDU.data['Flag'][swpIdx][j]>6:
                    j += 1
                p1 = pixelDataHDU.data['Pulse_Data'][swpIdx][j][:lenTrace]
                t1 = pixelDataHDU.data['Pulse_Time'][swpIdx][j] 
                s1 = 0
                s1 = self.findSegInData(p1,s1)
                if s1>0:
                    pulseFound = True
                else:
                    j += 1
                    attempts += 1
            if pulseFound:
                j = -1
                attempts = 0
                pulseFound = False
                while (not pulseFound) and (attempts<3):
                    while pixelDataHDU.data['Flag'][swpIdx][j]>6:
                        j -= 1
                    p2 = pixelDataHDU.data['Pulse_Data'][swpIdx][j][:lenTrace]
                    t2 = pixelDataHDU.data['Pulse_Time'][swpIdx][j]
                    s2 = int(s1 + (t2-t1)*fsGuess)
                    s2 = self.findSegInData(p2,s2)
                    if s2>0:
                        pulseFound = True
                    else:
                        j += 1
                        attempts += 1
                if pulseFound:
                    fs = (s2 - s1)/(t2 - t1)
                    T0Sample = s1 - (t1-T0)*fs
                    self.setSampleRate(fs)
                    self.setT0Sample(T0Sample) 
    def getTiming(self):
        useExisting = True
        haveTiming = bool(self.T0Sample())
        if haveTiming and self.interactive:
            print "Timing Information Already In .bsn File"
            useExisting = not raw_input("Use existing? (Y/n):")[:1].lower()=='n'
        if ((not haveTiming)or(not useExisting))and(self.haveData()):
            self.timingFromSwp()
        haveTiming = bool(self.T0Sample())
        if haveTiming:
            print "Timing Information Loaded. Sample Rate: %.1f" % self.sampleRate()
    def applyTiming(self):
        haveTiming = bool(self.T0Sample())
        if haveTiming:
            self._pulses['time'] = self.TtimeForSample(self._pulses['sample'])
        else:
            print "Timing information not available."
    def haveFilter(self):
        return (self._avgPulse is not None) and (self._noise2 is not None)
    def haveData(self):
        haveData = False
        if self._data is not None:
            haveData = bool(len(self._data))
            return haveData
    def goodFilter(self):
        return np.any(self._avgPulse) and np.any(self._noise2)
    def makeFilter(self,filterLength):
        minutesOfCalData = 30
        i0 = 0
        i1 = min(int(minutesOfCalData*60*self._sampleRate),len(self._data))
        calData = self.dataInVolts()[i0:i1]
        print "Triggering Pulses..."
        signalLoc,noiseLoc = signalAndNoiseLoc(calData,filterLength)
        print "Selecting Cal Pulses..."        
        goodCals = selectCals(calData,signalLoc,filterLength,self.sampleRate(),self.interactive)
        self.setAvgPulse(avgSignal(calData,goodCals,filterLength))
        self.setNoise2(avgNoise(calData,noiseLoc,filterLength,self.sampleRate()))
    def getFilter(self):
        useExisting = True
        defaultLength = True
        filterLength = 4096
        if self.haveFilter() and self.interactive:
            print "Filter Information Already In .bsn File"
            useExisting = not raw_input("Use existing? (Y/n):")[:1].lower()=='n'           
        if ((not self.haveFilter())or(not useExisting))and(self.haveData()):
            if self.interactive:
                defaultLength = not raw_input("Use default filter length of %d? (Y/n):" % filterLength)[:1].lower()=='n'
            if not defaultLength:
                while True:
                    try:
                        filterLength = int(raw_input("Enter filter length in samples:"))
                        break
                    except ValueError:
                        print "ERROR: INVALID INPUT. ENTER INTEGER"
            self.makeFilter(filterLength)
        if self.goodFilter():
            fo = filterObject(self._avgPulse,self._noise2,self.sampleRate())
            if self.interactive:
                plotSignalNoise(fo.avgSig,fo.avgNos,fo.fs)
                plotFilters(fo._phFilter,fo._ptFilter,fo.fs)
                plotFilteredPulses(fo._phTemplates[0],fo._ptTemplates[0],fo.fs)
            return fo
        else:
            return None
    def deadTime(self):
        return 100*(1-np.sum(self._mask)/float(len(self._mask)))
    def filterData(self,filterObj):
        if self.haveData():
            satPls = findSaturatedPulses(self._data,'minimum')
            normData = filterObj.normalizeData(self.dataInVolts())
            self._normData = normData
            print "Removing Saturated Pulses..."
            modData, mask = removeSaturatedPulses(normData,satPls)
            self._mask = mask  
            print  "%.1f%% Deadtime From Saturated Pulses"% self.deadTime()
            print "Convolving Data With Filters..."
            self._phData = filterObj.convolvePhFilter(modData)
            self._ptData = filterObj.convolvePtFilter(modData)
    def dtype(self):
        return self._pulses.dtype
    def dtypeLabels(self):
        return np.array_str(np.array(self.dtype().names))[1:-1]
    def resetPulses(self):
        dtype1 = [('sample', '<f4'), ('ph', '<f4'), ('pixel', '<i2'), 
                  ('shape', '<f4'), ('resolution', '<i2'), ('time', '<f4'),
                  ('gain', '<f4'), ('lin', '<f4'), ('quad', '<f4'), 
                  ('energy','<f4')]
        self._pulses = np.array([], dtype=dtype1)
    def addPulseRow(self):
        newRow = np.array([tuple(np.zeros(len(self.dtype())))],
                           dtype=self.dtype())
        self._pulses = np.append(self._pulses, newRow)
    def storeHiRes(self,fitterObj):
        params = fitterObj.storeParams()
        amps = params[::2]
        poss = params[1::2] + fitterObj.sampleForPulseArrival() + fitterObj.firstSample()
        shapes = fitterObj.shapes()
        for i in range(len(amps)):
            self.addPulseRow()
            self._pulses['ph'][-1] = amps[i]
            self._pulses['sample'][-1] = poss[i]
            self._pulses['pixel'][-1] = self.pixel
            self._pulses['shape'][-1] = shapes[i]
            self._pulses['resolution'][-1] = 0
            if bool(self.T0Sample()):
                self._pulses['time'][-1] = self.TtimeForSample(poss[i])
    def storeLoRes(self,fitterObj):
        params = fitterObj.lowResFit()
        amps = params[::2]
        poss = params[1::2] + fitterObj.firstSample()
        for i in range(len(amps)):
            self.addPulseRow()
            self._pulses['ph'][-1] = amps[i]
            self._pulses['sample'][-1] = poss[i]
            self._pulses['pixel'][-1] = self.pixel
            self._pulses['resolution'][-1] = 1
            if bool(self.T0Sample()):
                self._pulses['time'][-1] = self.TtimeForSample(poss[i])
    def storePulses(self,fitterObj):
        self.storeHiRes(fitterObj)
        self.storeLoRes(fitterObj)
    def segmentData(self,i1,i2):
        filteredSeg = np.append(self._phData[i1:i2],self._ptData[i1:i2])
        unfilteredSeg = self._normData[i1:i2]
        maskSeg = self._mask[i1:i2]
        segmentObj = fitSegements(filteredSeg,unfilteredSeg,maskSeg,i1)
        return segmentObj        
    def nextSegment(self,fitterObj,segmentLength=65536):
        i1 = fitterObj.nextSegmentStart()
        i2 = i1 + segmentLength
        fitterObj.resetOverlap()
        if i2 >= len(self._data):  # for last segment
            i2 = len(self._data)
            fitterObj.setOverlap(0)
        segmentObj = self.segmentData(i1,i2)
        fitterObj.nextSegment(segmentObj)
        return fitterObj    
    def badSegment(self,fitterObj):
        assert(fitterObj.numPulses()==0)
        i1 = fitterObj.firstSample()
        iGoal = fitterObj.nextSegmentStart()
        i2 = findNextFlat(self._normData,iGoal)
        # some good pulses may already have been subtracted off the beginning
        if fitterObj.previousFitSuccess():
            self._mask[i1+fitterObj.sampleForPulseArrival():i2] = 0
        else: 
            self._mask[i1:i2] = 0
        badSegmentObj = self.segmentData(i1,i2)
        fitterObj.setSegment(badSegmentObj)
        fitterObj.setOverlap(0)
        return fitterObj
    def fitFilteredData(self,filterObj):
        self.resetPulses()
        fitterObj = fitter(filterObj)
        i1 = 0
        i2 = i1+65536
        segmentObj = self.segmentData(i1,i2)
        fitterObj.setSegment(segmentObj)
        printWhen = 1
        t1 = time.time()
        print "Fitting Filtered Data..."
        while fitterObj.nextSegmentStart()<len(self._data):
            fitterObj.fit()
            if fitterObj.fitSuccess():
                if self.interactive:
                    stayInteractive = fitterObj.activeFit()
                    self.interactive = stayInteractive             
            else:
                fitterObj = self.badSegment(fitterObj)
            self.storePulses(fitterObj)
            fitterObj = self.nextSegment(fitterObj)
            fc = float(fitterObj.firstSample())/len(self._data)
            t2 = time.time()
            if fc//0.1==printWhen:
                msg = "%d%% Complete. Estimated " % (int(fc*100))
                msg += "%d Seconds Remaining" % (int((t2-t1)/fc-(t2-t1)))
                print msg
                printWhen += 1
        print  "%.1f%% Total Deadtime"% self.deadTime()
    def plotPulses(self):
        scatterPulses(self._pulses)
    def nonLinGain(self):
        pulseHeights = self._pulses['ph']
        lines = spectralLines(pulseHeights)  
        coefs = lines.nonLinGain()
        if np.any(coefs):
            self._pulses['lin'] = coefs[0]
            self._pulses['quad'] = coefs[1]
            self._pulses['energy'] = (self._pulses['lin']*self._pulses['ph']
                                  + self._pulses['quad']*np.square(self._pulses['ph']))
            print "Non-Linear Gain Correction Applied"
        else: 
            print "Unable To Fit Non-Linearities"
            
class spectralLines(object):
    def __init__(self, pulseHeights):
        self.ph = pulseHeights
        self.binWidth = 5.
        self.maxE = 4000.
        self.gl = {'B':183.3, 'C':277.0, 'N':392.4, 'O':524.9, 'F':676.8, 
                   'Ka':3313.8, 'Kb':3589.6}  
        self.realE = np.array([])
        self.fitE = np.array([])
        self.fitEErr = np.array([])
    def setBinWidth(self,binWidth):
        self.binWidth = binWidth
    def setMaxEnergy(self,maxEnergy):
        self.maxE = maxEnergy
    def numBins(self):
        return int(self.maxE/self.binWidth)
    def histogram(self):
        rng = [0,self.maxE]
        bins = self.numBins()
        hy,hx = np.histogram(self.ph,range=rng,bins=bins)
        hx = hx[:-1] + 0.5*self.binWidth
        return hx,hy
    def fitLine(self,lineID):
        guess = self.gl[lineID]
        halfWidth = 50.
        linePH = self.ph[np.abs(self.ph-guess)<halfWidth]
        if len(linePH)>0:
            mu = np.mean(linePH)
            sig = np.std(linePH)
            halfWidth = 7*sig
            hx,hy = self.histogram()
            yfit = hy[np.abs(hx-mu)<halfWidth]
            xfit = hx[np.abs(hx-mu)<halfWidth]
            if np.sum(yfit)>0:
                amp = np.max(yfit)
                guess = np.array([mu,sig,amp,0])
                try:
                    popt, pcov = curve_fit(fitGaussWBG, xfit, yfit, p0=guess, 
                                           bounds=(0,np.inf))
                    area = np.sqrt(2*np.pi)*popt[2]*popt[1]/self.binWidth
                    bkgdFlux = np.sqrt(popt[3]*(2*halfWidth)/self.binWidth)
                    if (area>bkgdFlux):
                        self.fitE = np.append(self.fitE,popt[0])
                        self.fitEErr = np.append(self.fitEErr,popt[1]/np.sqrt(area))
                        self.realE = np.append(self.realE,self.gl[lineID])
                except RuntimeError:
                    pass
    def fitLines(self):
        for line in self.gl:
            self.fitLine(line)
    def nonLinGain(self,order=2):
        self.fitLines()
        fit = np.zeros(order)
        if len(self.fitE)>order: # orthogonal distance regression
            guess = np.zeros(order)
            guess[0] = 1.
            myModel = Model(fitPolyBkwd)
            xErr = self.fitEErr
            yErr = .1*np.ones(len(self.realE))
            myData = RealData(self.fitE,self.realE,sx=xErr,sy=yErr)
            myODR = ODR(myData,myModel,beta0=guess)
            myOutput = myODR.run()
            fit = myOutput.beta 
        elif len(self.realE)==order: # analytic solution
            mat = np.zeros((order,order))
            vec = np.zeros(order)
            for i in range(order):
                vec[i] = self.realE[i]
                for j in range(order):
                    mat[i,j] = self.fitE[i]**(j+1)
            fit = np.linalg.solve(mat,vec)
        return fit

##############################################################################
# General Use Funciton Definitions 
# These are all needed for the above classes
# Be careful changing these, since they may be called by multiple subroutines 
##############################################################################
def fitPolyBkwd(B,x):
    func = np.zeros(len(x))
    for i in range(len(B)):
        func += B[i] * np.power(x,i+1)
    return func

def fitGaussWBG(x, *params):
    mu = params[0]
    sig = params[1]
    amp = params[2]
    yint = params[3]
    if len(params)>4:
        slp = params[4]
    else:
        slp = 0
    return yint + slp*x + amp*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.)))

def scatterPulses(pulses):
    goodPulses =  pulses[np.where(pulses['resolution']==0)]
    mpl.scatter(goodPulses['sample'],goodPulses['ph'],marker='+',lw=.2)
    mpl.xlabel('Sample')
    mpl.ylabel('Pulse Height (Aprox. eV)')
    mpl.title('Hi-Res Pulses')
    mpl.show(block=True)

def plotSignalNoise(signal,noise2,fs):
    filterLength = len(signal)
    dt = 1./fs
    df = fs/filterLength
    timeSeries = np.arange(filterLength)*dt
    frequencySeries = np.arange(filterLength/2 +1)*df
    mpl.figure()
    mpl.subplot(211)
    mpl.plot(timeSeries, signal, label = 'Avg. Pulse')
    mpl.xlabel('Time [s]')
    mpl.ylabel('Volts')
    mpl.legend(loc='best')
    mpl.subplot(212)
    mpl.loglog(frequencySeries, noise2, label = 'Noise Power Spectrum')
    mpl.xlabel('Frequency [Hz]')
    mpl.ylabel('Vrms^2 / Hz')
    mpl.legend(loc='best')
    mpl.show(block=True)
    
def plotFilters(amplitude,arrival,fs):
    filterLength = len(amplitude)
    dt = 1./fs
    timeSeries = np.arange(filterLength)*dt
    mpl.figure()
    mpl.subplot(211)
    mpl.plot(timeSeries, amplitude, label = 'Amplitude Filter')
    mpl.xlabel('Time [s]')
    mpl.ylabel('1/V')
    mpl.legend(loc='best')
    mpl.subplot(212)
    mpl.plot(timeSeries, arrival, label = 'Arrival Filter')
    mpl.xlabel('Time [s]')
    mpl.ylabel('1/V')
    mpl.legend(loc='best')
    mpl.show(block=True)
    
def plotFilteredPulses(amplitude,arrival,fs):
    templateLength = len(amplitude)
    dt = 1./fs
    timeSeries = np.arange(templateLength)*dt
    mpl.figure()
    mpl.subplot(211)
    mpl.plot(timeSeries, amplitude, label = 'Amplitude Filtered Pulse')
    mpl.ylabel('Arb. Units')
    mpl.xlabel('Time [s]')
    mpl.legend(loc='best')
    mpl.subplot(212)
    mpl.plot(timeSeries, arrival, label = 'Arrival Filtered Pulse')
    mpl.ylabel('Arb. Units')
    mpl.xlabel('Time [s]')
    mpl.legend(loc='best')
    mpl.show(block=True)

def findNextFlat(data,iGoal,length=8192,step=65536):
    flatFound = False
    j1 = max(0,iGoal-length)
    j2 = min(j1+step,len(data))
    while not flatFound:
        searchStream = data[j1:j2]
        firstFlat = findFlat(searchStream,iGoal-j1,length=length)[1]
        if (firstFlat<len(searchStream)) or (j2>=len(data)):
            flatFound = True
        else:
            j2 += step
    return firstFlat+j1 

def innerProduct(f1,f2):
    return np.inner(f1,f2).astype(float)/np.inner(f2,f2)       
    
def dataFunction(segmentLength, dictionaryList, *params):
    yValues = np.array([])
    for dictionary in dictionaryList:
        function = np.zeros(segmentLength)
        templateLength = len(dictionary[0])
        precision = len(dictionary)-1
        for i in range(0,len(params),2):
            amplitude = params[i]
            floatPosition = params[i+1]
            floorPosition = int(np.floor(floatPosition))
            if ((floorPosition>-1.0*templateLength)
                and(floorPosition<segmentLength)):
                templateStart = max(0, -1*floorPosition)
                templateStop = min(templateLength, segmentLength-floorPosition)
                fullDecimal = floatPosition - floorPosition
                equivalentDecimalLHS = ((fullDecimal // (1./precision)) 
                                        * 1./precision)
                equivalentDecimalLHS = np.around(equivalentDecimalLHS, 
                                                 decimals=2)
                equivalentDecimalRHS = np.around(equivalentDecimalLHS 
                                                 + 1./precision, decimals=2)
                templateWeight = ((equivalentDecimalRHS - fullDecimal)
                                  /(1./precision))
                templateToAdd = (templateWeight*dictionary[equivalentDecimalLHS] 
                                 + (1-templateWeight)*dictionary[equivalentDecimalRHS])
                templateToAdd = templateToAdd[templateStart:templateStop]
                segmentStart = max(0,floorPosition)
                segmentStop = min(segmentLength, floorPosition+templateLength)
                assert (templateStop-templateStart)==(segmentStop-segmentStart)
                function[segmentStart:segmentStop] += amplitude * templateToAdd
        yValues = np.append(yValues, function)
    return yValues
    
def fitFunction(x, segmentLength, dictionaryList, *params):
    totalFunction = dataFunction(segmentLength, dictionaryList, *params)    
    return totalFunction[x]

def findFlat(dataSeg,position,length=4096,calAmp=3313.8):    
    pls = trigger(dataSeg) 
    amp,rt = phRT(dataSeg,pls)
    minSep = length*np.ones(len(pls)).astype(int)
    # big pulses require more recovery time
    minSep[np.where(amp>calAmp)] = 2*length 
    posArg = len(np.where(pls<position)[0])
    pls = np.insert(pls,posArg,position)
    minSep = np.insert(minSep,posArg,0)
    pls = np.append(0,pls)
    minSep = np.append(2*length,minSep)
    cutLoc = pls+minSep
    for i in range(len(pls)):
        cutLoc[i] = np.max(cutLoc[:i+1])
    minSep = cutLoc-pls
    minSep = np.append(0,minSep)
    plsSep = np.diff(pls)
    plsSep = np.append(0,plsSep)
    plsSep = np.append(plsSep,len(dataSeg)-pls[-1])
    LHplsNum = np.where((plsSep[:-1]>=minSep[:-1])&(pls<=position))[0]
    if len(LHplsNum)>0:
        LH = pls[LHplsNum[-1]]
    else:
        LH = 0
    RHplsNum = np.where((plsSep[1:]>=minSep[1:])&(pls>=position))[0]
    if len(RHplsNum)>0:
        RH = pls[RHplsNum[0]]+minSep[RHplsNum[0]+1]
    else:
        RH = len(dataSeg)
    return LH, RH

def findSaturatedPulses(data, saturationType, threshold=None):    
    saturatedPulses = np.array([], dtype=int)
    saturationWidth = np.array([], dtype=int)
    if saturationType=='minimum':
        if threshold is None:
            saturationPoint = np.nanmin(data)
        else:
            saturationPoint = threshold
        saturatedSamples = np.where(data<=saturationPoint)[0] 
    elif saturationType == 'maximum':
        if threshold is None:
            saturationPoint = np.nanmax(data)
        else:
            saturationPoint = threshold
        saturatedSamples = np.where(data>=saturationPoint)[0]
    else:
        print 'Unrecognized Option'
        return saturatedPulses
    # if at least one region is found, continue. Else return empty list
    if len(saturatedSamples)!=0: 
        # Within indices find continuous regions so we can separate them
        saturatedSampleDiff = np.diff(saturatedSamples) 
        if saturationType == 'minimum':  
            saturatedRegions = np.where(saturatedSampleDiff > 1)[0]
        else:
            saturatedRegions = np.where(saturatedSampleDiff > 10)[0]
        i1 = saturatedSamples[0]
        if (len(saturatedRegions)==0)and(len(saturatedSampleDiff)>1):
            saturatedPulses = np.append(saturatedPulses, saturatedSamples[0])
            saturationWidth = np.append(saturationWidth, len(saturatedSamples))
        for region in saturatedRegions:        
            i2 = saturatedSamples[region]+1
            if i2-i1>2:
                saturatedPulses = np.append(saturatedPulses, i1)
                saturationLength = i2-i1 
                saturationWidth = np.append(saturationWidth, saturationLength)
            i1 = saturatedSamples[region+1]
    return saturatedPulses

def removeSaturatedPulses(data,positions,fade=512,searchWindow=65536,
                          minSeparation=8192):
    dataLen = len(data)
    mask = np.ones(dataLen)
    modData = data.copy()
    if len(positions)>0:
        j1s = positions - searchWindow//4
        j2s = j1s + searchWindow
        firstEligibleCut = positions + minSeparation
        # Select saturated pulses far enough apart to warrant separate sections
        keepPositions = np.where(positions[1:]>firstEligibleCut[:-1])[0] +  1
        goodPositions = positions[(keepPositions,)]
        goodPositions = np.append(positions[0],goodPositions)
        goodJ1s = j1s[(keepPositions,)]
        goodJ1s = np.append(j1s[0],goodJ1s)
        goodJ2s = j2s[(keepPositions-1,)]
        goodJ2s = np.append(goodJ2s,j2s[-1])
        keepPositions = np.append(0,keepPositions)
        smplComplete = 0
        for i in range(len(goodPositions)):
            if i%100==0:
                print i, "of", len(goodPositions), "Saturated Segments Removed."
            position = goodPositions[i]
            if position>smplComplete:
                LHCutFound = False
                RHCutFound = False
                j1 = max(0,goodJ1s[i])
                j2 = min(dataLen,goodJ2s[i])
                while not (LHCutFound and RHCutFound):
                    posInSeg = position-j1
                    segLen = j2-j1
                    LHcut,RHcut = findFlat(data[j1:j2],posInSeg)
                    if (LHcut>0 or j1==0) and not LHCutFound:
                        LH = LHcut + j1
                        LHCutFound = True
                    if (RHcut<segLen or j2==dataLen) and not RHCutFound:
                        RH = RHcut + j1
                        RHCutFound = True
                    if not LHCutFound:
                        j1 = max(0,j1-searchWindow//4)
                        j2 = j1 + searchWindow
                    elif not RHCutFound:
                        if j2<position:
                            j2 = min(dataLen,goodJ2s[i])
                        j2 = min(j2+3*searchWindow//4,dataLen)
                        j1 = j2 - searchWindow
                cut = RH-LH                        
                if LH-fade>0 and RH<dataLen:
                    fadein = fade
                    fadeout = min(fade,cut)
                    xSeries = np.append(np.arange(fadein),
                                        np.arange(fadein+cut-fadeout,fadein+cut))
                    ySeries = np.append(modData[LH-fadein:LH],modData[RH-fadeout:RH])
                    transition = (1-halfCosWindow(np.ones(fadein+fadeout),.5))
                    polyBSLN = np.polyfit(xSeries,ySeries,5,w=transition)
                    lineBSLN = np.polyfit(xSeries,ySeries,1,w=transition)
                    fadedY = ((ySeries - np.polyval(polyBSLN,xSeries))*transition 
                              + np.polyval(lineBSLN,xSeries)) 
                    if fadeout==fade:
                        patch =  np.hstack((fadedY[:fadein],
                                            np.polyval(lineBSLN,
                                                       np.arange(fadein,
                                                                 fadein+cut-fadeout)),
                                            fadedY[-1*fadeout:]))                          
                    else:
                        patch = fadedY
                    modData[LH-fadein:RH] = patch      
                    mask[LH:RH] = 0
                else:
                    k1 = max(0,LH-fade)
                    k2 = min(dataLen,RH)
                    modData[k1:k2] = 0
                    mask[k1:k2] = 0
                smplComplete = RH      
    return modData, mask.astype(bool)

def FWHM(xdata, ydata, xpos):
    try:
        xidx = np.where(xdata==xpos)[0][0]
    except IndexError:
        print "ERROR: XPOSITION NOT IN XVECTOR"
        return xdata[-1] - xdata[0]
    try:
        ymax = float(ydata[xidx])
        try:
            xLH = xdata[:xidx][ydata[:xidx]<ymax/2][-1]
            xRH = xdata[xidx:][ydata[xidx:]<ymax/2][0]
            return xRH - xLH
        except IndexError:
            print "ERROR: FWHM NOT CONTAINED IN GIVEN RANGE"
            return xdata[-1] - xdata[0]
    except ValueError:
        print "ERROR: NO MAXIMUM IN DATA"
        return xdata[-1] - xdata[0]

def avgNoise(data,locations,length,fs=1./96E-6):
    df = fs/length
    halfCosFraction = 0.5
    window = halfCosWindow(np.ones(length),halfCosFraction)    
    windowAdjustmentFactor = 1./(np.mean(np.square(window)))
    avgPowerSpec = np.zeros(length/2+1)
    numNoise = 0
    for nIdx in locations:
        sampledNoise = data[nIdx:nIdx+length]
        sig, med = trueSigma(sampledNoise)
        sampledNoise -= med
        if np.sum(np.abs(sampledNoise)>5*sig)==0:
            sampledNoise = halfCosWindow(sampledNoise,halfCosFraction)
            sampledFFT = np.fft.rfft(sampledNoise)
            sampledPS = np.square(np.abs(sampledFFT))
            normalizedPSD = 2*windowAdjustmentFactor*sampledPS/(df*length**2)
            avgPowerSpec += normalizedPSD
            numNoise += 1
    if numNoise<10:
        print "ERROR: Not Enough Noise Samples"
        avgPowerSpec = np.zeros(length/2+1)
    else:
        avgPowerSpec = avgPowerSpec/numNoise
        print numNoise, "Good Noise Samples" 
    return avgPowerSpec

def avgSignal(data,locations,length):
    numCal = len(locations)
    signal = np.zeros(length)
    if numCal<5:
        print "ERROR: Not Enough Cal Pulses (%d)" % numCal
    else:
        for i in range(numCal):
            j = int(locations[i])
            signal += data[j-length/4:j+3*length/4]
        signal = signal/float(numCal)
        print numCal, "Good Cal Pulses"    
    return signal

def boxcarSmooth(inArray, boxsize, algorithm=None):  
    box = np.ones(boxsize)/float(boxsize)
    if algorithm=='long':
        outArray = longConvolution(inArray,box,0)
    else:
        outArray = sig.fftconvolve(inArray,box,mode='full')
        outArray = outArray[:len(inArray)]
    return outArray
    
def pulseHeights(dataSegment,positions,discriminator=0):
    dat = triggerChannel(dataSegment,deriv=0)
    plsPeak = np.diff(dat)<=0
    plsPeak = np.append(plsPeak,True)
    pulseHeights = np.array([],dtype=float)
    for i in range(len(positions)):
        pulseStart = positions[i]
        preTrig = dat[pulseStart]
        startSearch = pulseStart + np.argmax((dat[pulseStart:]-preTrig)>
                                              discriminator)
        pulsePeak = startSearch + np.argmax(plsPeak[startSearch:])
        rawnp = dat[pulsePeak]
        ph = rawnp-preTrig
        pulseHeights = np.append(pulseHeights,ph)
    return pulseHeights
def phRT(dataSegment, positions, dt=96E-6, discriminator=0, doPlot=False): 
    dat = triggerChannel(dataSegment,deriv=0)
    plsPeak = np.diff(dat)<=0
    plsPeak = np.append(plsPeak,True)
    pulseHeights = np.array([],dtype=float)
    riseTimes = np.array([],dtype=float)
    fb = 0.2
    ft = 0.8
    for i in range(len(positions)):
        pulseStart = positions[i]
        preTrig = dat[pulseStart]
        startSearch = pulseStart + np.argmax((dat[pulseStart:]-preTrig)>
                                              discriminator)
        pulsePeak = startSearch + np.argmax(plsPeak[startSearch:])
        rawnp = dat[pulsePeak]
        ph = rawnp-preTrig
        rawfb = rawnp-(1.-fb)*(ph)
        rawft = rawnp-(1.-ft)*(ph)
        tfb = np.interp(rawfb,dat[pulseStart:pulsePeak+1],np.arange(pulseStart,
                        pulsePeak+1))
        tft = np.interp(rawft,dat[pulseStart:pulsePeak+1],np.arange(pulseStart,
                        pulsePeak+1))
        rt = dt*(tft-tfb)/np.log((1.-fb)/(1.-ft))
        pulseHeights = np.append(pulseHeights,ph)
        riseTimes = np.append(riseTimes,rt)
        if doPlot:
            mpl.plot(dat[pulseStart-248:pulseStart+776])
            mpl.vlines(248,rawnp-ph,rawnp)
            mpl.show(block=True)
    return pulseHeights, riseTimes

def selectCals(data,locations,length,fs=1./96E-6,interactive=False):
    amplitudes, RTs = phRT(data,locations,dt=1./fs)
    phPDFy,phPDFx = np.histogram(amplitudes,range=[2,6],bins=400)
    phPDFx = phPDFx[:-1]+.005
    phPDFy = boxcarSmooth(phPDFy,2)
    calPeak = phPDFx[phPDFy.argmax()]
    calFWHM = FWHM(phPDFx,phPDFy,calPeak)
    xmin,xmax = calPeak-calFWHM, calPeak+calFWHM    
    rtPDFy,rtPDFx = np.histogram(RTs[np.where((amplitudes>xmin)&(amplitudes<xmax))],
                                 range=[0,.0025],bins=250)
    rtPDFx = rtPDFx[:-1]+5e-6    
    rtPDFy = boxcarSmooth(rtPDFy,2)
    rtPeak = rtPDFx[rtPDFy.argmax()]
    rtFWHM = FWHM(rtPDFx,rtPDFy,rtPeak)
    ymin,ymax = rtPeak-rtFWHM, rtPeak+rtFWHM
    if interactive:
        mpl.scatter(amplitudes, RTs, marker='+')
        mpl.vlines([xmin,xmax],ymin,ymax)
        mpl.hlines([ymin,ymax],xmin,xmax)
        mpl.ion()
        mpl.show(block=False)
        acceptRange = not raw_input("Accept Pulse Selection? (Y/n)")[:1].lower()=='n'
        if not acceptRange:
            raw_input("Adjust Axis Limits To Only Include Cals, Then Press Enter.")
            ymin, ymax = mpl.ylim()
            xmin, xmax = mpl.xlim()
        mpl.close()
    goodCals = locations[(RTs>=ymin)&(RTs<=ymax)&(amplitudes>=xmin)&(amplitudes<=xmax)]
    return goodCals

def findMultiMax(level, data, minSep=0):
    maxlist = np.array([],dtype=int)
    # Find indices of all points that exceed threshold
    i = np.where(data > level)[0] 
    # if at least one region is found, continue. Else return empty list
    if len(i)!=0: 
        # Within indices find continuous regions so we can separate them
        z = np.diff(i) 
        j = np.where(z > 1)[0]
        i1 = i[0]
        # if there are not multiple regions, return the global max location
        if (len(j)==0): 
            maxloc = np.argmax(data)
            maxlist = np.append(maxlist, maxloc)
        else:
            for k in j:  # work to find max within each region         
                i2 = i[k]+1
                maxloc = np.argmax(data[i1:i2])+i1
                maxlist = np.append(maxlist, maxloc)
                i1 = i[k+1]
            maxloc = np.argmax(data[i1:])+i1 # needed to get last region
            maxlist = np.append(maxlist, maxloc)
            if len(maxlist)>1:
                maxlist = np.append(maxlist[0], 
                                    maxlist[1:][np.where(np.diff(maxlist)>minSep)])
    return maxlist    

def trueSigma(inArray):
    tolerance = .01
    maxIter = 10
    relChange = 1
    numIter = 0 
    currentPoints = inArray.copy()
    sig = np.std(currentPoints)
    med = np.median(currentPoints)
    while (relChange>tolerance)and(numIter<maxIter)and(len(currentPoints)>0):
        currentPoints = currentPoints[np.where(np.abs(currentPoints-med)<5*sig)]
        if len(currentPoints)>0:
            relChange = (sig - np.std(currentPoints)) / np.std(currentPoints)
            sig = np.std(currentPoints)
            med = np.median(currentPoints)
            numIter += 1
    return sig, med

def triggerChannel(dataSegment, deriv=1, invert=True):
    # filter parameters chosen for risetime of 1.5ms (~15 samples)
    # Calculate a smoothed 1st derivative
    savgol = sig.savgol_filter(dataSegment, window_length= 15, polyorder = 3, 
                               deriv=deriv)     
    if invert:    
        savgol = -1.0*savgol
    return savgol   

def trigger(dataSegment):
    trig = triggerChannel(dataSegment)
    triggerLevel = 5*trueSigma(trig)[0]    
    triggerList = findMultiMax(triggerLevel, trig)
    preTrig = (np.diff(trig)<=0)&(trig[1:]<=0)
    preTrig = np.append(True,preTrig)
    samples = np.array([t-np.argmax(preTrig[:t+1][::-1]) for t in triggerList],dtype=int)
    return samples
    
def signalAndNoiseLoc(inArray, filLen):
    pulseLocations = trigger(inArray) 
    pre = np.diff(pulseLocations)
    post = np.copy(pre)
    pre = np.append(pulseLocations[0],pre)
    post = np.append(post,len(inArray)-pulseLocations[-1]-1)
    signalLoc = pulseLocations[np.where((pre>5*filLen/4)&(post>3*filLen/4))]
    noiseLoc = pulseLocations[np.where(pre>9*filLen/4)]-5*filLen/4
    return signalLoc, noiseLoc

def fftShift(inArray, sampleShift, sampleForPulseArrival):
    inLength = len(inArray)
    inFft = np.fft.rfft(inArray)
    samples = np.arange(len(inFft))
    # +sampleForPulseArrival/templateLength needed to 
    # account for t=0 sample occuring part way through the trace
    shiftedFft = (np.exp(-2*np.pi*np.complex(0,1)*samples
                        *(sampleShift+float(sampleForPulseArrival)/float(inLength))
                        /inLength)*inFft) 
    outArray = np.fft.irfft(shiftedFft)
    if inLength%2 != 0:
        outArray = np.append(outArray, 0)
    assert(len(outArray)==inLength)
    return outArray
    
def longConvolution(data, filter_t, sampleForPulseArrival):    
    totalDataLength = len(data)
    convolutionLength = int(2**20) 
    filterLength = len(filter_t)
    i1 = 0
    i2 = min(int(i1+convolutionLength),totalDataLength)
    dataFiltered = np.array([], dtype=float)
    while i2 <= totalDataLength:
        dataConvolutionSegment = data[i1:i2]
        dataFilteredSeg = sig.fftconvolve(dataConvolutionSegment, 
                                          filter_t, mode='full') 
        # Some housekeeping needed to conserve sample-Ttime conversion 
        dataFilteredSeg = np.roll(dataFilteredSeg, -1*sampleForPulseArrival) 
        dataFilteredSeg = dataFilteredSeg[:-1*filterLength+1]
        assert(len(dataFilteredSeg)==len(dataConvolutionSegment))
        if i1==0:
            dataFiltered = np.append(dataFiltered, 
                                     dataFilteredSeg[:-1*filterLength])
        elif i2==totalDataLength:
            dataFiltered = np.append(dataFiltered, 
                                     dataFilteredSeg[filterLength+1:])
            break
        else:
            dataFiltered = np.append(dataFiltered, 
                                     dataFilteredSeg[filterLength+1
                                                     :-1*filterLength])
        i1 = i2 - 2*filterLength - 1
        i2 = min(int(i1+convolutionLength), totalDataLength)
    assert(len(dataFiltered)==totalDataLength)
    return dataFiltered

def halfCosWindow(inArray, halfCosFraction):
    window = np.ones(len(inArray))
    halfCosLength = int(halfCosFraction*len(inArray))
    xSeries =  np.arange(2*halfCosLength)/float(2*halfCosLength-1)
    cosFunction = 0.5*(1 - np.cos(2*np.pi*xSeries))
    cosFunction[-1] = 0
    window[:halfCosLength] = cosFunction[:halfCosLength]
    window[-1*halfCosLength:] = cosFunction[-1*halfCosLength:]
    outArray = inArray*window
    return outArray

    
    

    

    
    

    
    
                
     
    
    
        
