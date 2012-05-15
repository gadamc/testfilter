#!/usr/bin/env python
import couchdbkit
import operator
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
import json, sys
plt.ion()


def get_as_nparray(c_pointer, size):
  data = np.zeros(size)
  for i in range(size):
    data[i] = c_pointer[i]
  return data

def get_out(pta):
  return get_as_nparray(pta.GetOutputPulse(), pta.GetOutputPulseSize())

def get_in(pta):
  return get_as_nparray(pta.GetInputPulse(), pta.GetInputPulseSize())
  
def scout(f, kamp, channels):
  e = f.GetEvent()
  for i in range(f.GetEntries()):
    f.GetEntry(i)
    if operator.mod(i,100) == 0: print i
    for j in range(e.GetNumBoloPulses()):
      p = e.GetBoloPulse(j)
      if p.GetChannelName() not in channels:
        continue
      kamp.ScoutKampSite( e.GetBoloPulse(j) ,e)
      

#
gSystem.Load('libkds')  
gSystem.Load('libkpta')

myChannel = 'ionisA FID806'

lin = KLinearRemoval();
era = KEraPeakFinder();
r2hc = KRealToHalfComplexDFT()
hc2r = KHalfComplexToRealDFT()
hcp = KHalfComplexPower()
pat2 = KPatternRemoval()
pat2.SetPatternLength(6)
pat3 = KPatternRemoval()
pat3.SetPatternLength(100)
pat = KPatternRemoval()
pat.SetPatternLength(10)
window = KWindow()
windesign = KWindowDesign()
window.SetWindow( windesign.GetTukeyWindow(8192, 0.1), 8192)
whit = KNoiseWhitening();
polCalc = KPulsePolarityCalculator()

bbv2 = KBBv2TimeDomainFitKamper()
preProc = KPulseAnalysisChain()
preProc.AddProcessor(lin)
preProc.AddProcessor(pat3)
preProc.AddProcessor(pat2)
preProc.AddProcessor(pat)
preProc.AddProcessor(window)
preProc.AddProcessor(r2hc)
preProc.AddProcessor(whit)  
preProc.AddProcessor(hc2r)
bbv2.SetPreProcessor(preProc)

npChain = KPulseAnalysisChain() #processing chain to calculation the noise spectra
npChain.AddProcessor(window)
npChain.AddProcessor(r2hc)
npChain.AddProcessor(hcp)

#configure the ERA peak finder ???
era.SetOrder(3)
era.SetNumRms(3.5)

f = KDataReader(sys.argv[1])
e = f.GetEvent()

#now, we have to loop through the data and 
#determine the average noise spectrum for `myChannel'

aveSpec = std.vector("double")()  #will store the average spectrum here. 
numOfSpectra = 0.0  # the number of spectra used in calculation of the average


for i in range(f.GetEntries()):
  f.GetEntry(i)

  for ii in range(e.GetNumBoloPulses()):
    p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
  
    if p.GetPulseLength() == 0 or p.GetChannelName() != myChannel: continue  #skip empty pulses
              
    lin.SetInputPulse( p.GetTrace()  )
    lin.RunProcess()   
    pat3.SetInputPulse( lin )
    pat3.RunProcess()
    pat2.SetInputPulse( pat3 ) 
    pat2.RunProcess()
    pat.SetInputPulse( pat2 ) 
    pat.RunProcess()
    era.SetInputPulse( pat )
    era.SetPolarity(polCalc.GetExpectedPolarity(p)) 
    era.RunProcess()
    #print era.GetPeakBins().size()
    if era.GetPeakBins().size() > 0:  continue  #skip pulse if we didn't find a noise pulse
      
    #applying windowing and then calculate spectrum
    npChain.SetInputPulse( pat )
    if npChain.RunProcess() == False:
      print 'excusez-moi... ' #this shouldn't fail
      continue
           
    #i assume here that the pulse length never changes within a run
    if numOfSpectra == 0:
      aveSpec.resize( npChain.GetOutputPulseSize())
      for  k in range(aveSpec.size()):
        aveSpec[k] =  npChain.GetOutputPulse()[k]
    else: #calculate the running averge of the noise spectrum
      for k in range(aveSpec.size()):
        aveSpec[k] = aveSpec[k] * (numOfSpectra-1.0)/numOfSpectra + npChain.GetOutputPulse()[k] / numOfSpectra
  
    numOfSpectra += 1.0  

print numOfSpectra

plt.subplot(1,1,1)
plt.cla()
plt.loglog( get_out(hcp) )
plt.title('noise power')

whit.SetNoiseSpectrum(aveSpec)
raw_input('spin the data')


for i in range(f.GetEntries()):
  f.GetEntry(i)

  for ii in range(e.GetNumBoloPulses()):
    p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
  
    if p.GetPulseLength() == 0 or p.GetChannelName() != myChannel: continue  #skip empty pulses
    
    rec = KPulseAnalysisRecord()  
    bbv2.MakeKamp(p, rec)
    
    plt.subplot(4,1,1)
    plt.cla()
    plt.plot(np.array( p.GetTrace() ))
    plt.title('raw')
    
    plt.subplot(4,1,2)
    plt.cla()
    plt.plot( get_out(lin) )
    plt.title('lin')
    
    plt.subplot(4,1,3)
    plt.cla()
    plt.plot( get_out(pat) )
    plt.title('pat')
    
    
    plt.subplot(4,1,4)
    plt.cla()
    plt.plot( get_out(hc2r) )
    plt.title('hc2r')

    c1 = TCanvas()
    bbv2.GetGraph().Draw('A*')
    bbv2.GetFitFunction().Draw('same')
    
    raw_input('enter for next pulse')
    del c1

