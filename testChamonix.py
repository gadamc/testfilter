#!/usr/bin/env python

import couchdbkit
import operator
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np

plt.ion()



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

def kamp(f, cham, channels):  
  e = f.GetEvent()
  
  hc2p = KHalfComplexPower()
  r2hc = KRealToHalfComplexDFT()
  for ii in range(f.GetEntries()):
    f.GetEntry(ii)
    
    for j in range(e.GetNumBolos()):
      b = e.GetBolo(j)
      
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        if p.GetIsHeatPulse() == False: 
          continue
        if p.GetPulseLength() == 0: 
          continue
        if p.GetChannelName() not in channels:
          continue
          
        print 'Entry', ii, p.GetChannelName()
        
        pulse = p.GetTrace()
        
        tempDft = cham.GetTemplateSpectrum(p.GetChannelName())
        noisePower = cham.GetNoisePower(p.GetChannelName())
        #print noisePower.size()
        #print noisePower[noisePower.size()-1]
        print 'number of noise events found', cham.GetNumNoiseEventsFound(p.GetChannelName())
        optKamper = cham.GetOptimalKamper()
        optFilter = optKamper.GetOptimalFilter()
        optFilter.SetTemplateDFT(tempDft)
        optFilter.SetNoiseSpectrum(noisePower)
        optFilter.SetToRecalculate()
        #optFilter.BuildFilter()
        
        optKamper.SetWindow(cham.GetHeatWindow())
        optKamper.SetBaselineRemoval(cham.GetBaselineRemovalHeat())
        
        pRec= KPulseAnalysisRecord()
        optKamper.MakeKamp(p, pRec)
        
        kernel = std.vector("double")()
        kernel.resize(optFilter.GetOptimalFilterSize())
        for i in range(kernel.size()):
          kernel[i] = optFilter.GetOptimalFilter()[i]
          
        
        plt.subplot(7,1,1)
        plt.cla()
        plt.plot(np.array(pulse))
        plt.title('raw')
        
        win = cham.GetHeatWindow()
        windowfunction = std.vector("double")()
        windowfunction.resize(win.GetWindowSize())
        for i in range(windowfunction.size()):
          windowfunction[i] = win.GetWindow()[i]
        
        
        plt.subplot(7,1,2)
        plt.cla()
        plt.plot(np.array(windowfunction))
        plt.title('window')
        
        winpulse = std.vector("double")()
        winpulse.resize(win.GetOutputPulseSize())
        for i in range(winpulse.size()):
          winpulse[i] = win.GetOutputPulse()[i]
        
        plt.subplot(7,1,3)
        plt.cla()
        plt.plot(np.array(winpulse))
        plt.title('baseline removed and windowed')
        # plt.subplot(4,1,2)
        #         plt.cla()
        #         plt.plot(np.array(tempDft))
        #         plt.title('tempDft')
        
        hc2p.SetInputPulse(tempDft)
        print 'template power', hc2p.RunProcess()
        tempPower = std.vector("double")()
        tempPower.resize(hc2p.GetOutputPulseSize())
        for i in range(tempPower.size()):
          tempPower[i] = hc2p.GetOutputPulse()[i]
          
        hc2p.SetInputPulse(kernel)
        hc2p.RunProcess()
        kernelPower = std.vector("double")()
        kernelPower.resize(hc2p.GetOutputPulseSize())
        for i in range(kernelPower.size()):
          kernelPower[i] = hc2p.GetOutputPulse()[i]
        
       
        plt.subplot(7,1,4)
        plt.cla()
        # plt.loglog(np.array(noisePower[:len(noisePower)-1]))
        #         plt.loglog(np.array(tempPower[:len(tempPower)-1]))
        #         plt.loglog(np.array(kernelPower[:len(kernelPower)-1]))
        npp = np.array(noisePower)
        kpp = np.array(kernelPower)
        # print npp 
        #         print kpp
        #         print len(npp)
        #         print len(kpp)
        #         print 
        plt.loglog(npp)
        plt.loglog(np.array(tempPower))
        plt.loglog(kpp)
        plt.title('filter power')
        
        hc2p.SetInputPulse(optFilter.GetInputPulse(), optFilter.GetInputPulseSize())
        hc2p.RunProcess()
        pulsePower = std.vector("double")()
        pulsePower.resize(hc2p.GetOutputPulseSize())
        for i in range(pulsePower.size()):
          pulsePower[i] = hc2p.GetOutputPulse()[i]
            
            
        plt.subplot(7,1,5)
        plt.cla()
        plt.loglog(np.array(pulsePower))
        plt.title('pulse power')
        
        optOut = std.vector("double")()
        optOut.resize(optFilter.GetOutputPulseSize())
        for i in range(optOut.size()):
          optOut[i] = optFilter.GetOutputPulse()[i]
          
        plt.subplot(7,1,6)
        plt.cla()
        plt.plot(optOut)
        plt.title('opt out')
        
        #calculate chi^2 as a function of t and see what it looks like.
        chi2 = std.vector("double")()
        for time in range(pulse.size()):  #loop over time bins
          chi2.push_back(optFilter.GetChiSquared(time))
        
        plt.subplot(7,1,7)
        plt.cla()
        plt.plot(np.array(chi2))
        plt.title('chi squared')
        plt.show()
        raw_input()
          


#s = couchdbkit.Server('https://edwdbik.fzk.de:6984')
s = couchdbkit.Server('http://localhost:5984')
db = s['pulsetemplates']

cham = KChamonixKAmpSite()
f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_000.root')

#cham.GetHeatWindow().SetWindow(KWindowDesign.GetTukeyWindow(512,0.7), 512);
cham.GetHeatPeakDetector().SetNumRms(2.3)
chanList = ["chalA FID807", "chalB FID807", "chalA FID808", "chalB FID808"]
hc2p = KHalfComplexPower()
hc2r = KHalfComplexToRealDFT()


for chan in chanList:
  print chan
  vr = db.view('analytical/bychandate',reduce=False, descending=True, startkey=[chan, "2012-01-22 00:00:00.0"], limit=1, include_docs=True)
  doc = vr.first()['doc']
  vp = std.vector("double")()
  exec(doc['formula']['rise_double_decay']['python']) #defines 'template' function

  #doc['formula']['rise_double_decay']['par'][0]=300  #put it close to zero, but away from the windowing function
  for i in range( 512 ):
    vp.push_back( template(i, doc['formula']['rise_double_decay']['par']))
    
  #plot the pulse templates for documentation
  scaleFactor = 1./abs(min(np.array(vp)))
  print 'scaling by', scaleFactor
  print vp.size(), vp[i], vp[vp.size()-1], vp[i]*scaleFactor
  
  for i in range(vp.size()):
    vp[i] = scaleFactor * vp[i]
    
  plt.plot(np.array(vp))
  #raw_input('... continue')
  
  cham.GetPulseTemplateShifter().SetShift(-1 * doc['formula']['rise_double_decay']['par'][0] )
  print cham.SetTemplate(chan, vp)
  cham.SetTrapDecayConstant(chan, doc['formula']['rise_double_decay']['par'][3])
  
  #plot the power spectrum
  hc2p.SetInputPulse(cham.GetTemplateSpectrum(chan))
  hc2p.RunProcess()
  pw = std.vector("double")()
  pwh = TH1D('p', 'p', hc2p.GetOutputPulseSize(), 0, hc2p.GetOutputPulseSize())
  for i in range(hc2p.GetOutputPulseSize()):
    pw.push_back(hc2p.GetOutputPulse()[i])
    pwh.SetBinContent(i+1, hc2p.GetOutputPulse()[i])
       
  #print min(pw), len(pw), max(pw)
  #ppw = np.array(pw)
  #ppw = ppw*1e9
  #print ppw
  #plt.loglog(ppw)
  c1 = TCanvas()
  pwh.Draw()
  c1.SetLogx()
  c1.SetLogy()
  pwh.Draw()
     
     
  #and replot the inverse-fourier transform of the windowed template pulse
  hc2r.SetInputPulse(cham.GetTemplateSpectrum(chan))
  hc2r.RunProcess()
  newTemplate = std.vector("double")()
  for i in range(hc2r.GetOutputPulseSize()):
    newTemplate.push_back(hc2r.GetOutputPulse()[i])
       
  plt.plot(np.array(newTemplate))
  
  #raw_input('... continue')
  plt.cla()
  del pwh
  


#scout(f,cham, chanList)
#kamp(f, cham, chanList)
#f.Close()

k = KAmpKounselor()
k.AddKAmpSite(cham)
for i in range(10):
  k.RunKamp('/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_00%d.root' % i, '/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_00%d.amp.root' % i)



