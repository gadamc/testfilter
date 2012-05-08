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
        
        s = couchdbkit.Server('http://localhost:5984')
        db = s['pulsetemplates']
        
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
        
        plt.subplot(7,1,2)
        plt.cla()
        plt.plot( get_as_nparray(win.GetWindow(), win.GetWindowSize()) )
        plt.title('window')
        
        #winpulse = std.vector("double")()
        #winpulse.resize(win.GetOutputPulseSize())
        #for i in range(winpulse.size()):
        #  winpulse[i] = win.GetOutputPulse()[i]
        winpulse = get_as_nparray(win.GetOutputPulse(), win.GetOutputPulseSize())
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
        tempPower = get_out(hc2p)
          
        hc2p.SetInputPulse(kernel)
        hc2p.RunProcess()
        kernelPower = get_out(hc2p)
        
       
        plt.subplot(7,1,4)
        plt.cla()
        # plt.loglog(np.array(noisePower[:len(noisePower)-1]))
        #         plt.loglog(np.array(tempPower[:len(tempPower)-1]))
        #         plt.loglog(np.array(kernelPower[:len(kernelPower)-1]))
        npp = np.array(noisePower)
        # print npp 
        #         print kpp
        #         print len(npp)
        #         print len(kpp)
        #         print 
        plt.loglog(npp)
        plt.loglog(tempPower)
        plt.loglog(kernelPower)
        plt.title('filter power')
        
        hc2p.SetInputPulse(optFilter.GetInputPulse(), optFilter.GetInputPulseSize())
        hc2p.RunProcess()
        pulsePower = get_out(hc2p)
            
        plt.subplot(7,1,5)
        plt.cla()
        plt.loglog(pulsePower)
        plt.title('pulse power')
        
        optOut = get_out(optFilter)
        
        plt.subplot(7,1,6)
        plt.cla()
        plt.plot(optOut)
        plt.title('amplitude estimator')
        
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
        
        
        plt.subplot(1,1,1)
        plt.cla()
        
        bas = cham.GetBaselineRemovalHeat()
        baspulse = get_out(bas)
          
        plt.plot(baspulse)
        s = couchdbkit.Server('http://localhost:5984')
        #s = couchdbkit.Server('https://edwdbik.fzk.de:6984')
        db = s['pulsetemplates']
        vr = db.view('analytical/bychandate',reduce=False, descending=True, startkey=[p.GetChannelName(), "2012-01-22 00:00:00.0"], limit=1, include_docs=True)
        doc = vr.first()['doc']
        vp = std.vector("double")()
        exec(doc['formula']['rise_double_decay']['python']) #defines 'template' function
        #doc['formula']['rise_double_decay']['par'][2] = doc['formula']['rise_double_decay']['par'][2]/2.016
        #doc['formula']['rise_double_decay']['par'][3] = doc['formula']['rise_double_decay']['par'][3]/2.016
        #doc['formula']['rise_double_decay']['par'][5] = doc['formula']['rise_double_decay']['par'][5]/2.016
        #doc['formula']['rise_double_decay']['par'][0]=300  #put it close to zero, but away from the windowing function
        for i in range( 512 ):
          vp.push_back( template(i*2.016, doc['formula']['rise_double_decay']['par']))
          
        #plot the pulse templates for documentation
        scaleFactor = abs(min(np.array(winpulse)))/abs(min(np.array(vp)))
        print 'scaling by', scaleFactor
        print vp.size(), vp[i], vp[vp.size()-1], vp[i]*scaleFactor

        for i in range(vp.size()):
          vp[i] = scaleFactor * vp[i]

        plt.plot(np.array(vp))
        
        raw_input()

#s = couchdbkit.Server('https://edwdbik.fzk.de:6984')
s = couchdbkit.Server('http://localhost:5984')
db = s['pulsetemplates']

filenum = int(sys.argv[1])


#f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/ma14f004_000.root')
f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_00%d.root' % filenum)

#cham.GetHeatWindow().SetWindow(KWindowDesign.GetTukeyWindow(512,0.7), 512);

chanList = ["chalA FID807", "chalB FID807", "chalA FID808", "chalB FID808"]
#chanList = ["chalA FID803", "chalB FID803", "chalA FID806", "chalB FID806"]
#chanList = ["chalA FID802", "chalB FID802", "chalA FID804", "chalB FID804"]
#803 and 806 are on S5 (e) and 802 and 804 are on S6 (f)

hc2p = KHalfComplexPower()
hc2r = KHalfComplexToRealDFT()

cham = KChamonixKAmpSite()
cham.GetHeatPeakDetector().SetOrder(3)
cham.GetHeatPeakDetector().SetNumRms(2.9)


for chan in chanList:
  print chan
  
  vr = db.view('analytical/bychandate',reduce=False, descending=True, startkey=[chan, "2012-01-22 00:00:00.0"], limit=1, include_docs=True)
  doc = vr.first()['doc']
  vp = std.vector("double")()
  #print json.dumps(doc, indent=1)
  #doc['formula']['rise_double_decay']['par'][2] = doc['formula']['rise_double_decay']['par'][2]/2.016
  #doc['formula']['rise_double_decay']['par'][3] = doc['formula']['rise_double_decay']['par'][3]/2.016
  #doc['formula']['rise_double_decay']['par'][5] = doc['formula']['rise_double_decay']['par'][5]/2.016
  
  exec(doc['formula']['rise_double_decay']['python']) #defines 'template' function
  
  #doc['formula']['rise_double_decay']['par'][0]=300  #put it close to zero, but away from the windowing function
  for i in range( 512 ):
    #print i, template(i*2.016, doc['formula']['rise_double_decay']['par'])
    vp.push_back( template(i*2.016, doc['formula']['rise_double_decay']['par']))
    
  #plot the pulse templates for documentation
  #plt.plot(np.array(vp))
  #raw_input('... continue')
  
  scaleFactor = 1./abs(min(np.array(vp)))
  print 'scaling by', scaleFactor
  print vp.size(), vp[i], vp[vp.size()-1], vp[i]*scaleFactor
  
  for i in range(vp.size()):
    vp[i] = scaleFactor * vp[i]
    
  #plt.plot(np.array(vp))
  #raw_input('... continue')
  
  cham.GetPulseTemplateShifter().SetShift(-1 * int(doc['formula']['rise_double_decay']['par'][0]/2.016 + 0.5) )
  cham.SetTemplate(chan, vp)
  #cham.SetTrapDecayConstant(chan, doc['formula']['rise_double_decay']['par'][3])
  
  #plot the power spectrum
  #hc2p.SetInputPulse(cham.GetTemplateSpectrum(chan))
  #hc2p.RunProcess()
  #pw = std.vector("double")()
  #pwh = TH1D('p', 'p', hc2p.GetOutputPulseSize(), 0, hc2p.GetOutputPulseSize())
  #for i in range(hc2p.GetOutputPulseSize()):
  #  pw.push_back(hc2p.GetOutputPulse()[i])
  #  pwh.SetBinContent(i+1, hc2p.GetOutputPulse()[i])
       
  #print min(pw), len(pw), max(pw)
  #ppw = np.array(pw)
  #ppw = ppw*1e9
  #print ppw
  #plt.loglog(ppw)
  #c1 = TCanvas()
  #pwh.Draw()
  #c1.SetLogx()
  #c1.SetLogy()
  #pwh.Draw()
     
     
  #and replot the inverse-fourier transform of the windowed template pulse
  #hc2r.SetInputPulse(cham.GetTemplateSpectrum(chan))
  #hc2r.RunProcess()
  #newTemplate = std.vector("double")()
  #for i in range(hc2r.GetOutputPulseSize()):
  #  newTemplate.push_back(hc2r.GetOutputPulse()[i])
  #     
  #plt.plot(np.array(newTemplate))
  
  #raw_input('... continue')
  #plt.cla()
  #del pwh
  


#scout(f,cham, chanList)
#kamp(f, cham, chanList)
#f.Close()


k = KAmpKounselor()
k.AddKAmpSite(cham)
k.RunKamp('/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_00%d.root' % filenum, '/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_00%d.amp.root' % filenum)

  



