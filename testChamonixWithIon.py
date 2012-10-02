#!/usr/bin/env python

import operator
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
import json, sys
import KDataPy.util as kutil
import KDataPy.database as kdb

plt.ion()


def scout(f, kamp, channels, maxEvents=None):

  for e in f:
    if operator.mod(f.GetCurrentEntryNumber(),100) == 0: print f.GetCurrentEntryNumber()

    for p in e.boloPulseRecords():
      if p.GetChannelName() not in channels: continue
      kamp.ScoutKampSite( p ,e)   

    if maxEvents and f.GetCurrentEntryNumber() >= maxEvents: return   


def getPulseInfo(chan):
  pulseLength = 512
  binSize = 2.016
  min = -2
  max = 50
  type = 0

  if chan.startswith("ion"):
    pulseLength = 8192
    binSize = 1.0
    min = 0
    max = 600
    type  = 1

  return (pulseLength, binSize, min, max, type)


def kamp(p, pta, **kwargs):
  
  f  = kwargs["kdatafile"]
  channels = kwargs["channellist"]
  cham = kwargs["chamonix"]
  hc2p = kwargs["halfcomp2power"]
  r2hc = kwargs["real2halfcomp"]
  pulsePol = kwargs['pulsePol']
  pulseLength, binSize, minSearch, maxSearch, _pulsetype = getPulseInfo(p.GetChannelName())

  db = kwargs['templateDB']

  if p.GetChannelName() not in channels: 
    return None
          
  print 'Entry', f.GetCurrentEntryNumber(), p.GetChannelName()
        
        
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
  
  if p.GetIsHeatPulse():
    optKamper.SetWindow(cham.GetHeatWindow())
    optKamper.SetPreProcessor(cham.GetHeatPreProcessor())
    win = cham.GetHeatWindow()
    preproc= cham.GetHeatPreProcessor()

  elif p.GetBoloBoxVersion() > 1.9:
    optKamper.SetWindow(cham.GetIonWindow())
    optKamper.SetPreProcessor(cham.GetBBv2IonPreProcessor())
    win = cham.GetIonWindow()
    preproc= cham.GetBBv2IonPreProcessor()
  else:
    optKamper.SetWindow(cham.GetIonWindow())
    optKamper.SetPreProcessor(cham.GetBBv1IonPreProcessor())
    win = cham.GetIonWindow()
    preproc= cham.GetBBv1IonPreProcessor()

  optKamper.SetPulseTemplateShiftFromPreTrigger( cham.GetTemplateShift(p.GetChannelName()) );
  optKamper.SetAmplitudeEstimatorSearchRangeMax(maxSearch)
  optKamper.SetAmplitudeEstimatorSearchRangeMin(minSearch)

  plt.figure(1)
  plt.subplot(7,1,1)
  plt.cla()
  plt.plot(np.array(pulse))
  plt.title('raw')
  
          
  
  plt.subplot(7,1,2)
  plt.cla()
  plt.plot( kutil.get_as_nparray(win.GetWindow(), win.GetWindowSize()) )
  plt.title('window')
  
  #winpulse = std.vector("double")()
  #winpulse.reserve(win.GetOutputPulseSize())
  #for i in range(winpulse.size()):
  #  winpulse[i] = win.GetOutputPulse()[i]
  
  preproc.SetInputPulse(pulse)
  preproc.RunProcess()
  win.SetInputPulse(preproc)
  win.RunProcess()
  winpulse = kutil.get_out(win)
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
  tempPower = kutil.get_out(hc2p)


  plt.subplot(7,1,4)
  plt.cla()
  npp = np.array(noisePower)

  plt.loglog(npp)
  plt.loglog(tempPower)

  #raw_input('hit enter to run filter....')

  #wait, doesn't the optKamper do this?
  #r2hc.SetInputPulse(win)
  #r2hc.RunProcess()
  #optFilter.SetInputPulse(r2hc)
  #optFilter.BuildFilter()

  optFilterResults = optKamper.MakeKamp(p)
  
  kernel = kutil.get_as_nparray(optFilter.GetOptimalFilter(), optFilter.GetOptimalFilterSize())

  hc2p.SetInputPulse(kernel, len(kernel))
  hc2p.RunProcess()
  kernelPower = kutil.get_out(hc2p)
  
  plt.loglog(kernelPower)
  plt.title('filter power')
  
  hc2p.SetInputPulse(optFilter.GetInputPulse(), optFilter.GetInputPulseSize())
  hc2p.RunProcess()
  pulsePower = kutil.get_out(hc2p)
      
  plt.subplot(7,1,5)
  plt.cla()
  plt.loglog(pulsePower)
  plt.title('pulse power')
  
  optOut = kutil.get_out(optFilter)
  
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
  #raw_input()
  
  plt.figure(2)
  #plt.subplot(1,1,1)
  plt.cla()
    
  plt.plot(winpulse)

  vr = db.view('analytical/bychandate',reduce=False, descending=True, startkey=[p.GetChannelName(), "2012-01-22 00:00:00.0"], limit=1, include_docs=True)
  doc = vr.first()['doc']
  vp = std.vector("double")()
  vp.reserve(pulseLength)
  exec(doc['formula']['python']) #defines 'template' function
  #doc['formula']['par'][2] = doc['formula']['par'][2]/2.016
  #doc['formula']['par'][3] = doc['formula']['par'][3]/2.016
  #doc['formula']['par'][5] = doc['formula']['par'][5]/2.016
  #doc['formula']['par'][0]=300  #put it close to zero, but away from the windowing function
  for i in range( pulseLength ):
    vp.push_back( template(i*binSize, doc['formula']['par']))
    
  #plot the pulse templates for documentation
  if p.GetChannelName().startswith('chal'): thePolarity = 1 #template heat pulses are already negative polarity
  else: thePolarity = pulsePol.GetExpectedPolarity(p)
  scaleFactor = thePolarity*max(abs(np.array(winpulse)))/max(abs(np.array(vp)))
  optScaleFactor = thePolarity*optFilterResults['amp'].fValue/max(abs(np.array(vp)))
  print 'simple scaling by', scaleFactor
  print 'opt filter scaling by at amp', optScaleFactor
  print vp.size(), vp[i], vp[vp.size()-1], vp[i]*scaleFactor
  optPulse = std.vector("double")()
  optPulse.resize(vp.size())
  diffPulse = std.vector("double")()
  diffPulse.resize(vp.size())

  for i in range(vp.size()):
    optPulse[i] = optScaleFactor * vp[i]
    vp[i] = scaleFactor * vp[i]
    diffPulse[i] = optPulse[i] - vp[i]

  plt.plot(np.array(vp))
  plt.plot(np.array(optPulse))
  #plt.plot(np.array(diffPulse))
  
  print 'optimal filter results'
  for key, val in optFilterResults:
    print key, val.fValue, val.fUnit
  raw_input()


db = kdb.pulsetemplates()

#cham.GetHeatWindow().SetWindow(KWindowDesign.GetTukeyWindow(512,0.7), 512);

#chanList = ["chalA FID807", "chalB FID807", "chalA FID808", "chalB FID808"]
chanList = ["chalA FID807", "chalB FID807", "ionisA FID807", "ionisB FID807", "ionisC FID807", "ionisD FID807"]
#chanList = ["chalA FID803", "chalB FID803", "chalA FID806", "chalB FID806"]
#chanList = ["chalA FID802", "chalB FID802", "chalA FID804", "chalB FID804"]
#803 and 806 are on S5 (e) and 802 and 804 are on S6 (f)

hc2p = KHalfComplexPower()
hc2r = KHalfComplexToRealDFT()
r2hc = KRealToHalfComplexDFT()
pulsePol = KPulsePolarityCalculator()

cham = KChamonixKAmpSite()
cham.GetBBv1IonPeakDetector().SetOrder(2)  
cham.GetBBv1IonPeakDetector().SetNumRms(6.5) 
cham.NeedScout(True)

#
runname = sys.argv[1]
filenum = int(sys.argv[2])
display = False
try:
  if sys.argv[3] == 'display':
    display = True
except: pass

f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/%s_%03d.root' % (runname, filenum))

def getPulseRecord(kdatareader, channelName):
  for p in kdatareader.GetEvent().boloPulseRecords():
    if p.GetChannelName() == channelName: return p
  return None


for chan in chanList:
  print chan
  pulseLength, binSize, _blah, _blah2, pulseType = getPulseInfo(chan)

  vr = db.view('analytical/bychandate',reduce=False, descending=True, startkey=[chan, "2012-01-22 00:00:00.0"], limit=1, include_docs=True)
  doc = vr.first()['doc']
  vp = std.vector("double")()

  exec(doc['formula']['python']) #defines 'template' function
  
  #doc['formula']['par'][0]=300  #put it close to zero, but away from the windowing function
  if chan.startswith("chal") == False: aPol = pulsePol.GetExpectedPolarity( getPulseRecord(f,chan))
  else: aPol = 1  #heat pulses are already negative polarity in the database. should i change this??
  print chan, aPol
  for i in range( pulseLength ):
    vp.push_back( aPol*template(i*binSize, doc['formula']['par']))

  scaleFactor = 1./max(abs(np.array(vp)))

  print 'scaling by', scaleFactor
  
  for i in range(vp.size()):
    vp[i] = scaleFactor * vp[i]
    
  plt.plot(np.array(vp))
  #raw_input('... continue')
  
  cham.SetTemplate(chan, vp, -1 * int(doc['formula']['par'][0]/binSize+ 0.5), pulseType)
  
  #plot the power spectrum
  hc2p.SetInputPulse(cham.GetTemplateSpectrum(chan))
  hc2p.RunProcess()
  #pw = kutil.get_out(hc2p)
  pwh = TH1D('p', 'p', hc2p.GetOutputPulseSize(), 0, hc2p.GetOutputPulseSize())
  for i in range(hc2p.GetOutputPulseSize()):
    pwh.SetBinContent(i+1, hc2p.GetOutputPulse()[i])
       
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
  
  raw_input('... continue')
  plt.cla()
  del pwh
  



if display:

  scout(f,cham, chanList,maxEvents=1000)

  kutil.looppulse(f, name=None, pta=None, analysisFunction = kamp, maxEvents=1000, kdatafile=f, channellist=chanList, 
    chamonix=cham, halfcomp2power=hc2p, real2halfcomp=r2hc, templateDB=db, pulsePol=pulsePol)

  #kamp(f, cham, chanList)
  f.Close()

else:
  k = KAmpKounselor()
  k.AddKAmpSite(cham)
  k.RunKamp('/Users/adam/analysis/edelweiss/data/kdata/raw/%s_%03d.root' % (runname, filenum), '/Users/adam/analysis/edelweiss/data/kdata/raw/%s_%03d.amp.root' % (runname, filenum))
  for chan in chanList:
    print chan
    print cham.GetNumNoiseEventsFound(chan), 'noise events found'
  



