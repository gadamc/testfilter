#!/usr/bin/env python
import operator
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

#make sure you $> source /sps/edelweis/kdata/code/dev/config/setup_kdata.sh  (or .csh if you use tcsh or csh)

gSystem.Load('libkds')  
gSystem.Load('libkpta')
gSystem.Load('libkamping')

f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/ma22a000_000.root')  # all data in /sps/edelweis/kdata/data/raw. ex: /sps/edelweis/kdata/data/raw/ma22a000_000.root
myChannel = 'chalA FID808'
outputHistFile = 'myOutfile.root'

cham = KChamonixKAmpSite()

#these are the parameters of the KEraPeakFinder. See the webpage: https://edwdev-ik.fzk.de/kdata_dev_ref/KEraPeakFinder.html
cham.GetHeatPeakDetector().SetNumRms(2.3)  
cham.GetHeatPeakDetector().SetOrder(5)

go = ''
e = f.GetEvent()
for i in range(f.GetEntries()):
  f.GetEntry(i)
  if operator.mod(i,100) == 0: print i
  for j in range(e.GetNumBoloPulses()):
    pulse = e.GetBoloPulse(j)
    if pulse.GetChannelName() != myChannel:
      continue
    if cham.ScoutKampSite( pulse ,e) == True and go != 'skip': #KChamonixKAmpSite::ScoutKampSite() returns True when it finds a noise pulse
      plt.cla()
      plt.plot( np.array(pulse.GetTrace()) )
      go = raw_input('hit enter to continue, or type "skip" to continue without stopping on noise pulses')
      
print 'number of noise events found', cham.GetNumNoiseEventsFound(myChannel)
powerVec = cham.GetNoisePower(myChannel)
hist = TH1D('noisepower', 'noisepower', powerVec.size(), 0, powerVec.size())
for i in range(powerVec.size()):
  hist.SetBinContent(i+1, powerVec[i])

plt.cla()
plt.loglog( np.array(powerVec) )
raw_input('hit enter to close')

ff = TFile.Open(outputHistFile, 'recreate')
hist.Write()
ff.Close()

#open root file by starting ROOT
# $> root
# root [0] TFile f("myOutfile.root")
# root [1] f.ls()
# root [2] noisepower.Draw()
# root [3] c1.SetLogx()
# root [4] c1.SetLogy()
# root [5] noisepower.Draw()