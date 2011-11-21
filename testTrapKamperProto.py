#!/usr/bin/env python
from ROOT import *
import sys
import matplotlib.pyplot as plt
import numpy as np

global figObj
global trapKamp = KTrapKamperProto
#
def getAmplitudeTrap(p):
  global trapKamp
  if p.GetIsHeatPulse():  return trapKamp.GetTrapHeatAmplitude()
  else:   return trapKamp.GetTrapIonAmplitude()
    
#
def getTimingTrapVector(p):
  global trapKamp
  if p.GetIsHeatPulse():  return trapKamp.GetTrapHeatTime()
  else:  return trapKamp.GetTrapIonTime()
    
#
def changeTrapParameters(trap, message = None):
  d_orig_t = trap.GetDecayTimeConstant()
  ri_orig_t = trap.GetRiseTime()
  ff_orig_t = trap.GetFlatTopWidth()
  
  if message:
    print message
    
  vals = raw_input('Set trap parameters (decay %f, rise %d, flat %d): ' % (d_orig_t, ri_orig_t, ff_orig_t)) 
  if vals != '':
    d = float(vals.split(' ')[0])
    ri = int(vals.split(' ')[1])
    ff = int(vals.split(' ')[2])
    print str(d), ri, ff
    trap.SetParams(d,ri,ff)
    # return (d, ri, ff)
    #   else:
    #     return (d_orig_t, ri_orig_t, ff_orig_t)

#
def changePulseFinderTrapParameters():
  
  for trap in getTimingTrapVector()
    changeTrapParameters(trap)
       
#
def plot(p, r):

  global trapKamp
  
  print p.GetChannelName()
  polCalc = KPulsePolarityCalculator()
  print 'bias', p.GetPolarity()
  print 'polarity', polCalc.GetExpectedPolarity(p)
  
  pulse = np.array(p.GetTrace())
  
  plt.subplot(5,1,1)
  plt.cla()
  plt.plot(pulse)
  plt.title('raw')
  
  trapAmp = getAmplitudeTrap(p)
    
  cleanPulse = np.array([])
  for i in range(trapAmp.GetInputPulseSize()):
    cleanPulse = np.append(cleanPulse, trapAmp.GetInputPulse()[i])
  
  plt.subplot(5,1,2)
  plt.cla()
  plt.plot(cleanPulse)
  plt.title('cleaned')
  
  peakPosResult = trapKamp.GetPeakPositionResult()
  
  plt.subplot(5,1,3)
  plt.cla()
  plt.plot(np.array(peakPosResult))
  plt.title('peakPos')
  
  numNonZero = 0
  for val in peakPosResult:
    if val != 0:
      numNonZero += 1
  print 'num peaks', numNonZero
  
  ampPulse = np.array([])
  for i in range(trapAmp.GetOutputPulseSize()):
    ampPulse = np.append(ampPulse, trapAmp.GetOutputPulse()[i])
      
  plt.subplot(5,1,4)
  plt.cla()
  plt.plot(ampPulse)
  plt.title('amplitude trapezoid')
  
  print 'results:'
  print 'amp / peak position', r.GetAmp, '/', r.GetPeakPosition()
  
#
def analyzePulse(p):

  global trapKamp
  r = KPulseAnalysisRecord()
  trapKamp.MakeKamp(p,r)
  plot(p,r)
      
  go = raw_input('Change Parameters?(y/N/s/quit)')
  while(go == 'y'):
    
    changeTrapParameters(getAmplitudeTrap(p), 'amplitude trap filter' )
    changePeakAmplitudeFinder()
    gogo = raw_input('Change Pulse Finder Trap Parameters? y/N')
    if gogo == 'y':
      changePulseFinderTrapParameters()
      
    #run analysis and plot the results
    trapKamp.MakeKamp(p,r)
    plot(p, r)
    
    go = raw_input('Change Trap Parameters?(y/N/s/quit)')
    
  return go
  
#
def runApp(*argv):
  global figObj
 
  f = KDataReader(argv[0])
  e = f.GetEvent()
  figObj = plt.figure(1)
  plt.ion()
  
  for entry in range(f.GetEntries()):
    f.GetEntry(entry)

    for j in range(e.GetNumBolos()):
      b = e.GetBolo(j)
            
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        
        answer = plotPulse(p)
        if answer == 'quit':
          sys.exit(0)
        if answer == 's':
          break  #skip to the next bolometer 
       
if __name__ == '__main__':
  runApp(*sys.argv[1:])