#!/usr/bin/env python
from ROOT import *
import sys
import matplotlib.pyplot as plt
import numpy as np

global figObj
global trapKamp
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
def changePulseFinderTrapParameters(p):
  
  ### why does this break!±!!
  traps = getTimingTrapVector(p)
  for i in range(len(traps)):
    changeTrapParameters(traps[i])
       
#
def changePeakAmplitudeFinder():
  global trapKamp
  val = raw_input('Set peak finder amplitude (%f): ' % trapKamp.GetPeakPositionSearchAmplifier()) 
  if val != '':
    trapKamp.SetPeakPositionSearchAmplifier(float(val))
#
def plotTrapResults(p, trap, o2):

  global trapKamp
  global figObj
  
  polCalc = KPulsePolarityCalculator()
  print 'bias', p.GetPolarity()
  print 'polarity', polCalc.GetExpectedPolarity(p)
  
  cleanPulse = np.array([])
  for i in range(trap.GetInputPulseSize()):
    cleanPulse = np.append(cleanPulse, trap.GetInputPulse()[i])
  
  plt.subplot(4,1,1)
  plt.cla()
  plt.plot(cleanPulse)
  plt.title('cleaned')
  
  trapOut = np.array([])
  for i in range(trap.GetOutputPulseSize()):
    trapOut = np.append(trapOut, trap.GetOutputPulse()[i])
  
  plt.subplot(4,1,2)
  plt.cla()
  plt.plot(trapOut)
  plt.title('timing trap')
  
  
  orderOut = np.array([])
  for i in range(o2.GetOutputPulseSize()):
    orderOut = np.append(orderOut, o2.GetOutputPulse()[i])
        
  plt.subplot(4,1,3)
  plt.cla()
  plt.plot(orderOut)
  plt.title('second order')
  
  peakPosResult = trapKamp.GetPeakPositionResult()
  
  plt.subplot(4,1,4)
  plt.cla()
  plt.plot(np.array(peakPosResult))
  plt.title('peakPos')
  
  numNonZero = 0
  for val in peakPosResult:
    if val != 0:
      numNonZero += 1
  print 'num peaks', numNonZero
  
  
  
  plt.show()
#

def analyzeWithNewTrap(p):
  
  global trapKamp
  
  trap = KTrapezoidalFilter()
  if p.GetIsHeatPulse():
    trap.SetParams(20, 7, 10)
  else:
    trap.SetParams(400, 7, 40)
  
  changeTrapParameters(trap, 'testing this trap')
  
  trap.SetInputPulse(getAmplitudeTrap(p).GetInputPulse(), getAmplitudeTrap(p).GetInputPulseSize())
  trap.RunProcess()
  o1 = KOrderFilter()
  o2 = KOrderFilter()
  o1.SetOrder(1)
  o2.SetOrder(1)
  o1.SetInputPulse(trap.GetOutputPulse(), trap.GetOutputPulseSize())
  o1.RunProcess()
  o2.SetInputPulse(o1.GetOutputPulse(), o1.GetOutputPulseSize())
  o2.RunProcess()
  
  polCalc = KPulsePolarityCalculator()
  trapKamp.ClearPeakPositionResult()
  trapKamp.ResizePeakPositionResult( getAmplitudeTrap(p).GetInputPulseSize() )
  trapKamp.FillPeakPositionResult(o2, trap, polCalc.GetExpectedPolarity(p))
  
  plotTrapResults(p, trap, o2)
  
  

#
def plot(p, r):

  global trapKamp
  global figObj
  
  polCalc = KPulsePolarityCalculator()
  print 'bias', p.GetPolarity()
  print 'polarity', polCalc.GetExpectedPolarity(p)
  
  pulse = np.array(p.GetTrace())
  
  plt.subplot(4,1,1)
  plt.cla()
  plt.plot(pulse)
  plt.title('raw')
  
  trapAmp = getAmplitudeTrap(p)
    
  cleanPulse = np.array([])
  for i in range(trapAmp.GetInputPulseSize()):
    cleanPulse = np.append(cleanPulse, trapAmp.GetInputPulse()[i])
  
  plt.subplot(4,1,2)
  plt.cla()
  plt.plot(cleanPulse)
  plt.title('cleaned')
  
  peakPosResult = trapKamp.GetPeakPositionResult()
  
  plt.subplot(4,1,3)
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
      
  plt.subplot(4,1,4)
  plt.cla()
  plt.plot(ampPulse)
  plt.title('amplitude trapezoid')
  
  print 'results:'
  print 'amp / peak position', r.GetAmp(), '/', r.GetPeakPosition()
  
  plt.show()
#
def analyzePulse(p):

  global trapKamp
  r = KPulseAnalysisRecord()
  print ''
  print p.GetChannelName()
  print 'make kamp returned', trapKamp.MakeKamp(p,r)
  plot(p,r)
      
  go = raw_input('Change Parameters?(y/N/s/quit)')
  while(go == 'y'):
    
    changeTrapParameters(getAmplitudeTrap(p), 'amplitude trap filter' )
    changePeakAmplitudeFinder()
    gogo = raw_input('Change Pulse Finder Trap Parameters? y/N')
    if gogo == 'y':
      changePulseFinderTrapParameters(p)
      
    #run analysis and plot the results
    print 'make kamp returned', trapKamp.MakeKamp(p,r)
    plot(p, r)
    
    go = raw_input('Change Parameters?(y/N/s/quit)')
  
  if go != 's' and go != 'quit':
    gogo = raw_input('Analyze clean pulse with new trap filter? y/N')
    while(gogo == 'y'):
      analyzeWithNewTrap(p)
      gogo = raw_input('Analyze clean pulse with new trap filter? y/N')
    
  return go
  
#
def runApp(*argv):
  global figObj
  global trapKamp
  
  figObj = plt.figure(1)
  trapKamp = KTrapKamperProto()
  plt.ion()
  
  f = KDataReader(argv[0])
  e = f.GetEvent()
  t = f.GetTTree()
  
  for entry in t:

    for j in range(e.GetNumBolos()):
      b = e.GetBolo(j)
            
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        
        answer = analyzePulse(p)
        if answer == 'quit':
          sys.exit(0)
        if answer == 's':
          break  #skip to the next bolometer 
       
if __name__ == '__main__':
  runApp(*sys.argv[1:])