#!/usr/bin/env python

from ROOT import *
import sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig


global figObj

global heat_decay_fast, heat_rt_fast, heat_ft_fast
global ion_decay_fast, ion_rt_fast, ion_ft_fast
global heat_decay_slow, heat_rt_slow, heat_ft_slow
global ion_decay_slow, ion_rt_slow, ion_ft_slow

global trapKamp

global mostRecentIonPulseStartTime

def plotTheTraps(p, trapTop, toptitle, trapBottom, bottomtitle):
  
  global trapKamp
  
  trapBottomOut = []
  trapTopOut = []
      
  for i in range(trapTop.GetOutputPulseSize()):
    trapTopOut.append(trapTop.GetOutputPulse()[i])
      
  for i in range(trapBottom.GetOutputPulseSize()):
    trapBottomOut.append(trapBottom.GetOutputPulse()[i])
  
  plt.subplot(6,1,3)
  plt.cla()
  plt.plot(np.array(trapTopOut))
  plt.title(toptitle)
  
  diffone = np.zeros(len(trapTopOut))
  for i in range(1, len(trapTopOut)):
    diffone[i] = trapTopOut[i] - trapTopOut[i-1]
  
  difftwo = np.zeros(len(diffone))
  for i in range(1, len(diffone)):
    difftwo[i] = diffone[i] - diffone[i-1]
  
  plt.subplot(6,1,4)
  plt.cla()
  plt.plot(difftwo)
  plt.title(toptitle + ' second derivative')
    
  
  polCalc = KPulsePolarityCalculator()
  trapKamp.ClearPeakPositionResult()
  trapKamp.ResizePeakPositionResult( trapTop.GetOutputPulseSize() )
  o1 = KOrderFilter()
  o2 = KOrderFilter()
  o1.SetOrder(1)
  o2.SetOrder(1)
  o1.SetInputPulse(trapTop.GetOutputPulse(), trapTop.GetOutputPulseSize())
  o1.RunProcess()
  o2.SetInputPulse(o1.GetOutputPulse(), o1.GetOutputPulseSize())
  o2.RunProcess()
  
  trapKamp.FillPeakPositionResult(o2, trapTop, polCalc.GetExpectedPolarity(p))  
  
  peakPosResult = trapKamp.GetPeakPositionResult()
  
  plt.subplot(6,1,5)
  plt.cla()
  plt.plot(np.array(peakPosResult))
  plt.title('peakPos')
  
  plt.subplot(6,1,6)
  plt.cla()
  plt.plot(np.array(trapBottomOut))
  plt.title(bottomtitle)
  plt.show()
  
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
    return (d, ri, ff)
  else:
    return (d_orig_t, ri_orig_t, ff_orig_t)
      
def changePeakAmplitudeFinder():
  global trapKamp
  val = raw_input('Set peak finder amplitude (%f): ' % trapKamp.GetPeakPositionSearchAmplifier()) 
  if val != '':
    trapKamp.SetPeakPositionSearchAmplifier(float(val))
#
     
def plotPulse(p):
  
  global figObj
  
  global heat_decay_fast, heat_rt_fast, heat_ft_fast
  global ion_decay_fast, ion_rt_fast, ion_ft_fast
  global heat_decay_slow, heat_rt_slow, heat_ft_slow
  global ion_decay_slow, ion_rt_slow, ion_ft_slow
  
  global mostRecentIonPulseStartTime
  
  r = KPulseAnalysisRecord()
  bas = KBaselineRemoval()
  lin = KLinearRemoval()
  lin.SetBaselineStop(0.20)
  pat = KPatternRemoval()
  pat.SetPatternLength(200)
  (b,a) = sig.iirfilter(4,0.05, btype='lowpass')
  
  iir = KIIRFourthOrder(a[1], a[2], a[3], a[4], b[0], b[1], b[2], b[3], b[4])
  trapFast = KTrapezoidalFilter()
  trapSlow = KTrapezoidalFilter()
  graphFit = TGraph(p.GetPulseLength())
  polCalc = KPulsePolarityCalculator()
  theSign = "+"
  if polCalc.GetExpectedPolarity(p) == -1:
    theSign = "-"
  theForm = '[0] %s [1]/(1 + exp(-2*[2]*(x - [3])))' % theSign
  
  f1 = TF1('bbv2',theForm, 2000, 6000)
  f1.SetParameter(0, 0)
  f1.SetParameter(1, 1)
  f1.SetParameter(2, 0.149)
  #f1.FixParameter(2, 0.149)
  #f1.SetParLimits(2, 0, 0.3)
  f1.SetParameter(3, 4100)
  f1.SetParLimits(3, 2000, 6000)
  
  ptaToTraps = pat
  
  if p.GetIsHeatPulse() == False and p.GetBoloBoxVersion() < 2.0: #why is there no pattern on bbv2
    print ' bbv1 '
    bas.SetInputPulse(p.GetTrace())
    bas.RunProcess()
    pat.SetPatternLength( 200 )
    pat.SetInputPulse(bas.GetOutputPulse(), bas.GetOutputPulseSize())
    #print 'here'
    pat.RunProcess()
    pat.SetInputPulse(pat.GetOutputPulse(), pat.GetOutputPulseSize())
    pat.SetPatternLength(400)
    #print 'here'
    pat.RunProcess()
    ptaToTraps = pat
    #print 'here'
    
  elif p.GetIsHeatPulse() == False and p.GetBoloBoxVersion >= 2.0:
    print 'lowpass filter'
    lin.SetInputPulse(p.GetTrace())
    lin.RunProcess()
    iir.SetInputPulse(lin.GetOutputPulse(), lin.GetOutputPulseSize())
    iir.RunProcess()
    ptaToTraps = iir
    print 'b', b
    print 'a', a
    
    #pat.SetPatternLength( 10 )
    #pat.SetInputPulse(bas.GetOutputPulse(), bas.GetOutputPulseSize())
    #pat.RunProcess()
    
    for i in range(ptaToTraps.GetOutputPulseSize()):
      graphFit.SetPoint(i, i, ptaToTraps.GetOutputPulse()[i])
    
    res = graphFit.Fit(f1, 'SR')
    rep = res.Get()
    #print rep.Print('V')
    print ''
    print rep.Status()
    if rep.Status == 0:
      mostRecentIonPulseStartTime = f1.GetParmater(3)
    
    print theForm
    
  else:
    bas.SetInputPulse(p.GetTrace())
    bas.RunProcess()
    ptaToTraps = bas
    
  if p.GetIsHeatPulse():
    trapFast.SetParams(heat_decay_fast, heat_rt_fast, heat_ft_fast)
    trapSlow.SetParams(heat_decay_slow, heat_rt_slow, heat_ft_slow)
    
  else:
    trapFast.SetParams(ion_decay_fast, ion_rt_fast, ion_ft_fast)
    trapSlow.SetParams(ion_decay_slow, ion_rt_slow, ion_ft_slow)
    
  print 'here'
  trapFast.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
  trapFast.RunProcess()
  trapSlow.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
  trapSlow.RunProcess()
  print 'here'
  
  #plot the results and output information to terminal
  print p.GetChannelName()
  print 'fast trap (decay, rise, flat):', trapFast.GetDecayTimeConstant(), trapFast.GetRiseTime(), trapFast.GetFlatTopWidth()
  print 'slow trap (decay, rise, flat):', trapSlow.GetDecayTimeConstant(), trapSlow.GetRiseTime(), trapSlow.GetFlatTopWidth()
  pulse = np.array(p.GetTrace())
  
  plt.subplot(6,1,1)
  plt.cla()
  plt.plot(pulse)
  plt.title('raw')
  
  plt.subplot(6,1,2)
  plt.cla()
  clean = []
  for i in range(ptaToTraps.GetOutputPulseSize()):
    clean.append(ptaToTraps.GetOutputPulse()[i])
    
  plt.plot(np.array(clean))
  plt.title('clean')
  
  if p.GetIsHeatPulse() == False and p.GetBoloBoxVersion >= 2.0:
    'plotting fit for bbv2'
    fitArr = np.zeros(p.GetPulseLength())
    for i in range(p.GetPulseLength()):
      fitArr[i] = f1.Eval(i)
    plt.plot(fitArr)
  
  
  plotTheTraps(p, trapFast, 'fast trap', trapSlow, 'slow trap')
  
  
  #change trap parameters and plot as desired.
  go = raw_input('Change Trap Parameters?(y/n/s/quit)')
  while(go == 'y'):
    
    if p.GetIsHeatPulse():
      heat_decay_fast, heat_rt_fast, heat_ft_fast = changeTrapParameters(trapFast, 'fast trap' )
      heat_decay_slow, heat_rt_slow, heat_ft_slow = changeTrapParameters(trapSlow, 'slow trap')
    else:
      ion_decay_fast, ion_rt_fast, ion_ft_fast = changeTrapParameters(trapFast, 'fast trap')
      ion_decay_slow, ion_rt_slow, ion_ft_slow = changeTrapParameters(trapSlow, 'slow trap')
      
    trapSlow.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
    trapSlow.RunProcess()
    trapFast.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
    trapFast.RunProcess()
    
    changePeakAmplitudeFinder()
    #plot the results
    plotTheTraps(p, trapFast, 'fast trap', trapSlow, 'slow trap')
    
    go = raw_input('Change Trap Parameters?(y/n/s/quit)')
    
  return go
  
def runApp(*argv):
  global figObj
  
  global heat_decay_fast, heat_rt_fast, heat_ft_fast
  global ion_decay_fast, ion_rt_fast, ion_ft_fast
  global heat_decay_slow, heat_rt_slow, heat_ft_slow
  global ion_decay_slow, ion_rt_slow, ion_ft_slow
  
  global trapKamp
  
  global mostRecentIonPulseStartTime
  mostRecentIonPulseStartTime =  4096  
  trapKamp = KTrapKamperProto()
  
  heat_decay_fast = 20.0
  heat_rt_fast = 15
  heat_ft_fast = 60
  ion_decay_fast = 1000.0
  ion_rt_fast = 15
  ion_ft_fast = 70
  
  heat_decay_slow = 20.0
  heat_rt_slow = 10
  heat_ft_slow = 30
  ion_decay_slow = 1000.0
  ion_rt_slow = 15
  ion_ft_slow = 1000
  
  print argv
  f = KDataReader(argv[0])
  e = f.GetEvent()
  figObj = plt.figure(1)
  plt.ion()
  
  for entry in range(f.GetEntries()):
    f.GetEntry(entry)
    print "event number in file", entry
    for j in range(e.GetNumBolos()):
      b = e.GetBolo(j)
      
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        if p.GetPulseLength() == 0:
          continue
          
        answer = plotPulse(p)
        if answer == 'quit':
          sys.exit(0)
        if answer == 's':
          break  #skip to the next bolometer 

       
if __name__ == '__main__':
  runApp(*sys.argv[1:])