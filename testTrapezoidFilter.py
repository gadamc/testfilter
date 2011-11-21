#!/usr/bin/env python

from ROOT import *
import sys
import matplotlib.pyplot as plt
import numpy as np

global figObj

global heat_decay_fast, heat_rt_fast, heat_ft_fast
global ion_decay_fast, ion_rt_fast, ion_ft_fast
global heat_decay_slow, heat_rt_slow, heat_ft_slow
global ion_decay_slow, ion_rt_slow, ion_ft_slow

def plotTheTraps(trapTop, toptitle, trapBottom, bottomtitle):
  
  trapBottomOut = []
  trapTopOut = []
      
  for i in range(trapTop.GetOutputPulseSize()):
    trapTopOut.append(trapTop.GetOutputPulse()[i])
      
  for i in range(trapBottom.GetOutputPulseSize()):
    trapBottomOut.append(trapBottom.GetOutputPulse()[i])
  
  plt.subplot(5,1,3)
  plt.cla()
  plt.plot(np.array(trapTopOut))
  plt.title(toptitle)
  
  diffone = np.zeros(len(trapTopOut))
  for i in range(1, len(trapTopOut)):
    diffone[i] = trapTopOut[i] - trapTopOut[i-1]
  
  difftwo = np.zeros(len(diffone))
  for i in range(1, len(diffone)):
    difftwo[i] = diffone[i] - diffone[i-1]
  
  plt.subplot(5,1,4)
  plt.cla()
  plt.plot(difftwo)
  plt.title(toptitle + ' second derivative')
    
    
  plt.subplot(5,1,5)
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
      
    
def plotPulse(p):
  
  global figObj
  
  global heat_decay_fast, heat_rt_fast, heat_ft_fast
  global ion_decay_fast, ion_rt_fast, ion_ft_fast
  global heat_decay_slow, heat_rt_slow, heat_ft_slow
  global ion_decay_slow, ion_rt_slow, ion_ft_slow
  
  r = KPulseAnalysisRecord()
  bas = KBaselineRemoval()
  pat = KPatternRemoval()
  trapFast = KTrapezoidalFilter()
  trapSlow = KTrapezoidalFilter()
  
  bas.SetInputPulse(p.GetTrace())
  bas.RunProcess()
  
  ptaToTraps = pat
  
  if p.GetIsHeatPulse() == False:
    pat.SetPatternLength( int(p.GetHeatPulseStampWidth()) )
    pat.SetInputPulse(bas.GetOutputPulse(), bas.GetOutputPulseSize())
    pat.RunProcess()
    pat.SetInputPulse(pat.GetOutputPulse(), pat.GetOutputPulseSize())
    pat.SetPatternLength(int(2*p.GetHeatPulseStampWidth()))
    pat.RunProcess()
  else:
    ptaToTraps = bas
    
  if p.GetIsHeatPulse():
    trapFast.SetParams(heat_decay_fast, heat_rt_fast, heat_ft_fast)
    trapSlow.SetParams(heat_decay_slow, heat_rt_slow, heat_ft_slow)
    
  else:
    trapFast.SetParams(ion_decay_fast, ion_rt_fast, ion_ft_fast)
    trapSlow.SetParams(ion_decay_slow, ion_rt_slow, ion_ft_slow)
    
  trapFast.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
  trapFast.RunProcess()
  trapSlow.SetInputPulse(ptaToTraps.GetOutputPulse(), ptaToTraps.GetOutputPulseSize())
  trapSlow.RunProcess()
  
  
  #plot the results and output information to terminal
  print p.GetChannelName()
  print 'fast trap (decay, rise, flat):', trapFast.GetDecayTimeConstant(), trapFast.GetRiseTime(), trapFast.GetFlatTopWidth()
  print 'slow trap (decay, rise, flat):', trapSlow.GetDecayTimeConstant(), trapSlow.GetRiseTime(), trapSlow.GetFlatTopWidth()
  pulse = np.array(p.GetTrace())
  
  plt.subplot(5,1,1)
  plt.cla()
  plt.plot(pulse)
  plt.title('raw')
  
  plt.subplot(5,1,2)
  plt.cla()
  clean = []
  for i in range(ptaToTraps.GetOutputPulseSize()):
    clean.append(ptaToTraps.GetOutputPulse()[i])
    
  plt.plot(np.array(clean))
  plt.title('clean')
  
  plotTheTraps(trapFast, 'fast trap', trapSlow, 'slow trap')
  
  
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
    
    #plot the results
    plotTheTraps(trapFast, 'fast trap', trapSlow, 'slow trap')
    
    go = raw_input('Change Trap Parameters?(y/n/s/quit)')
    
  return go
  
def runApp(*argv):
  global figObj
  
  global heat_decay_fast, heat_rt_fast, heat_ft_fast
  global ion_decay_fast, ion_rt_fast, ion_ft_fast
  global heat_decay_slow, heat_rt_slow, heat_ft_slow
  global ion_decay_slow, ion_rt_slow, ion_ft_slow
  
  heat_decay_fast = 50.0
  heat_rt_fast = 3
  heat_ft_fast = 10
  ion_decay_fast = 1100.0
  ion_rt_fast = 3
  ion_ft_fast = 10
  
  heat_decay_slow = 50.0
  heat_rt_slow = 10
  heat_ft_slow = 20
  ion_decay_slow = 1100.0
  ion_rt_slow = 50
  ion_ft_slow = 200
  
  print argv
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