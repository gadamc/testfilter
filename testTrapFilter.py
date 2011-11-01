#!/usr/bin/env python

from ROOT import TSystem, KDataReader, KTrapezoidalKamper, KPulseAnalysisRecord, KTrapKamperProto
import sys
import matplotlib.pyplot as plt
import numpy as np

global trap, p, r, entry, j, k, b

def runPlotTrap(clear = True):
  global trap, p, r, entry, j, k, b
  
  print p, r
  kampResult = trap.MakeKamp(p,r)
  print 'kampresult', kampResult
  print 'event, bolo, pulse'
  print  entry, j, k
  print 'bolo', b.GetDetectorName()
  print 'pulse', p.GetChannelName()
  print 'amplitude', r.GetAmp()
  print 'peak position', r.GetPeakPosition()
  for xx in range(6):
    print 'extra', xx, r.GetExtra(xx)

    #figObj.clf()
  
  
  pulse = np.array(p.GetTrace())
  
  plt.subplot(3,1,1)
  if clear: plt.cla()
  plt.plot(pulse)
  plt.title('raw')
    
  trapAmp = []
  trapTime = []
  tt = trap.GetTrapIonTime()
  ta = trap.GetTrapIonAmplitude()   
            
  if p.GetIsHeatPulse():
    tt = trap.GetTrapHeatTime()
    ta = trap.GetTrapHeatAmplitude()
      
  for i in range(tt.GetOutputPulseSize()):
    trapTime.append(tt.GetOutputPulse()[i])
      
  for i in range(ta.GetOutputPulseSize()):
    trapAmp.append(ta.GetOutputPulse()[i])
  
  plt.subplot(3,1,2)
  if clear: plt.cla()
  plt.plot(np.array(trapTime))
  plt.title('pulse time trap')
  
  plt.subplot(3,1,3)
  if clear: plt.cla()
  plt.plot(np.array(trapAmp))
  plt.title('pulse amp trap')
  plt.show()
  
  

if __name__ == '__main__':
  
  global trap, p, r, entry, j, k, b
  
  
  f = KDataReader(sys.argv[1])
  e = f.GetEvent()
  figObj = plt.figure(1)
  plt.ion()
  trap = KTrapKamperProto()

  for entry in range(f.GetEntries()):

    f.GetEntry(entry)

    for j in range(e.GetNumBolos()):

      b = e.GetBolo(j)
   
      if b.GetDetectorName() == 'ID3':
        for k in range(b.GetNumPulseRecords()):
          p = b.GetPulseRecord(k)
      
          r = KPulseAnalysisRecord()
      
          #runPlotTrap(trap, p, r, entry, j, k, b)
          runPlotTrap()
          
          go = raw_input('continue? (y/n/s/t): ')
          print 'you chose: ', str(go)
          
          #figObj.close()
          
          
          tt = ta = None
            
          while go == 't':
            if p.GetIsHeatPulse():
              tt = trap.GetTrapHeatTime()
              ta = trap.GetTrapHeatAmplitude()
              
            else:
              tt = trap.GetTrapIonTime()
              ta = trap.GetTrapIonAmplitude()
             
            
            d_orig_t = tt.GetDecayTimeConstant()
            ri_orig_t = tt.GetRiseTime()
            ff_orig_t = tt.GetFlatTopWidth()
            
            vals = ''
            vals = raw_input('Set TrapFilter for Time (decay %f, rise %d, flat %d): ' % (d_orig_t, ri_orig_t, ff_orig_t)) 
            if vals != '':
              d = float(vals.split(' ')[0])
              ri = int(vals.split(' ')[1])
              ff = int(vals.split(' ')[2])
              print str(d), ri, ff
              tt.SetParams(d,ri,ff)
             
            d_orig_a = ta.GetDecayTimeConstant()
            ri_orig_a = ta.GetRiseTime()
            ff_orig_a = ta.GetFlatTopWidth()
            
            vals = ''
            vals = raw_input('Set TrapFilter for Amp (decay %f, rise %d, flat %d): ' % (d_orig_a, ri_orig_a, ff_orig_a)) 
            if vals != '':
              d = float(vals.split(' ')[0])
              ri = int(vals.split(' ')[1])
              ff = int(vals.split(' ')[2])
              print str(d), ri, ff
              ta.SetParams(d,ri,ff)
             
            vals = raw_input('Clear? (y/n default = n)')
            if vals == 'y':
              runPlotTrap()
            else:
              runPlotTrap(False)
            go = ''
            go = raw_input('continue? (y/n/s/t): ')
            print 'you chose: ', str(go)
          
          if tt != None:
            tt.SetParams(d_orig_t,ri_orig_t,ff_orig_t) 
          if ta != None:
            ta.SetParams(d_orig_a,ri_orig_a,ff_orig_a)
            
          if go == 's':
            break
                
          if go == 'n':
            print 'quiting'
            sys.exit(0)
            
