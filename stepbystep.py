#!/usr/bin/env python

from ROOT import TSystem, KDataReader, KTrapezoidalKamper, KPulseAnalysisRecord
import sys
import matplotlib.pyplot as plt
import numpy as np


global e, f, r

f = KDataReader(sys.argv[1])
e = f.GetEvent()


figObj = plt.figure(1)

plt.ion()

trap = KTrapezoidalKamper()
trap.SetDebugMode()

global entry, j, k


for entry in range(f.GetEntries()):

  f.GetEntry(entry)

  for j in range(e.GetNumBolos()):

    b = e.GetBolo(j)
   
    if b.GetDetectorName() == 'ID3':
        for k in range(b.GetNumPulseRecords()):
          p = b.GetPulseRecord(k)
      
          r = KPulseAnalysisRecord()
      
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
          
          plt.subplot(3,3,1)
          plt.cla()
          plt.plot(pulse)
          plt.title('raw')
        
          vout = trap.GetDebugResults()
          vname = trap.GetDebugSteps()
          print 'size of debug results', vname.size()
      
          for i in range(8):
            plt.subplot(3,3,i+2)
            plt.cla()
            
          for i in range(vname.size()):
            plt.subplot(3,3,i+2)  
            plt.cla()
            plt.plot(np.array(vout[i]))
            plt.title(vname[i])
            print i, vname[i]
                               
          plt.show()
          
          go = raw_input('continue? (y/n/s): ')
          print 'you chose: ', str(go)
          
          #figObj.close()
          
          if go == 's':
            break
          
          if go == 'n':
            print 'quiting'
            sys.exit(0)
            
              