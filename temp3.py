#!/usr/bin/env python
# -*- coding: utf-8 -*-
from couchdbkit import Server, Database
from ROOT import *
import matplotlib.pyplot as plt
import pywt
import numpy as np
import string, os
import KDataPy.util as kutil
plt.ion()


myChannel = 'chalB Gc1'

mydatadir = '/kalinka/home/luo/data'
filelist = ['ma23a003_013.root']

try:
  lvl=int(raw_input('enter decomposition level: '))
except ValueError:
    print "Not a number"

for aFile in filelist:
  #open the file
  f = KDataReader( os.path.join(mydatadir,aFile) )
  e = f.GetEvent()  #get the KRawEvent object
  bas = KBaselineRemoval()
  bas.SetBaselineStart(0)	
  bas.SetBaselineStop(0.4)
  
  #now loop through the data again, applying the optimal filter that was just built.
  #for i in range(f.GetEntries()):
    
   # f.GetEntry(i) 
  
  for ii in range(e.GetNumBoloPulses()):
      
    p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
    
    if p.GetPulseLength() == 0: continue  #skip empty pulses
    
    if p.GetChannelName() == myChannel:    
      bas.SetInputPulse( p.GetTrace() )
      y=bas.RunProcess()
     #plt.subplot(3,3,ii+1)
      x =  kutil.get_out(bas)
      
      coeffs = []
      for i in range(lvl):
	
        coeffs.append(pywt.wavedec(x, 'db2', level = i+1))
      
      
        if max(abs(coeffs[i][0])) > 2*max(abs(coeffs[i][1])):
	  print 'for level', i+1, ':true'
        else:
	  print 'for level', i+1, ':false'
	
        plt.figure(2)
        plt.subplot(lvl+1,2,1)
        plt.plot(x)
        plt.subplot(lvl+1,2,2*i+3)
        plt.plot(coeffs[i][0])
        plt.subplot(lvl+1,2,2*i+4)
        plt.plot(coeffs[i][1])
         
      #plt.figure(1)
      #plt.subplot(5,2,1)
      #plt.plot(x)
      #plt.subplot(5,2,3)
      #plt.plot(cA1)
      #plt.subplot(5,2,4)
      #plt.plot(cD11)
      #plt.subplot(5,2,5)
      #plt.plot(cA2)
      #plt.subplot(5,2,6)
      #plt.plot(cD22)
      #plt.subplot(5,2,7)
      #plt.plot(cA3)
      #plt.subplot(5,2,8)
      #plt.plot(cD33)
      #plt.subplot(5,2,9)
      #plt.plot(cA4)
      #plt.subplot(5,2,10)
      #plt.plot(cD44)
      
      
      plt.show()
      
      
	
      
#      plt.plot( np.array(p.GetTrace()) )
#      raw_input('hit enter to close')
     # plt.cla() #clears the plot