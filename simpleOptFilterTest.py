#!/usr/bin/env python
import sys
from ROOT import *
from couchdbkit import Server, Database
import matplotlib.pyplot as plt
import random
import scipy.signal as sig
import math
import numpy as np

def getgauss(a, b, c, x):
  return a * exp(-(x - b)**2 / (2*c**2))
  
def main(*args):

  tr = {}
  
  whitenoise = std.vector("double")()
  template = std.vector("double")()
  python_template = []
  for i in range(8196):
    whitenoise.push_back(random.gauss(0,1))
    template.push_back(getgauss(1./(50*sqrt(2*math.pi)),4096., 50., float(i))) #the template needs to have a maximum amplitude of 1
    #template.push_back(getgauss(10.,4096., 50., float(i))) #the template needs to have a maximum amplitude of 1
    python_template.append(getgauss(1./(50*sqrt(2*math.pi)),4096., 50., float(i)))
    #python_template.append(getgauss(10.,4096., 50., float(i)))
  
  signal = std.vector("double")()
  for i in range(8196):
    signal.push_back( 50*sqrt(2*math.pi)*1.*python_template[i-2000] + whitenoise[i])
    #signal.push_back( python_template[i-2000] + whitenoise[i])
    
  filter = KOptimalFilter()
  
  r2hc = KRealToHalfComplexDFT()
  r2hc.SetInputPulse(whitenoise)
  print 'real to half complex white noise', r2hc.RunProcess()
  hcp = KHalfComplexPower()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'power white noise', hcp.RunProcess()
  wp = []
  for i in range(hcp.GetOutputPulseSize()):
    wp.append(hcp.GetOutputPulse()[i])
    
  tr['noisepower'] = wp
  aveNoisePower = np.mean(wp)
  avewhitepower = std.vector("double")()
  wavepower= []
  for i in range(len(wp)):
    avewhitepower.push_back(wp[0])
    wavepower.append(aveNoisePower)
  
  tr['avenoisepower'] = wavepower
  
  #use ave white noise power
  filter.SetNoiseSpectrum(avewhitepower)
  #or use the white noise power at this instance
  #filter.SetNoiseSpectrum(hcp.GetOutputPulse(), hcp.GetOutputPulseSize())
  
  r2hc.SetInputPulse(template)
  print 'real to half complex template', r2hc.RunProcess()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'power template', hcp.RunProcess()
  tmppower = []
  for i in range(hcp.GetOutputPulseSize()):
    tmppower.append(hcp.GetOutputPulse()[i])
    
  tr['templatepower'] = tmppower
  
  filter.SetTemplateDFT(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  
  r2hc.SetInputPulse(signal)
  print 'real to half complex signal', r2hc.RunProcess()
  
  filter.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  
  print 'building filter', filter.BuildFilter()
  print 'applying filter', filter.RunProcess()
  
  output = []
  for i in range(filter.GetOutputPulseSize()):
    output.append(filter.GetOutputPulse()[i])
  
  tr['amplitude'] = output
  
  wn = []
  for i in range(8196):
    wn.append(whitenoise[i])
    
  tr['noise'] = wn
  
  tmp = []
  for i in range(8196):
    tmp.append(template[i])
  
  tr['template'] = tmp
  
  sig = []
  for i in range(8196):
    sig.append(signal[i])
  
  tr['signal'] = sig
  
  optfil = []
  for i in range(filter.GetOptimalFilterSize()):
    optfil.append(filter.GetOptimalFilter()[i])
  
  tr['optfilter'] = optfil
  
  hcp.SetInputPulse(filter.GetOptimalFilter(), filter.GetOptimalFilterSize())
  print 'get optimal filter power', hcp.RunProcess()
  
  optfilpower = []
  for i in range(hcp.GetOutputPulseSize()):
    optfilpower.append(hcp.GetOutputPulse()[i])
    
  tr['optfilpower'] = optfilpower
  

  return tr

if __name__ == '__main__':
  main(*sys.argv[1:])