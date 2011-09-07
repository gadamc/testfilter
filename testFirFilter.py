#!/usr/bin/env python
import sys
from ROOT import *
from couchdbkit import Server, Database
import matplotlib.pyplot as plt
import random
import scipy.signal as sig


  
def main(*args):

  tr = {}
  
  whitenoise = std.vector("double")()
  for i in range(8196):
    whitenoise.push_back(random.gauss(0,1))
  
  wn = []  #white noise
  for i in range(whitenoise.size()):
    wn.append(whitenoise[i])
  
  tr['whitenoise'] = wn
  
  r2hc = KRealToHalfComplexDFT()
  r2hc.SetInputPulse(whitenoise)
  print 'real to half complex white noise', r2hc.RunProcess()
  
  hcp = KHalfComplexPower()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'power white noise', hcp.RunProcess()
  
  wnp = []  #white noise power
  for i in range(hcp.GetOutputPulseSize()):
    wnp.append(hcp.GetOutputPulse()[i])
  
  tr['whitenoisepower'] = wnp
  
  coef = sig.firwin(50, 0.0001) #10-tap fir filter with a cut-off at 0.1 nyquist
  tr['fircoefficients'] = coef
  
  filter = KFIRFilter()
  coeff = std.vector("double")()
  for i in range(len(coef)):
    coeff.push_back(coef[i])
  
  filter.SetInputPulse(whitenoise)
  filter.SetCoefficients(coeff)
  print 'run filter', filter.RunProcess()
  
  filteredpulse = []
  for i in range(filter.GetOutputPulseSize()):
    filteredpulse.append(filter.GetOutputPulse()[i])
  
  tr['filteredpulse'] = filteredpulse
  
  r2hc.SetInputPulse(filter.GetOutputPulse(), filter.GetOutputPulseSize())
  r2hc.RunProcess()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'calculate power of filtered noise', hcp.RunProcess()
  
  fp = []
  for i in range(hcp.GetOutputPulseSize()):
    fp.append(hcp.GetOutputPulse()[i])
  
  tr['filteredpower'] = fp
  
  
  return tr

if __name__ == '__main__':
  main(*sys.argv[1:])