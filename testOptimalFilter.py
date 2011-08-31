#!/usr/bin/env python
import sys
from ROOT import *
from couchdbkit import Server, Database
import matplotlib.pyplot as plt
import random, math

class testresults(object):
  def __init__(self):

def getNoise():
  tr = {}
  random.seed()
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

  order = 4
  (b,a) = sig.iirfilter(order,0.1, btype='lowpass')
  tr['b'] = b
  tr['a'] = a
  tr['filter order'] =  order

  filter = KIIRFourthOrder(-a[1], -a[2], -a[3], -a[4], b[0], b[1], b[2], b[3], b[4])
  filter.SetInputPulse(whitenoise)
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
  
def main(*args):
  s = Server('https://edwdbik.fzk.de:6984')
  db = s['analysis']
  doc = db['run13_templatepulse_centre_FID804AB']
  
  
  templatepulse = std.vector("double")()
  maxvalue = abs(min(doc['pulse']))
  for i in doc['pulse']:
    templatepulse.push_back(i/maxvalue)
  
  
  tempNoise = getNoise()
  
  maxfilter = max(tempNoise['filteredpulse'])
  #now we add the filteredpulse, scaled to smaller values, to the template pulse to 
  
  return tr

if __name__ == '__main__':
  main(*sys.argv[1:])