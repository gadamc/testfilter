#!/usr/bin/env python

from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import random


def makePulse(slope, offset):
  pulse = std.vector("double")()
  for i in range(8196):
    pulse.push_back(slope * i + offset)
    
  return pulse
  
  
def main(*args):
  
  pulse = makePulse(0.5, 3000.)
  
  for i in range(len(pulse)):
    pulse[i] += random.gauss(0,50)
    
  #print 'offset ', pulse[0]
  
  lin = KLinearRemoval()
  lin.SetInputPulse(pulse)
  #print 'input pulse', lin.GetInputPulse()
  #print 'output pulse', lin.GetOutputPulse()
  lin.RunProcess()
  #print 'running linear removal', lin.RunProcess()
  #print 'slope', lin.CalculateSlope()
  
  output = np.zeros(len(pulse))
   
  for i in range(len(pulse)):
    output[i] = lin.GetOutputPulse()[i]
    
  plt.figure(1)
  plt.plot(pulse)
  plt.figure(2)
  plt.plot(output)

  
if __name__ == '__main__':
  main(*sys.argv[1:])