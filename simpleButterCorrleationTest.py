#!/usr/bin/env python
import sys
from ROOT import *
from couchdbkit import Server, Database
import matplotlib.pyplot as plt
import random
import scipy.signal as sig
import math
import numpy as np
import copy

def getTemplate():
  '''
  returns a normalized template pulse packed in a tuple. 
  the first element of the tuple is a std.vector and
  the second element is a regular python arrays

  '''
  s = Server('https://edwdbik.fzk.de:6984')
  db = s['analysis']
  doc = db['run13_templatepulse_centre_FID804AB']
  
  integral = 0
  for i in range(len(doc['pulse'])):
    integral += doc['pulse'][i]
    
  t = std.vector("double")()
  pt = copy.deepcopy(doc['pulse'])
  for i in range(len(pt)):
    pt[i] = pt[i]/abs(float(integral))
    t.push_back(pt[i])
    
  return (t, pt)

def createSignal(pl, template, noisepulse, amplitude = 1.):
  signal = std.vector("double")()
  for i in range(pl):
    signal.push_back( amplitude*template[i-100] + noisepulse[i])
    
  return signal
  

def main(*argv):
  
  tr = {}
  ### create the pulse template and save it to the return
  (template, python_template) = getTemplate()
  tr['template'] = python_template

  pulseLength = len(python_template) #get the length from the template
  
  # get a noise pulse 
  # save it to the return
  random.seed()
  noisepulse = getNoisePulse(pulseLength)
  wn = []
  for i in range(pulseLength):
    wn.append(noisepulse[i])
    
  tr['noise'] = wn
  
  #create the signal by adding the template to the noise
  #save it to the return
  signal = createSignal(pulseLength, python_template, noisepulse, 10000.)
  signal2 = createSignal(pulseLength, python_template, noisepulse, 20000.)
  sig = []
  for i in range(pulseLength):
    sig.append(signal[i])

  tr['signal'] = sig
  
  
if __name__ == '__main__':
  main(*sys.argv[1:])