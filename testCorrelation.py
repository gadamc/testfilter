#!/usr/bin/env python

from couchdbkit import Server, Database
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import random
import scipy.signal as sig
import copy
import pickle

def getNoisePulse(pl, amp):
  whitenoise = ROOT.std.vector("double")()

  for i in range(pl):
    whitenoise.push_back(random.gauss(0,amp))
  
  #now filter the white noise with an iir filter to make it more realistic
  order = 4
  (b,a) = sig.iirfilter(order,0.5, btype='lowpass')
  iir4 = ROOT.KIIRFourthOrder(a[1], a[2], a[3], a[4], b[0], b[1], b[2], b[3], b[4])
  iir4.SetInputPulse(whitenoise)
  print 'filter white noise', iir4.RunProcess()

  noise = ROOT.std.vector("double")()
  for i in range(iir4.GetOutputPulseSize()):
    noise.push_back(iir4.GetOutputPulse()[i])

  return noise

####
def chalA_FID803_template(x,par):
  import math
  xx = x
  if xx < par[0]: return 0
  else: return par[1]*(1 - math.exp(-(xx-par[0])/par[2]))*(math.exp(-(xx-par[0])/par[3]) + par[4]*math.exp(-(xx-par[0])/par[5]))
  

def getTemplate():
  '''
  returns a normalized template pulse packed in a tuple. 
  the first element of the tuple is a std.vector and
  the second element is a regular python arrays

  '''
  s = Server('https://edwdbik.fzk.de:6984')
  #db = s['analysis']
  #doc = db['run13_templatepulse_centre_FID804AB']
  #pulse = doc['pulse']
  db = s['pulsetemplates']
  doc = db['8b20ed91c4dee6782b47a899d89e6ba5']
  
  pulse = []
  for i in range(512):
    pulse.append(chalA_FID803_template(i, doc['formula']['double exponential heat template']['par']))
    
  integral = 0
  for i in range(len(pulse)):
    integral += pulse[i]
  
  
  t = ROOT.std.vector("double")()
  pt = copy.deepcopy(pulse)
  for i in range(len(pt)):
    pt[i] = pt[i]/abs(float(integral))
    t.push_back(pt[i])
    
  return (t, pt)
  
#####
random.seed()

(t, pt) = getTemplate()
noise = np.array(getNoisePulse(len(pt), max(pt)/3.))
signal = np.array(pt) + noise
#signal = noise 
ps = ROOT.KPulseShifter()
ps.SetShift(-20)
ps.SetInputPulse(signal, len(signal))
ps.RunProcess()

plt.subplot(2,1,1)
plt.ion()
#plt.plot(signal)
#raw_input()

plt.plot(np.array(pt))
plt.plot(np.array(pt[250:350]))
#raw_input()

for i in range(ps.GetOutputPulseSize()):
  signal[i] = ps.GetOutputPulse()[i]
  

plt.plot(signal)

cor = ROOT.KCorrelation()
cor.SetResponse(t[250:350])
cor.SetInputPulse(signal, len(signal))
if cor.RunProcess():
  out= np.zeros(cor.GetOutputPulseSize())
  for i in range(cor.GetOutputPulseSize()):
    out[i] = cor.GetOutputPulse()[i]
    
  plt.subplot(2,1,2)    
  plt.cla()
  plt.plot(out)
#  raw_input()
  
  #out = sig.convolve(signal, np.array(pt), 'same')
  #   #plt.cla()
  #plt.plot(out)

  #   
  out = sig.convolve(signal, np.array(pt[250:350]), 'full')
  #   #plt.cla()
  plt.plot(out)
  # raw_input()
  
  out = sig.correlate(signal, np.array(pt[250:350]), 'valid')
  #plt.cla()
  plt.plot(out)
  
  #out = sig.correlate(signal, np.array(pt), 'full')
  #plt.cla()
  #plt.plot(out)
  
  #out = sig.correlate(signal, np.array(pt), 'valid')
  #plt.cla()
  #plt.plot(out)
  plt.show()
  raw_input()
  

else:
  print 'bahhh'
    
