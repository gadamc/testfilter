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


    
def getNoisePulse(pl):
  whitenoise = std.vector("double")()
  
  for i in range(pl):
    whitenoise.push_back(random.gauss(0,1))
  
  #now filter the white noise with an iir filter to make it more realistic
  order = 4
  (b,a) = sig.iirfilter(order,0.5, btype='lowpass')
  iir4 = KIIRFourthOrder(a[1], a[2], a[3], a[4], b[0], b[1], b[2], b[3], b[4])
  iir4.SetInputPulse(whitenoise)
  print 'filter white noise', iir4.RunProcess()

  noise = std.vector("double")()
  for i in range(iir4.GetOutputPulseSize()):
    noise.push_back(iir4.GetOutputPulse()[i])

  return noise
  
  
def calculatePower(signal):
  r2hc = KRealToHalfComplexDFT()
  r2hc.SetInputPulse(signal)
  r2hc.RunProcess()
  hcp = KHalfComplexPower()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  hcp.RunProcess()
  #save it for the return
  sigp = []
  for i in range(hcp.GetOutputPulseSize()):
    sigp.append(hcp.GetOutputPulse()[i])

  return sigp
  
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

def createSignal(pl, template, noisepulse, amplitude = 1., delay = 2095):
  signal = std.vector("double")()
  for i in range(pl):
    if i >= delay:
      signal.push_back( amplitude*template[i-delay] + noisepulse[i])
    else:
      signal.push_back(noisepulse[i])
    
  return signal
  

def shiftSignalLeft(signal, shift):
  newsignal = copy.copy(signal)
  for i in range(len(newsignal)):
    if i + shift < len(signal):
      newsignal[i] = signal[i+shift]
    else:
      newsignal[i] = 0
  return newsignal
  
def shiftSignalRight(signal, shift):
  newsignal = copy.copy(signal)
  for i in range(len(newsignal)):
    if i > shift:
      newsignal[i] = signal[i-shift]
    else:
      newsignal[i] = 0
  
  return newsignal
  
def addWindowFunction(signal, window, center):
  #this assumes that outside of the range of the window function, the window = 0
  #the middle position of the window function is aligned with the value 'center'
  firstPos = center - len(window)/2
  j = 0
  for i in range(len(signal)):
    if i < firstPos:
      signal[i] = 0
    elif i >= firstPos and i < firstPos + len(window):
      signal[i] = signal[i] * window[j]
      j+=1
    elif i >= firstPos + len(window):
      signal[i] = 0
      

def main(*args):

  tr = {}  #the object that is returned by this function - used to analyze the results
  
  tr['nyquist'] = 50000  #50 kHz is the nyquist frequncy since we have a sample rate of 100 kHz
  
  ### create the pulse template and save it to the return
  (template, python_template) = getTemplate()
  tr['template'] = python_template
  tr['templatepower'] = calculatePower(template) ### calculate the power 
   
  pulseLength = len(python_template) #get the length from the template
  
  # get a noise pulse 
  # save it to the return
  random.seed()
  noisepulse = getNoisePulse(pulseLength)
  wn = []
  for i in range(pulseLength):
    wn.append(noisepulse[i])
    
  tr['noise'] = wn
  tr['noisepower'] = calculatePower(noisepulse) ### calculate the power 
  
  
    
  #create the signal by adding the template to the noise
  #save it to the return
  signal = createSignal(pulseLength, python_template, noisepulse, 500.)
  signalpy = []
  for i in range(pulseLength):
    signalpy.append(signal[i])

  tr['signal'] = signalpy
  tr['signalpower'] = calculatePower(signal) ### calculate the power 
  
  
  #add a windowing function
  x = sig.blackman(3200)
  xx = KWindowDesign.GetBlackmanWindow(3200)
  x = np.zeros(3200)
  for i in range(3200):
    x[i] = xx[i]
  tr['window'] = x.tolist()
  addWindowFunction(signal, x, 6556)
  addWindowFunction(template, x, np.argmin(template))
  
  tr['window_signal'] = []
  tr['window_template'] = []
  for i in range(len(signal)):
    tr['window_signal'].append(signal[i])
  
  for i in range(len(template)):
    tr['window_template'].append(template[i])
  #perform the filter
  
  (b,a) = sig.iirfilter(2,[0.001, 0.01])
  filter = KIIRFourthOrder(a[1], a[2], a[3], a[4], b[0], b[1], b[2], b[3], b[4])
  
  #here's the filter's frequency response function
  tr['b'] = b
  tr['a'] = a
  
  
  
  #for fun, pass the noise through the bandpass
  filter.SetInputPulse(noisepulse)
  filter.RunProcess()
  bp_noise = std.vector("double")()
  bp_noisepy = []
  for i in range(pulseLength):
    bp_noise.push_back(filter.GetOutputPulse()[i])
    bp_noisepy.append(filter.GetOutputPulse()[i])
  
  tr['bp_noise'] = bp_noisepy
  tr['bp_noisepower'] = calculatePower(bp_noise)
  
  #pass the template through the bandpass
  filter.SetInputPulse(template)
  filter.RunProcess()
  bp_template = std.vector("double")()
  bp_templatepy = []
  for i in range(pulseLength):
    bp_template.push_back(filter.GetOutputPulse()[i])
    bp_templatepy.append(filter.GetOutputPulse()[i])
  
  tr['bp_template'] = bp_templatepy
  tr['bp_templatepower'] = calculatePower(bp_template)
  
  #pass the signal through the bandpass
  filter.SetInputPulse(signal)
  filter.RunProcess()
  bp_signal = std.vector("double")()
  bp_signalpy = []
  for i in range(pulseLength):
    bp_signal.push_back(filter.GetOutputPulse()[i])
    bp_signalpy.append(filter.GetOutputPulse()[i])
    
  tr['bp_signal'] = bp_signalpy  
  tr['bp_signalpower'] = calculatePower(bp_signal)


  correlation = KCorrelation()

  #set the template as the response to the correlation function
  bp_tempLeft = shiftSignalRight(bp_template, 2090)
  bp_tempMiddle = shiftSignalRight(bp_template, 2095)
  bp_tempRight = shiftSignalRight(bp_template, 2100)
  
  bp_py = []
  for i in range(len(bp_tempLeft)):
    bp_py.append(bp_tempLeft[i])
  tr['bp_tempLeft'] = bp_py
  
  bp_py = []
  for i in range(len(bp_tempMiddle)):
    bp_py.append(bp_tempMiddle[i])
  tr['bp_tempMiddle'] = bp_py
  
  bp_py = []
  for i in range(len(bp_tempRight)):
    bp_py.append(bp_tempRight[i])
  tr['bp_tempRight'] = bp_py
  
  
  correlation.SetResponse(bp_tempLeft)
  correlation.SetInputPulse(bp_signal)
  correlation.RunProcess()
  corrOut = std.vector("double")()
  corrOutpy = []
  for i in range(correlation.GetOutputPulseSize()):
    corrOut.push_back(correlation.GetOutputPulse()[i])
    corrOutpy.append(correlation.GetOutputPulse()[i])
    
  tr['corrLeft'] = corrOutpy
  
  
  correlation.SetResponse(bp_tempMiddle)
  correlation.SetInputPulse(bp_signal)
  correlation.RunProcess()
  corrOut = std.vector("double")()
  corrOutpy = []
  for i in range(correlation.GetOutputPulseSize()):
    corrOut.push_back(correlation.GetOutputPulse()[i])
    corrOutpy.append(correlation.GetOutputPulse()[i])
    
  tr['corrMiddle'] = corrOutpy
  
  
  correlation.SetResponse(bp_tempRight)
  correlation.SetInputPulse(bp_signal)
  correlation.RunProcess()
  corrOut = std.vector("double")()
  corrOutpy = []
  for i in range(correlation.GetOutputPulseSize()):
    corrOut.push_back(correlation.GetOutputPulse()[i])
    corrOutpy.append(correlation.GetOutputPulse()[i])
    
  tr['corrRight'] = corrOutpy
  
  
    
  return tr

if __name__ == '__main__':
  main(*sys.argv[1:])