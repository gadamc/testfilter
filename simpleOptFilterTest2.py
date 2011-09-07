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
    
def getNoisePulse(pl):
  whitenoise = std.vector("double")()

  for i in range(pl):
    whitenoise.push_back(random.gauss(0,1))
  
  #now filter the white noise with an iir filter to make it more realistic
  order = 4
  (b,a) = sig.iirfilter(order,0.1, btype='lowpass')
  iir4 = KIIRFourthOrder(-a[1], -a[2], -a[3], -a[4], b[0], b[1], b[2], b[3], b[4])
  iir4.SetInputPulse(whitenoise)
  print 'filter white noise', iir4.RunProcess()

  noise = std.vector("double")()
  for i in range(iir4.GetOutputPulseSize()):
    noise.push_back(iir4.GetOutputPulse()[i])

  return noise
  
def getNoiseSpectrum(pl):
  npulse = getNoisePulse(pl)
  r2hc = KRealToHalfComplexDFT()
  r2hc.SetInputPulse(npulse)
  r2hc.RunProcess()
  hcp = KHalfComplexPower()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  hcp.RunProcess()
  spectrum = std.vector("double")()
  for i in range(hcp.GetOutputPulseSize()):
    spectrum.push_back(hcp.GetOutputPulse()[i])
    
  return spectrum
  
def getTemplate(pl):
  t = std.vector("double")()
  pt = []
  for i in range(pl):
    t.push_back(getgauss(1./(50*sqrt(2*math.pi)), pl/2., 50., float(i))) #the template needs to have a maximum amplitude of 1
    #t.push_back(getgauss(10., pl/2., 50., float(i))) #the template needs to have a maximum amplitude of 1
    pt.append(getgauss(1./(50*sqrt(2*math.pi)), pl/2., 50., float(i)))
    #pt.append(getgauss(10., pl/2., 50., float(i)))
    
  return (t, pt)

def createSignal(pl, template, noisepulse, amplitude = 1.):
  signal = std.vector("double")()
  for i in range(pl):
    signal.push_back( 50*sqrt(2*math.pi)*amplitude*template[i-2000] + noisepulse[i])
    #signal.push_back( template[i-2000] + noisepulse[i])
    
  return signal
  
def getAverageNoisePower(pl):
  
  avePower = getNoiseSpectrum(pl)
  
  num = 50
  
  for i in range(1,num):
    aNoise = getNoiseSpectrum(pl)
    for ii in range(avePower.size()):
      avePower[ii] += aNoise[ii]
    
  for j in range(avePower.size()):
    avePower[j] = avePower[j]/float(num)

  
  return avePower

def main(*args):

  tr = {}
  
  pulseLength = 8196
  
  ### create the pulse template and save it to the return
  (template, python_template) = getTemplate(pulseLength)
  tr['template'] = python_template
  
  # get a noise pulse 
  # save it to the return
  noisepulse = getNoisePulse(pulseLength)
  wn = []
  for i in range(pulseLength):
    wn.append(noisepulse[i])
    
  tr['noise'] = wn
  
  #create the signal by adding the template to the noise
  #save it to the return
  signal = createSignal(pulseLength, python_template, noisepulse, 1.)
  sig = []
  for i in range(pulseLength):
    sig.append(signal[i])

  tr['signal'] = sig
  
  
  
  ### calculate the power of the noise pulse for this instance and save it for the return
  r2hc = KRealToHalfComplexDFT()
  r2hc.SetInputPulse(noisepulse)
  print 'real to half complex noisepulse', r2hc.RunProcess()
  hcp = KHalfComplexPower()
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'power noise', hcp.RunProcess()
  #save it for the return
  wp = []
  for i in range(hcp.GetOutputPulseSize()):
    wp.append(hcp.GetOutputPulse()[i])
    
  tr['noisepower'] = wp
  #######
  

  filter = KOptimalFilter()
  
  #### get the average noise power in order to pass it to the filter
  avePower = getAverageNoisePower(pulseLength)
  filter.SetNoiseSpectrum(avePower)
  
  #save it for the return
  ap= []
  for i in range(len(wp)):
    ap.append(avePower[i])
  
  tr['avenoisepower'] = ap
  
  ####  calculate the fourier transform and the power of the pulse template.
  #pass the DFT of the transform to the filter and then save the template power in a
  #python arry for the return
  r2hc.SetInputPulse(template)
  print 'real to half complex template', r2hc.RunProcess()
  filter.SetTemplateDFT(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  
  hcp.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  print 'power template', hcp.RunProcess()
  tmppower = []
  for i in range(hcp.GetOutputPulseSize()):
    tmppower.append(hcp.GetOutputPulse()[i])
    
  tr['templatepower'] = tmppower
  
  
  # now calculate the fourier transform of the event signal and pass it to the filter
  r2hc.SetInputPulse(signal)
  print 'real to half complex signal', r2hc.RunProcess()
  filter.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
  
  #build the filter and run
  print 'building filter', filter.BuildFilter()
  print 'applying filter', filter.RunProcess()
  
  #save the output pulse amplitude estimates to the return
  output = []
  for i in range(filter.GetOutputPulseSize()):
    output.append(filter.GetOutputPulse()[i])
  
  tr['amplitude'] = output
  
  #save the optimal filter for the return
  optfil = []
  for i in range(filter.GetOptimalFilterSize()):
    optfil.append(filter.GetOptimalFilter()[i])
  
  tr['optfilter'] = optfil
  
  #save the power of the optimal filter for the return
  hcp.SetInputPulse(filter.GetOptimalFilter(), filter.GetOptimalFilterSize())
  print 'get optimal filter power', hcp.RunProcess()
  optfilpower = []
  for i in range(hcp.GetOutputPulseSize()):
    optfilpower.append(hcp.GetOutputPulse()[i])
    
  tr['optfilpower'] = optfilpower
  

  return tr

if __name__ == '__main__':
  main(*sys.argv[1:])