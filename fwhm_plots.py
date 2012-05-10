#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

freq = [0.5,0.75,1.,1.5,2,3,5,7.5,10,12.5,15]
fwhm = {}
fwhm['chalA FID807'] = [3.64,3.05,3.08,3.12,3.01,3.02,2.9,2.98,2.97,3.01,3.12]
fwhm['chalB FID807'] = [4.84,3.18,3,2.59,2.28,2.27,1.66,1.72,1.79,1.88,1.98]
fwhm['chalA FID808'] = [1.7,1.39,1.28,1.19,1.12,1.11,1.03,1.08,1.13,1.19,1.24]
fwhm['chalB FID808'] = [2.55,1.99,1.86,1.63,1.51,1.51,2.01,2.91,4.03,5.35,6.84]

plt.ion()

# chalA FID807
plt.figure(1)
plt.clf()
plt.loglog(freq,fwhm['chalA FID807'],marker='o',ms=15)
plt.axes().set_xlim([0.25,80])
plt.axes().set_ylim([0.5,7.5])
plt.axes().set_title('chalA FID807, lowpass 2nd order at 50Hz',fontsize=16)
plt.axes().set_xlabel('2nd order highpass cut-off frequency [Hz]',fontsize=16)
plt.axes().set_ylabel('FWHM [keV]',fontsize=16)
for tick in plt.axes().xaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().xaxis.get_ticklines(minor=False):
  line.set_markersize(30)

for line in plt.axes().xaxis.get_ticklines(minor=True):
  line.set_markersize(15)

for tick in plt.axes().yaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().yaxis.get_ticklines(minor=False):
  line.set_markersize(20)

plt.text(0.75,6.0,'5.0 Hz, 2.90 keV',horizontalalignment = 'center',verticalalignment='top',fontsize=15)

raw_input()
plt.savefig('chalAfid807.eps',transparent=True)
#chalB FID807
plt.figure(2)
plt.clf()
plt.loglog(freq,fwhm['chalB FID807'],marker='o',ms=15)
plt.axes().set_xlim([0.25,80])
plt.axes().set_ylim([0.5,7.5])
plt.axes().set_title('chalB FID807, lowpass 2nd order at 50Hz',fontsize=16)
plt.axes().set_xlabel('2nd order highpass cut-off frequency [Hz]',fontsize=16)
plt.axes().set_ylabel('FWHM [keV]',fontsize=16)
for tick in plt.axes().xaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().xaxis.get_ticklines(minor=False):
  line.set_markersize(30)

for line in plt.axes().xaxis.get_ticklines(minor=True):
  line.set_markersize(15)

for tick in plt.axes().yaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().yaxis.get_ticklines(minor=False):
  line.set_markersize(20)

plt.text(0.75,6.0,'5.0 Hz, 1.66 keV',horizontalalignment = 'center',verticalalignment='top',fontsize=15)


raw_input()
plt.savefig('chalBfid807.eps',transparent=True)


#chalA FID808
plt.figure(3)
plt.clf()
plt.loglog(freq,fwhm['chalA FID808'],marker='o',ms=15)
plt.axes().set_xlim([0.25,80])
plt.axes().set_ylim([0.5,7.5])
plt.axes().set_title('chalA FID808, lowpass 2nd order at 50Hz',fontsize=16)
plt.axes().set_xlabel('2nd order highpass cut-off frequency [Hz]',fontsize=16)
plt.axes().set_ylabel('FWHM [keV]',fontsize=16)
for tick in plt.axes().xaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().xaxis.get_ticklines(minor=False):
  line.set_markersize(30)

for line in plt.axes().xaxis.get_ticklines(minor=True):
  line.set_markersize(15)

for tick in plt.axes().yaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().yaxis.get_ticklines(minor=False):
  line.set_markersize(20)

plt.text(0.75,6.0,'5.0 Hz, 1.03 keV',horizontalalignment = 'center',verticalalignment='top',fontsize=15)

raw_input()
plt.savefig('chalAfid808.eps', transparent=True)


#chalB FID808
plt.figure(4)
plt.clf()
plt.loglog(freq,fwhm['chalB FID808'],marker='o',ms=15)
plt.axes().set_xlim([0.25,80])
plt.axes().set_ylim([0.5,7.5])
plt.axes().set_title('chalB FID808, lowpass 2nd order at 50Hz',fontsize=16)
plt.axes().set_xlabel('2nd order highpass cut-off frequency [Hz]',fontsize=16)
plt.axes().set_ylabel('FWHM [keV]',fontsize=16)
for tick in plt.axes().xaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().xaxis.get_ticklines(minor=False):
  line.set_markersize(30)

for line in plt.axes().xaxis.get_ticklines(minor=True):
  line.set_markersize(15)

for tick in plt.axes().yaxis.get_major_ticks():
  tick.label.set_fontsize(18)

for line in plt.axes().yaxis.get_ticklines(minor=False):
  line.set_markersize(20)

plt.text(0.75,6.0,'2.0 Hz, 1.51 keV',horizontalalignment = 'center',verticalalignment='top',fontsize=15)

raw_input()
plt.savefig('chalBfid808.eps',transparent=True)
