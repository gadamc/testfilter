#!/usr/bin/env python

import testBandPassFilterAndCorrelation as t
import scipy.signal as sig
r = t.main()

import matplotlib.pyplot as plt
import numpy as np
noise = np.array(r['noise'])
noisepower = np.array(r['noisepower'])
template = np.array(r['template'])
templatepower = np.array(r['templatepower'])
signal = np.array(r['signal'])
signalpower = np.array(r['signalpower'])
bp_noise = np.array(r['bp_noise'])
bp_noisepower = np.array(r['bp_noisepower'])
bp_template = np.array(r['bp_template'])
bp_templatepower = np.array(r['bp_templatepower'])
bp_signal = np.array(r['bp_signal'])
bp_signalpower = np.array(r['bp_signalpower'])

h, w = sig.freqz(r['b'],r['a'])

corr = np.array(r['corr'])

f, axarr = plt.subplots(5, 3)
axarr[0,0].plot(noise)
axarr[0,1].plot(template)
axarr[0,2].plot(signal)
axarr[1,0].semilogy(noisepower)
axarr[1,1].semilogy(templatepower)
axarr[1,2].semilogy(signalpower)
axarr[2,0].plot(bp_noise)
axarr[2,1].plot(bp_template)
axarr[2,2].plot(bp_signal)
axarr[3,0].semilogy(bp_noisepower)
axarr[3,1].semilogy(bp_templatepower)
axarr[3,2].semilogy(bp_signalpower)
axarr[4,0].title('Digital filter frequency response')
axarr[4,0].title('Digital filter frequency response')
axarr[4,0].semilogy(h, np.abs(w), 'b')
axarr[4,0].ylabel('Amplitude (dB)', color='b')
#axarr[4,0].xlabel('Frequency (rad/sample)')
axarr[4,0].grid()
axarr[4,0].legend()

ax2 = ax1.twinx()
angles = np.unwrap(np.angle(w))
axarr[4,0].ylabel('Angle (radians)', color='g')
axarr[4,0].plot(h, angles, 'g')
axarr[4,2].plot(corr)
