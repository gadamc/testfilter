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
windowsignal = np.array(r['window_signal'])
windowtemplate = np.array(r['window_template'])
bp_noise = np.array(r['bp_noise'])
bp_noisepower = np.array(r['bp_noisepower'])
bp_template = np.array(r['bp_template'])
bp_templatepower = np.array(r['bp_templatepower'])
bp_signal = np.array(r['bp_signal'])
bp_signalpower = np.array(r['bp_signalpower'])

#lets make the horizontal axis of the power spectrums show the correct frequency
hzaxis = np.array(range(0,len(signalpower)))
hzperpoint = r['nyquist']/len(signalpower)
hzaxis *= hzperpoint

h, w = sig.freqz(r['b'],r['a'])

corrL = np.array(r['corrLeft'])
corrM = np.array(r['corrMiddle'])
corrR = np.array(r['corrRight'])

bp_tempLeft = np.array(r['bp_tempLeft'])
bp_tempMiddle = np.array(r['bp_tempMiddle'])
bp_tempRight = np.array(r['bp_tempRight'])


f, axarr = plt.subplots(7, 3)
axarr[0,0].plot(noise)
axarr[0,1].plot(template)
axarr[0,2].plot(signal)
axarr[1,0].loglog(hzaxis, noisepower)
axarr[1,1].loglog(hzaxis, templatepower)
axarr[1,2].loglog(hzaxis, signalpower)
#axarr[2,0].title('Digital filter frequency response')
#axarr[2,0].title('Digital filter frequency response')

#axarr[2,0].ylabel('Amplitude (dB)', color='b')
#axarr[2,0].xlabel('Frequency (rad/sample)')
#axarr[2,0].grid()
#axarr[2,0].legend()
#ax1 = axarr[2,0].add_subplot(111)
#ax2 = ax1.twinx()
# hzperpoint = tr['nyquist']/len(signalpower)
hzperpoint = r['nyquist']/h[len(h)-1]
h *= hzperpoint
axarr[2,0].loglog(h, np.abs(w), 'b')
angles = np.unwrap(np.angle(w))
#axarr[2,1].ylabel('Angle (radians)', color='g')
#axarr[2,1].plot(h, angles, 'g')
axarr[2,1].plot(windowtemplate)
axarr[2,2].plot(windowsignal)
axarr[3,0].plot(bp_noise)
axarr[3,1].plot(bp_template)
axarr[3,2].plot(bp_signal)
axarr[4,0].loglog(hzaxis, bp_noisepower)
axarr[4,1].loglog(hzaxis, bp_templatepower)
axarr[4,2].loglog(hzaxis, bp_signalpower)

axarr[5,0].plot(corrL, 'r')
axarr[5,1].plot(corrM, 'r')
axarr[5,2].plot(corrR, 'r')

axarr[6,0].plot(bp_tempLeft, 'r')
axarr[6,1].plot(bp_tempMiddle, 'r')
axarr[6,2].plot(bp_tempRight, 'r')
