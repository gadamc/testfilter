#!/usr/bin/env python

import testBandPassFilterAndCorrelation as t
r = t.main()

import matplotlib.pyplot as plt
import numpy as np
noise = np.array(r['noise'])
noisepower = np.array(r['noisepower'])
template = np.array(r['template'])
templatepower = np.array(r['templatepower'])
signal = np.array(r['signal'])
signalpower = np.array(r['signalpower'])
bp_template = np.array(r['bp_template'])
bp_templatepower = np.array(r['bp_templatepower'])
bp_signal = np.array(r['bp_signal'])
bp_signalpower = np.array(r['bp_signalpower'])

corr = np.array(r['corr'])

f, axarr = plt.subplots(5, 3)
axarr[0,0].plot(noise)
axarr[0,1].plot(template)
axarr[0,2].plot(signal)
axarr[1,0].loglog(noisepower)
axarr[1,1].loglog(templatepower)
axarr[1,2].loglog(signalpower)
axarr[2,1].plot(bp_template)
axarr[2,2].plot(bp_signal)
axarr[3,1].loglog(bp_templatepower)
axarr[3,2].loglog(bp_signalpower)
axarr[4,0].plot(corr)
