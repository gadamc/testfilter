#!/usr/bin/env python
import math
import simpleOptFilterTest3 as t
r = t.main()

import matplotlib.pyplot as plt
import numpy as np
noise = np.array(r['noise'])
template = np.array(r['template'])
signal = np.array(r['signal'])
noisepower = np.array(r['noisepower'])
temppower = np.array(r['templatepower'])
optfilpower = np.array(r['optfilpower'])
optfilter = np.array(r['optfilter'])
avenoisepower = np.array(r['avenoisepower'])


amp = np.array(r['amplitude'])
amp2 = np.array(r['amplitude2'])

f, axarr = plt.subplots(3, 3)
axarr[0,0].plot(noise)
axarr[0,1].plot(template)
axarr[0,2].plot(signal)
hzaxis = np.array(range(0,len(temppower)))
hzperpoint = r['nyquist']/len(temppower)
hzaxis *= hzperpoint
axarr[1,0].loglog(hzaxis, noisepower)
axarr[1,1].loglog(hzaxis, temppower)
axarr[1,2].loglog(hzaxis, optfilpower)
axarr[2,0].loglog(hzaxis, avenoisepower)
axarr[2,1].plot(abs(amp))
axarr[2,2].plot(abs(amp2))
