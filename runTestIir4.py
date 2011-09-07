#!/usr/bin/env python

import testiir4 as t
r = t.main()

import matplotlib.pyplot as plt
import numpy as np
wn = np.array(r['whitenoise'])
wnp = np.array(r['whitenoisepower'])
r['a']
r['b']
r['filter order']
fp  = np.array(r['filteredpulse'])
fpp  = np.array(r['filteredpower'])


f, axarr = plt.subplots(2, 2)
axarr[0,0].plot(wn)
axarr[0,1].loglog(wnp)
axarr[1,0].plot(fp)
axarr[1,1].loglog(fpp)