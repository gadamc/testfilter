from opFilTest import *
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt

def findPeaks(pl,  pulse , r,  t):
  peaks = []
  amp = 1*np.std(np.array(pulse))
  rmax = 5
  for i in range(0, pl-2*r-t, 1):
    if ((pulse[i]<-amp) and (pulse[i+r]>amp) and (pulse[i+r+t]>amp) and (pulse[i+2*r+t]<-amp)):
#      r1 = abs(pulse[i]/pulse[i+2*r+f])
#      r2 = abs(pulse[i+r]/pulse[i+r+f])
#      r3 = abs(pulse[i]/pulse[i+r])
#      r4 = abs(pulse[i+r+f]/pulse[i+2*r+f])
#      if r1<1:
#        r1 = 1/r1
#      if r2<1:
#        r2 = 1/r2
#      if r3<1:
#        r3 = 1/r3
#      if r4<1:
#        r4 = 1/r4
      peaks.append(pulse[i+r])
     # peaks.append(1.0)
#      if r1<rmax and r2<rmax and r3<rmax and r4<rmax:
#        peaks.append(1.)
#      else:
#        peaks.append(0.)
    else:
      peaks.append(0.)
  return peaks
def main(entry):
  f = KDataReader('/kalinka/home/edelweiss/uspace/gadamc/le31b007_007.root')
  f.GetEntry(entry)
  e = f.GetEvent()
  b = e.GetBolo(2)
  p = b.GetPulseRecord(1)
  raw = p.GetTrace()
  pattrem = KPatternRemoval()
  basrem = KBaselineRemoval()
  pattrem.SetInputPulse(p.GetTrace())
  pattrem.RunProcess()
  basrem.SetInputPulse(pattrem.GetOutputPulse(), pattrem.GetOutputPulseSize())
  basrem.RunProcess()
  noise = std.vector("double")()
  for i in range(basrem.GetOutputPulseSize()):
    noise.push_back(basrem.GetOutputPulse()[i])
  rms = np.std(np.array(noise))
  
  risetime = 100
  flat = 200
  template, pytemplate = getTemplate()
  tempcoeff = abs(1/float(pytemplate[4464]))
  zeronoise = np.zeros(len(pytemplate))
  signal,  pysignal = createSignal(len(pytemplate), template, noise, 5.*tempcoeff*rms, 0)
  #signal = lowpass(signal1)
  trapfilter = KTrapezoidalFilter()
  trapfilter.SetInputPulse(signal)
  list = []
  for i in range(len(pytemplate)):
    list.append(0.0)
  
  for r in range(3, 5, 2):
    for f in [20]:
      trapfilter.SetParams(1100,r,  f)
      trapfilter.RunProcess()
      out = std.vector("double")()
      for i in range(trapfilter.GetOutputPulseSize()):
        out.push_back(trapfilter.GetOutputPulse()[i])
      diff = std.vector("double")()
      diff.push_back(0.0)
      for i in range(1, trapfilter.GetOutputPulseSize(), 1):
        diff.push_back(trapfilter.GetOutputPulse()[i]-trapfilter.GetOutputPulse()[i-1])
      diff2 = std.vector("double")()
      diff2.push_back(0.0)
      for i in range(1, diff.size(), 1):
        diff2.push_back(diff[i]-diff[i-1])
      peaks = findPeaks(diff2.size(), diff2, r, f)
      for i in range(len(peaks)):
        list[i]=list[i]+out[i+r+f]*peaks[i]
  return np.argmax(np.absolute(list))
#  trapfilter.SetParams(1100,100,  100)
#  trapfilter.RunProcess()
#  out = std.vector("double")()
#  for i in range(trapfilter.GetOutputPulseSize()):
#    out.push_back(trapfilter.GetOutputPulse()[i])
#    
#  flat = np.mean([out[x] for x in range(4564, 4664)])
#  amp = np.zeros(len(pytemplate))
#  for i in range(4564, 4664):
#    amp[i] = flat
#  list1 = sorted(np.absolute(list), reverse=True)
#  print list1[:5]
#  print "SNR: ", list1[0]/float(list1[1])
#  
  
#  plt.figure(1)
#  plt.subplot(311)
#  plt.plot(np.array(signal[4400:4800]))
#  plt.subplot(312)
#  plt.plot(np.array(out[4400:4800]))
#  plt.plot(amp[4400:4800])
#  plt.subplot(313)
#  #plt.plot(np.array(noise))
#  #plt.subplot(414)
#  plt.plot(np.array(list[4400:4800]))
  
  
#  plt.show()

if __name__ == '__main__':
  pos = main()
