from opFilTest import *
from ROOT import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pickle

def main():
  # Parameter
  amplitude = 15.
  decay = 1100.
  rise = 100
  top = 100
  #--------------
  file = open('/kalinka/home/unrau/noiseevents-le31b007_007-all.pkl', 'rb')
  noiseentries = pickle.load(file)
  rms = pickle.load(file) 
  file.close()
  noiselist = []
  rmslist = []
  # RMS cut
  for i in range(len(noiseentries)):
    if rms[i] >0 and rms[i]<10:
      noiselist.append(noiseentries[i])
      rmslist.append(rms[i])
  print len(noiselist)
  AmpEstimate = []
  template,  pytemplate = getTemplate()
  tempcoeff = 1/float(pytemplate[4464])
  zeronoise = np.zeros(len(template))
  f = KDataReader('/kalinka/home/edelweiss/uspace/gadamc/le31b007_007.root')
  patternrv = KPatternRemoval()
  baselinerv = KBaselineRemoval()
  trapfilter = KTrapezoidalFilter()
  trapfilter.SetParams(decay, rise, top)
  event = f.GetEvent()
  bolo = event.GetBolo(2)
  for entry in noiselist:
    f.GetEntry(entry)
    if noiselist.index(entry)%10==0:
      print "Entry: ", entry
    patternrv.SetInputPulse(bolo.GetPulseRecord(1).GetTrace())
    patternrv.RunProcess()
    baselinerv.SetInputPulse(patternrv.GetOutputPulse(), patternrv.GetOutputPulseSize())
    baselinerv.RunProcess()
    noise = std.vector("double")()
    for i in range(baselinerv.GetOutputPulseSize()):
      noise.push_back(baselinerv.GetOutputPulse()[i])
  
    # Amplitude estimation for the template
    signal,  pysignal = createSignal(len(pytemplate), template, zeronoise, amplitude*tempcoeff*rmslist[noiselist.index(entry)], 0)
    trapfilter.SetInputPulse(signal)
    trapfilter.RunProcess()
    out = std.vector("double")()
    for i in range(trapfilter.GetOutputPulseSize()):
      out.push_back(trapfilter.GetOutputPulse()[i])
    TemplateAmp = np.mean([out[x] for x in range(4464+rise, 4464+rise+top)])
    
    # Amplitude estimation for template+noise
    signal,  pysignal = createSignal(len(pytemplate), template, noise, amplitude*tempcoeff*rmslist[noiselist.index(entry)], 0)
    trapfilter.SetInputPulse(signal)
    trapfilter.RunProcess()
    out = std.vector("double")()
    for i in range(trapfilter.GetOutputPulseSize()):
      out.push_back(trapfilter.GetOutputPulse()[i])
    Amp = np.mean([out[x] for x in range(4464+rise, 4464+rise+top)])
    AmpEstimate.append(Amp/float(TemplateAmp))
  
  #print AmpEstimate
  return noiselist,  AmpEstimate

if __name__ == '__main__':
  noiselist,  AmpEstimate = main()
  mean = np.mean(AmpEstimate)
  std = np.std(AmpEstimate)
  print mean,  std
  plt.figure(1)
  n,  bins,  patches = plt.hist(AmpEstimate, 50, normed = 1, facecolor='green')
  y = mlab.normpdf(bins, mean, std)
  plt.plot(bins, y, 'r--', linewidth=2)
  plt.title('accuracy of amplitude estimation')
  plt.xlabel('estimated amplitude/true amplitude')
  plt.ylabel('relative frequency')
  plt.show()
