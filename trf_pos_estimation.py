from opFilTest import *
from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import RNT2 as rnt
import pickle

def main():
  # Parameter
  amplitude = 15.
  decay = 1100.
  rise = 3
  top = 20
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
  pos = []
  for entry in noiselist:
    if noiselist.index(entry)%10==0:
      print "Entry:", entry
    pos.append(rnt.main(entry))
  return pos
    


if __name__ == '__main__':
  pos = main()
  plt.figure(1)
  n,  bins,  patches = plt.hist(pos, 20, normed = 1, facecolor='green')
  plt.title('accuracy of peakposition estimation')
  plt.xlabel('estimated peak position')
  plt.ylabel('relative frequency')
  plt.show()
  print pos

    
    
