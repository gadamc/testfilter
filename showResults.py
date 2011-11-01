#!/usr/bin/env python

from ROOT import *
import sys

gSystem.Load("libkamping")
f = KDataReader('jf07b002/jf07b002_000.local.simpleamp.root')
e = f.GetEvent()

smin = smax = 0
bmin = bmax = 0


for i in range(f.GetEntries()):
  f.GetEntry(i)
  for j in range(e.GetNumBolos()):
    b = e.GetBolo(j)
   
    if b.GetDetectorName() == 'ID3':
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        
        if p.GetChannelName() == 'centre ID3CD':
          for n in range(p.GetNumPulseAnalysisRecords()):
            r = p.GetPulseAnalysisRecord(n)
          
            if r.GetName() == "KTrapKamperProto":
              if r.IsBaseline()==0:
                if r.GetAmp() < smin: smin = r.GetAmp()
                if r.GetAmp() > smax: smax = r.GetAmp()
              else:
                if r.GetAmp() < bmin: bmin = r.GetAmp()
                if r.GetAmp() > bmax: bmax = r.GetAmp()

                
h2 = TH1D('h2','h2',1000,bmin, bmax)
h = TH1D('h','h',1000,smin, smax)

for i in range(f.GetEntries()):
  f.GetEntry(i)
  for j in range(e.GetNumBolos()):
    b = e.GetBolo(j)
   
    if b.GetDetectorName() == 'ID3':
      for k in range(b.GetNumPulseRecords()):
        p = b.GetPulseRecord(k)
        
        if p.GetChannelName() == 'centre ID3CD':
          for n in range(p.GetNumPulseAnalysisRecords()):
            r = p.GetPulseAnalysisRecord(n)
          
            if r.GetName() == "KTrapKamperProto":
              if r.IsBaseline()==0:
                print 'baseline amp', r.GetAmp()
                print 'baseline pos', r.GetPeakPosition()
                h.Fill(r.GetAmp())
              else:
                print 'pulse amp', r.GetAmp()
                print 'pulse pos', r.GetPeakPosition()
                h2.Fill(r.GetAmp())
                
              