#!/usr/bin/env python
from couchdbkit import Server, Database
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
plt.ion()

#small utility functions
def get_as_nparray(c_pointer, size):
  data = np.zeros(size)
  for i in range(size):
    data[i] = c_pointer[i]
  return data
          
def get_out(pta):
  return get_as_nparray(pta.GetOutputPulse(), pta.GetOutputPulseSize())

def get_in(pta):
  return get_as_nparray(pta.GetInputPulse(), pta.GetInputPulseSize())


gSystem.Load('libkds')  
gSystem.Load('libkpta')

myChannel = 'chalA FID807'

s = Server('http://edwdbik.fzk.de:5984')
db = s['datadb']  
viewResults = db.view('proc/raw')  #proc/raw returns all data files that have a raw .root file

dbtemplates = s['pulsetemplates']   #holds the pulse templates
vr = dbtemplates.view('analytical/bychandate', startkey=[myChannel, '2012-01-22 20:38:00'], limit=1, reduce=False, descending=True, include_docs=True)  #this is the view call to get the template
template_doc = vr.first()['doc']

#have to 'massage' the template a bit here. the template parameters are stored in units of ms. 
#so make sure to convert properly between time and digitized sample number.

print 'grabbing template from database'
vp = std.vector("double")()
exec( template_doc['formula']['rise_double_decay']['python'])  #defines 'template' function used in the following for-loop
for i in range( 512 ):
  vp.push_back( template(i * 2.016, template_doc['formula']['rise_double_decay']['par']))  #2.016 ms / sample
          
scaleFactor = 1./abs(min(np.array(vp)))
for i in range(vp.size()):
  vp[i] = scaleFactor * vp[i]
   
print 'creating pulse processors for analysis'
#create the pulse processors that are needed for this analysis
bas = KBaselineRemoval() 
r2hc = KRealToHalfComplexDFT()
era = KEraPeakFinder()
hc2p = KHalfComplexPower()
optFilter = KOptimalFilter()

#configure the ERA peak finder
era.SetOrder(5)
era.SetNumRms(2.3)
era.SetPolarity(-1)  #tell the ERA peak finder to only look for negative polarity 

#put the pulse template fourier transform into the optimal filter.
r2hc.SetInputPulse(vp)
r2hc.RunProcess()
optFilter.SetTemplateDFT( r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize() )

#calculate the template Power here just to plot it later...
hc2p.SetInputPulse( r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize() )
hc2p.RunProcess()
templatePower = get_out(hc2p)

#before starting, create a dictionary to store results.
results = dict()  
results['amp'] = []
results['time'] = []

print 'starting loop over data files'
for row in viewResults: 
  
  if row['key'] != 'ma22a000': continue #skip all runs but ma22a000
  doc = db[row['id']] #grab the document from the database.
  
  print 'searching for noise events in', doc['file']
  #open the file
  f = KDataReader( doc['file'] )
  e = f.GetEvent()  #get the KRawEvent object

  #now, we have to loop through the data and 
  #determine the average noise power spectrum for `myChannel'

  power = std.vector("double")()  #will store the average power spectrum here. 
  numOfSpectra = 0  # the number of power spectra used in calculation of the average
  
  for i in range(f.GetEntries()):
    f.GetEntry(i)
  
    for ii in range(e.GetNumBolosPulses()):
      p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
    
      if p.GetPulseLength() == 0: continue  #skip empty pulses
      
      if p.GetChannelName() == myChannel:        
            
        bas.SetInputPulse( p.GetTrace()  )
        bas.RunProcess()    
        era.SetInputPulse( bas.GetOutputPulse(), bas.GetOuputPulseSize() )
        era.RunProcess()
        if era.GetPeakBins().size() > 0:  continue  #skip pulse if we didn't find a noise pulse
        
        r2hc.SetInputPulse( bas.GetOutputPulse(), bas.GetOutputPulseSize() )
        r2hc.RunProcess()
        hc2p.SetInputPulse ( r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize() )
        hc2p.RunProcess()
      
        #i assume here that the pulse length never changes within a run
        if numOfSpectra == 0:
          power.resize(hc2p.GetOutputPulseSize())
          for  k in range(power.size()):
            power[k] = hc2p.GetOutputPulse()[k]
        else:
          for k in range(power.size()):
            power[k] = power[k] * (numOfSpectra-1.0)/numOfSpectra + hc2p.GetOutputPulse()[k] / numOfSpectra
      
        numOfSpectra += 1.0  
  
  #done looping through the data, now set the power spectrum in the optimal filter      
  optFilter.SetNoiseSpectrum(power)
  
  #take a quick look at the noise power spectra and the pulse template power
  plt.loglog(templatePower)
  plt.loglog( np.array(power) )
  raw_input('hit Enter to continue...')
  
  #now loop through the data again, applying the optimal filter that was just built.
  for i in range(f.GetEntries()):
    f.GetEntry(i)
  
    for ii in range(e.GetNumBolosPulses()):
      p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
    
      if p.GetPulseLength() == 0: continue  #skip empty pulses
    
      if p.GetChannelName() == myChannel:    
        
        bas.SetInputPulse( p.GetTrace() )
        bas.RunProcess()
        optFilter.SetInputPulse( bas.GetOutputPulse(), bas.GetOutputPulseSize() )
        optFilter.RunProcess()
      
        amp = get_out( optFilter )
        subamp = amp[:20]  #only look at the first 20 bins of the output amplitude estimator
      
        results['amp'].append( subamp[subamp.argmax()] )
        results['peaktime'].append( subamp.argmax() )
      
    
hall = TH1D('hall', 'hall', 15000, -1000, 14000)      
hpeak =   TH1D('hpeak', 'hpeak', 15000, -1000, 14000)
  
for i in range( len(results['amp']) ):
  hall.Fill( results['amp'][i] )
  if results['peaktime'] > 4 and results['peaktime'] < 8:
    hpeak.Fill( results['amp'][i] )
  
hpeak.SetLineColor(kRed)
hall.Draw()
hpeak.Draw('same')

saveme = raw_input('type "yes" if you want to save this')
fresults = TFile('output.root', 'recreate')
hall.Write()
hpeak.Write()
fresults.Close()
