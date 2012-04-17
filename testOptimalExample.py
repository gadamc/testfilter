#!/usr/bin/env python
from couchdbkit import Server, Database
from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
import string
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

myChannel = 'chalB FID808'

s = Server('http://edwdbik.fzk.de:5984')
db = s['datadb']  
viewResults = db.view('proc/raw')  #proc/raw returns all data files that have a raw .root file

dbtemplates = s['pulsetemplates']   #holds the pulse templates
vr = dbtemplates.view('analytical/bychandate', startkey=[myChannel, '2012-01-22 20:38:00'], limit=1, reduce=False, descending=True, include_docs=True)  #this is the view call to get the template
template_doc = vr.first()['doc']  #the template_doc stores all of the information of the pulse template

   
print 'creating pulse processors for analysis'
#create the pulse processors that are needed for this analysis
bas = KBaselineRemoval() 
r2hc = KRealToHalfComplexDFT()
era = KEraPeakFinder()
hc2p = KHalfComplexPower()
optFilter = KOptimalFilter()
window = KWindow()
windesign = KWindowDesign()
window.SetWindow( windesign.GetTukeyWindow(512, 0.7), 512)
pulseshift = KPulseShifter()


#configure the ERA peak finder
era.SetOrder(5)
era.SetNumRms(2.3)
era.SetPolarity(-1)  #tell the ERA peak finder to only look for negative polarity 

#configure the pulsehsifter. this makes it easier to find the pulse peak position
shiftVal = int(template_doc['formula']['rise_double_decay']['par'][0]/2.016 + 0.5) 
pulseshift.SetShift(-1 * shiftVal )
pulseshift.SetMode(2)

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


#put the pulse template fourier transform into the optimal filter.
tempChain = KPulseAnalysisChain()
tempChain.AddProcessor(window)
tempChain.AddProcessor(pulseshift)
tempChain.AddProcessor(r2hc)

tempChain.SetInputPulse(vp)
tempChain.RunProcess()
plt.plot( np.array(vp) )
plt.plot( get_out(window) )
plt.plot( get_out(pulseshift) )
raw_input('enter to continue')
plt.cla()

optFilter.SetTemplateDFT( tempChain.GetOutputPulse(), tempChain.GetOutputPulseSize() )

hc2p.SetInputPulse(r2hc.GetOutputPulse(), r2hc.GetOutputPulseSize())
hc2p.RunProcess()
plt.plot( get_out(hc2p) )
raw_input('enter to continue')

npChain = KPulseAnalysisChain() #processing chain to calculation the noise power
npChain.AddProcessor(window)
npChain.AddProcessor(r2hc)
npChain.AddProcessor(hc2p)

optChain = KPulseAnalysisChain() #processing chain to apply the optimal filter to a raw pulse
optChain.AddProcessor(bas)
optChain.AddProcessor(window)
optChain.AddProcessor(r2hc)
optChain.AddProcessor(optFilter)

#before starting, create a dictionary to store results.
results = dict()  
results['amp'] = []
results['time'] = []

print 'starting loop over data files'
for row in viewResults: 
  
  if row['key'] != 'ma22a000': continue #skip all runs but ma22a000
  doc = db[row['id']] #grab the document from the database.
  if doc['file_number'] > 4: continue #just do the first 5 hours of data for this example
  
  print 'searching for noise events on', myChannel, 'in', doc['proc1']['file']
  #open the file
  f = KDataReader( doc['proc1']['file'] )
  e = f.GetEvent()  #get the KRawEvent object

  #now, we have to loop through the data and 
  #determine the average noise power spectrum for `myChannel'

  power = std.vector("double")()  #will store the average power spectrum here. 
  numOfSpectra = 0.0  # the number of power spectra used in calculation of the average
  
  for i in range(f.GetEntries()):
    f.GetEntry(i)
  
    for ii in range(e.GetNumBoloPulses()):
      p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
    
      if p.GetPulseLength() == 0: continue  #skip empty pulses
      
      if p.GetChannelName() == myChannel:        
            
        bas.SetInputPulse( p.GetTrace()  )
        bas.RunProcess()    
        era.SetInputPulse( bas.GetOutputPulse(), bas.GetOutputPulseSize() )
        era.RunProcess()
        if era.GetPeakBins().size() > 0:  continue  #skip pulse if we didn't find a noise pulse
        
        #applying windowing and then calculate power spectrum
        npChain.SetInputPulse( bas.GetOutputPulse(), bas.GetOutputPulseSize() )
        #print npChain.GetInputPulseSize(), npChain.GetOutputPulseSize()
        #npChain.SetInputPulse( p.GetTrace() )
        if npChain.RunProcess() == False:
          print 'excusez-moi... ' #this shouldn't fail
          continue
             
        #i assume here that the pulse length never changes within a run
        if numOfSpectra == 0:
          power.resize( npChain.GetOutputPulseSize())
          for  k in range(power.size()):
            power[k] =  npChain.GetOutputPulse()[k]
        else: #calculate the running averge of the noise power spectrum
          for k in range(power.size()):
            power[k] = power[k] * (numOfSpectra-1.0)/numOfSpectra + npChain.GetOutputPulse()[k] / numOfSpectra
      
        numOfSpectra += 1.0  
  
  #done looping through the data, now set the power spectrum in the optimal filter      
  optFilter.SetNoiseSpectrum(power)
  
  #take a quick look at the noise power spectra and the pulse template power
  #print 'plotting average noise power as calculated from %d noise pulses' % numOfSpectra
  plt.loglog( np.array(power) )
  plt.draw()
  
  #now loop through the data again, applying the optimal filter that was just built.
  for i in range(f.GetEntries()):
    f.GetEntry(i)
  
    for ii in range(e.GetNumBoloPulses()):
      p = e.GetBoloPulse(ii)  #returns a KRawBoloPulseRecord
    
      if p.GetPulseLength() == 0: continue  #skip empty pulses
    
      if p.GetChannelName() == myChannel:    
        
        optChain.SetInputPulse( p.GetTrace() )
        optChain.RunProcess()
        
        amp = get_out( optChain )
        subamp = amp[ shiftVal - 50: shiftVal + 50]  #only look around the pretrigger
       
        results['amp'].append( subamp[subamp.argmax()] )
        results['time'].append( subamp.argmax() + shiftVal - 50 )
      
    
hall = TH1I('hall', 'hall', 1500, -100, 1400)      
hpeak =   TH1I('hpeak', 'hpeak', 1500, -100, 1400)
peakdist = TH1I('peak','peak relative to trigger time',512, 0, 512)
for peaktime in results['time']:
  peakdist.Fill(peaktime)

peakdist.Draw()

peakmin = float(raw_input('set minimum value ... '))
peakmax = float(raw_input('set maximum value ... '))

for i in range( len(results['amp']) ):
  hall.Fill( results['amp'][i] )
  if results['time'][i] < peakmax  and results['time'][i] > peakmin:
    hpeak.Fill( results['amp'][i] )
  
c2 = TCanvas('c2')
hpeak.SetLineColor(kRed)
hall.Draw()
hpeak.Draw('same')

saveme = raw_input('type "yes" if you want to save this: ').strip()
if saveme == 'yes':
  fresults = TFile('output.root', 'recreate')
  hall.Write()
  hpeak.Write()
  fresults.Close()
