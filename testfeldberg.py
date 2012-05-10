#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
import json
from couchdbkit import Server, Database
import sys,  datetime,  os, fnmatch, string, time, subprocess, signal, getpass
import numpy as np

def SetParameters(kampsite, serverurl, config, dataFile):

  f = ROOT.KDataReader(dataFile)
  e = f.GetEvent()
  chanList = {}
  fileDate = ''
  for entry in range(f.GetEntries()):
    f.GetEntry(entry)
    if (fileDate == '') and (e.GetNumSambas()>0):
      timestamp = e.GetSamba(0).GetNtpDateSec()
      fileDate = str(datetime.datetime.utcfromtimestamp(timestamp))
    for pulse in range(e.GetNumBoloPulses()):
      prec = e.GetBoloPulse(pulse)
      if prec.GetChannelName() not in chanList.keys():
        chanList[prec.GetChannelName()] = prec.GetPulseLength()
  f.Close()
  
  if os.path.isfile(config):
    configfile = open(config,'ro')
    json_data = json.load(configfile)
    configfile.close()
    
  server = Server(serverurl)
  db_template = server['pulsetemplates']
  for chan in chanList.keys():
    print "Set settings for channel:", chan, " pulse length:", chanList[chan]
    vr_template = db_template.view('analytical/bychandate', descending=True, reduce=False, startkey=[chan,'2013'], limit=1, include_docs=True)
    try:
      configpath = config.split('/')
      if configpath[0]=='db':
        db = server['multiprocessanadb']
        vr = db.view('channel/bydate', descending=True, reduce=False,startkey=[chan,'2013'], limit=1, include_docs=True)
        doc = vr.first()['doc']['ana']['multiFilterAndCorr'][0]
      elif os.path.isfile(config):
        if chan in json_data.keys():
          db = json_data
          doc = db[chan]['ana']['multiFilterAndCorr'][0]
        else:
          db = server['multiprocessanadb']
          vr = db.view('channel/bydate', descending=True, reduce=False,startkey=[chan,'2013'], limit=1, include_docs=True)
          doc = vr.first()['doc']['ana']['multiFilterAndCorr'][0] 
      delta_t = db[doc['templateparams']]['delta_t']
      #print "HERE"
      width = db[doc['templateparams']]['templatewidth']
      pretrigger = db[doc['templateparams']]['pretrigger']
      filters = doc['filters']
      
      for f in filters:
        kampsite.AddIIRFilter(chan,np.array(db[f]['coef']['a'][1:]).astype(float),len(db[f]['coef']['a'][1:]),np.array(db[f]['coef']['b']).astype(float),len(db[f]['coef']['b']))

      #DEBUG:
      #print delta_t, width, pretrigger
      #print filters
      
      doc_template = vr_template.first()['doc']
      if 'rise_double_decay' in doc_template['formula'].keys():
        funct = 'rise_double_decay'
      elif 'single_decay' in doc_template['formula'].keys():
        funct = 'single_decay'
      else:
        print "No known template function found!"
        
      exec(doc_template['formula'][funct]['python'])
       
      #!!! workaround for heat templates in ms scale, needs to be removed in the future
      templatePulse = ROOT.std.vector("double")() 
      par = doc_template['formula'][funct]['par']
      if funct == 'rise_double_decay':          
        par[0] /= 2.016
        par[2] /= 2.016
        par[3] /= 2.016
        par[5] /= 2.016
      for i in range(chanList[chan]):
        templatePulse.push_back(template(i,par))
      #print "Template for ",chan,":",np.array(templatePulse)
      
      kampsite.SetNormalizeTemplate(True)
      kampsite.SetTemplate(chan,templatePulse,pretrigger,delta_t,width)
      #print "test"
      if chan.find('chal') > -1:
        kampsite.SetPeakPositionSearchRange(chan,240,270)
      elif chan.find('ionis') > -1:
        kampsite.SetPeakPositionSearchRange(chan,2191,6191)
      kampsite.SetDoFit(doc['do_fit']==1)
      
      print "Settings for the channel ", chan, "were loaded"
    except:
      print "Something went wrong!"
      pass



if __name__ == '__main__':
  #serverurl = 'http://edwdbik.fzk.de:5984'
  serverurl = 'http://localhost:5984'
  
  config = 'db'

  filenum = int(sys.argv[2])
  run = sys.argv[1]
  starttime = datetime.datetime.now()
  
  #f = KDataReader('/Users/adam/analysis/edelweiss/data/kdata/raw/ma14f004_000.root')
  inputFile = '/Users/adam/analysis/edelweiss/data/kdata/raw/%s_00%d.root' % (run, filenum)
  outputFile = '/Users/adam/analysis/edelweiss/data/kdata/raw/%s_00%d.amp.feld.root' % (run ,filenum)
  kounselor = ROOT.KAmpKounselor()
  fbk = ROOT.KFeldbergKAmpSite()
  
  SetParameters(fbk, serverurl, config,  inputFile)
  print "KAmpSite is set up"
    
  kounselor.AddKAmpSite(fbk)
  print 'starting process'
  theRet = kounselor.RunKamp(inputFile,  outputFile)
  print "Kounselor return value: ",  theRet
  
  endtime = datetime.datetime.now()
  print 'time elapsed', endtime - starttime  
  
  
