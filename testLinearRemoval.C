void go(TH1D &input, TH1D &output){
  
  double slope = 0.005;
  double offset = 3000.;
  vector<double> pulse;
  for(unsigned int i = 0; i < 8196; i++) {
    pulse.push_back(slope * i + offset);
  }
  KLinearRemoval lin  ;
  lin.SetInputPulse(pulse);
  cout <<  "input pulse " <<  lin.GetInputPulse() << endl;
  cout << "output pulse " <<  lin.GetOutputPulse() << endl;
  unsigned int theRet = lin.RunProcess();
  cout << "running linear removal " <<  theRet << endl;
  //cout << "slope "  << lin.CalculateSlope() << endl;

  for(unsigned int i = 0; i < 8196; i++){
    input.SetBinContent(i+1, lin.GetInputPulse()[i]);
    output.SetBinContent(i+1, lin.GetOutputPulse()[i]);
  }
  
}