#!/usr/bin/python

import sys
sys.argv.append('-b')
import ROOT

ff = ROOT.TFile(sys.argv[1])
f1 = ROOT.TF1('f1', 'gaus(0)+pol0(3)', 0,1)

h1 = ff.Get("ebeam_hebth").ProjectionY()
f1.SetParameter(1, h1.GetBinCenter(h1.GetMaximumBin()))
f1.SetParameter(2, h1.GetRMS())
mu,sig = [f1.GetParameter(i) for i in range(1,3)]

f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
h1.Fit(f1)
mu,sig = [f1.GetParameter(i) for i in range(1,3)]

f1.SetRange(mu-2.5*abs(sig), mu+2.5*abs(sig))
h1.Fit(f1, 'R')

c1 = ROOT.TCanvas('c1','c1',1100,800)

h1.Draw()
c1.Print("ebeam.pdf")
