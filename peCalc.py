# pyROOT script to calculate the photoelectron (pe) yeild of a tile/SiPM cosmicRay setup


import ROOT as rt
import sys, os.path, argparse

def truncMean(histo, acc = 0.001, maxIter = 100):
	oldMean = histo.GetMean()
	newMean = oldMean + 1.0
	i = 0
	while abs(oldMean - newMean) > acc:
		histo.GetXaxis().SetRangeUser(newMean*0.2, newMean*2.0)
		oldMean = newMean
		newMean = histo.GetMean()
		histo.GetXaxis().SetRange()
		i += 1
		if i > maxIter:
			sys.exit("maximum number of iterations reached for truncated mean")
	else:
		histo.GetXaxis().SetRangeUser(newMean*0.2, newMean*2.0)
		meanPE = histo.GetMean()
		sigma = histo.GetStdDev()
		meanErr = histo.GetMeanError()
		sigmaErr = histo.GetStdDevError()
		print("TruncMean " + str(meanPE) + " " + str(meanErr) + " " +str(sigma) + " " + str(sigmaErr))


parser = argparse.ArgumentParser(description="Calculate the photoelectron (p.e.) yield of a cosmicRay data set.")

parser.add_argument('-f','--file', dest='inFileName', action='store', default='test.root', help='The input file\'s name. Full path if not in running directory. .root extension')
parser.add_argument('-s','--startPulse',dest='pulseStart', action='store',type=int, help='The beginning of the pulse in data. Use -d flag to display plot.')
parser.add_argument('-e','--endPulse',dest='pulseEnd', action='store',type=int, help='The end of the pulse in data. Use -d flag to display plot.')
parser.add_argument('-g','--gainScale', dest='gainScale', action='store', default=1, type=float, help='The gain scale. Default = 1')
parser.add_argument('-d','--display', dest='showPlot', action='store_true', default=False, help='Creates image of all pulses overlayed on one another and quits the program.')
parser.add_argument('-peaks', '--peaks', dest='doPEconversionScaleCalculation',action='store_true',default=False, help="Uses the first five peaks in the RAW output to calculate the Volts-to-PE scalar. You may need to inspect and adjust the limits on each PE gaus function.")
parser.add_argument('-ch2', dest='dataInCh2', action='store_true', default=False, help="Use this flag if the data was stored in channel 2 of the DRS4. (post may 14, 2018. Channel 1 died ?)")

args = parser.parse_args()

inputFile = rt.TFile(args.inFileName, "READ") # saving only the parts of the root file we're intrested in. For later use and easier access
oldtree = inputFile.Get("T")
oldtree.SetBranchStatus("*",0)
oldtree.SetBranchStatus("event",1) # we may not even need this branch...
oldtree.SetBranchStatus("c1",1)
oldtree.SetBranchStatus("c2",1)
oldtree.SetBranchStatus("c3",1)
oldtree.SetBranchStatus("c4",1)

#auto detect what channel data is in
oldtree.Draw('c1>>ch1Hist')
ch1Hist = rt.gDirectory.Get("ch1Hist")
oldtree.Draw('c2>>ch2Hist')
ch2Hist = rt.gDirectory.Get("ch2Hist")
oldtree.Draw('c3>>ch3Hist')
ch3Hist = rt.gDirectory.Get("ch3Hist")
oldtree.Draw('c4>>ch4Hist')
ch4Hist = rt.gDirectory.Get("ch4Hist")
chMeans = [ch1Hist.GetMean(),ch2Hist.GetMean(),ch3Hist.GetMean(),ch4Hist.GetMean()]

dataCh = chMeans.index(min(chMeans))+1

print("data in channel: " + str(dataCh))

oldtree.SetBranchStatus("c1",0)
oldtree.SetBranchStatus("c2",0)
oldtree.SetBranchStatus("c3",0)
oldtree.SetBranchStatus("c4",0)


oldtree.SetBranchStatus('c{}'.format(dataCh),1)
oldtree.SetBranchStatus('t{}'.format(dataCh),1)

if args.showPlot:
	c = rt.TCanvas('c','c',2000,1000)
	oldtree.Draw("-c{}:t{}".format(dataCh,dataCh))
	c.SaveAs(args.inFileName[:-5]+".png")
	sys.exit("Examine the plot and find the pulse edges. Then rerun this script with the correct arguments. (-s [pulseStartBin] -e [pulseEndBin])")

outputFile = rt.TFile(args.inFileName[:-5]+"_analysed.root","RECREATE")
tree = oldtree.CloneTree()
inputFile.Close()
#rename data leaf for easy access later
tree.GetLeaf("c{}".format(dataCh)).SetTitle("data")
tree.GetLeaf("c{}".format(dataCh)).SetName("data")

tree.GetLeaf("t{}".format(dataCh)).SetTitle("time")
tree.GetLeaf("t{}".format(dataCh)).SetName("time")

nEntries = tree.GetEntries()
print("Total Events:" +str(nEntries))

totalEventsOver0p5 = 0
totalBinsOver0p5 = 0



pulseDelta  = args.pulseEnd - args.pulseStart

if args.pulseStart - pulseDelta <= 30:
	sys.exit("You dumbo! Pedestal calculation will be in negative time! Make sure 2*pS-pE > 30!")

#Values for signal->pe conversion taken from Ping
# 0.96 is Ping's GainScale.
# I wanted to set the default to 1, but not messup these numbers, so I hardcoded it in.
# these values are depreciated anyway.

gain_at_1p8V = 8.98e5 * 0.96
gain_at_3V = 1.7e6 * 0.96 * 0.978543
gain_at_3V_2050VE = 1.76e6 * 0.96
conversion_factor_old = 1e9/gain_at_3V_2050VE*6.24/50.0/13.0

conversion_factor = 5.1783804406 # from D1 Low Light test

hist_pe_All = rt.TH1F("pe","Calculated photoelectron count, all events;p.e.;count",400,0,100)
hist_pe_Used = rt.TH1F('hist_pe_Used','Calculated photoelectron count;p.e.;count, no OV',400,0,100)
hist_RAW = rt.TH1F('hist_RAW','Raw Output;ADC/Integrated Voltage;Count',200,-1,2)
hist_Ped = rt.TH1F("hist_Ped","Pedestal Output;ADC/Voltage;Count",500,0.5,1.1)

rt.gStyle.SetOptStat("MRen")


#dataVector = []
#tree.SetBranchAddress('c{}'.format(dataCh),dataVector)

for iEvent in range(nEntries):
	eventOver0p5flag = False
	tree.GetEvent(iEvent)
	dataVector = tree.data
	ped = 0 # Calculate the pedestal value
	for pedBin in range(args.pulseStart-30-pulseDelta,args.pulseStart-30):
		ped -= dataVector[pedBin] # minus because signal is negative
	hist_Ped.Fill(ped)
	sig = -ped # Calculate the integral of the pulse, corrected for the pedestal
	#NOTE: This assumes the dt of the pedestal is the same as the dt of the signal region
	# this dt is the variable 'deltaPulse'
	# if dt is different for the pedestal and pulse, you need to scale the pedestal subtraction corretly. You're smart, you can do it!
	for sigBin in range(args.pulseStart, args.pulseEnd):
		sig -= dataVector[sigBin] # minus because signal is negative
		if -dataVector[sigBin] >= 0.499908: # check to make sure that the event didnt produce a signal greater than 0.5 volts. 
			totalBinsOver0p5 += 1
			eventOver0p5flag = True
	if eventOver0p5flag:
		totalEventsOver0p5 += 1
	
	# Convert integrated pulse into # p.e.
	pe = sig*conversion_factor
	hist_pe_All.Fill(pe)
	hist_RAW.Fill(sig)
	if (pe > 0.5 and not eventOver0p5flag):
		hist_pe_Used.Fill(pe)

print("Total number of Events (Bins) over voltage: "+str(totalEventsOver0p5)+ " ("+str(totalBinsOver0p5)+")")

if args.doPEconversionScaleCalculation:
	pe1 = rt.TF1("pe1",'gaus',0.13,0.29)
	pe1.SetLineColor(rt.kGreen)
	pe2 = rt.TF1("pe2",'gaus',0.29,0.48)
	pe2.SetLineColor(rt.kGreen)
	pe3 = rt.TF1("pe3",'gaus',0.48,0.67)
	pe3.SetLineColor(rt.kGreen)
	pe4 = rt.TF1("pe4",'gaus',0.67,0.85)
	pe4.SetLineColor(rt.kGreen)
	pe5 = rt.TF1("pe5",'gaus',0.85,1.05)
	pe5.SetLineColor(rt.kGreen)
	pe6 = rt.TF1("pe6",'gaus',1.05,1.27)
	pe6.SetLineColor(rt.kGreen)
	total = rt.TF1("sixPeaks","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",0.13,1.27)
	
	hist_RAW.Fit(pe1,"R")
	hist_RAW.Fit(pe2,"R+")
	hist_RAW.Fit(pe3,"R+")
	hist_RAW.Fit(pe4,"R+")
	hist_RAW.Fit(pe5,"R+")
	hist_RAW.Fit(pe6,"R+")
	
	total.SetParameter(0, pe1.GetParameter(0))
	total.SetParameter(1, pe1.GetParameter(1))
	total.SetParameter(2, pe1.GetParameter(2))
	total.SetParameter(3, pe2.GetParameter(0))
	total.SetParameter(4, pe2.GetParameter(1))
	total.SetParameter(5, pe2.GetParameter(2))
	total.SetParameter(6, pe3.GetParameter(0))
	total.SetParameter(7, pe3.GetParameter(1))
	total.SetParameter(8, pe3.GetParameter(2))
	total.SetParameter(9, pe4.GetParameter(0))
	total.SetParameter(10, pe4.GetParameter(1))
	total.SetParameter(11, pe4.GetParameter(2))
	total.SetParameter(12, pe5.GetParameter(0))
	total.SetParameter(13, pe5.GetParameter(1))
	total.SetParameter(14, pe5.GetParameter(2))
	total.SetParameter(15, pe6.GetParameter(0))
	total.SetParameter(16, pe6.GetParameter(1))
	total.SetParameter(17, pe6.GetParameter(2))
	
	hist_RAW.Fit(total,"R+")
	
	print("Current Conversion Factor = " +str(conversion_factor))
	newCF = (1.0/total.GetParameter(1)+2.0/total.GetParameter(4)+3.0/total.GetParameter(7)+4.0/total.GetParameter(10)+5.0/total.GetParameter(13)+6.0/total.GetParameter(16))/6.0
	print("New Conversion Factor = " + str(newCF))

outputFile.Write()

meanPE = hist_pe_Used.GetMean()
sigma = hist_pe_Used.GetStdDev()
meanErr = hist_pe_Used.GetMeanError()
sigmaErr = hist_pe_Used.GetStdDevError()

print("Events Counted: " +str(int(hist_pe_Used.GetEntries())) + " (p.e. > 0.5 and not overvoltage)")
print("                    p.e.           err        stdDev            err")
print("Full Mean " + str(meanPE) + " " + str(meanErr) + " " +str(sigma) + " " + str(sigmaErr))
truncMean(hist_pe_Used)
print("Pulse Range: " + str(args.pulseStart) + " - " + str(args.pulseEnd))


	
outputFile.Close()

