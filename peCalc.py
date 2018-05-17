# pyROOT script to calculate the photoelectron (pe) yeild of a tile/SiPM cosmicRay setup


import ROOT as rt
import sys, os.path, argparse

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
if args.dataInCh2:
	oldtree.SetBranchStatus("c2",1) # pulse height for channel 2
	oldtree.SetBranchStatus("t2",1) # time values for channel 2
else:
	oldtree.SetBranchStatus("c1",1) # pulse height for channel 1
	oldtree.SetBranchStatus("t1",1) # time values for channel 1

if args.showPlot:
	c = rt.TCanvas('c','c',2000,1000)
	if args.dataInCh2:
		oldtree.Draw("-c2:t2")
	else:
		oldtree.Draw("-c1:t1")		
	c.SaveAs(args.inFileName[:-5]+".png")
	sys.exit("Examine the plot and find the pulse edges. Then rerun this script with the correct arguments. (-s [pulseStartBin] -e [pulseEndBin])")

outputFile = rt.TFile(args.inFileName[:-5]+"_analysed.root","RECREATE")
tree = oldtree.CloneTree()
inputFile.Close()

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

hist_pe_All = rt.TH1F("pe","Calculated photoelectron count, all events",400,0,100)
hist_pe_All.SetXTitle("# p.e.")
hist_pe_All.SetYTitle("Count")

hist_pe_Used = rt.TH1F('hist_pe_Used','Calculated photoelectron count, no OV',400,0,100)
hist_pe_Used.SetXTitle("# p.e.")
hist_pe_Used.SetYTitle("Count")

hist_RAW = rt.TH1F('hist_RAW','Raw Output',200,-1,2)
hist_RAW.SetXTitle("ADC/Integrated Voltage")
hist_RAW.SetYTitle("Count")

hist_Ped = rt.TH1F("hist_Ped","Pedestal Output;ADC/Voltage;Count",500,0.5,1.1)

rt.gStyle.SetOptStat("MRen")

for iEvent in range(nEntries):
	eventOver0p5flag = False
	tree.GetEvent(iEvent)
	if args.dataInCh2:
		dataVector = tree.c2
	else:
		dataVector = tree.c1
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
print("         p.e.           err        stdDev            err")
print(str(meanPE) + " " + str(meanErr) + " " +str(sigma) + " " + str(sigmaErr))
print("Pulse Range: " + str(args.pulseStart) + " - " + str(args.pulseEnd))


	
outputFile.Close()

