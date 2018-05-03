# pyROOT script to calculate the photoelectron (pe) yeild of a tile/SiPM cosmicRay setup

# Porgram flow:
# Convert integrated Voltage into p.e.
# calculate mean and RMS of all events' p.e. yeild
# report

import ROOT as rt

import sys, os.path, argparse

parser = argparse.ArgumentParser(description="Calculate the photoelectron (p.e.) yield of a cosmicRay data set.")

parser.add_argument('-f','--file', dest='inFileName', action='store', default='test.root', help='The input file\'s name. Full path if not in script directory. .root extension')
parser.add_argument('-s','--startPulse',dest='pulseStart', action='store',type=int, help='The beginning of the pulse in data. Use -d flag to display plot.')
parser.add_argument('-e','--endPulse',dest='pulseEnd', action='store',type=int, help='The end of the pulse in data. Use -d flag to display plot.')
parser.add_argument('-g','--gainScale', dest='gainScale', action='store', default=0.96, type=float, help='The gain scale. Default = 0.96')
parser.add_argument('-d','--display', dest='showPlot', action='store_true', default=False, help='Creates image of all pulses overlayed on one another and quits the program.')

args = parser.parse_args()

inputFile = rt.TFile(args.inFileName, "READ") # saving only the parts of the root file we're intrested in. For later use and easier access
oldtree = inputFile.Get("T")
oldtree.SetBranchStatus("*",0)
oldtree.SetBranchStatus("event",1) # we may not even need this branch...
oldtree.SetBranchStatus("c1",1) # pulse height for channel 1
oldtree.SetBranchStatus("t1",1) # time values for channel 1

outputFile = rt.TFile(args.inFileName[:-5]+"_analysed.root","RECREATE")
tree = oldtree.CloneTree()
inputFile.Close()

nEntries = tree.GetEntries()
print("Total Events:" +str(nEntries))

totalEventsOver0p5 = 0
totalBinsOver0p5 = 0

if args.showPlot:
	c = rt.TCanvas('c','c',2000,1000)
	tree.Draw("-c1:t1")
	c.SaveAs(args.inFileName[:-5]+".png")
	sys.exit("Examine the plot and find the pulse edges. Then rerun this script with the correct arguments.")

pulseDelta  = args.pulseEnd - args.pulseStart

if args.pulseStart - pulseDelta <= 30:
	sys.exit("You dumbo! Pedestal calculation will be in negative time! Make sure 2*pS-pE > 30!")
peList = []

#Values for signal->pe conversion taken from Ping
gain_at_1p8V = 8.98e5 * args.gainScale
converstion_factor = 1e9/gain_at_1p8V*6.24/50.0/13.0
gain_at_3V = 1.7e6 * args.gainScale * 0.978543
gain_at_3V_2050VE = 1.76e6 * args.gainScale
converstion_factor = 1e9/gain_at_3V_2050VE*6.24/50.0/13.0

for iEvent in range(nEntries):
	eventOver0p5flag = False
	tree.GetEvent(iEvent)
	ped = 0 # Calculate the pedestal value
	for pedBin in range(args.pulseStart-30-pulseDelta,args.pulseStart-30):
		ped -= tree.c1[pedBin] # minus because signal is negative

	sig = -ped # Calculate the integral of the pulse, corrected for the pedestal
	#NOTE: This assumes the dt of the pedestal is the same as the dt of the signal region
	# this dt is the variable 'deltaPulse'
	# if dt is different for the pedestal and pulse, you need to scale the pedestal subtraction corretly. You're smart, you can do it!
	for sigBin in range(args.pulseStart, args.pulseEnd):
		sig -= tree.c1[sigBin] # minus because signal is negative
		if -tree.c1[sigBin] >= 0.499908: # check to make sure that the event didnt produce a signal greater than 0.5 volts. 
			totalBinsOver0p5 += 1
			print(str(iEvent) + " " + str(sigBin))
			eventOver0p5flag = True
	if eventOver0p5flag:
		totalEventsOver0p5 += 1
	
	# Convert integrated pulse into # p.e.
	pe = sig*converstion_factor
	if (pe>3 and not eventOver0p5flag):
		peList.append(sig*converstion_factor)

print("Total number of Events over voltage: "+str(totalEventsOver0p5))
print("Total number of Bins over voltage: "+str(totalBinsOver0p5))

sumPE = 0

peHist = rt.TH1F('peHist','Calculated photoelectron count',500,min(peList)-1,max(peList)+1)
peHist.SetXTitle("Number of photoelectrons")
peHist.SetYTitle("Count")
rt.gStyle.SetOptStat("MRen")

for value in peList:
	sumPE += value
	peHist.Fill(value)

c = rt.TCanvas('c','c',2000,1000)
peHist.Draw()
c.SaveAs(args.inFileName[:-5]+"_pe.png")
outputFile.Write()
meanPE = sumPE/len(peList) # arithmetic average

sigma = 0
for value in peList: # sigma^2 = (1/N)*Sigma{n=1->N}((x_n-mu)^2) definition of standard deviation
	sigma += (value-meanPE)**2
sigma = rt.TMath.Sqrt(sigma/len(peList))

meanErr = sigma/rt.TMath.Sqrt(len(peList)) # standard error of the mean, defined as standard deviation of the sample divided by sqrt of sample size


print("Events Counted: " +str(len(peList)) + " (p.e. > 3 and not overvoltage)")
print("p.e. = "+ str(meanPE) + " +- " + str(meanErr) + " +- " +str(sigma) + " (mu +- err_mu +- stdDev)")

