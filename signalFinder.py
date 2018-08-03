import ROOT as rt
import sys, os.path, argparse


# function to find the start (and end) of a signal in a (large) set of samples


def findSignal(tree):
	averageSignal = [0 for x in range(len(tree.data))]
	nEntries = tree.GetEntries()
	for iEvent in range(nEntries):
		tree.GetEvent(iEvent)
		for x in range(len(tree.data)):
			averageSignal[x] -= tree.data[x] #negative because singal is negative
	for x in range(len(averageSignal)-6):
		test1 = int(abs(averageSignal[x]/float(nEntries)-averageSignal[x+1]/float(nEntries)) > 0.0001)
		test2 = int(abs(averageSignal[x+1]/float(nEntries)-averageSignal[x+2]/float(nEntries)) > 0.0001)
		test3 = int(abs(averageSignal[x+2]/float(nEntries)-averageSignal[x+3]/float(nEntries)) > 0.0001)
		test4 = int(abs(averageSignal[x+3]/float(nEntries)-averageSignal[x+4]/float(nEntries)) > 0.0001)
		test5 = int(abs(averageSignal[x+4]/float(nEntries)-averageSignal[x+5]/float(nEntries)) > 0.0001)
		test6 = int(abs(averageSignal[x+5]/float(nEntries)-averageSignal[x+6]/float(nEntries)) > 0.0001)
		if test1+test2+test3+test4+test5+test6 >= 4:
			return int(x/10)*10, int(x/10)*10+200
	sys.exit("Auto-Detect Failed, please manually find signal region")
