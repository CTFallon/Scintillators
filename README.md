# Scintillators
Code for processing Cosmic Ray scintillator data for HGCAL Upgrade to CMS

## peCalc2 options:

* **-f, --file** : input file path. File should be produced using https://github.com/ChristopherRSettles/DRS4-Data-Analysis-Interface.
* **-s, --startPulse**: for manual detection of the beginning of the signal pulse
* **-e --endPulse**: for manual detenction of the end of the signal pulse
* **-g, --gainScale**: multiplier on the converstion factor. Mostly unused. Its easier to manually update the conversion factor that is hard-coded.
* **-d, --display**: creates plot to manually detect the pulse area (use alone, will stop program after producing the plot)
* **-peaks, --peaks**: used with Low-Light (i.e., calibration) files. You may need to manually adjust the limits of the fit functions.
* **-A, --AutoPulse**: automatically sets -s and -e based on the average pulse shape.

##### The output of the program is any one or more of the following:
###### A ROOT file named "<inputFileName>_analysed.root" that contains various histograms created (some may not be included based on which options were used above:
* A histogram of the pedestal.
* A histogram of the pedestal subtracted signal.
* A histogram of the PE count for all events.
* A histogram of the PE count for events that did not saturate the DRS4 and have more than half a PE.

###### In addition to the ROOT file, several .PNGs will be made:
* A scatter plot of time vs. signal, for all events.
* A histogram plot of PE counts for all events, with vertical lines showing the limits of the truncated mean calculation.

###### Finally, the most valuable output is written directly to terminal:
* The number of events recorded.
* The number of events used (i.e, events that did not saturate the DRS4 and have more than half a PE)
* The range that the signal was found to be in (repeats input if -s and -e were used, or tells us the values found if -A was used.)
* The mean and standard deviation (and their errors) of the PE count (Used) histogram.
* The truncated mean and standard deviation (and their errors) of the PE count (Used) histogram.
* Repeats this information in a single line for copying into a spreadsheet program.
