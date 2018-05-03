#include <iostream>
#include <fstream>



#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

#include "TFile.h"
#include "TTree.h"


void setTDRStyle() {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(480); //Height of canvas
  //  tdrStyle->SetCanvasDefH(450); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //  tdrStyle->SetErrorMarker(20);
  //  tdrStyle->SetErrorX(0.);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
  // tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.03);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.055, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  //tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.00);
  tdrStyle->SetTitleYOffset(1.35);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.04, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(508, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
 cout << endl << "    using TDR format: gROOT->SetStyle(\"tdrStyle\");" << endl;
  gROOT->SetStyle("tdrStyle");


}




Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0];

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]) + par[4]; // add linear term
}





TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());
   //   std::cout << FunName << std::endl;

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],5);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma", "Pedstal");

   for (i=0; i<5; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   // his->Draw("HIST");
   his->Fit(ffit,"R", "", fitrange[0], fitrange[1] );   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<5; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);

      if (l < lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);

      if (l > lold)
         step = -step/10;

      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}



void convertToCSV(const char *filename, int ientry = -1, const char *csv_file_tag = ".csv", int channel=3) {


  TFile *myfile= new TFile(filename,"READ");
  TTree *tree= (TTree* )myfile->Get("T");
  Int_t nentries= tree->GetEntries();

  std::string rootfile_tag(".root");
  std::string output_filename(filename);
  rootfile_tag_pos = output_filename.find(rootfile_tag);

  if (int(rootfile_tag_pos)>=0) {

    output_filename.replace(rootfile_tag_pos, rootfile_tag.size(), "");
  }
  output_filename += csv_file_tag;

  std::fstream output_file;
  //  try {

  output_filename.resize(0);
  output_filename += csv_file_tag;
  //  output_filename += "testme.csv";
  output_file.open(output_filename, ios::out|ios::trunc) ;
    //  } catch (const ifstream::failure& e){
    //   std::cout << "exception opening file " << std::endl;
    // break;
    // }
  if (!output_file.is_open() ) {
    std::cout << "output file is not properly opened" << std::endl;
    return;
  }
  std::cout << "convert root file into .csv file ---> " 
	    << output_filename 
	    << std::endl;

  float c1[1024];
  float c2[1024];
  float c3[1024];
  float c4[1024];
  
  int t1[1024];
  int t2[1024];
  int t3[1024];
  int t4[1024];



  tree->SetBranchAddress( (std::string("c")+= std::to_string(channel)).data(), c1);
  tree->SetBranchAddress( (std::string("t")+= std::to_string(channel)).data(), t1);


  std::cout << " total events = " << nentries << std::endl;
  for (int ii=00; ii<nentries; ii++) {
    //    for (int ii=11000; ii<14000; ii++) {

    if (ientry>=0) {
      if (ii==ientry)      tree->GetEntry(ientry);
      else continue;
    } else {

      tree->GetEntry(ii);

      if (ii%1000 ==0) std::cout << "processing " << (ii*1.0/nentries)*100 << "%" << std::endl;
    }


    for (int jj =0; jj < 1024; jj++) {

      float bin_content  = c1[jj];

      bin_content = -bin_content;

      if (jj != (1024-1) ) output_file << bin_content << ",";
	else output_file << bin_content;
    }
    output_file << std::endl;
  }
  output_file.close();
}

void analyzePulse(const char *filename, int ientry=12, const char *output_filename="test.root", float start_pulse=480, float end_pulse=580, bool calibration=false, float gain_scale = 1.0, bool invert_signal=true) {


  setTDRStyle();

  TFile *myfile= new TFile(filename,"READ");
  TTree *tree= (TTree* )myfile->Get("T");
  Int_t nentries= tree->GetEntries();


  float c1[1024];
  float c2[1024];
  float c3[1024];
  float c4[1024];
  
  int t1[1024];
  int t2[1024];
  int t3[1024];
  int t4[1024];

  tree->SetBranchAddress("c1", c1);
  tree->SetBranchAddress("t1", t1);

  TH1F *apulse = 0;
  TH1F *aindex = 0;


  // fitting function
  TF1 *func=0;
  //  TF1 *func = new TF1("func","[0]*TMath::Landau(x,[1],[2]) + [3]",800,1024);
  // func->SetParameters(0.01,850,10, 0.001); //for example

  TF1 *func_exp = new TF1("func_exp", "[0]*(exp(-(x-[1])/[2]) - exp(-(x-[1])/[3]) )+[4]", 1, 1024);
  func_exp->SetLineColor(kRed);


  TF1 *func_pol0 = new TF1("func_pol0", "[0]", 1, 1024);

  std::vector<float> values;
  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[5], pllo[5], plhi[5], fp[5], fpe[5];
  // default 120 bins
  fr[0]=start_pulse;
  fr[1]=end_pulse;

  //  fr[0]=640;
  // fr[1]=740;
  // fr[1]=720;

  double pulse_range =  fr[1] - fr[0];


  Double_t chisqr;
  Int_t    ndf;

  Double_t SNRPeak, SNRFWHM;


  //  TH1F *hist_noise = new TH1F("hist_noise", "", 100, -0.004*pulse_range , 0.004 * pulse_range);
  TH1F *hist_noise = new TH1F("hist_noise", "", 100, 0,20);
  hist_noise->GetXaxis()->SetTitle("p.e. equivalent");
  hist_noise->GetYaxis()->SetTitle("Counts");

  // the binning was 8 *250; now try 10*250; Jan 23, 2018
  TH1F *hist_charge = new TH1F("hist_charge", "", 10*250, -5, 245);
  TH1F *hist_charge_zoomin = new TH1F("hist_charge_zoomin", "", 30, -5, 175);
  TH1F *hist_charge_exclzero = new TH1F("hist_charge_exclzero", "", 6*242, 3, 245);
  hist_charge->GetXaxis()->SetTitle("Charge (arbitary unit)");
  hist_charge->GetXaxis()->SetTitle("p.e.");
  hist_charge->GetYaxis()->SetTitle("Counts");


  hist_charge_zoomin->GetXaxis()->SetTitle("p.e.");
  hist_charge_zoomin->GetYaxis()->SetTitle("Counts");

  hist_charge_exclzero->GetXaxis()->SetTitle("p.e.");
  hist_charge_exclzero->GetYaxis()->SetTitle("Counts");

  TH1F *hist_charge_calibration = (TH1F *) hist_charge->Clone("hist_charge_calibration");


  int navg = 10;
  int navg_counter =0;
  double avg =0;

  TH1F *hist_noise_entry = new TH1F("hist_noise_entry", "", nentries, 0, nentries);
  TH1F *hist_noise_entry_avg = new TH1F("hist_noise_entry_avg", "", (int)(nentries/navg), 0, (int)(nentries/navg));
  hist_noise_entry->GetXaxis()->SetTitle("Event number");
  hist_noise_entry->GetYaxis()->SetTitle("Pedstal (arbitary unit)");

  TCanvas *cc1 = new TCanvas("cc1", "", 1200, 600);
  cc1->Divide(2, 1);
  char temp[1024];
  bool first_processed_event = true;
  int  passed_events    = 0;

  //  nentries = 4000;

  apulse = new TH1F("base_pulse", "", 1024, 0, 1024);
  apulse->SetMarkerSize(0.005);
  apulse->GetXaxis()->SetTitle("Time (ns)");
  apulse->GetYaxis()->SetTitle("Negative of Pulse Voltage (V)");
  TH1F *mean_pulse = (TH1F *) apulse->Clone("mean_pulse");

  double gain_at_1p8V = 8.98e5 * gain_scale;
  double conversion_factor = 1e9/gain_at_1p8V*6.24/50.0/13.0;

  double gain_at_3V = 1.7e6 * gain_scale * 0.978543; // final calibration from gain curve of DCR.
 
  double gain_at_3V_2050VE = 1.7e6 * gain_scale;

  conversion_factor = 1e9/gain_at_3V_2050VE*6.24/50.0/13.0;

  // old note
  //  double conversion_factor = 1.05*1.0/50*6.24/13/4.0*(1e4)/1.23/gain_scale; //V/50*pulse_range*1e-9*6.24 * 10e18/13/4e5/1.23; // 1.23 from the actual single p.e. counting
  int num_above_1pe = 0, num_above_3pe = 0;
  std::cout << " total events = " << nentries << std::endl;
  for (int ii=00; ii<nentries; ii++) {
    //    for (int ii=11000; ii<14000; ii++) {

    if (ientry>=0) {
      if (ii==ientry)      tree->GetEntry(ientry);
      else continue;
    } else {

      tree->GetEntry(ii);

      if (ii%1000 ==0) std::cout << "processing " << (ii*1.0/nentries)*100 << "%" << std::endl;
    }

    sprintf( temp, "apulse_%d", ii);
    apulse = (TH1F *) mean_pulse->Clone(temp);

    if (first_processed_event== true)  {
      aindex = (TH1F *)apulse->Clone("index");
      first_processed_event = false;
    }
    bool is_noise_event = false;
    int noise_bins =0;
    double maximum_signal_pulse=-9999;
    for (int jj =0; jj < 1024; jj++) {

      float bin_content  = c1[jj];

      if (invert_signal) bin_content = -bin_content;
      apulse->SetBinContent( jj+1, bin_content);

      if (jj>= fr[0] && jj<fr[1]) {

	if (bin_content>maximum_signal_pulse) maximum_signal_pulse = bin_content;
      }
      //      if (jj>= (fr[1]+50) && (maximum_signal_pulse > 0.01) && (fabs(c1[jj]) > maximum_signal_pulse)) {
      if (jj>= (fr[1]+20) && fabs(bin_content) >0.01) {

	noise_bins ++;
	//	is_noise_event = true; 
	//	break;
      }
    }
    //        if (noise_bins>=3) continue;  // filter out noise event

    passed_events ++;


    
    //    return;

    apulse->SetMinimum(0.0);
    apulse->SetMaximum(0.25);


   
    pllo[0]=1;  pllo[1]=800;  pllo[2]=0.1;  pllo[3]=0.01;  pllo[4]=0; 
    plhi[0]=30;  plhi[1]=900;  plhi[2]=20;  plhi[3]=20;  plhi[4]=0.1;
    sv[0]=10; sv[1] = 870; sv[2]=0.27; sv[3]=1; sv[4]=0.002; 


    double pedastal = 0;
    double charge =0;

   
    pedastal = apulse->Integral( aindex->Fill( fr[0] - 30 - pulse_range ), 
				 aindex->Fill( fr[0] - 30));


    // filter
    // if ( fabs(fabs(pedastal) - 0.3547)>0.2*0.0164) continue;

    ///( aindex->Fill( fr[0]) -  aindex->Fill( fr[0] - pulse_range ) );    
    charge = apulse->Integral( aindex->Fill(fr[0]), aindex->Fill(fr[1] ) );///(aindex->Fill( fr[1]) - aindex->Fill( fr[0])  );

    //    std::cout << charge-pedastal << std::endl;

    if (pedastal/( aindex->Fill( fr[0]) -  aindex->Fill( fr[0] - pulse_range ) )>0.01) continue; // skip events with unstable baseline



    hist_noise->Fill( pedastal * conversion_factor);
    // covert to p.e.
    // 1.05 is the gain drift on DRS

    hist_charge->Fill( (charge-pedastal) * conversion_factor);
    hist_charge_exclzero->Fill( (charge-pedastal) * conversion_factor);
    hist_charge_zoomin->Fill( (charge-pedastal) * conversion_factor);
    if ( (charge-pedastal) * conversion_factor > 1)  num_above_1pe ++;
    if ( (charge-pedastal) * conversion_factor > 3)  num_above_3pe ++;

    for (int kk =0; kk < values.size(); kk++) {
      hist_charge_calibration->Fill( values[kk] * conversion_factor);
    }

    if ((ii+1)%navg !=0 ) {

      avg +=   pedastal;
      navg_counter ++;

    } else {

      if (navg_counter)      hist_noise_entry_avg->SetBinContent( (ii+1)/navg + 1, avg/navg_counter);
      navg_counter = 0;
      avg =0;

    }

    hist_noise_entry->SetBinContent( ii+1, pedastal);

    cc1->cd(1); 
    
    apulse->SetMarkerColor(kRed);
    apulse->SetLineColor(kRed);
    //    apulse->SetLineWidth(0.01);
    if (passed_events ==1) apulse->Draw("p"); else apulse->Draw("samep");
    mean_pulse->Add( apulse);
      
    //    apulse->Draw("hist");
    // return;

    /*
    if (ientry>=0) {
      func = langaufit(apulse,fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
      func->SetLineColor(kRed);
      //apulse->Fit(func, "R", "", 800, 1024);
    //std::cout << c1[0] << std::endl;


    //    maximum_height->Fill( apulse->GetMaximum() - func->GetParameter(4) );
      apulse->Draw("HIST");
    //    apulse->Fit(func, "", "", fr[0], fr[1]);
      func->Draw("SAME");
    //langaupro(fp,SNRPeak,SNRFWHM);
    }
    */
  }


  std::cout << "numbers above 1 p.e. = " << num_above_1pe << std::endl;
  std::cout << "numbers above 3 p.e. = " << num_above_3pe << std::endl;


  cc1->cd(2); 
  mean_pulse->Scale( 1.0/passed_events);
  mean_pulse->Draw("HIST");


  TCanvas *cc2 = new TCanvas("cc2", "", 1200, 600);
  cc2->Divide(2, 1);
  cc2->cd(1);   hist_noise->Draw("HIST");
  cc2->cd(2);   hist_charge->Draw("HIST");
  
  std::cout << "total area" << hist_charge->Integral(1, hist_charge->GetNbinsX(), "width") << std::endl;

  TCanvas *cc3 = new TCanvas("cc3", "", 1200, 600);
  hist_noise_entry->Draw("hist");

  TCanvas *cc4 = new TCanvas("cc4", "", 1200, 600);
  hist_noise_entry_avg->Draw("hist");



  TCanvas *cc5 = new TCanvas("cc5", "");
  hist_charge_calibration->Draw("hist");

  TCanvas *cc6 = new TCanvas("cc6", "");
  hist_charge_exclzero->Draw("hist");



  TCanvas *cc7 = new TCanvas("cc7", "");
  hist_charge_zoomin->Draw("PE");
  cc7->Print( (TString(hist_charge_zoomin->GetName()) += ".pdf").Data() );




  TFile *output_file = new TFile(output_filename, "RECREATE");

  output_file->cd();
  mean_pulse->Write();
  hist_noise->Write();
  hist_charge->Write();
  hist_charge_zoomin->Write();
  hist_charge_exclzero->Write();

  std::cout << hist_charge_exclzero->GetMean() << "+/-" << hist_charge_exclzero->GetMeanError() << "; RMS = " << hist_charge_exclzero->GetRMS()  <<"+/-"<< hist_charge_exclzero->GetRMSError()  << std::endl;

  hist_noise_entry->Write();
  hist_noise_entry_avg->Write();
  hist_charge_calibration->Write();

  cc1->Write();
  cc2->Write();
  cc3->Write();
  cc4->Write();
  cc5->Write();
  cc6->Write();
  cc7->Write();

  cc2->Print("hist_signal_pedstal.pdf");
  cc6->Print("hist_charge_exclzero.pdf");

  //  output_file->Write();
}


Double_t signal_shape(Double_t *x, Double_t *par) {

  Double_t N0 = par[0];
  Double_t mean0= par[1];
  Double_t sigma_ele=par[2];
  Double_t sigma0 = sigma_ele *sqrt(2);

  Double_t N1 = par[3]*1;
  Double_t mean1= par[4]+mean0;
  Double_t sigma_pho=par[5];
  Double_t sigma1 = sqrt( pow(sigma0, 2) +pow(sigma_pho, 2) );

  Double_t N2 = par[6];
  Double_t mean2= par[7]*(mean1-mean0);
  Double_t sigma2 = sqrt( pow(sigma0, 2) +2*pow(sigma_pho, 2) );

  Double_t N3 = par[8];
  Double_t mean3= par[9]*(mean1-mean0);
  Double_t sigma3 = sqrt( pow(sigma0, 2) +3*pow(sigma_pho, 2) );


  Double_t c0= par[10];
  Double_t c1= par[11];
  Double_t c2= par[12];
  //  Double_t C3= par[13];

  Double_t xx = x[0];
  Double_t value =0;
  //  value = c0+c1*x+c2*x*x;
  value = c0*TMath::Landau(xx, c1, c2);

  value += N0 * TMath::Gaus(xx, mean0, sigma0);
  value += N1 * TMath::Gaus(xx, mean1, sigma1);
  value += N2 * TMath::Gaus(xx, mean2, sigma2);
  value += N3 * TMath::Gaus(xx, mean3, sigma3);


  return value;

}


void testfit(const char *filename, const char *histname="") {


  setTDRStyle();
  TCanvas *c1 = new TCanvas("c1", "");


  TFile *file = new TFile(filename, "READ");

  TH1F *hist = (TH1F *) file->Get(histname);

  TF1 *myfunc = new TF1("myfunc", signal_shape, -0.5, 3.5, 13);

  /*
  Double_t parInits[]={ 
			2301, -0.02, 0.1,
			238, 1.01, 0.15,
			100, 2.0, 
			10, 3.0, 
			67, -0.1, 1};
  */
  Double_t parInits[]={ 
			230.1, -0.02, 0.1,
			23.8, 1.01, 0.15,
			10.0, 2.0, 
			10, 3.0, 
			67, -0.1, 1};
  myfunc->SetParameters(parInits); 
  myfunc->SetParName(0, "N_{0}"); 
  myfunc->SetParName(1, "m_{0}"); 
  myfunc->SetParName(2, "#sigma_{ele}"); 
  myfunc->SetParName(3, "N_{1}"); 
  myfunc->SetParName(4, "m_{1}"); 
  myfunc->SetParName(5, "#sigma_{p.e.}"); 
  myfunc->SetParName(6, "N_{2}"); 
  myfunc->SetParName(7, "m_{2}"); 
  myfunc->SetParName(8, "N_{3}"); 
  myfunc->SetParName(9, "m_{3}"); 

  myfunc->SetParName(10, "N_{bkg}"); 
  myfunc->SetParName(11, "c1_{bkg}"); 
  myfunc->SetParName(12, "c2_{bkg}"); 

  /*
  // nn setting
  //myfunc->FixParameter(0, parInits[0]);
  //myfunc->FixParameter(1, parInits[1]);
  //myfunc->FixParameter(2, parInits[2]);
  //myfunc->FixParameter(3, parInits[3]);
  myfunc->FixParameter(4, parInits[4]);
  myfunc->FixParameter(5, parInits[5]);
  myfunc->FixParameter(6, 0);
  myfunc->FixParameter(7, 0);
  myfunc->FixParameter(8, 0); 
  myfunc->FixParameter(9, 0);
  */

  /*
  myfunc->FixParameter(0, parInits[0]);
  myfunc->FixParameter(1, parInits[1]);
  myfunc->FixParameter(2, parInits[2]);
  myfunc->FixParameter(3, parInits[3]);
  myfunc->FixParameter(4, parInits[4]);
  myfunc->FixParameter(5, parInits[5]);
  myfunc->FixParameter(6, parInits[6]);
  myfunc->FixParameter(7, parInits[7]);
  myfunc->FixParameter(8, parInits[8]);

  
  */
  // sn setting
  //myfunc->FixParameter(8, 0);
  //  myfunc->FixParameter(9, parInits[9]);

  // the following 3 are set for S13360
  myfunc->SetParameter(5, 0.05);
  myfunc->FixParameter(6, 0);
  myfunc->FixParameter(7, 0);
  myfunc->FixParameter(8, 0);
  myfunc->FixParameter(9, 0);
    myfunc->FixParameter(10, 0);


    myfunc->FixParameter(11, 0.2);
   myfunc->FixParameter(12, 0.18);

   //   myfunc->FixParameter(12, 0.25);// default
  



  hist->Fit(myfunc, "", "", -1, 5);

  gStyle->SetOptFit(1);

  gPad->SetLogy();
}


void draw_func(void) {

  TF1 *f1 = new TF1("f1", "TMath::Poisson(x, [0])", 0, 50);
  f1->SetParameter(0, 20);

  f1->Draw();

}
