//===========Program to plot TDC and ADC data from root trees by Davison, Hesse, Losasda, Kerver, Porter, Bonilla, Horn, Thompson, Watkins, Brash==========
//edited by Laurence Carlucci
#include <iostream>
#include <fstream>
#include <TRandom.h>

//Global Variables
bool remakepedfile= false;
const int bin = 100; 
const Int_t Nadc = 16;
const Int_t Ntdc = 16;
const Int_t pedrun = 172;
const Int_t bl=Ntdc-16+0;//Start at channel 0
const Int_t br=Ntdc-16+1;
const Int_t tl=Ntdc-16+2;
const Int_t tr=Ntdc-16+3;
const Double_t adjadcto=1500.0;//value to ADJust ADC TO
Double_t t_fullscale = 140.0E-09; // full scale TDC range in seconds
Double_t t_convert=t_fullscale/4096.0;

const Int_t nthetabins = 51;
const Double_t thetalow = -100.5;
const Double_t thetahigh = 100.5;
const int num = 25;
//=============================================Method starts here for plotting===================================

void cosmic_chi2(Int_t nrun) {

  Double_t vn, resolution, granularity, xpos_range, dscint;  
  int xposbin;
        Double_t nscint_max=1.6;  //upper bound for index of refraction
	Double_t nscint_min=1.45; //lower bound for index of refraction
	Double_t chi2_temp, chi2_lo, chi2_prev, nscint_temp, nscint_lo, mult;
	chi2_lo = -1; //uninitialized chi2 value to determine
	nscint_temp = nscint_lo = nscint_min;
	mult = 4.0;
	Double_t nscintvals[num][2];	
	Double_t chi2vals[num][2];
	TH1F *htdcraw[Ntdc], *htdcadjusted[Ntdc], *hadcraw[Nadc], *hadccut[Nadc], *hadcadjusted[Nadc]; //Initialize histograms for drawing after loop for finding lowest chi2
	TH1F *htpos, *hbpos, *htheta, *htheta2;	

	//TCanvas *c1;

	for (int j = 0; j<num; j++){ 
	  cout<<"iteration "<<j<<endl; //test loop
	  vn = 2.997E08/nscint_temp;
	  resolution = 0.0232*nscint_temp*nscint_temp-0.1061*nscint_temp+0.1617;
	  granularity = t_convert*vn/2.0;
	  xpos_range = 0.30;
	  dscint = 0.105; // distance between scintillators in metres
	  xposbin = 2.0*xpos_range/granularity;

	  TRandom r;
	  Double_t rnd;

	  //create new pedestal values
	  Char_t tdcnames[][Ntdc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};	
	  Char_t adcnames[][Nadc]={"Bottom Left","Bottom Right","Top Left","Top Right","4","5","6","7","8","9","10","11","12","13","14","15"};

	  //Set correction values
	  Double_t tdccorrect[Ntdc];
	  Double_t ped[Nadc];
	  Double_t gain[Nadc];
	
	  //Read files and tree branches 
	  gStyle->SetOptStat(1);
	  TFile *froot =  new TFile(Form("./rootfiles/test%d.root",nrun));
	  TTree *troot = (TTree*)froot->Get("tdata");
	  Int_t tdc[Ntdc]; 
	  Int_t adc[Nadc]; 
	  troot->SetBranchAddress("tdc",&tdc);
	  troot->SetBranchAddress("adc",&adc);
	  const int nevents_in_file = (int)troot->GetEntries();

	  //cout << "Opened rootfile and read in " << nevents_in_file << " events." << endl;

	  //Fill Histograms
	  htpos = new TH1F("htpos","Top Position",xposbin,-1.0*xpos_range,xpos_range); //Histogram for top scintillator
	  hbpos = new TH1F("hbpos","Bottom Position",xposbin,-1.0*xpos_range,xpos_range); //Histogram for bottom scintillator 
	  //htheta = new TH1F("htheta","Angle (Degrees)",10*bin,-100.5,100.5); //Histogram for incidence angle
	  htheta = new TH1F("htheta","Angle (Degrees)",nthetabins,thetalow,thetahigh); //Histogram for incidence angle
	  htheta2 = new TH1F("htheta2","Angle (Degrees)",nthetabins,thetalow,thetahigh); //Histogram for simulated incidence angle
	  //cout << "Defined first set of histograms ... " << endl;

	  for ( int i = 0; i < Ntdc ; i++) {
	    htdcraw[i] = new TH1F(Form("htdcraw%02d", i),Form("%s   raw tdc",tdcnames[i]),bin,1700,2100);
	    htdcadjusted[i] = new TH1F(Form("htdcadjusted%02d", i),Form("%s   adjusted tdc",tdcnames[i]),bin,1500,2500);
	  }
    
	  for ( int i = 0; i < Nadc ; i++) {
	    hadcraw[i] = new TH1F(Form("hadcraw%02d", i),Form("%s   raw adc",adcnames[i]),bin,0,5500);
	    hadccut[i] = new TH1F(Form("hadccut%02d", i),Form("%s  cut adc",adcnames[i]),bin,0,5500);
	    hadcadjusted[i] = new TH1F(Form("hadcadjusted%02d", i),Form("%s   adjusted adc",adcnames[i]),bin,0,5500);
	  }
    
	  //cout << "Defined second set of histograms ... " << endl;
	
	  //=====================================GET PED for ADC==========================
	  if(remakepedfile){cout << "Executing luterpedestals.C" << endl; gROOT->ProcessLine(Form(".x luterpedestals.C(%d)",pedrun));cout<<"pedestal file made"<<endl;}
	  //cout << "Opening pedestal values file ... " << endl;
	  FILE *adcpeds = fopen(Form("./pedestalfiles/pedestalrun%d.dat",pedrun),"r");
	  //cout << "Filling ADC pedestal array ..." << endl;
	  for( int i = 0; i < Nadc ; i++) {//Start ADC filling loop
	    fscanf(adcpeds,"%lf\n",&ped[i]);
	    //printf("%lf\n",ped[i]);
	  }

	  //cout << "Read in pedestal file ... " << endl;

	  //=====================================GET ADJ FACTORS for ADC&TDC==========================
	  //cout << "Starting adc and tdc calibration loop ... " << endl;
	  for (int ie = 0; ie < nevents_in_file; ie++) {
	    troot->GetEntry(ie);
	    if(ie%1000==0.0&&ie!=0){cout<<"Progress: "<<((double)ie)/((double)nevents_in_file)*100<<"%"<<endl;
	    }
	    for( int i = 0; i < Ntdc ; i++) {//Start TDC filling loop
	      if(adc[i] < 3600){
		htdcraw[i]->Fill(tdc[i]);
	      }
	    }
	    for( int i = 0; i < Nadc ; i++) {//Start ADC-ADJUSTED filling loop
	      if (adc[i] < 3600){
		hadcraw[i]->Fill(adc[i]);
		hadccut[i]->Fill(adc[i]-ped[i]);
	      }
	    }//End ADC for loop
	  }

	  //c1 = new TCanvas("c1","tdc raw fitter",75,75,600,600);
	  //c1->Divide(2,2);
	  for ( int i = bl; i <= tr ; i++) {
	    //c1->cd(i-bl+1);		
	    //htdcraw[i]->Draw();
	    //the next 7 lines allow us to automatically set the fitting range
	    int highbin=5;		
	    for(int ii=0;ii<=bin;ii++){
	      if((htdcraw[i]->GetBinContent(ii))>(htdcraw[i]->GetBinContent(highbin))&&ii>=5){
		highbin=ii;
	      }
	    }
	    int median = htdcraw[i]->GetBinCenter(highbin);
	    htdcraw[i]->Fit("gaus","Q","",(median-200),(median+200));
	    tdccorrect[i] = 2000 - (htdcraw[i]->GetFunction("gaus")->GetParameter(1));
	  }  
	
	//c1->~TCanvas(); 
	//closes default canvas built within histo_>Fit

	  //==============LOOP THROUGH EVENTS IN FILE AND DO MATHS=============================
	  for (int ie = 0; ie < nevents_in_file; ie++) {
	    troot->GetEntry(ie);  
	    //=========================================ADC===============================
	    if(ie%1000==0.0&&ie!=0){cout<<"Progress: "<<((double)ie)/((double)nevents_in_file)*100<<"%"<<endl;
	    }
	    for( int i = 0; i < Nadc ; i++) {//Start ADC-ADJUSTED filling loop
	      gain[i]=adjadcto/(hadccut[i]->GetMean());
	      if(gain[i]*(adc[i]-ped[i])<4000){	 
		hadcadjusted[i]->Fill(gain[i]*(adc[i]-ped[i]));
	      }
	    }//End ADC for loop
		    
	    Double_t ebl = gain[bl]*(adc[bl]-ped[bl]);
	    Double_t ebr = gain[br]*(adc[br]-ped[br]);
	    Double_t etl = gain[tl]*(adc[tl]-ped[tl]);
	    Double_t etr = gain[tr]*(adc[tr]-ped[tr]);
		
	    //==========================================TDC===============================
		    
	    for( int i = 0; i < Ntdc ; i++) {//Start TDC filling loop
	      if(i<4){
		htdcadjusted[i]->Fill(tdc[i]+tdccorrect[i]);		
	      }
	      else{	
		htdcadjusted[i]->Fill(tdc[i]);
	      }
	    }//End TDC for loop
			
	    Bool_t good_bl = (abs(tdc[bl]+tdccorrect[bl]-2000)<200.0);
	    Bool_t good_br = (abs(tdc[br]+tdccorrect[br]-2000)<200.0);
	    Bool_t good_tl = (abs(tdc[tl]+tdccorrect[tl]-2000)<200.0);
	    Bool_t good_tr = (abs(tdc[tr]+tdccorrect[tr]-2000)<200.0);
	    Bool_t good_event = good_bl&&good_br&&good_tl&&good_tr;

	    if (good_event){
	      Double_t tbl = tdc[bl]+tdccorrect[bl];
	      Double_t tbr = tdc[br]+tdccorrect[br];
	      Double_t ttl = tdc[tl]+tdccorrect[tl];
	      Double_t ttr = tdc[tr]+tdccorrect[tr];
	      Double_t rtod = 180.0/3.14159265;
	      Double_t tdiff = (ttl-ttr)/2.0;
	      Double_t bdiff = (tbl-tbr)/2.0;
	      Double_t xmeantime = ((ttl+ttr)/2.0-(tbl+tbr)/2.0)*t_convert*vn;	
	      Double_t xtop = tdiff*t_convert*vn;
	      Double_t xbottom = bdiff*t_convert*vn;
	      rnd = r.Gaus(0.0,1.5);
	      Double_t rnd_top_pos = r.Gaus(0.0,resolution);
	      Double_t rnd_bottom_pos = r.Gaus(0.0,resolution);
	      Double_t rndxt = r.Uniform(-0.10,0.10);
	      Double_t rndyt = r.Uniform(-0.15,0.15);
	      Double_t rndrt = sqrt(rndxt*rndxt+rndyt*rndyt)+rnd_top_pos;
	      Double_t rndxb = r.Uniform(-0.10,0.10);
	      Double_t rndyb = r.Uniform(-0.15,0.15);
	      Double_t rndrb = sqrt(rndxb*rndxb+rndyb*rndyb)+rnd_bottom_pos;
	      Double_t theta = rtod*atan((xbottom-xtop)/dscint)+rnd;
	      Double_t theta2 = rtod*atan((rndrb-rndrt)/dscint);

	      if(xtop > -0.4 && xtop < 0.4 && xbottom > -0.4 && xbottom < 0.4) {
		htpos->Fill(xtop);
		hbpos->Fill(xbottom);
		if(abs(theta)<=85.0) htheta->Fill(theta);
		if(abs(theta2)<=85.0) htheta2->Fill(theta2);
	      }
	    }

	  } //End of loop over events
	  //const Int_t nxbins = htheta->GetNbinsX();
	  //Double_t res[nxbins],x[nxbins];
	  chi2_temp = htheta->Chi2Test(htheta2,"UW P CHI2");
	  cout << chi2_temp<< endl;
	  if (chi2_lo < 0){		//initially setting chi2_lo to the first chi2 value calculated (chi2_temp)	
	    chi2_lo = chi2_temp;
	  }	
	  else if (chi2_temp < chi2_lo){	//if a lower chi2 value is found, set chi2_lo & nscint_lo to the current values (chi2_temp & nscint_temp)
	    chi2_lo = chi2_temp;
	    nscint_lo = nscint_temp;
	  }
	  else if (chi2_temp > chi2_prev) {				//reverse the sign of the multiplier if the current chi2 value is greater than chi2_lo
	    mult*=-1;
	  }				//increment loop
	chi2_prev = chi2_temp;
	nscintvals[j][0] = nscint_temp;
	nscintvals[j][1] = nscint_lo;  
	chi2vals[j][0] = chi2_temp;
	chi2vals[j][1] = chi2_lo;
	nscint_temp+=(nscint_min+nscint_max)/(mult*(j+1));	//add/subtract [(nscint_min+nscint_max)/(4*j)] to/from  nscint_temp 
	}	
	//print run vals
	for (int j = 0; j<num; j++){
		cout << "nscint_temp: " << nscintvals[j][0] << " nscint_lo: " << nscintvals[j][1] << " chi2_temp: "<<chi2vals[j][0] << " chi2_lo: " <<chi2vals[j][1] << endl; 
	}
//=======================================Canvases=============================
//Draw X Positions Canvas
    
 	TCanvas *tdcpos = new TCanvas("tdcpos","X Positions",200,200,600,600);
 
	tdcpos->Divide(2,3);
	gStyle->SetOptFit(1);
	tdcpos->cd(1);
	//htpos->Draw();
	tdcpos->cd(2);
	//hbpos->Draw();
	tdcpos->cd(3);
	//htheta->Draw();
	htheta2->SetLineColor(kRed);
	//htheta2->Draw("SAME");
	tdcpos->cd(4);
	//htheta2->Draw();

        for (unsigned int jj = 0; jj < htheta->GetNbinsX();jj++) {
            htheta->SetBinError(jj,sqrt(htheta->GetBinContent(jj)));
        }
        for (unsigned int jj = 0; jj < htheta2->GetNbinsX();jj++) {
            htheta2->SetBinError(jj,sqrt(htheta2->GetBinContent(jj)));
        }

        const Int_t nxbins = htheta->GetNbinsX();
        Double_t res[nxbins],x[nxbins];
        htheta->Chi2Test(htheta2,"UW P",res);
        Double_t thetabinsize = (thetahigh-thetalow)/nthetabins;
        for (Int_t iii=0; iii<nxbins; iii++) x[iii]=thetalow+iii*thetabinsize;
        TGraph *resgr = new TGraph(nxbins,x,res);
        resgr->GetXaxis()->SetRangeUser(thetalow,thetahigh);
        resgr->SetMarkerStyle(22);
        resgr->SetMarkerColor(2);
        resgr->SetMarkerSize(.2);
	tdcpos->cd(5);
        //resgr->Draw("APL");

}
