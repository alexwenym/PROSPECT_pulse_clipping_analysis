#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TAttMarker.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TRandom.h"  

#include <dirent.h>     //for locating directories
#include <math.h> 
#include <iostream>


// helper function to get max y value of TGraph
Double_t findMaxPosition(TGraph *G) {
    Double_t x, y;
    G->GetPoint(0, x, y);
    Double_t max_x = x, max_y = y;
    for(int i = 1; i < G->GetN(); i++)
    {
        G->GetPoint(i, x, y);
        if(y > max_y) {
           max_x = x;
           max_y = y;
        }
    }
    return max_y;
}


TGraph *g_copy;
double myfunc(double *xx, double *) {return g_copy->Eval(xx[0]);}

Double_t Get_Time(TGraph* g){

  //find the maximum y value 
  // This is erroneous, TGraph->GetHistogram does't return substantial content, it only draws the axis. 
  //double half_max = g->GetHistogram()->GetMaximum()/2.0;
  // this is a slower but correct way
  double half_max = findMaxPosition(g)/2.0; 
      
  // fill g_copy
  Int_t n; n = g->GetN(); //get plotted array dimention
  Double_t * ax = new Double_t[n];
  Double_t * ay = new Double_t[n];

  for(Int_t i = 0; i < n; i++) g->GetPoint(i,ax[i],ay[i]);
 
  double min = ax[0]; double max = ax[n-1];
  g_copy = new TGraph(n,ax,ay);
  delete [] ax;   delete [] ay; 

  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
  TF1 *f1 = new TF1("f1",myfunc,min,max);
  
  //find coresponding leading x value
  return f1->GetX(half_max, -10, 20); 
}



//as the waveform grows, the axis will shift
pair<double, double> find_zero(TGraph* g, double range[]){

  double new_zero = Get_Time(g);
  cout << "new zero: " << new_zero << "\n\n"<< endl;
  
  //shift axis
 return std::make_pair(range[0]-new_zero, range[1]-new_zero);

 
}

Double_t Get_Energy(TGraph* g){

  // the area of the total pulse is scaled to match an energy
  // the scaling is arbitrary (seen in calibration file)

  // fill g_copy
  Int_t n; n = g->GetN(); //get plotted array dimention
  Double_t * ax = new Double_t[n];
  Double_t * ay = new Double_t[n];

  for(Int_t i = 0; i < n; i++) g->GetPoint(i,ax[i],ay[i]);
  
  double min = ax[0]; double max = ax[n-1];
  g_copy = new TGraph(n,ax,ay);
  delete [] ax;   delete [] ay; 

  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
  TF1 *f1 = new TF1("f1",myfunc,min,max);
  return  (f1->Integral(min,max));
}


double Get_PSD(TGraph* g){
  
  // full integral from -12 to 100ns
  // tail from 44 to 100ns
  
  //fill g_copy
  Int_t n; n = g->GetN(); //get plotted array dimention
  Double_t * ax = new Double_t[n];
  Double_t * ay = new Double_t[n];

  for(Int_t i = 0; i < n; i++){
    g->GetPoint(i,ax[i],ay[i]);
  }
  
  g_copy = new TGraph(n,ax,ay);  
  ROOT::Math::IntegratorOneDimOptions::SetDefaultRelTolerance(1.E-6);
  TF1 *f1 = new TF1("f1",myfunc,ax[0],ax[n-1]);
  delete [] ax;   delete [] ay;

  double new_zero = Get_Time(g);
  double range[2] = {44+new_zero, 200+new_zero};
  double full_range[2] = {-12+new_zero, 200+new_zero};
 
  return (f1->Integral(range[0], range[1])) /
    (f1->Integral(full_range[0], full_range[1]));
}

void AddPlots(TString datafile, TString name, TH1* hist){

  TH1* h;
  TFile *f = TFile::Open(datafile);
  f->cd("MuonTracking");
  gDirectory->GetObject(name,h);
  hist->Add(h);
  f->Close();
}






void Clipped_Wave(){

  
  //template directory
  TString template_file = "./templates/AD_Templates.root";

  TFile *f = TFile::Open(template_file);
  TGraph* graph = (TGraph*)f->Get("gPSD_Template_0_0");
  //TGraph* graph = (TGraph*)f->Get("gPSD_Template_108_1");
  f->Close();

  //TH2F* PSDvE = new TH2F("PSD verses Energy", "", 200,0,200, 1000,0,1);
  //TString PSDvE_file = "/home/austin/Documents/test/data/simulated/combined/AD1_Wet_MuonTracking.root";
  //AddPlots(PSDvE_file, "hPSDvE", PSDvE);

  
  // from graph, get maximum and gradually scale it
  Int_t n; n = graph->GetN(); //get plotted array dimention
  Double_t max_y = findMaxPosition(graph);
  Double_t cutoff= max_y*1.5; // when clipping takes place
  Double_t dh = 0.0;

  // create a list of events to contrast as the wave scales by dh
  struct compare { Double_t dh, PSD, area, unclipped_area,time;};
  vector<compare> list;

  
  int counter = 0;
  while (dh < cutoff) {

    Double_t * ax = new Double_t[n];
    Double_t * ay = new Double_t[n];
    Double_t * ay_new = new Double_t[n];
    Double_t * ay_scaled = new Double_t[n];

    // scaling fraction
    double s = (max_y + dh) / max_y;
    for(Int_t i = 0; i < n; i++){
      graph->GetPoint(i,ax[i],ay[i]);
      
      ay_new[i] = ay[i]*s;
      ay_scaled[i] = ay[i]*s;
      if(ay_new[i] >= cutoff) ay_new[i] = cutoff;
    }

    TGraph* graph_new = new TGraph(n,ax,ay_new);
    TGraph* graph_scaled = new TGraph(n,ax,ay_scaled);

    delete [] ax;      delete [] ay;
    delete [] ay_new;  delete [] ay_scaled;

    counter++; 
    
    compare C;
    C.dh = dh; C.PSD = Get_PSD(graph_new); C.time = Get_Time(graph_new);
    C.area = Get_Energy(graph_new); C.unclipped_area = Get_Energy(graph_scaled);
    list.push_back(C);
    dh += 1.0;
  }


  
  Double_t * ah = new Double_t[(int)list.size()];
  Double_t * a_time = new Double_t[(int)list.size()];
  Double_t * aPSD = new Double_t[(int)list.size()];
  Double_t * aE = new Double_t[(int)list.size()];
  Double_t * aE_unclipped = new Double_t[(int)list.size()];

  cout << "__________________________"<< endl;
  cout << "\n\n"<< endl;

  //double offset = list[0].area / 2000.0;
  for(int i = 0; i < (int)list.size(); ++i){
    ah[i] = list[i].dh;
    aPSD[i] = list[i].PSD;
    aE[i] = list[i].area; /// 2000.0 -offset;
    a_time[i] = list[i].time;
    aE_unclipped[i] = list[i].unclipped_area;
  }

  TGraph* g_PSD_vs_dh = new TGraph((int)list.size(),ah,aPSD);
  TGraph* g_area_vs_dh = new TGraph((int)list.size(),ah,aE);
  TGraph* g_time_vs_dh = new TGraph((int)list.size(),ah,a_time);
  TGraph* g_PSD_vs_area = new TGraph((int)list.size(),aE,aPSD);
  TGraph* g_PSD_vs_time = new TGraph((int)list.size(), a_time, aPSD);
  //TGraph* Calculated_PSDvE = new TGraph((int)list.size, aE,aPSD);
  g_copy = new TGraph((int)list.size(), aE,aPSD);
  TF1 *f1 = new TF1("f1",myfunc,aE[0],aE[(int)list.size()-1]);
  
  delete [] ah; delete [] aPSD; delete [] a_time;
  delete [] aE; delete [] aE_unclipped;
  
  
 
  TCanvas *c1 = new TCanvas();
  graph->SetTitle("Waveform Template Example");
  graph->GetXaxis()->SetTitle("time [ns]");
  graph->SetMarkerColor(4);
  graph->Draw("");
  
  TCanvas *c2 = new TCanvas();
  g_PSD_vs_dh->GetXaxis()->SetTitle("#Deltah [arbitrary]");
  g_PSD_vs_dh->GetYaxis()->SetTitle("PSD");
  g_PSD_vs_dh->SetTitle("PSD vs Change in Height");
  g_PSD_vs_dh->Draw();

  TCanvas *c3 = new TCanvas();
  g_area_vs_dh->GetXaxis()->SetTitle("#Deltah [arbitrary]");
  g_area_vs_dh->GetYaxis()->SetTitle("Area");
  g_area_vs_dh->SetTitle("Area vs Change in Height");
  g_area_vs_dh->Draw();

  TCanvas *c4 = new TCanvas();
  g_PSD_vs_area->GetXaxis()->SetTitle("Area");
  g_PSD_vs_area->GetYaxis()->SetTitle("PSD");
  g_PSD_vs_area->SetTitle("PSD vs Clipped Waveform Area");
  g_PSD_vs_area->Draw();

  TCanvas *c5 = new TCanvas();
  g_time_vs_dh->GetXaxis()->SetTitle("#Deltah [arbitrary]");
  g_time_vs_dh->GetYaxis()->SetTitle("time");
  g_time_vs_dh->SetTitle("time vs Change in Height");
  g_time_vs_dh->Draw();

  TCanvas *c6 = new TCanvas();
  g_PSD_vs_time->GetXaxis()->SetTitle("#Deltat [ns]");
  g_PSD_vs_time->GetYaxis()->SetTitle("PSD");
  g_PSD_vs_time->SetTitle("time vs PSD");
  g_PSD_vs_time->Draw();
  
//   gStyle->SetOptStat(0);
//   TCanvas *c7 = new TCanvas();
//   PSDvE->SetTitle("PSD vs Energy");
//   PSDvE->GetYaxis()->SetTitle("PSD");
//   PSDvE->GetXaxis()->SetTitle("Enegy");
//   PSDvE->Draw();
//   g_PSD_vs_area->Draw("SAME");
 
  
  TCanvas *c8 = new TCanvas();
  f1->Draw();
 
}
