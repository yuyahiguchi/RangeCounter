//#define GAUSS
//#define SCALE

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <string.h>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TColor.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TColor.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "../include/ParticleID.h"
#include "../include/functions.h"
#include "../include/DetectorClass.h"
#include "../include/VirtualDetectorClass.h"
#include "../include/NewParticleClass.h"
#include "../include/BeamLossClass.h"


using namespace std;
using std::filesystem::directory_iterator;


int main(int argc, char** argv){
  int select_hour=1;
  int select_deg=5;
  string absorberlength="5.0";
  string material="Cu";
  bool TestMode = false;
  bool NewVirtual = false;
  string Path = "";
  for(int i=1; i<argc; ++i){
    string sargv(argv[i]);
    if((sargv == "-a" || sargv == "--absorberlength") && argc>i){
      absorberlength=argv[i+1];
      Path = "../../data/Simulation/Absorber"+absorberlength;
    }
    else if((sargv == "-m" || sargv == "--material") && argc>i){
      material=argv[i+1];
    }
    else if((sargv == "-d" || sargv == "--degraderlength") && argc>i){
      TestMode = true;
      select_deg=stoi(argv[i+1]);
    }
    else if((sargv == "-n" || sargv == "--newvirtual") && argc>i){
      NewVirtual = true;
      Path = "../../data/Simulation/NewGeo";
    }
    else if((sargv == "-h" || sargv == "--help") && argc>i){
      cout << "-a : absorber length\n-m : material\n-d : degrader length. This option is used to analyze for a particular length.\n-n : use for additional virtual detector." << endl;
      return 0;
    }
  }
  
  if(Path == ""){
    cout << "ディレクトリパスが指定されていません" << endl;
    return 0;
  }
  
  gROOT -> SetBatch(true);
  gStyle->SetOptStat(11);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = {0.0,0.34,0.61,0.84,1.00};
  Double_t Red[NRGBs] = {0.0,0.0,0.87,1.0,0.51};
  Double_t Green[NRGBs] = {0.0,0.81,1.0,0.2,0.0};
  Double_t Blue[NRGBs] = {0.51,1.0,0.12,0.0,0.0};
  TColor::CreateGradientColorTable(NRGBs,stops,Red,Green,Blue,NCont);
  gStyle -> SetNumberContours(NCont);
  gStyle -> SetPalette(1);

  
  map<int, float> F_scale;
  ifstream fin("ScalingFactor.dat");
  if(!fin.is_open()){
    cout << "ScalingFactor.dat does't exist!!" << endl;
    return EXIT_FAILURE;
  }
  int length;
  float f;
  while(fin >> length >> f){
    //cout << length << " " << f << endl;
    F_scale[length] = f;
  }
    
  map<int,int> degraderID;
  int degradercount=0;
  map <string,int> kindofparticle;
 
  vector<string> fname;
  for(const auto & file : directory_iterator(Path)){
    string buf=file.path();
    string strbuf = buf.substr(buf.find_last_of("/")+1);
    if(strncmp(strbuf.data(), ".",1)==0) continue;
    string ini_particle = strbuf.substr(strbuf.length()-13,3);
    int hour = stoi(strbuf.substr(18,1));
    int degraderlength = stoi(strbuf.substr(28,strbuf.length()-56));
    string sub_material = strbuf.substr(strbuf.length()-7,2);
    if(sub_material!=material) continue;
    //if(ini_particle!="mu-") continue;
    if(TestMode){
      if(degraderlength!=select_deg) continue;
    }
    if(degraderID.find(degraderlength)==degraderID.end()){
      degraderID[degraderlength]=degradercount;
      degradercount+=1;
    }
    kindofparticle[ini_particle] = 0;
    if(hour != select_hour) continue;
    fname.push_back(file.path());
    cout << file.path() << endl;
  }

  const int nofdeg = (int)fname.size()/(int)kindofparticle.size();

  int n;
  TFile *file;
  TFile *OUTFILE;

  map<string, int> particleNo;
  particleNo["mu-"] = 0;
  particleNo["mu+"] = 1;
  particleNo["pi+"] = 2;
  particleNo["pi-"] = 3;
  
  map<string, int> ini_pid;
  ini_pid["mu+"]=ParticleID::anti_muon;
  ini_pid["mu-"]=ParticleID::muon;
  ini_pid["pi+"]=ParticleID::positive_pion;
  ini_pid["pi-"]=ParticleID::negative_pion;
    
  TH1D *momentum_distribution_before_deg[4*nofdeg];
  TH1D *stopped_absorber_momentum[4*nofdeg];
  TH1D *MomDist_afterdeg[4*nofdeg];
  for(auto itr=degraderID.begin(); itr!=degraderID.end(); itr++){
    n=itr->second;
    for(int i=0; i<4; i++){
      momentum_distribution_before_deg[4*n+i] = new TH1D(Form("momentum_distribution_before_deg%d",4*n+i),"momentum distribution;momentum[MeV/c];entry/[MeV/c]",300,0,300);
      momentum_distribution_before_deg[4*n+i] -> SetLineColor(2*(i+1));
      momentum_distribution_before_deg[4*n+i] -> SetMinimum(0.1);
      momentum_distribution_before_deg[4*n+i] -> SetMaximum(1e6);
      stopped_absorber_momentum[4*n+i] = new TH1D(Form("stopped_absorber_momentum%d",4*n+i),Form("initial momentum stopped at absorber ( degrader %dmm);momentum[MeV/c];entry/[MeV/c]",itr->first),300,0,300);
      stopped_absorber_momentum[4*n+i] -> SetLineColor(2*(i+1));
      stopped_absorber_momentum[4*n+i] -> SetMinimum(0.1);
      MomDist_afterdeg[4*n+i] = new TH1D(Form("MomDist_afterdeg_%d",4*n+i),Form("momentum distribution after deg(Absorber %smm,degrader %dmm);[MeV/c];entry/[MeV/c]",absorberlength.data(),itr->first),300,0,300);
      MomDist_afterdeg[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_afterdeg[4*n+i] -> SetMinimum(0.1);
    }
  }
  TH1D *decay_spectrum[nofdeg];
  TH1D *PionStopMomentumLayer3[nofdeg];
  for(auto itr=degraderID.begin(); itr!=degraderID.end(); itr++){
    n=itr->second;
    decay_spectrum[n] = new TH1D(Form("decay_spectrum%d",n),Form("decay spectrum in %s(#mu+,#mu-,#pi+, degrader %dmm);time[ns];entry",material.data(),itr->first),5000,0,5000);
    decay_spectrum[n] -> SetMinimum(0.1);
    PionStopMomentumLayer3[n] = new TH1D(Form("PionStopMomentumLayer3%d",n),Form("initial #pi+ momentum stopped at layer3 (degrader %dmm, absorber %s)",itr->first,absorberlength.data()),300,0,300);
    PionStopMomentumLayer3[n] -> SetMinimum(0.1);
  }
  TH1D *layer1_and_not2[4*nofdeg];
  for(auto ITR=degraderID.begin();ITR!=degraderID.end();ITR++){
    for(auto itr=particleNo.begin();itr!=particleNo.end();itr++){
      layer1_and_not2[4*ITR->second+itr->second] = new TH1D(Form("layer1_and_not2_%d",4*ITR->second+itr->second),Form("stopped at absorber(#%s, degarder %dmm);time[ns];entry",itr->first.data(),ITR->first),5000,0,5000);
      layer1_and_not2[4*ITR->second+itr->second] -> SetMinimum(0.1);
    }
  }
    
  bool pene_event;
  bool IsFill;
  
  TString Name;
  Name=Form("DegradeAna_%dhour_absorber%smm_%s.pdf",select_hour,absorberlength.data(),material.data());
  
  vector<float> remove;
  float P0[nofdeg];//muon
  float P0_error[nofdeg];  
  float P1[nofdeg];//anti muon
  float P1_error[nofdeg];
  float P2[nofdeg];//positive pion
  float P2_error[nofdeg];
  float length_array[nofdeg];
  float len_error[nofdeg];
  float length_array_posimu[nofdeg];
  float len_error_posimu[nofdeg];
  float length_array_posipi[nofdeg];
  float len_error_posipi[nofdeg];
  for(int i=0; i<nofdeg; ++i){
    P0[i]=0;
    P0_error[i]=0;
    P1[i]=0;
    P1_error[i]=0;
    length_array[i]=0;
    len_error[i]=0;
    length_array_posimu[i]=0;
    len_error_posimu[i]=0;
  }


  //**********************************************************************************************
  //**********************************************************************************************
  for(int No=0; No<fname.size(); ++No){
    int numofPion=0;
    string strbuf = fname[No].substr(fname[No].find_last_of("/")+1);
    string ini_particle = strbuf.substr(strbuf.length()-13,3);
    int hour = stoi(strbuf.substr(18,1));
    int degraderlength = stoi(strbuf.substr(28,strbuf.length()-56));
    file = new TFile(fname[No].data());
    
    //###############################################
    // NewParticleNTuple、BeamLossNTupleのtreeをセット
    //###############################################
    NewParticle *NP = new NewParticle(file,"NP");
    BeamLoss *blnt = new BeamLoss(file,"blnt");
    
    
    //#############################        
    // VirtualDetectorのtreeをセット
    //#############################        
    VirtualDetector *Vdet1 = new VirtualDetector(file,"bm0");
    VirtualDetector *Vdet2 = new VirtualDetector(file,"bm1");
    VirtualDetector *Vdet3 = new VirtualDetector(file,"bm2");
    ULong64_t totalEvent = Vdet1 -> GetTotalEvent();
    
    // <---必要があれば運動量の絶対値を求める式を書く
    
    //#######################     
    // Detectorのtreeをセット
    //#######################  
    Detector *Layer1 = new Detector(file,"rc0");
    Detector *Layer2 = new Detector(file,"rc1");
    Detector *Layer3 = new Detector(file,"rc2");
    
    cout << "totalEvent=" << totalEvent << endl;

    //###########
    // Event loop
    //###########
    for(ULong64_t ev=0; ev<(ULong64_t)totalEvent; ++ev){
      if(ev%100000==0)
	cout << "\r"  << fname[No].substr(fname[No].find_last_of("/")+1) << " : " /*<< fixed << setprecision(2) */<< (100*ev)/(ULong64_t)totalEvent << "%" << flush;
      float t1 = 0;
      float t2 = 0;

      //######################
      //NewParticle & BeamLoss
      //######################
      NP -> SetEventInfo(ev);
      float ini_momentum = NP -> SetInitialInfo() -> Get_Pz();
      momentum_distribution_before_deg[4*degraderID[degraderlength]+particleNo[ini_particle]] -> Fill(ini_momentum);
      
      //###############
      //VirtualDetector
      //###############


      //#######
      // layer1
      //#######

      Layer1 -> SetEventInfo(ev);
      if(Layer1->GetEntries() >= 2){
	continue;
      }

      if(Layer1->GetEntries() == 0 ){
	continue;
      }
      t1 = Layer1 -> Get_t(0);

      //#######
      // layer2
      //#######

      Layer2 -> SetEventInfo(ev);
      if(Layer2 -> GetEntries() != 0){
	t2 = Layer2 -> Get_t(0);
      }
      if(Layer2 -> GetEntries() >= 2){
	continue;
      }
      
      //#######
      // layer3
      //#######
      
      Layer3 -> SetEventInfo(ev);
      if(ini_particle=="pi+" && Layer3 -> GetEntries() >= 3){
	PionStopMomentumLayer3[degraderID[degraderlength]] -> Fill(ini_momentum);
      }
      if(Layer3 -> GetEntries() >= 3){
	numofPion++;
      }
      
      //////////////////////
      //崩壊スペクトルを作る
      /////////////////////
      if( t1!=0 && t2!=0 && t1<5 && (t2-t1>1.5) ){
	IsFill = false;
	for(size_t j=0; j<Layer2->GetEntries();++j){
	  if(Layer2->Get_t(j)>t1 && Layer2->Get_Edep(j)>0.05 /*&& !InContainer(remove,Layer2->GetNthTime(j))*/){
	    layer1_and_not2[4*degraderID[degraderlength]+particleNo[ini_particle]]->Fill(Layer2->Get_t(j)-t1);
	    decay_spectrum[degraderID[degraderlength]]->Fill(Layer2->Get_t(j)-t1);
	    IsFill = true;
	  }
	}
	if(IsFill)
	  stopped_absorber_momentum[4*degraderID[degraderlength]+particleNo[ini_particle]]->Fill(ini_momentum);
      }
      //////////////////////
      //イベントごとの処理
      //////////////////////
    }
    file->Close();
    cout << "\r"  << fname[No].substr(fname[No].find_last_of("/")+1) << " : " << 100 << "%"  << flush;
    cout << endl;
    P2[degraderID[degraderlength]] = numofPion;
    P2_error[degraderID[degraderlength]] = sqrt(numofPion);
    delete NP;
    delete blnt;
    delete Vdet1;
    delete Vdet2;
    delete Vdet3;
    delete Layer1;
    delete Layer2;
    delete Layer3;
  }

  TCanvas *c1 = new TCanvas(Name.Data(),Name.Data());
  c1 -> Print(Name+"[","pdf");
  OUTFILE = new TFile("DegraderAna.root","RECREATE");
  cout << "Open DegraderAna.root" << endl;
  OUTFILE->cd();
  c1->cd();
  c1 -> SetGridy();
  int page=1;
  TPaveStats *st1;
  TPaveStats *st2;
  TPaveStats *st3;
  TF1 *f1 = new TF1("f1","[0]*TMath::Gaus(x,[1],[2])");
  int maxEntry;
  float min_p[4];
  float max_p[4];
  float width[4];
  float mean_momentum[4];

  for(auto itr=degraderID.begin(); itr!=degraderID.end(); ++itr){
    n=itr->second;
    //#############
    //運動量分布
    //#############
    gPad->SetLogy();
    momentum_distribution_before_deg[4*n+0]->Draw();
    for(int i=1; i<3; i++){
      momentum_distribution_before_deg[4*n+i]->Draw("sames");
    }
    c1->Update();
    st1 = (TPaveStats*)momentum_distribution_before_deg[4*n+0]->FindObject("stats");
    st1->SetTextColor(2*(0+1));
    st2 = (TPaveStats*)momentum_distribution_before_deg[4*n+1]->FindObject("stats");
    st2->SetTextColor(2*(1+1));
    st3 = (TPaveStats*)momentum_distribution_before_deg[4*n+2]->FindObject("stats");
    st3->SetTextColor(2*(2+1));
    st1->SetX1NDC(0.78);
    st1->SetX2NDC(0.98);
    st1->SetY1NDC(0.83);
    st1->SetY2NDC(0.91);
    st2->SetX1NDC(0.78);
    st2->SetX2NDC(0.98);
    st2->SetY1NDC(0.72);
    st2->SetY2NDC(0.80);
    st3->SetX1NDC(0.78);
    st3->SetX2NDC(0.98);
    st3->SetY1NDC(0.61);
    st3->SetY2NDC(0.69);
    TLegend *lg1 = new TLegend(0.8,0.5,0.90,0.58);
    lg1->AddEntry(momentum_distribution_before_deg[4*n+0],"#mu-","l");
    lg1->AddEntry(momentum_distribution_before_deg[4*n+1],"#mu+","l");
    lg1->AddEntry(momentum_distribution_before_deg[4*n+2],"#pi+","l");
    lg1->Draw();
    c1->Print(Name,"pdf");
    c1->Write(Form("pageNo%d",page));
    page+=1;
    
    if(NewVirtual){
      MomDist_afterdeg[4*n+0] -> Draw();
      for(int i=1; i<3; ++i){
	MomDist_afterdeg[4*n+i] -> Draw("sames");
      }
      c1->Print(Name,"pdf");
      c1->Write(Form("pageNo%d",page));
      page+=1;
    }
      

#ifdef GAUSS
    //for mu-
    mean_momentum[0] = stopped_absorber_momentum[4*n+0] -> GetMaximumBin()-1;
    maxEntry = stopped_absorber_momentum[4*n+0]->GetBinContent(stopped_absorber_momentum[4*n+0]->GetMaximumBin());
    f1 -> SetParameters(maxEntry,mean_momentum[0],1);
    stopped_absorber_momentum[4*n+0] -> Fit(f1,"L","",mean_momentum[0]-20,mean_momentum[0]+20);
    width[0] = sqrt(f1 -> GetParameter(2));
    cout << "mu- : degrader " << itr->first << "mm width=" << width[0] << endl;
    //for mu+
    mean_momentum[1] = stopped_absorber_momentum[4*n+1] -> GetMaximumBin()-1;
    maxEntry = stopped_absorber_momentum[4*n+1]->GetBinContent(stopped_absorber_momentum[4*n+1]->GetMaximumBin());
    f1 -> SetParameters(maxEntry,mean_momentum[1],1);
    stopped_absorber_momentum[4*n+1] -> Fit("f1","L","",mean_momentum[1]-20,mean_momentum[1]+20);
    width[1] = sqrt(f1 -> GetParameter(2));
    cout << "mu+ : degrader " << itr->first << "mm width=" << width[1] << endl;
    //for pi+
    mean_momentum[2] = stopped_absorber_momentum[4*n+2] -> GetMaximumBin()-1;
    maxEntry = stopped_absorber_momentum[4*n+2]->GetBinContent(stopped_absorber_momentum[4*n+2]->GetMaximumBin());
    f1 -> SetParameters(maxEntry,mean_momentum[2],1);
    stopped_absorber_momentum[4*n+2] -> Fit("f1","L","",mean_momentum[2]-20,mean_momentum[2]+20);
    width[2] = sqrt(f1 -> GetParameter(2));
    cout << "pi+ : degrader " << itr->first << "mm width=" << width[2] << endl;
#else
    //for mu-
    bool th=false;
    maxEntry = stopped_absorber_momentum[4*n+0]->GetBinContent(stopped_absorber_momentum[4*n+0]->GetMaximumBin());
    for(int t=0; t<250; ++t){
      if(!th && stopped_absorber_momentum[4*n+0]->GetBinContent(t+1)>maxEntry*0.5){
	min_p[0]=t;
	th=true;
      }
      if(th && stopped_absorber_momentum[4*n+0]->GetBinContent(t+1)<maxEntry*0.5){
	max_p[0]=t;
	th=false;
	break;
      }
    }
    width[0] = max_p[0]-min_p[0];
    mean_momentum[0] = (min_p[0]+max_p[0])/2;
    cout << "mu- : degrader " << itr->first << "mm width=" << width[0] << endl;
    
    //for mu+
    maxEntry = stopped_absorber_momentum[4*n+1]->GetBinContent(stopped_absorber_momentum[4*n+1]->GetMaximumBin());
    for(int t=0; t<250; ++t){
      if(!th && stopped_absorber_momentum[4*n+1]->GetBinContent(t+1)>maxEntry*0.5){
	min_p[1]=t;
	th=true;
      }
      if(th && stopped_absorber_momentum[4*n+1]->GetBinContent(t+1)<maxEntry*0.5){
	max_p[1]=t;
	th=false;
	break;
      }
    }
    width[1] = max_p[1]-min_p[1];
    mean_momentum[1] = (min_p[1]+max_p[1])/2;
    cout << "mu+ : degrader " << itr->first << "mm width=" << width[1] << endl;
    
    //for pi+
    maxEntry = stopped_absorber_momentum[4*n+2]->GetBinContent(stopped_absorber_momentum[4*n+2]->GetMaximumBin());
    for(int t=0; t<250; ++t){
      if(!th && stopped_absorber_momentum[4*n+2]->GetBinContent(t+1)>maxEntry*0.5){
	min_p[2]=t;
	th=true;
      }
      if(th && stopped_absorber_momentum[4*n+2]->GetBinContent(t+1)<maxEntry*0.5){
	max_p[2]=t;
	th=false;
	break;
      }
    }
    width[2] = max_p[2]-min_p[2];
    mean_momentum[2] = (min_p[2]+max_p[2])/2;
    cout << "pi+ : degrader " << itr->first << "mm width=" << width[2] << endl;
#endif

    
    gPad -> SetLogy(0);
    for(int i=0; i<3; ++i){
      stopped_absorber_momentum[4*n+i]->Draw();
      c1->Print(Name,"pdf");
      c1->Write(Form("pageNo%d",page));
      page+=1;
    }
    /* for(int i=1; i<3; i++){ */
    /*   stopped_absorber_momentum[4*n+i]->Draw("sames"); */
    /* } */
    /* c1->Update(); */
    /* st1 = (TPaveStats*)stopped_absorber_momentum[4*n+0]->FindObject("stats"); */
    /* st1->SetTextColor(2*(0+1)); */
    /* st2 = (TPaveStats*)stopped_absorber_momentum[4*n+1]->FindObject("stats"); */
    /* st2->SetTextColor(2*(1+1)); */
    /* st3 = (TPaveStats*)stopped_absorber_momentum[4*n+2]->FindObject("stats"); */
    /* st3->SetTextColor(2*(2+1)); */
    /* st1->SetX1NDC(0.78); */
    /* st1->SetX2NDC(0.98); */
    /* st1->SetY1NDC(0.83); */
    /* st1->SetY2NDC(0.91); */
    /* st2->SetX1NDC(0.78); */
    /* st2->SetX2NDC(0.98); */
    /* st2->SetY1NDC(0.72); */
    /* st2->SetY2NDC(0.80); */
    /* st3->SetX1NDC(0.78); */
    /* st3->SetX2NDC(0.98); */
    /* st3->SetY1NDC(0.61); */
    /* st3->SetY2NDC(0.69); */
    /* TLegend *l1 = new TLegend(0.8,0.5,0.90,0.58); */
    /* l1->AddEntry(stopped_absorber_momentum[4*n+0],"#mu-","l"); */
    /* l1->AddEntry(stopped_absorber_momentum[4*n+1],"#mu+","l"); */
    /* l1->AddEntry(stopped_absorber_momentum[4*n+2],"#pi+","l"); */
    /* l1->Draw(); */
    /* c1->Print(Name,"pdf"); */
    /* c1->Write(Form("pageNo%d",page)); */
    /* page+=1; */
    /* delete l1; */

    PionStopMomentumLayer3[n] -> Draw();
    c1->Print(Name,"pdf");
    c1->Write(Form("pageNo%d",page));
    page+=1;
    

    //#############
    //時間分布
    //#############

    //#############
    //崩壊スペクト
    //#############
    c1->Clear();
    gPad->SetLogy();
    for(int i=0; i<3; i++){
      layer1_and_not2[4*n+i]->Draw();
      c1->Print(Name,"pdf");
      c1->Write(Form("pageNo%d",page));
      page+=1;
    }
    auto[p0,p0_error,p1,p1_error,LT] = GetFit(decay_spectrum[n],material);
    cout << "mu- : " << itr->first << " " << itr->second << " " << mean_momentum[0] << " " << width[0] << " " << p1 << " " << p1_error << endl;
    cout << "mu+ : " << itr->first << " " << itr->second << " " << mean_momentum[1] << " " << width[1] << " " << p0 << " " << p0_error << endl;
    P0[n] = p0*2200/width[0];//LT=2200;
    P0_error[n] = p0_error*2200/width[0];
    if(material=="Cu"){
#ifdef SCALE
      P1[n] = p1*LT*F_scale[itr->first]/width[0];//LT=163;
#else
      P1[n] = p1*LT/width[0];//LT=163;
#endif
      P1_error[n] = p1_error*LT/width[0];
    }
    else if(material=="Al"){
      P1[n] = p1*859/width[0];
      P1_error[n] = p1_error*859/width[0];
    }
    length_array[n] = mean_momentum[0];;
    len_error[n] = width[0];;
    length_array_posimu[n] = mean_momentum[1];
    len_error_posimu[n] = mean_momentum[1];
    length_array_posipi[n] = width[1];
    len_error_posipi[n] = mean_momentum[2];
    decay_spectrum[n]->Draw();
    c1->Print(Name,"pdf");
    c1->Write(Form("pageNo%d",page));
    page+=1;
  }

  //運動量分布の図にreconstructした点を載せる
  gPad->SetLogy();
  c1 -> SetGridy();
  double decay_rate;
  if(material=="Al"){
    decay_rate=0.40;
  }
  else if(material=="Cu"){
    decay_rate=0.074;
  }
  double geo_acceptance=0.45;
  //double acceptance = decay_rate * geo_acceptance;
  double acceptance = 0.16;
  for(auto &x : P0){
    x = x/geo_acceptance;
  }
  for(auto &x : P0_error){
    x = x/geo_acceptance;
  }
  for(auto &x : P1){
    x = x/acceptance;
  }
  for(auto &x : P1_error){
    x = x/acceptance;
  }
  TGraphErrors *gr1 = new TGraphErrors(nofdeg,length_array,P1,len_error,P1_error);//mu-
  gr1->SetMarkerStyle(20);
  gr1->SetMarkerSize(1);
  gr1->SetTitle("# of stopping mu-");
  gr1->GetXaxis()->SetTitle("momentum (MeV/c)");
  gr1->GetYaxis()->SetTitle("# of decay");
  TGraphErrors *gr2 = new TGraphErrors(nofdeg,length_array_posimu,P0,len_error_posimu,P0_error);//mu+
  gr2->SetMarkerStyle(20);
  gr2->SetMarkerSize(1);
  gr2->SetTitle("# of stopping mu+");
  gr2->GetXaxis()->SetTitle("momentum (MeV/c)");
  gr2->GetYaxis()->SetTitle("# of decay");
  TGraphErrors *gr3 = new TGraphErrors(nofdeg,length_array_posipi,P2,len_error_posipi,P2_error);//pi+
  
  momentum_distribution_before_deg[4*n+0]->SetMinimum(1e3);
  momentum_distribution_before_deg[4*n+0]->Draw();
  gr1->Draw("P");
  c1->Print(Name,"pdf");
  c1->Write(Form("pageNo%d",page));
  page+=1;

  momentum_distribution_before_deg[4*n+1]->SetMinimum(1e2);
  momentum_distribution_before_deg[4*n+1]->Draw();
  gr2->Draw("P");
  c1->Print(Name,"pdf");
  c1->Write(Form("pageNo%d",page));
  page+=1;

  momentum_distribution_before_deg[4*n+2]->SetMinimum(1);
  momentum_distribution_before_deg[4*n+2]->Draw();
  gr3->Draw("P");
  c1->Print(Name,"pdf");
  c1->Write(Form("pageNo%d",page));
  page+=1;

  c1->Print(Name+"]","pdf");
  OUTFILE->Close();
  cout << "Close DegraderAna.root" << endl;
  return 0;
}
