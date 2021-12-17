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
#include "TLine.h"
#include "TPaveStats.h"
#include "TColor.h"
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
  string absorberlength = "5.0";
  string material="Cu";
  bool TestMode = false;
  bool NewVirtual = false;
  string Path = "";

  
  for(int i=1; i<argc; ++i){
    string sargv(argv[i]);
    if((sargv == "-a" || sargv == "--absorberlength") && argc>i){
      absorberlength=argv[i+1];
      Path = "/Users/higuchiyuya/range-counter/higuchi/data/Simulation/Absorber"+absorberlength;
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
      Path = "/Users/higuchiyuya/range-counter/higuchi/data/Simulation/NewGeo";
      //Path = "/Users/higuchiyuya/range-counter/higuchi/data/Simulation/IncludeDecayInfo";
    }
    else if((sargv == "-h" || sargv == "--help") && argc>i){
      cout << "-a : absorber length\n-m : material\n-d : degrader length. This option is used to analyze for a particular length.\n-n : use for additional virtual detector." << endl;
      return 0;
    }
  }

  if(Path == ""){
    cout << "No directory path!!" << endl;
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
    //if(ini_particle!="mu+") continue;
    if(TestMode){
      if(degraderlength!=select_deg) continue;
    }
    if(degraderID.find(degraderlength)==degraderID.end()){
      degraderID[degraderlength]=degradercount++;
    }
    kindofparticle[ini_particle] = 0;
    if(hour != select_hour) continue;
    cout << file.path() << endl;
    fname.push_back(file.path());
  }

  const int nofParticle=(int)kindofparticle.size();
  const int nofdeg = (int)fname.size()/(int)kindofparticle.size();
  int n;
  TFile *file;
  
  map<string,int> particleNo;
  particleNo["mu-"]=0;
  particleNo["mu+"]=1;
  particleNo["pi+"]=2;
  particleNo["pi-"]=3;
  map<string,int> ini_pid;
  ini_pid["mu+"]=ParticleID::anti_muon;
  ini_pid["mu-"]=ParticleID::muon;
  ini_pid["pi+"]=ParticleID::positive_pion;
  ini_pid["pi-"]=ParticleID::negative_pion;

  
  TH1D *momentum_distribution[4*nofdeg];
  TH2D *momentum_direction[4*nofdeg];
  TH1D *MomDist_Vdet1[4*nofdeg];
  TH1D *MomDist_Vdet2[4*nofdeg];
  TH1D *MomDist_Vdet3[4*nofdeg];
  TH1D *MomDist_InDeg[4*nofdeg];
  TH1D *EdepWhenStopAbso[4*nofdeg];
  TH1D *MomDist_InAbso[4*nofdeg];
  TH1D *PID[4*nofdeg];
  TH1D *PionEdep[nofdeg];
  for(auto itr=degraderID.begin(); itr!=degraderID.end(); itr++){
    n=itr->second;
    PionEdep[n] = new TH1D(Form("PionEdep_%d",n),Form("pion edep in layer2 (Absoeber %smm, degrader %dmm);Edep[MeV];entry/[MeV]",absorberlength.data(),itr->first),200,0,20);
    for(auto it=particleNo.begin();it!=particleNo.end();it++){
      int i=it->second;
      momentum_distribution[4*n+i] = new TH1D(Form("momentum_distribution_%d",4*n+i),Form("momentum distribution (Absorber %smm);momentum [MeV/c];entry/[MeV/c]",absorberlength.data()),300,0,300);
      momentum_distribution[4*n+i] -> SetLineColor(2*(i+1));
      momentum_distribution[4*n+i] -> SetMinimum(0.1);
      
      momentum_direction[4*n+i] = new TH2D(Form("momentum_direction_%d",4*n+i),Form("momentum direction from absorber(Absorber %smm, #%s);pz/p;py/p;entry",absorberlength.data(),it->first.data()),40,-1,1,40,-1,1);
      momentum_direction[4*n+i] -> SetLineColor(2*(i+1));
      momentum_direction[4*n+i] -> SetMinimum(0.1);
      
      MomDist_Vdet1[4*n+i] = new TH1D(Form("MomDist_Vdet1_%d",4*n+i),Form("momentum distribution at Vdet1 (Absorber %smm, degrader %dmm);momentum [MeV/c];entry/[MeV/c]",absorberlength.data(),itr->first),300,0,300);
      MomDist_Vdet1[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_Vdet1[4*n+i] -> SetMinimum(0.1);

      MomDist_Vdet2[4*n+i] = new TH1D(Form("MomDist_Vdet2_%d",4*n+i),Form("momentum distribution at Vdet2 (Absorber %smm, degrader %dmm);momentum [MeV/c];entry/[MeV/c]",absorberlength.data(),itr->first),300,0,300);
      MomDist_Vdet2[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_Vdet2[4*n+i] -> SetMinimum(0.1);
      
      MomDist_Vdet3[4*n+i] = new TH1D(Form("MomDist_Vdet3_%d",4*n+i),Form("momentum distribution at Vdet3 (Absorber %smm, degrader %dmm);momentum [MeV/c];entry/[MeV/c]",absorberlength.data(),itr->first),300,0,300);
      MomDist_Vdet3[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_Vdet3[4*n+i] -> SetMinimum(0.1);
      
      MomDist_InDeg[4*n+i] = new TH1D(Form("MomDist_InDeg_%d",4*n+i),Form("momentum distribution in degrader (Absorber %smm, degrader %dmm);momentum [MeV/c];entry/[MeV/c]",absorberlength.data(),itr->first),300,0,300);
      MomDist_InDeg[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_InDeg[4*n+i] -> SetMinimum(0.1);

      EdepWhenStopAbso[4*n+i] = new TH1D(Form("VisibleEdepWhenStopAbso_%d",4*n+i),Form("Energy deposit when stopped at absorber (Absorber %smm, degrader %dmm, #%s);Edep[MeV];count [/MeV]",absorberlength.data(),itr->first,it->first.data()),600,0,30);
      MomDist_InAbso[4*n+i] = new TH1D(Form("MomDist_InAbso_%d",4*n+i),Form("Particle when stopped at absorber (Absorber %smm, degrader %dmm, #%s);momentum [MeV/c];entry / [MeV/c]",absorberlength.data(),itr->first,it->first.data()),300,0,300);
      MomDist_InAbso[4*n+i] -> SetLineColor(2*(i+1));
      MomDist_InAbso[4*n+i] -> SetMinimum(0.1);
      PID[4*n+i] = new TH1D(Form("PID_%d",4*n+i),Form("Particle identification (Absorber %smm, degrader %dmm, #%s);PID;entry",absorberlength.data(),itr->first,it->first.data()),200,-100,100);
      PID[4*n+i] -> SetLineColor(2*(i+1));
      PID[4*n+i] -> SetMinimum(0.1);
    }
  }

  
  TH1D *decay_spectrum[4*nofdeg];
  TH2D *DecayPosition[4*nofdeg];
  for(auto ITR=degraderID.begin();ITR!=degraderID.end();ITR++){
    int n = ITR->second;
    for(auto itr=particleNo.begin();itr!=particleNo.end();itr++){
      int m = itr->second;
      decay_spectrum[4*n+m] = new TH1D(Form("decay_spectrum_%d",4*n+m),Form("decay spectrum(#%s, degarder %dmm);time [ns];entry",itr->first.data(),ITR->first),5000,0,5000);
      decay_spectrum[4*n+m] -> SetMinimum(0.1);
      
      DecayPosition[4*n+m] = new TH2D(Form("decay_posi_%d",4*n+m),Form("decay position VS initial momentum (#%s,degrader %dmm);z [mm];initial momentum [MeV/c]",itr->first.data(),ITR->first),600,-100,200,300,0,300);
    }
  }
  
  
  TGraph *MaximumStopBin[3];
  TGraph *ratio[3];
  for(int i=0; i<3; ++i){
    MaximumStopBin[i] = new TGraph;
    MaximumStopBin[i] -> SetMarkerStyle(20);
    MaximumStopBin[i] -> SetMarkerSize(1);
    MaximumStopBin[i] -> SetMarkerColor(40);
    ratio[i] = new TGraph;
    ratio[i] -> SetMarkerStyle(20);
    ratio[i] -> SetMarkerSize(1);
    ratio[i] -> SetTitle("input VS stop momentum");
    ratio[i] -> GetXaxis() -> SetTitle("stop momentum [MeV]");
    ratio[i] -> GetYaxis() -> SetTitle("input momentum / stop momentum");
  }

  
  float degraderArray[3][nofdeg];
  float degraderArrayError[3][nofdeg];
  float nofStop[3][nofdeg];
  float nofStopError[3][nofdeg];
  float acceptance[3][nofdeg];
  float acceptanceError[3][nofdeg];
  for(int i=0;i<3;++i){
    for(int j=0;j<nofdeg;++j){
      degraderArray[i][j]=0;
      degraderArrayError[i][j]=0;
      nofStop[i][j]=0;
      nofStopError[i][j]=0;
      acceptance[i][j]=0;
      acceptanceError[i][j]=0;      
    }
  }
  
  vector<float> remove;

  //**********************************************************************************************
  // File loop
  //**********************************************************************************************

  ofstream ScaleFile;
  ScaleFile.open("ScalingFactor.dat",ios::out);
    
  for(int No=0, fsize=fname.size(); No<fsize; ++No){
    string strbuf = fname[No].substr(fname[No].find_last_of("/")+1);
    string ini_particle = strbuf.substr(strbuf.length()-13,3);
    int hour = stoi(strbuf.substr(18,1));
    int degraderlength = stoi(strbuf.substr(28,strbuf.length()-56));

    int ini_pdgid = ini_pid[ini_particle];
    int pdgNo = particleNo[ini_particle];
    int degNo = degraderID[degraderlength];
    
    file = new TFile(fname[No].data());


    //#############################    
    // Setting NewParticle and BeamLoss tree
    //#############################
    NewParticle *NP = new NewParticle(file,"NP");
    BeamLoss *blnt = new BeamLoss(file,"blnt");
    
    //#############################        
    // Setting VirtualDetector tree
    //#############################        

    VirtualDetector *Vdet1 = new VirtualDetector(file,"bm0");
    ULong64_t totalEvent = Vdet1 -> GetTotalEvent();
    VirtualDetector *Vdet2 = new VirtualDetector(file,"bm1");
    VirtualDetector *Vdet3 = new VirtualDetector(file,"bm2");
    
    //#######################     
    // Setting Detector tree
    //#######################
    
    Detector *Layer1 = new Detector(file,"rc0");
    Detector *Layer2 = new Detector(file,"rc1");
    //Detector *Layer3 = new Detector(file,"rc2");

    cout << "totalEvent=" << totalEvent << endl;

    //###########
    // Event loop
    //###########
    
    for(ULong64_t ev=0; ev</*10000*/totalEvent; ++ev){
      if(ev%100000==0)
	cout << "\r"  << fname[No].substr(fname[No].find_last_of("/")+1) << " : " /*<< fixed << setprecision(2) */<< (100*ev)/(int)totalEvent << "%" << flush;
      float t1=0;
      float t2=0;

      //######################
      //NewParticle & BeamLoss
      //######################
      
      NP -> SetEventInfo(ev);
      float ini_momentum = NP -> SetInitialInfo() -> Get_Pz();
      momentum_distribution[4*degNo+pdgNo] -> Fill(ini_momentum);

      blnt -> SetEventInfo(ev);
      float stop_x = blnt -> SetInitialInfo() -> Get_x();
      float stop_y = blnt -> SetInitialInfo() -> Get_y();
      float stop_z = blnt -> SetInitialInfo() -> Get_z();
      DecayPosition[4*degNo+pdgNo] -> Fill(stop_z,ini_momentum);
      if(stop_z >= 0 && stop_z <= degraderlength){
	MomDist_InDeg[4*degNo+pdgNo] -> Fill(ini_momentum);
      }
      
      //################
      // VirtualDetector
      //################
      
      /* Vdet1 -> SetEventInfo(ev); */
      /* if(Vdet1->InParentID(0)){ */
      /* 	MomDist_Vdet1[4*degNo+pdgNo] -> Fill(ini_momentum); */
      /* } */
      
      /* Vdet2 -> SetEventInfo(ev); */
      /* if(Vdet2->InParentID(0)){ */
      /* 	MomDist_Vdet2[4*degNo+pdgNo] -> Fill(ini_momentum); */
      /* }       */
      
      /* Vdet3 -> SetEventInfo(ev); */
      /* if(Vdet3->InParentID(0)){ */
      /* 	MomDist_Vdet3[4*degNo+pdgNo] -> Fill(ini_momentum); */
      /* } */
            
      //#######
      // layer1
      //#######

      Layer1 -> SetEventInfo(ev);
      if(Layer1 -> GetEntries() != 0 && Layer1 -> InPDGid(ini_pdgid)){
	t1 = Layer1 -> Get_t(0);
      }
  
      //#######
      // layer2
      //#######

      Layer2 -> SetEventInfo(ev);
      if(Layer2 -> GetEntries() != 0 && Layer2 -> InPDGid(ini_pdgid)){
	t2 = Layer2 -> Get_t(0);
      }

      if(ini_pdgid==ParticleID::positive_pion && Layer2 -> GetEntries() != 0){
	for(int i=0, entries=Layer2->GetEntries(); i<entries;++i){
	  PionEdep[degNo] -> Fill(Layer2 -> Get_VisibleEdep(i));
	}
      }	    
      
      
      //#######
      // layer3
      //#######
      //Layer3 -> SetEventInfo(ev);


      //#####################
      // stop in absorber
      //#####################
      //if(Vdet2->InParentID(0) && !Vdet3->InParentID(0)){
      if(stop_z >= degraderlength+40 && stop_z <= degraderlength+40+stod(absorberlength)){
	EdepWhenStopAbso[4*degNo+pdgNo] -> Fill(Layer1->SetInitialInfo()->Get_VisibleEdep());
	MomDist_InAbso[4*degNo+pdgNo] -> Fill(ini_momentum);
	for(int i=0, entries=Layer2->GetEntries(); i<entries; ++i){
	  if(Layer2->Get_ParentID(i)==1 && Layer2->Get_PDGid(i)<2000){
	    decay_spectrum[4*degNo+pdgNo] -> Fill(Layer2->Get_t(i)-t1);
	  }
	}
	for(int i=0, entries=NP->GetEntries(); i<entries; ++i){
	  if(NP->Get_x(i)==stop_x && NP->Get_y(i)==stop_y && NP->Get_z(i)==stop_z){
	    if(NP->Get_ParentID(i)==1 && NP->Get_PDGid(i)<2000){
	      PID[4*degNo+pdgNo] -> Fill(NP->Get_PDGid(i));
	      float p_x = NP -> Get_Px(i);
	      float p_y = NP -> Get_Py(i);
	      float p_z = NP -> Get_Pz(i);
	      float p = sqrt(p_x*p_x+p_y*p_y+p_z*p_z);
	      momentum_direction[4*degNo+pdgNo] -> Fill(p_z/p,p_y/p);
	    }
	  }
	}
      }
      
      //###################
      //Decay spectrum
      //###################
      
      /*   if( t1!=0 && t1<5 && (t2-t1>1.5) ){ */
      /* 	for(size_t j=0; j<Layer2->GetEntries();++j){ */
      /* 	  if(Layer2->Get_t(j)>t1 && Layer2->Get_Edep(j)>0.05 /\*&& !InContainer(remove,Layer2->GetNthTime(j))*\/){ */
      /* 	    decay_spectrum[4*degNo+pdgNo]->Fill(Layer2->Get_t(j)-t1); */
      /* 	  } */
      /* 	} */
      /*   } */
    } 
    if(ini_particle=="mu-"){
      int mom_mode = MomDist_InAbso[4*degNo+0]->GetMaximumBin()-1;
      MaximumStopBin[pdgNo] -> SetPoint(degNo,mom_mode,MomDist_InAbso[4*degNo+0]->GetBinContent(mom_mode+1));
      ratio[pdgNo] -> SetPoint(degNo,mom_mode,momentum_distribution[0]->GetBinContent(mom_mode+1)/(float)MomDist_InAbso[4*degNo+0]->GetBinContent(mom_mode+1));
      ScaleFile << degraderlength << "\t" << momentum_distribution[0]->GetBinContent(mom_mode+1)/(float)MomDist_InAbso[4*degNo+0]->GetBinContent(mom_mode+1) << endl;
    }
   
    int DecayEntry = decay_spectrum[4*degNo+pdgNo]->GetEntries();
    int StopNumber = MomDist_InAbso[4*degNo+pdgNo]->GetEntries();
    //int IniParticle_inVdet1 = MomDist_Vdet1[4*degNo+pdgNo]->GetEntries();
    //int IniParticle_inVdet2 = MomDist_Vdet2[4*degNo+pdgNo]->GetEntries();
    //int IniParticle_inVdet3 = MomDist_Vdet3[4*degNo+pdgNo]->GetEntries();
      
    degraderArray[pdgNo][degNo]=(float)degraderlength;
    degraderArrayError[pdgNo][degNo]=0;
    nofStop[pdgNo][degNo]=(float)StopNumber;
    nofStopError[pdgNo][degNo]=sqrt(StopNumber);
    acceptance[pdgNo][degNo]=(float)DecayEntry/StopNumber;
    acceptanceError[pdgNo][degNo]=0/*sqrt((float)DecayEntry/(IniParticle_inVdet2-IniParticle_inVdet3))*/;
    
    file->Close();
    cout << "\r"  << fname[No].substr(fname[No].find_last_of("/")+1) << " : " << 100 << "%\n"  << flush;
    /* cout << ini_particle << " number:" << endl; */
    /* cout << IniParticle_inVdet1 << " in Vdet1." << endl; */
    /* cout << IniParticle_inVdet2 << " in Vdet2." << endl; */
    /* cout << IniParticle_inVdet3 << " in Vdet3." << endl; */

    delete NP;
    delete blnt;
    delete Vdet1;
    delete Vdet2;
    delete Vdet3;
    delete Layer1;
    delete Layer2;
    //delete Layer3;
  }
  ScaleFile.close();
  
  TString Name;
  Name=Form("TrueAna_absorber%s_%s.pdf",absorberlength.data(),material.data());
  TCanvas *c1 = new TCanvas(Name.Data(),Name.Data());
  c1 -> Print(Name+"[","pdf");
  TFile *OUTFILE = new TFile("TrueAna.root","RECREATE");
  OUTFILE -> cd();
  int page=1;
    
  TGraphErrors *StoppedParticle[3];
  TGraphErrors *Acceptance[3];
  for(int i=0; i<3;++i){
    StoppedParticle[i] = new TGraphErrors(nofdeg,degraderArray[i],nofStop[i],degraderArrayError[i],nofStopError[i]);
    StoppedParticle[i] -> SetMarkerStyle(20);
    StoppedParticle[i] -> SetMarkerSize(1);
    StoppedParticle[i] -> GetXaxis() -> SetTitle("degraderlength[mm]");
    StoppedParticle[i] -> GetYaxis() -> SetTitle("# of stopped in absorber");
    StoppedParticle[i] -> SetMinimum(0);
    Acceptance[i] = new TGraphErrors(nofdeg,degraderArray[i],acceptance[i],degraderArrayError[i],acceptanceError[i]);
    Acceptance[i] -> SetMarkerStyle(20);
    Acceptance[i] -> SetMarkerSize(1);
    Acceptance[i] -> GetXaxis() -> SetTitle("degraderlength[mm]");
    Acceptance[i] -> GetYaxis() -> SetTitle("Acceptance");
    Acceptance[i] -> SetMaximum(0);
    Acceptance[i] -> SetMinimum(0);
  }
  for(auto itr : particleNo){
    if(itr.second == 3) continue;
    StoppedParticle[itr.second] -> SetTitle(Form("# of #%s stopped at absorber%smm",itr.first.data(),absorberlength.data()));
    Acceptance[itr.second] -> SetTitle(Form("Acceptance of %s",itr.first.data()));
  }


  TPaveStats *st1;
  TPaveStats *st2;
  TPaveStats *st3;
  TLegend *lg1 = new TLegend(0.8,0.40,0.87,0.55);

  gPad -> SetLogy();
  momentum_distribution[0] -> Draw();
  for(int i=0; i<3; ++i){
    momentum_distribution[i] -> Draw("sames");
  }
  c1 -> Update();
  st1 = (TPaveStats*)momentum_distribution[0]->FindObject("stats");
  st1->SetTextColor(2*(0+1));
  st2 = (TPaveStats*)momentum_distribution[1]->FindObject("stats");
  st2->SetTextColor(2*(1+1));
  st3 = (TPaveStats*)momentum_distribution[2]->FindObject("stats");
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
  lg1->AddEntry(momentum_distribution[0]," #mu-","l");
  lg1->AddEntry(momentum_distribution[1]," #mu+","l");
  lg1->AddEntry(momentum_distribution[2]," #pi+","l");
  lg1->Draw();
  c1 -> Print(Name,"pdf");
  c1 -> Write(Form("pageNo%d",page++));
    
  for(int i=0; i<3; ++i){
    gPad -> SetLogy(0);
    c1 -> SetGridy();
    //c1 -> Clear();
    Acceptance[i] -> Draw("AP");
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));
    //c1 -> Clear();
    gPad -> SetLogy();
    StoppedParticle[i] -> Draw("AP");
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));
    if(i==0){
      momentum_distribution[0] -> Draw();
      MaximumStopBin[0] -> Draw("P");
      c1 -> Print(Name,"pdf");
      ratio[0] -> Draw("AP");
      c1 -> Print(Name,"pdf");
      c1 -> Write(Form("pageNo%d",page++));
    }
  }

  for(auto itr=degraderID.begin(); itr!=degraderID.end(); ++itr){
    n=itr->second;
    int deglen = itr->first;

    MomDist_InDeg[4*n+0] -> Draw();
    for(int i=1; i<3; ++i){
      MomDist_InDeg[4*n+i] -> Draw("sames");
    }
    c1 -> Update();
    st1 = (TPaveStats*)MomDist_InDeg[4*n+0]->FindObject("stats");
    st1->SetTextColor(2*(0+1));
    st2 = (TPaveStats*)MomDist_InDeg[4*n+1]->FindObject("stats");
    st2->SetTextColor(2*(1+1));
    st3 = (TPaveStats*)MomDist_InDeg[4*n+2]->FindObject("stats");
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
    lg1->Clear();
    lg1->AddEntry(MomDist_InDeg[4*n+0]," #mu-","l");
    lg1->AddEntry(MomDist_InDeg[4*n+1]," #mu+","l");
    lg1->AddEntry(MomDist_InDeg[4*n+2]," #pi+","l");
    lg1->Draw();
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));

    
    MomDist_InAbso[4*n+0] -> Draw();
    for(int i=1; i<3; ++i){
      MomDist_InAbso[4*n+i] -> Draw("sames");
    }
    c1 -> Update();
    st1 = (TPaveStats*)MomDist_InAbso[4*n+0]->FindObject("stats");
    st1->SetTextColor(2*(0+1));
    st2 = (TPaveStats*)MomDist_InAbso[4*n+1]->FindObject("stats");
    st2->SetTextColor(2*(1+1));
    st3 = (TPaveStats*)MomDist_InAbso[4*n+2]->FindObject("stats");
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
    lg1->Clear();
    lg1->AddEntry(MomDist_InAbso[4*n+0]," #mu-","l");
    lg1->AddEntry(MomDist_InAbso[4*n+1]," #mu+","l");
    lg1->AddEntry(MomDist_InAbso[4*n+2]," #pi+","l");
    lg1->Draw();
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));

    
    MomDist_InDeg[4*n+0] -> Draw();
    MomDist_InDeg[4*n+0] -> SetLineColor(1);
    MomDist_InAbso[4*n+0] -> Draw("sames");
    MomDist_InAbso[4*n+0] -> SetLineColor(2);
    c1 -> Update();
    st1 = (TPaveStats*)MomDist_InDeg[4*n+0]->FindObject("stats");
    st1->SetTextColor(1);
    st2 = (TPaveStats*)MomDist_InAbso[4*n+0]->FindObject("stats");
    st2->SetTextColor(2);
    st1->SetX1NDC(0.78);
    st1->SetX2NDC(0.98);
    st1->SetY1NDC(0.83);
    st1->SetY2NDC(0.91);
    st2->SetX1NDC(0.78);
    st2->SetX2NDC(0.98);
    st2->SetY1NDC(0.72);
    st2->SetY2NDC(0.80);
    lg1->Clear();
    lg1->AddEntry(MomDist_InDeg[4*n+0]," in degrader","l");
    lg1->AddEntry(MomDist_InAbso[4*n+0]," in absorber","l");
    lg1->Draw();
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));
    
    
    /* MomDist_Vdet2[4*n+0] -> Draw(); */
    /* for(int i=1; i<3; ++i){ */
    /*   MomDist_Vdet2[4*n+i] -> Draw("sames"); */
    /* } */
    /* c1 -> Print(Name,"pdf"); */
    /* c1->Write(Form("pageNo%d",page++)); */
    
    /* MomDist_Vdet3[4*n+0] -> Draw(); */
    /* for(int i=1; i<3; ++i){ */
    /*   MomDist_Vdet3[4*n+i] -> Draw("sames"); */
    /* } */
    /* c1 -> Print(Name,"pdf"); */
    /* c1->Write(Form("pageNo%d",page++)); */

    for(int i=0; i<3; ++i){
      EdepWhenStopAbso[4*n+i] -> Draw();
      c1 -> Print(Name,"pdf");
      c1 -> Write(Form("pageNo%d",page++));

      gPad -> SetLogy(0);
      DecayPosition[4*n+i] -> Draw("colz");
      TLine *l1 = new TLine(0,0,0,300);
      l1 -> SetLineWidth(1);
      l1 -> SetLineStyle(3);
      l1 -> Draw();
      TLine *l2 = new TLine(deglen,0,deglen,300);
      l2 -> SetLineWidth(1);
      l2 -> SetLineStyle(3);
      l2 -> Draw();
      TLine *l3 = new TLine(deglen+10,0,deglen+10,300);
      l3 -> SetLineWidth(1);
      l3 -> SetLineStyle(3);
      l3 -> Draw();
      TLine *l4 = new TLine(deglen+10.5,0,deglen+10.5,300);
      l4 -> SetLineWidth(1);
      l4 -> SetLineStyle(3);
      l4 -> Draw();
      TLine *l5 = new TLine(deglen+40,0,deglen+40,300);
      l5 -> SetLineWidth(1);
      l5 -> SetLineStyle(3);
      l5 -> Draw();
      TLine *l6 = new TLine(deglen+40+stod(absorberlength),0,deglen+40+stod(absorberlength),300);
      l6 -> SetLineWidth(1);
      l6 -> SetLineStyle(3);
      l6 -> Draw();
      TLine *l7 = new TLine(deglen+70,0,deglen+70,300);
      l7 -> SetLineWidth(1);
      l7 -> SetLineStyle(3);
      l7 -> Draw();
      TLine *l8 = new TLine(deglen+75,0,deglen+75,300);
      l8 -> SetLineWidth(1);
      l8 -> SetLineStyle(3);
      l8 -> Draw();
      c1 -> Print(Name,"pdf");
      c1->Write(Form("pageNo%d",page++));
      c1 -> Clear();

      PID[4*n+i] -> Draw();
      c1 -> Print(Name,"pdf");
      c1 -> Write(Form("pageNo%d",page++));
      
      momentum_direction[4*n+i] -> Draw("colz");
      c1 -> Print(Name,"pdf");
      c1 -> Write(Form("pageNo%d",page++));
      
      gPad -> SetLogy();
      if(i==2) continue;
      decay_spectrum[4*n+i] -> Draw();
      c1 -> Print(Name,"pdf");
      c1 -> Write(Form("pageNo%d",page++));
    }
    PionEdep[n] -> Draw();
    c1 -> Print(Name,"pdf");
    c1 -> Write(Form("pageNo%d",page++));
  }
  
  c1 -> Print(Name+"]","pdf");
  OUTFILE -> Close();
  
  return 0;
}

