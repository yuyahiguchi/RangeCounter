#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>
#include "TApplication.h"
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
#include "TGraphErrors.h"

#include "../include/ParticleID.h"
#include "../include/functions.h"
//vectorのサイズが0やった場合どうするか？

using namespace std;
using std::filesystem::directory_iterator;


int main(int argc, char** argv){
  const int vdet_count = 5;
  const int  det_count = 5;
  vector<int> nuetrino;
  nuetrino.push_back(ParticleID::nu_e);
  nuetrino.push_back(ParticleID::nu_mu);
  nuetrino.push_back(ParticleID::nu_tau);
  nuetrino.push_back(ParticleID::anti_nu_e);
  nuetrino.push_back(ParticleID::anti_nu_mu);
  nuetrino.push_back(ParticleID::anti_nu_tau);
  
  string fname = argv[1];
  string fout = fname.substr(0,fname.length()-5) + "_compressed.root";
  TFile *InputFile;
  TFile *OutputFile;
  
  TTree *VirtualTree[vdet_count];
  float vir_x[vdet_count]={};
  float vir_y[vdet_count]={};
  float vir_z[vdet_count]={};	
  float vir_t[vdet_count]={};
  float vir_Px[vdet_count]={};
  float vir_Py[vdet_count]={};
  float vir_Pz[vdet_count]={};
  float vir_EventID[vdet_count]={};
  float vir_PDGid[vdet_count]={};
  float vir_TrackID[vdet_count]={};
  float vir_ParentID[vdet_count]={};
  int vir_nEntry[vdet_count]={};
  
  TTree *DetectorTree[det_count];
  float det_x[det_count]={};
  float det_y[det_count]={};
  float det_z[det_count]={};
  float det_t[det_count]={};
  float det_Px[det_count]={};
  float det_Py[det_count]={};
  float det_Pz[det_count]={};
  float det_Edep[det_count]={};
  float det_VisibleEdep[det_count]={};
  float det_EventID[det_count]={};
  float det_PDGid[det_count]={};
  float det_TrackID[det_count]={};
  float det_ParentID[det_count]={};
  int det_nEntry[det_count]={};

  TTree *NewParticleTree;
  float NP_x;
  float NP_y;
  float NP_z;
  float NP_t;
  float NP_Px;
  float NP_Py;
  float NP_Pz;
  float NP_EventID;
  float NP_PDGid;
  float NP_TrackID;
  float NP_ParentID;
  int NP_nEntry;

  TTree *BeamLossTree;
  float BL_x;
  float BL_y;
  float BL_z;
  float BL_t;
  float BL_Px;
  float BL_Py;
  float BL_Pz;
  float BL_EventID;
  float BL_PDGid;
  float BL_TrackID;
  float BL_ParentID;
  int BL_nEntry;


  float x;
  float y;
  float z;
  float t;
  float Px;
  float Py;
  float Pz;
  float Edep;
  float VisibleEdep;
  float EventID;
  float PDGid;
  float TrackID;
  float ParentID;

  //TFile *OutputFile = new TFile("compressed.root","recreate");
  TTree *virtualtree[vdet_count];
  for(int i=0; i<vdet_count; ++i){
    virtualtree[i] = new TTree(Form("bm%d",i),"virtual detector");
    virtualtree[i]->Branch("x",&x);
    virtualtree[i]->Branch("y",&y);
    virtualtree[i]->Branch("z",&z);
    virtualtree[i]->Branch("t",&t);
    virtualtree[i]->Branch("Px",&Px);
    virtualtree[i]->Branch("Py",&Py);
    virtualtree[i]->Branch("Pz",&Pz);
    virtualtree[i]->Branch("EventID",&EventID);
    virtualtree[i]->Branch("PDGid",&PDGid);
    virtualtree[i]->Branch("TrackID",&TrackID);
    virtualtree[i]->Branch("ParentID",&ParentID);
  }
  
  TTree *detectortree[det_count];
  for(int i=0; i<det_count; ++i){
    detectortree[i] = new TTree(Form("rc%d",i),"range counter");
    detectortree[i]->Branch("x",&x);
    detectortree[i]->Branch("y",&y);
    detectortree[i]->Branch("z",&z);
    detectortree[i]->Branch("t",&t);
    detectortree[i]->Branch("Px",&Px);
    detectortree[i]->Branch("Py",&Py);
    detectortree[i]->Branch("Pz",&Pz);
    detectortree[i]->Branch("Edep",&Edep);
    detectortree[i]->Branch("VisibleEdep",&VisibleEdep);
    detectortree[i]->Branch("EventID",&EventID);
    detectortree[i]->Branch("PDGid",&PDGid);
    detectortree[i]->Branch("TrackID",&TrackID);
    detectortree[i]->Branch("ParentID",&ParentID);
  }

  TTree *newparticlentupletree;
  newparticlentupletree = new TTree("NP","new particle Ntuple");
  newparticlentupletree -> Branch("x",&x);
  newparticlentupletree -> Branch("y",&y);
  newparticlentupletree -> Branch("z",&z);
  newparticlentupletree -> Branch("t",&t);
  newparticlentupletree -> Branch("Px",&Px);
  newparticlentupletree -> Branch("Py",&Py);
  newparticlentupletree -> Branch("Pz",&Pz);
  newparticlentupletree -> Branch("Edep",&Edep);
  newparticlentupletree -> Branch("VisibleEdep",&VisibleEdep);
  newparticlentupletree -> Branch("EventID",&EventID);
  newparticlentupletree -> Branch("PDGid",&PDGid);
  newparticlentupletree -> Branch("TrackID",&TrackID);
  newparticlentupletree -> Branch("ParentID",&ParentID);

  TTree *beamlossntupletree;
  beamlossntupletree = new TTree("blnt","beam loss Ntuple");
  beamlossntupletree -> Branch("x",&x);
  beamlossntupletree -> Branch("y",&y);
  beamlossntupletree -> Branch("z",&z);
  beamlossntupletree -> Branch("t",&t);
  beamlossntupletree -> Branch("Px",&Px);
  beamlossntupletree -> Branch("Py",&Py);
  beamlossntupletree -> Branch("Pz",&Pz);
  beamlossntupletree -> Branch("Edep",&Edep);
  beamlossntupletree -> Branch("VisibleEdep",&VisibleEdep);
  beamlossntupletree -> Branch("EventID",&EventID);
  beamlossntupletree -> Branch("PDGid",&PDGid);
  beamlossntupletree -> Branch("TrackID",&TrackID);
  beamlossntupletree -> Branch("ParentID",&ParentID);


  //**********************************************************************************************
  //**********************************************************************************************

    
  InputFile = new TFile(fname.data());
  OutputFile = new TFile(fout.data(),"RECREATE");
  cout << "Open " << fname << "."<< endl;
  cout << "Remove nuetrinos and anti-nietrinos from tree." << endl;
  // Detectorのtreeをセット
  for(int i=0; i<det_count; ++i){
    DetectorTree[i] = (TTree*)InputFile->Get(Form("Detector/rc%d",i));
    if(DetectorTree[i] == nullptr){
      continue;
      //DetectorTree[i] = (TTree*)InputFile->Get(Form("rc%d",i));
    }
    DetectorTree[i]->SetBranchAddress("x",&det_x[i]);
    DetectorTree[i]->SetBranchAddress("y",&det_y[i]);
    DetectorTree[i]->SetBranchAddress("z",&det_z[i]);
    DetectorTree[i]->SetBranchAddress("t",&det_t[i]);
    DetectorTree[i]->SetBranchAddress("Px",&det_Px[i]);
    DetectorTree[i]->SetBranchAddress("Py",&det_Py[i]);
    DetectorTree[i]->SetBranchAddress("Pz",&det_Pz[i]);
    DetectorTree[i]->SetBranchAddress("Edep",&det_Edep[i]);
    DetectorTree[i]->SetBranchAddress("VisibleEdep",&det_VisibleEdep[i]);
    DetectorTree[i]->SetBranchAddress("EventID",&det_EventID[i]);
    DetectorTree[i]->SetBranchAddress("PDGid",&det_PDGid[i]);
    DetectorTree[i]->SetBranchAddress("TrackID",&det_TrackID[i]);
    DetectorTree[i]->SetBranchAddress("ParentID",&det_ParentID[i]);
    det_nEntry[i]=DetectorTree[i]->GetEntries();
    for(int j=0; j<det_nEntry[i]; ++j){
      DetectorTree[i] -> GetEntry(j);
      if(!InContainer(nuetrino,det_PDGid[i])){
	x = det_x[i];
	y = det_y[i];
	z = det_z[i];
	t = det_t[i];
	Px = det_Px[i];
	Py = det_Py[i];
	Pz = det_Pz[i];
	Edep = det_Edep[i];
	VisibleEdep = det_VisibleEdep[i];
	EventID = det_EventID[i];
	PDGid = det_PDGid[i];
	TrackID = det_TrackID[i];
	ParentID = det_ParentID[i];
	detectortree[i] -> Fill();
      }
    }
    detectortree[i] -> Write();
  }

  // VirtualDetectorのtreeをセット
  for(int i=0; i<vdet_count; ++i){
    VirtualTree[i] = (TTree*)InputFile->Get(Form("VirtualDetector/bm%d",i+1));
    if(VirtualTree[i] == nullptr){
      continue;
      //VirtualTree[i] = (TTree*)InputFile->Get(Form("bm%d",i));
    }
    VirtualTree[i]->SetBranchAddress("x",&vir_x[i]);
    VirtualTree[i]->SetBranchAddress("y",&vir_y[i]);
    VirtualTree[i]->SetBranchAddress("z",&vir_z[i]);
    VirtualTree[i]->SetBranchAddress("t",&vir_t[i]);
    VirtualTree[i]->SetBranchAddress("Px",&vir_Px[i]);
    VirtualTree[i]->SetBranchAddress("Py",&vir_Py[i]);
    VirtualTree[i]->SetBranchAddress("Pz",&vir_Pz[i]);
    VirtualTree[i]->SetBranchAddress("EventID",&vir_EventID[i]);
    VirtualTree[i]->SetBranchAddress("PDGid",&vir_PDGid[i]);
    VirtualTree[i]->SetBranchAddress("TrackID",&vir_TrackID[i]);
    VirtualTree[i]->SetBranchAddress("ParentID",&vir_ParentID[i]);
    vir_nEntry[i]=VirtualTree[i]->GetEntries();
    for(int j=0; j<vir_nEntry[i]; ++j){
      VirtualTree[i] -> GetEntry(j);
      if(!InContainer(nuetrino,vir_PDGid[i])){
	x = vir_x[i];
	y = vir_y[i];
	z = vir_z[i];
	t = vir_t[i];
	Px = vir_Px[i];
	Py = vir_Py[i];
	Pz = vir_Pz[i];
	EventID = vir_EventID[i];
	PDGid = vir_PDGid[i];
	TrackID = vir_TrackID[i];
	ParentID = vir_ParentID[i];
	virtualtree[i] -> Fill();
      }
    }
    virtualtree[i] -> Write();
  }  

  //newparticlentupleのtreeセット
  NewParticleTree = (TTree*)InputFile->Get("NTuple/NP");
  if(NewParticleTree == nullptr){
    goto next1;
  }
  NewParticleTree -> SetBranchAddress("x",&NP_x);
  NewParticleTree -> SetBranchAddress("y",&NP_y);
  NewParticleTree -> SetBranchAddress("z",&NP_z);
  NewParticleTree -> SetBranchAddress("t",&NP_t);
  NewParticleTree -> SetBranchAddress("Px",&NP_Px);
  NewParticleTree -> SetBranchAddress("Py",&NP_Py);
  NewParticleTree -> SetBranchAddress("Pz",&NP_Pz);
  NewParticleTree -> SetBranchAddress("EventID",&NP_EventID);
  NewParticleTree -> SetBranchAddress("PDGid",&NP_PDGid);
  NewParticleTree -> SetBranchAddress("TrackID",&NP_TrackID);
  NewParticleTree -> SetBranchAddress("ParentID",&NP_ParentID);
  NP_nEntry = NewParticleTree->GetEntries();
  for(int j=0; j<NP_nEntry; ++j){
    NewParticleTree -> GetEntry(j);
    if(!InContainer(nuetrino,NP_PDGid)){
      x = NP_x;
      y = NP_y;
      z = NP_z;
      t = NP_t;
      Px = NP_Px;
      Py = NP_Py;
      Pz = NP_Pz;
      EventID = NP_EventID;
      PDGid = NP_PDGid;
      TrackID = NP_TrackID;
      ParentID = NP_ParentID;
      newparticlentupletree -> Fill();
    }
  }
  newparticlentupletree -> Write();
 next1:
  cout << "NTuple/NP doesn't exist." << endl;

  //beamlossntupleのtreeセット
  BeamLossTree = (TTree*)InputFile->Get("NTuple/blnt");
  if(BeamLossTree == nullptr){
    goto next2;
  }
  BeamLossTree -> SetBranchAddress("x",&BL_x);
  BeamLossTree -> SetBranchAddress("y",&BL_y);
  BeamLossTree -> SetBranchAddress("z",&BL_z);
  BeamLossTree -> SetBranchAddress("t",&BL_t);
  BeamLossTree -> SetBranchAddress("Px",&BL_Px);
  BeamLossTree -> SetBranchAddress("Py",&BL_Py);
  BeamLossTree -> SetBranchAddress("Pz",&BL_Pz);
  BeamLossTree -> SetBranchAddress("EventID",&BL_EventID);
  BeamLossTree -> SetBranchAddress("PDGid",&BL_PDGid);
  BeamLossTree -> SetBranchAddress("TrackID",&BL_TrackID);
  BeamLossTree -> SetBranchAddress("ParentID",&BL_ParentID);
  BL_nEntry = BeamLossTree->GetEntries();
  for(int j=0; j<BL_nEntry; ++j){
    BeamLossTree -> GetEntry(j);
    if(!InContainer(nuetrino,BL_PDGid)){
      x = BL_x;
      y = BL_y;
      z = BL_z;
      t = BL_t;
      Px = BL_Px;
      Py = BL_Py;
      Pz = BL_Pz;
      EventID = BL_EventID;
      PDGid = BL_PDGid;
      TrackID = BL_TrackID;
      ParentID = BL_ParentID;
      beamlossntupletree -> Fill();
    }
  }
  beamlossntupletree -> Write();
 next2:
  cout << "NTuple/blnt doesn't exist." << endl;

  
  InputFile -> Close();
  OutputFile -> Close();
  cout << "Close " << fout << "."<< endl;
  return 0;
}

