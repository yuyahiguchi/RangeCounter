// #####################################################################################################################################
// Virtual detector & Base class
// #####################################################################################################################################

#ifndef _VirtualDetectorClass_h_
#define _VirtualDetectorClass_h_


class VirtualDetector{
 public:
  VirtualDetector(TFile *file, const char *treename);
  virtual void SetEventInfo(ULong64_t ev);
  ULong64_t GetTotalEvent();
  virtual size_t GetEntries();
  virtual float Get_x(int n);
  virtual float Get_y(int n);
  virtual float Get_z(int n);
  virtual float Get_Px(int n);
  virtual float Get_Py(int n);
  virtual float Get_Pz(int n);
  virtual float Get_t(int n);
  virtual float Get_PDGid(int n);
  virtual float Get_ParentID(int n);
  virtual float Get_EventID(int n);
  virtual float Get_TrackID(int n);
  bool InPDGid(float x);
  bool InParentID(float x);
  virtual void Clear();
  virtual VirtualDetector* SetInitialInfo();
  virtual VirtualDetector* SetNthEntry(int n);
  //virtual void Print();
 
 protected:
  ULong64_t nRow;
  ULong64_t iRow;
  size_t entry;
  TTree *tree;
  float x_1,y_1,z_1,Px_1,Py_1,Pz_1,t_1,PDGid_1,EventID_1,ParentID_1,TrackID_1;
  vector<float> x,y,z,Px,Py,Pz,t,PDGid,EventID,ParentID,TrackID;
};


//--------------------------------------------------------------------------------------------------------------------------------------


VirtualDetector::VirtualDetector(TFile *file, const char *treename){
  entry=0;
  iRow=0;
  tree = (TTree*)file -> Get(treename);
  if(tree == nullptr){
    cout << treename << " doesn't exist in this file." << endl;
    exit(0);
  }
  tree -> SetBranchAddress("x",&x_1);
  tree -> SetBranchAddress("y",&y_1);
  tree -> SetBranchAddress("z",&z_1);
  tree -> SetBranchAddress("Px",&Px_1);
  tree -> SetBranchAddress("Py",&Py_1);
  tree -> SetBranchAddress("Pz",&Pz_1);
  tree -> SetBranchAddress("t",&t_1);
  tree -> SetBranchAddress("PDGid",&PDGid_1);
  tree -> SetBranchAddress("EventID",&EventID_1);
  tree -> SetBranchAddress("ParentID",&ParentID_1);
  tree -> SetBranchAddress("TrackID",&TrackID_1);
  nRow = tree -> GetEntries();
}


inline void VirtualDetector::SetEventInfo(ULong64_t ev){
  this -> Clear();
  for(; iRow<nRow; ++iRow){
    tree -> GetEntry(iRow);
    if(PDGid_1==ParticleID::gamma) continue;
    if(EventID_1==ev){
      x.push_back(x_1);
      y.push_back(y_1);
      z.push_back(z_1);
      Px.push_back(Px_1);
      Py.push_back(Py_1);
      Pz.push_back(Pz_1);
      t.push_back(t_1);
      PDGid.push_back(PDGid_1);
      EventID.push_back(EventID_1);
      ParentID.push_back(ParentID_1);
      TrackID.push_back(TrackID_1);
    }
    else if(EventID_1 < ev){
      continue;
    }
    else if(EventID_1 > ev){
      break;
    }
  }
}


inline ULong64_t VirtualDetector::GetTotalEvent(){
  tree -> GetEntry(nRow-1);
  return EventID_1+1;
}


inline size_t VirtualDetector::GetEntries(){
  return PDGid.size();
}


inline float VirtualDetector::Get_x(int n=-1){
  int m = (n==-1) ? entry : n;
  return x[m];
}

inline float VirtualDetector::Get_y(int n=-1){
  int m = (n==-1) ? entry : n;
  return y[m];
}

inline float VirtualDetector::Get_z(int n=-1){
  int m = (n==-1) ? entry : n;
  return z[m];
}

inline float VirtualDetector::Get_Px(int n=-1){
  int m = (n==-1) ? entry : n;
  return Px[m];
}

inline float VirtualDetector::Get_Py(int n=-1){
  int m = (n==-1) ? entry : n;
  return Py[m];
}

inline float VirtualDetector::Get_Pz(int n=-1){
  int m = (n==-1) ? entry : n;
  return Pz[m];
}

inline float VirtualDetector::Get_t(int n=-1){
  int m = (n==-1) ? entry : n;
  return t[m];
}

inline float VirtualDetector::Get_PDGid(int n=-1){
  int m = (n==-1) ? entry : n;
  return PDGid[m];
}

inline float VirtualDetector::Get_EventID(int n=-1){
  int m = (n==-1) ? entry : n;
  return EventID[m];
}

inline float VirtualDetector::Get_ParentID(int n=-1){
  int m = (n==-1) ? entry : n;
  return ParentID[m];
}

inline float VirtualDetector::Get_TrackID(int n=-1){
  int m = (n==-1) ? entry : n;
  return TrackID[m];
}


inline bool VirtualDetector::InPDGid(float x){
  return InContainer(PDGid,x);
}


inline bool VirtualDetector::InParentID(float x){
  return InContainer(ParentID,x);
}



inline void VirtualDetector::Clear(){
  x.clear();
  y.clear();
  z.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  t.clear();
  PDGid.clear();
  EventID.clear();
  ParentID.clear();
  TrackID.clear();
}



inline VirtualDetector* VirtualDetector::SetInitialInfo(){
  auto itr = find(ParentID.begin(), ParentID.end(),0);
  entry = distance(ParentID.begin(), itr);
  return this;
}


inline VirtualDetector* VirtualDetector::SetNthEntry(int n){
  entry = n;
  return this;
}


// void VirtualDetector::Print(){
//   cout << "Current EventID = " << EventID_1-1 << " **********" << endl;
//   for(int i=0; i<this->GetEntries()(); ++i){
//     cout << "entry " << i << "  :  t=" << t[i] << " | PDGid=" << PDGid[i] << " | ParentID=" << ParentID[i] << " | \n" << flush;
//   }
// }

#endif
