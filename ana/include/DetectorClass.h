// #####################################################################################################################################
// Detector & derived class
// #####################################################################################################################################
#ifndef _DetectorClass_h_
#define _DetectorClass_h_

#include "VirtualDetectorClass.h"

class Detector : public VirtualDetector{
 public:
  Detector(TFile *file, const char *treename);
  void SetEventInfo(ULong64_t ev);
  //size_t GetEntries()(bool collected);
  float Get_Edep(int n);
  float Get_VisibleEdep(int n);
  
  // float Get_x(int n, bool collected);
  // float Get_y(int n, bool collected);
  // float Get_z(int n, bool collected);
  // float Get_Px(int n, bool collected);
  // float Get_Py(int n, bool collected);
  // float Get_Pz(int n, bool collected);
  // float Get_t(int n, bool collected);
  // float Get_PDGid(int n, bool collected);
  // float Get_ParentID(int n, bool collected);
  // float Get_EventID(int n, bool collected);
  // float Get_TrackID(int n, bool collected);
  // float Get_Edep(int n, bool collected);
  // float Get_VisibleEdep(int n, bool collected);

  Detector* SetInitialInfo();
  Detector* SetNthEntry(int n);
  void Clear();
  //void Print();

 private:
  float Edep_1, VisibleEdep_1;
  vector<float> Edep, VisibleEdep;
  vector<float> collected_x, collected_y, collected_z, collected_Px, collected_Py, collected_Pz, collected_t, collected_PDGid, collected_EventID, collected_ParentID, collected_TrackID, collected_Edep, collected_VisibleEdep;
  vector<vector<float> > tempvec;
};

//--------------------------------------------------------------------------------------------------------------------------------------

Detector::Detector(TFile *file, const char *treename):VirtualDetector(file,treename){
  tree -> SetBranchAddress("Edep",&Edep_1);
  tree -> SetBranchAddress("VisibleEdep",&VisibleEdep_1);
}
  
// Detector::Detector(TFile *file, const char *treename){
//   iRow=0;
//   entry=0;
//   tree = (TTree*)file -> Get(treename);
//   if(tree == nullptr){
//     cout << treename << " doesn't exist in this file." << endl;
//     exit(0);
//   }
//   tree -> SetBranchAddress("x",&x_1);
//   tree -> SetBranchAddress("y",&y_1);
//   tree -> SetBranchAddress("z",&z_1);
//   tree -> SetBranchAddress("Px",&Px_1);
//   tree -> SetBranchAddress("Py",&Py_1);
//   tree -> SetBranchAddress("Pz",&Pz_1);
//   tree -> SetBranchAddress("t",&t_1);
//   tree -> SetBranchAddress("PDGid",&PDGid_1);
//   tree -> SetBranchAddress("EventID",&EventID_1);
//   tree -> SetBranchAddress("ParentID",&ParentID_1);
//   tree -> SetBranchAddress("TrackID",&TrackID_1);
//   tree -> SetBranchAddress("Edep",&Edep_1);
//   tree -> SetBranchAddress("VisibleEdep",&VisibleEdep_1);

//   nRow = tree -> GetEntries();
// }


inline void Detector::SetEventInfo(ULong64_t ev){
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
      Edep.push_back(Edep_1);
      VisibleEdep.push_back(VisibleEdep_1);
    }
    else if(EventID_1 < ev){
      continue;
    }
    else if(EventID_1 > ev){
      break;
    }
  }


  tempvec.resize((int)t.size());
  for(int i=0; i<(int)t.size();++i){
    tempvec[i].push_back(t[i]);
    tempvec[i].push_back(Edep[i]);
    tempvec[i].push_back(VisibleEdep[i]);
    tempvec[i].push_back(PDGid[i]);
    tempvec[i].push_back(ParentID[i]);
  }
  t.clear();
  Edep.clear();
  VisibleEdep.clear();
  PDGid.clear();
  ParentID.clear();
  sort(tempvec.begin(),tempvec.end(),[](const vector<float> a,const vector<float> b){return a[0]<b[0];});
  for(int i=0; i<tempvec.size();++i){
    t.push_back(tempvec[i][0]);
    Edep.push_back(tempvec[i][1]);
    VisibleEdep.push_back(tempvec[i][2]);
    PDGid.push_back(tempvec[i][3]);
    ParentID.push_back(tempvec[i][4]);

    collected_t.push_back(tempvec[i][0]);
    collected_Edep.push_back(tempvec[i][1]);
    collected_VisibleEdep.push_back(tempvec[i][2]);
    collected_PDGid.push_back(tempvec[i][3]);
    collected_ParentID.push_back(tempvec[i][4]);
  }
  tempvec.clear();
  tempvec.shrink_to_fit();
  
  float val1;
  float val2;
  for(int i=collected_t.size()-1; i>0; --i){
    if(collected_t[i]-collected_t[i-1] < 5){
      val1=collected_Edep[i]+collected_Edep[i-1];
      val2=collected_VisibleEdep[i]+collected_VisibleEdep[i-1];
      collected_t.erase(collected_t.begin()+i);
      collected_Edep.erase(collected_Edep.begin()+i);
      collected_VisibleEdep.erase(collected_VisibleEdep.begin()+i);
      collected_PDGid.erase(collected_PDGid.begin()+i);
      collected_ParentID.erase(collected_ParentID.begin()+i);
      
      collected_Edep[i-1]=val1;
      collected_VisibleEdep[i-1]=val2;   
    }
  }
}



// size_t Detector::GetEntries()(bool collected=true){
//   if(collected==true){
//     return collected_t.size();
//   }
//   else if(collected==false){
//     return t.size();
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }


inline float Detector::Get_Edep(int n=-1){
  int m = (n==-1) ? entry : n;
  return Edep[m];
}

inline float Detector::Get_VisibleEdep(int n=-1){
  int m = (n==-1) ? entry : n;
  return VisibleEdep[m];
}



inline Detector* Detector::SetInitialInfo(){
  auto itr = find(ParentID.begin(), ParentID.end(),0);
  entry = distance(ParentID.begin(), itr);
  return this;
}



inline Detector* Detector::SetNthEntry(int n){
  entry = n;
  return this;
}



inline void Detector::Clear(){
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
  Edep.clear();
  VisibleEdep.clear();
  
  collected_t.clear();
  collected_Edep.clear();
  collected_VisibleEdep.clear();
  collected_PDGid.clear();
  collected_ParentID.clear();
}



// float Detector::Get_t(int n=-1, bool collected=true){
//   int m = (n==-1) ? entry : n;
//   if(collected==true){
//     return collected_t[m];
//   }
//   else if(collected==false){
//     return t[m];
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }



// float Detector::Get_Edep(int n=-1, bool collected=true){
//   int m = (n==-1) ? entry : n;
//   if(collected==true){
//     return collected_Edep[m];
//   }
//   else if(collected==false){
//     return Edep[m];
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }



// float Detector::Get_VisibleEdep(int n=-1, bool collected=true){
//   int m = (n==-1) ? entry : n;
//   if(collected==true){
//     return collected_VisibleEdep[m];
//   }
//   else if(collected==false){
//     return VisibleEdep[m];
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }



// float Detector::Get_PDGid(int n=-1, bool collected=true){
//   int m = (n==-1) ? entry : n;
//   if(collected==true){
//     return collected_PDGid[m];
//   }
//   else if(collected==false){
//     return PDGid[m];
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }



// float Detector::Get_ParentID(int n=-1, bool collected=true){
//   int m = (n==-1) ? entry : n;
//   if(collected==true){
//     return collected_ParentID[m];
//   }
//   else if(collected==false){
//     return ParentID[m];
//   }
//   else{
//     cout << "unexpected things happen" << endl;
//     exit(1);
//   }
// }


#endif
