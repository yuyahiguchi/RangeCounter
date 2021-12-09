#ifndef functions_h_
#define functions_h_

#include <iostream>
#include <tuple>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;

std::tuple<float,float,float,float,float> GetFit(TH1D *hist,string material="Cu"){
  double tau_mu=2193.5;
  double tau_mu_Cu=180.8;
  double tau_mu_Al=880;
  double tau_pi=22.2;
  
  float P0;
  float P0_error;
  float P1;
  float P1_error;
  float LT;

  TF1 *fit_func;
  if(material=="Al"){
    fit_func = new TF1("fit_func","[0]*TMath::Exp(-x/2200)+[1]*TMath::Exp(-x/859)",0,5000);
  }
  else if(material=="Cu"){
    //fit_func = new TF1("fit_func","[0]*TMath::Exp(-x/2200)+[1]*TMath::Exp(-x/160)",0,5000);
    fit_func = new TF1("fit_func","[0]*TMath::Exp(-x/[3])+[1]*TMath::Exp(-x/[2])",0,5000);
    fit_func->SetParameter(3,2200);
    fit_func->SetParameter(2,163);
  }
  //fit_func = new TF1("fit_func","[0]*TMath::Exp(-x/[3])+[1]*TMath::Exp(-x/[4])+[2]*((1/[5])/(1/[5]-1/[3]))*(TMath::Exp(-x/[3])-TMath::Exp(-x/[5]))",0,1000);
  // fit_func->SetParameter(3,2200);
  //fit_func->SetParameter(2,163);
  // fit_func->SetParameter(5,26);

  hist->Fit("fit_func","L","");
  cout << "Life time=" << fit_func->GetParameter(2) << endl;
  P0 = fit_func->GetParameter(0);
  P0_error = fit_func->GetParError(0);
  P1 = fit_func->GetParameter(1);
  P1_error = fit_func->GetParError(1);
  LT = fit_func->GetParameter(2);
  return {P0,P0_error,P1,P1_error,LT};
}

//-------------------------------------------------------------------------

bool InContainer(vector<float> vec, float value){
  if(find(vec.begin(),vec.end(),value)==vec.end())
    return false;
  else
    return true;
}

bool InContainer(vector<int> vec, int value){
  if(find(vec.begin(),vec.end(),value)==vec.end())
    return false;
  else
    return true;
}

#endif





//-------------------------------------------------------------------------

// void GetNewVector(vector<float>&vec1, vector<float>&vec2){ //time,Edep
//   vector<vector<float> > tempvec;
//   tempvec.resize((int)vec1.size());
//   for(int i=0; i<(int)vec1.size();++i){
//     tempvec[i].push_back(vec1[i]);
//     tempvec[i].push_back(vec2[i]);
//   }
//   sort(tempvec.begin(),tempvec.end(),[](const vector<float> a,const vector<float> b){return a[0]<b[0];});
//   vec1.clear();
//   vec1.shrink_to_fit();
//   vec2.clear();
//   vec2.shrink_to_fit();
//   for(int i=0; i<tempvec.size();++i){
//     vec1.push_back(tempvec[i][0]);
//     vec2.push_back(tempvec[i][1]);
//   }
//   tempvec.clear();
//   tempvec.shrink_to_fit();
//   float val;
//   for(int i=vec1.size()-1; i>0; --i){
//     if(vec1[i]-vec1[i-1]<5){
//       val=vec2[i]+vec2[i-1];
//       vec1.erase(vec1.begin()+i);
//       vec2.erase(vec2.begin()+i);
//       vec2[i-1]=val;
//     }
//   }
// }


