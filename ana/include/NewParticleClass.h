// #####################################################################################################################################
// New Particle NTuple & derived class
// #####################################################################################################################################
#ifndef _NewParticleClass_h_
#define _NewParticleClass_h_

#include "VirtualDetectorClass.h"

class NewParticle : public VirtualDetector{
private:
  
public:
  NewParticle(TFile *file, const char *treename);
  //NewParticle* SetInitialInfo();
  
};

//------------------------------------------------------------------------------------------------------------------------------------

NewParticle::NewParticle(TFile *file, const char *treename):VirtualDetector(file,treename){
}

// NewParticle* NewParticle::SetInitialInfo(){
//   auto itr = find(ParentID.begin(), ParentID.end(),0);
//   entry = distance(ParentID.begin(), itr);
//   return this;
// }

#endif
