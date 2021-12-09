// #####################################################################################################################################
// Beam Loss NTuple & derived class
// #####################################################################################################################################
#ifndef _BeamLossClass_h_
#define _BeamLossClass_h_

#include "VirtualDetectorClass.h"

class BeamLoss : public VirtualDetector{
private:
  
public:
  BeamLoss(TFile *file, const char *treename);
  
};

//------------------------------------------------------------------------------------------------------------------------------------

BeamLoss::BeamLoss(TFile *file, const char *treename):VirtualDetector(file,treename){
}


#endif
