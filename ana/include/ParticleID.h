#include <map>

#ifndef ParticleID_h_
#define ParticleID_h_

using namespace std;

struct ParticleID{
  enum{
    muon = 13,
    anti_muon = -13,
    negative_pion = -211,
    positive_pion = 211,
    gamma = 22,
    electron = 11,
    positron = -11,
    proton = 2212,
    neutron = 2112,
    nu_e = 12,
    anti_nu_e = -12,
    nu_mu = 14,
    anti_nu_mu = -14,
    nu_tau = 16,
    anti_nu_tau = -16,
  };
};

#endif
