/*
 * utility.h
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */

#ifndef UTILITY_H_
#define UTILITY_H_

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "correction.h"
#include <string>

using floats =  ROOT::VecOps::RVec<float>;
using doubles =  ROOT::VecOps::RVec<double>;
using ints =  ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using uchars = ROOT::VecOps::RVec<unsigned char>;

using FourVector = ROOT::Math::PtEtaPhiMVector;
using FourVectorVec = std::vector<FourVector>;
using FourVectorRVec = ROOT::VecOps::RVec<FourVector>;


struct hist1dinfo
{
	ROOT::RDF::TH1DModel hmodel;
	std::string varname;
	std::string weightname;
	std::string mincutstep;
} ;


struct varinfo
{
	std::string varname;
	std::string vardefinition;
	std::string mincutstep;
};

struct cutinfo
{
	std::string cutdefinition;
	std::string idx;
};


// generates vectors of 4 vectors given vectors of pt, eta, phi, mass
FourVectorVec generate_4vec(floats &pt, floats &eta, floats &phi, floats &mass);

// return a vector size equal to length of x all filled with evWeight value
floats weightv(floats &x, float evWeight);

floats sphericity(FourVectorVec &p);

double foxwolframmoment(int l, FourVectorVec &p, int minj=0, int maxj=-1);

floats btvcorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, floats &pts, floats &etas, ints &hadflav, floats &btags);

float pucorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, float ntruepileup);

ints good_idx(ints good); 
float calculate_deltaEta( FourVector &p1, FourVector &p2);
float calculate_deltaPhi( FourVector &p1, FourVector &p2);
float calculate_deltaR( FourVector &p1, FourVector &p2);
float calculate_invMass( FourVector &p1, FourVector &p2);
FourVector sum_4vec( FourVector &p1, FourVector &p2);
floats sort_discriminant( floats discr, floats obj );
FourVector select_leadingvec( FourVectorVec &v );
floats PrintVector(floats myvector);
floats w_reconstruction (FourVectorVec &jets);
floats compute_DR (FourVectorVec &muons, ints goodMuons_charge);
floats dimuon (FourVectorVec &muons, ints Muon_charge);
ints nearest (floats &eta1, floats &eta2, floats &phi1, floats &phi2, ints charge1, ints charge2);
ints nearestSameFlavor (floats &eta1, floats &eta2, floats &phi1, floats &phi2, ints charge1, ints charge2, ints pdgId);

#endif /* UTILITY_H_ */
