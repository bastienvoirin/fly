/*
 * utility.cpp
 *
 *  Created on: Dec 4, 2018
 *      Author: suyong
 */
#include "utility.h"
#include "TMatrixDSym.h"
#include "TVectorT.h"
#include "Math/SpecFuncMathMore.h"
#include "correction.h"
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/Rotation3D.h"
#include "Math/Math.h"

// Utility function to generate fourvector objects for thigs that pass selections

using namespace std;

FourVectorVec generate_4vec(floats &pt, floats &eta, floats &phi, floats &mass)
{
	const int nsize = pt.size();
	FourVectorVec fourvecs;
	fourvecs.reserve(nsize);
	for (auto i=0; i<nsize; i++)
	{
		fourvecs.emplace_back(pt[i], eta[i], phi[i], fabs(mass[i]));
	}

	return fourvecs;
}

floats weightv(floats &x, float evWeight)
{
	const int nsize = x.size();
	floats weightvector(nsize, evWeight);
	return weightvector;
}

floats sphericity(FourVectorVec &p)
{
	TMatrixDSym NormMomTensor(3);

	NormMomTensor = 0.0;
	double p2sum = 0.0;
	for (auto x: p)
	{
		p2sum += x.P2();
		double mom[3] = {x.Px(), x.Py(), x.Pz()};
		for (int irow=0; irow<3; irow++)
		{
			for (int icol=irow; icol<3; icol++)
			{
				NormMomTensor(irow, icol) += mom[irow] * mom[icol];
			}
		}
	}
	NormMomTensor *= (1.0/p2sum);
	TVectorT<double> Qrev;
	NormMomTensor.EigenVectors(Qrev);
	floats Q(3);
	for (auto i=0; i<3; i++) Q[i] = Qrev[2-i];

	return Q;
}


double foxwolframmoment(int l, FourVectorVec &p, int minj, int maxj)
{   // PRD 87, 073014 (2013)
	double answer = 0.0;

	double ptsum=0.0;

	if (maxj==-1) // process everything
	{
		maxj = p.size();
	}
	//for (auto x: p)
	for (auto i=minj; i<maxj; i++)
	{
		auto x = p[i];
		ptsum += x.Pt();
		//for (auto y: p)
		for (auto j=minj; j<maxj; j++)
		{
			auto y = p[j];
			double wij = x.Pt() * y.Pt();
			double cosdOmega = x.Vect().Dot(y.Vect()) / (x.P() * y.P());
			if (cosdOmega>1.0) cosdOmega=1.0;
			if (cosdOmega<-1.0) cosdOmega=-1.0;
			answer += wij * ROOT::Math::legendre(l, cosdOmega);
		}
	}
	answer /= ptsum*ptsum;
	if (fabs(answer)>1.0) std::cout << "FW>1 " << answer << std::endl;
	return answer;
}

floats btvcorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string type, std::string sys, floats &pts, floats &etas, ints &hadflav, floats &btags)
{
	floats scalefactors;
	auto nvecs = pts.size();
	scalefactors.reserve(nvecs);
	for (auto i=0; i<nvecs; i++)
	{
		float sfi = cset->at(type)->evaluate({sys, int(hadflav[i]), fabs(float(etas[i])), float(pts[i]), float(btags[i])});
		scalefactors.emplace_back(sfi);
	}
	return scalefactors;
}

float pucorrection(std::unique_ptr<correction::CorrectionSet> &cset, std::string name, std::string syst, float ntruepileup)
{
	return float(cset->at(name)->evaluate({ntruepileup, syst.c_str()}));
}

ints good_idx(ints good)
{
	ints out;
	for(unsigned int i = 0; i < good.size(); i++){
		if( good[i] ){
			out.emplace_back(i);
		}
	}
	return out;

}

float calculate_deltaEta( FourVector &p1, FourVector &p2){
	return p1.Eta() - p2.Eta();
}
float calculate_deltaPhi( FourVector &p1, FourVector &p2){
	return ROOT::Math::VectorUtil::DeltaPhi(p1, p2);
}
float calculate_deltaR( FourVector &p1, FourVector &p2){
	return ROOT::Math::VectorUtil::DeltaR(p1, p2);
}
float calculate_invMass( FourVector &p1, FourVector &p2){
	return ROOT::Math::VectorUtil::InvariantMass(p1, p2);
}
FourVector sum_4vec( FourVector &p1, FourVector &p2){
	return p1+p2;
}
//Get indices that sort the object vectors in descending order
floats sort_discriminant( floats discr, floats obj ){
	auto sorted_discr = Reverse(Argsort(discr));
	floats out;
	for (auto idx : sorted_discr){
		out.emplace_back(obj[idx]);
	}
	return out;
}
FourVector select_leadingvec( FourVectorVec &v ){
	FourVector vout;
	if(v.size() > 0) return v[0];
	else return vout;
}

floats PrintVector(floats myvector){
	for (size_t i = 0; i < myvector.size(); i++){
		std::cout<<myvector[i]<<"\n";
	}

}

//=====================build pair example=============================================//
//W reconstruction : 2 jets
//====================================================================================//
floats w_reconstruction (FourVectorVec &jets){
floats out;
float dijetMass;
float dijetMass_out;
float dijetDR;
const float Mass_W = 80.9; //W mass 
const float Width_W = 10.8; //W width
float X_Min = 99999;
float X_Recwidth_W; 

//Loop on all selected jets
	for(unsigned int j1 = 0; j1<jets.size()-1; j1++){
		for(unsigned int j2 = j1+1; j2<jets.size(); j2++){
			 //select 2 jets
			dijetMass = (jets[j1]+jets[j2]).M();
			  //select a best W candidate in min width in a event
			X_Recwidth_W = std::pow((Mass_W-dijetMass)/Width_W,2);
			if (X_Recwidth_W<X_Min){
				dijetDR = ROOT::Math::VectorUtil::DeltaR(jets[j1],jets[j2]);
				X_Min = X_Recwidth_W;
				dijetMass_out = dijetMass;
			}

		}
	}
	out.push_back(dijetMass_out);   //0:w_mass
	out.push_back(dijetDR);  
	return out;
}

floats compute_DR (FourVectorVec &muons, ints goodMuons_charge){
	floats out;
	float mu_ss_DR;
	float mu_os_DR;
	//std::cout<<"Muonsize: " << muons.size()<<std::endl;
	if(muons.size()>0)
	//Loop on all selected muons
	for(unsigned int mu1 = 0; mu1<muons.size()-1; mu1++){
		for(unsigned int mu2 = mu1+1; mu2<muons.size(); mu2++){
			 //select 2 muons with same sign
			 if (goodMuons_charge[mu1]!=goodMuons_charge[mu2]) continue; //check charge of muons
			mu_ss_DR = ROOT::Math::VectorUtil::DeltaR(muons[mu1],muons[mu2]);
			 //select 2 muons with same sign
			if (goodMuons_charge[mu1]==goodMuons_charge[mu2]) continue;
			mu_os_DR = ROOT::Math::VectorUtil::DeltaR(muons[mu1],muons[mu2]);

		}
	}
	out.push_back(mu_ss_DR);		//0: same sign dimuon DR
	out.push_back(mu_os_DR);        //1: opposite sign dimuon dR
	return out;

}


floats dimuon (FourVectorVec &muons, ints Muon_charge)
{
      floats out;
      float dimuonDphi;
      float dimuonDR;
      float minDimuonDphi = 5.;
      float minDimuonDR = 5.;
        
      // Loop over all selected muons
      for (unsigned int j1 = 0; j1 < muons.size(); j1++)
      {
        for (unsigned int j2 = j1 + 1; j2 < muons.size(); j2++)
        {
          // Select 2 muons of opposite charge
          if (Muon_charge[j1] != Muon_charge[j2])
          {
            dimuonDphi = ROOT::Math::VectorUtil::DeltaPhi(muons[j1],muons[j2]);
            dimuonDR = ROOT::Math::VectorUtil::DeltaR(muons[j1],muons[j2]);
            if (dimuonDR < minDimuonDR)
            {
              minDimuonDphi = dimuonDphi;
              minDimuonDR = dimuonDR;            
            }
          }
        }
      }
      out.push_back(minDimuonDphi);
      out.push_back(minDimuonDR);
      return out;
}

ints nearest(floats &eta1, floats &eta2, floats &phi1, floats &phi2, ints charge1, ints charge2)
{
  ints out = {-1, -1, -1, -1};
  float mindr = FLT_MAX;
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = i+1; j < 3; j++)
    {
      auto dr = ROOT::VecOps::DeltaR(eta1[i], eta2[j], phi1[i], phi2[j]);
      if (dr < mindr and charge1[i] != charge2[j])//dr > FLT_EPSILON and 
      {
        // nearest muons of opposite charge: {i, j}
        // let k = 3-i-j be the index of the third muon
        if (charge1[i] != charge2[3-i-j])
        {
            out = {i, j, i, 3-i-j};
        }
        else
        {
            out = {i, j, j, 3-i-j};
        }
        //out = {i, j};
        mindr = dr;
      }
    }
  }
  return out;
}

ints nearestSameFlavor(floats &eta1, floats &eta2, floats &phi1, floats &phi2, ints charge1, ints charge2, ints pdgId1, ints pdgId2)
{
  ints out = {-1, -1};
  float mindr = FLT_MAX;
  for (unsigned int i = 0; i < 3; i++)
  {
    for (unsigned int j = i+1; j < 3; j++)
    {
      auto dr = ROOT::VecOps::DeltaR(eta1[i], eta2[j], phi1[i], phi2[j]);
      if (dr < mindr and charge1[i] != charge2[j] and abs(pdgId1[i]) == abs(pdgId2[j]))
      {
        out = {i, j};
        mindr = dr;
      }
    }
  }
  return out;
}

