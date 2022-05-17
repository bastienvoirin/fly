/*
 * NanoAODAnalyzerrdframe.cpp
 *
 *  Created on: Sep 30, 2018
 *      Author: suyong
 */

#include "NanoAODAnalyzerrdframe.h"
#include <iostream>
#include <algorithm>
#include <typeinfo>

#include "TCanvas.h"
#include "Math/GenVector/VectorUtil.h"
#include <vector>
#include <fstream>
#include "utility.h"
#include <regex>
#include "ROOT/RDFHelpers.hxx"

using namespace std;

NanoAODAnalyzerrdframe::NanoAODAnalyzerrdframe(TTree *atree, std::string outfilename)
:_rd(*atree), _jsonOK(false),_outfilename(outfilename)
	, _outrootfile(0), _rlm(_rd)
	, _rnt(&_rlm)
{
	_atree=atree;
	//cout<< " run year=====" << _year <<endl;



}

NanoAODAnalyzerrdframe::~NanoAODAnalyzerrdframe() {
	// TODO Auto-generated destructor stub
}

bool NanoAODAnalyzerrdframe::isDefined(string v)
{
	auto result = std::find(_originalvars.begin(), _originalvars.end(), v);
	if (result != _originalvars.end()) return true;
	else return false;
}

void NanoAODAnalyzerrdframe::setTree(TTree *t, std::string outfilename)
{
	_rd = ROOT::RDataFrame(*t);
	_rlm = RNode(_rd);
	_outfilename = outfilename;
	_hist1dinfovector.clear();
	_th1dhistos.clear();
	_varstostore.clear();
	_hist1dinfovector.clear();
	_selections.clear();

	this->setupAnalysis();
}


void NanoAODAnalyzerrdframe::setupCorrections(string goodjsonfname, string pufname, string putag, string btvfname, string btvtype, string jercfname, string jerctag, string jercunctag)
{
	if (_isData) _jsonOK = readgoodjson(goodjsonfname); // read golden json file

	if (!_isData) {
		// using correctionlib
		// btag corrections
		_correction_btag1 = correction::CorrectionSet::from_file(btvfname);
		_btvtype = btvtype;
		assert(_correction_btag1->validate());

		// pile up weights
		_correction_pu = correction::CorrectionSet::from_file(pufname);
		assert(_correction_pu->validate());
		_putag = putag;
		auto punominal = [this](float x) { return pucorrection(_correction_pu, _putag, "nominal", x); };
		auto puplus = [this](float x) { return pucorrection(_correction_pu, _putag, "up", x); };
		auto puminus = [this](float x) { return pucorrection(_correction_pu, _putag, "down", x); };

		if (!isDefined("puWeight")) _rlm = _rlm.Define("puWeight", punominal, {"Pileup_nTrueInt"});
		if (!isDefined("puWeight_plus")) _rlm = _rlm.Define("puWeight_plus", puplus, {"Pileup_nTrueInt"});
		if (!isDefined("puWeight_minus")) _rlm = _rlm.Define("puWeight_minus", puminus, {"Pileup_nTrueInt"});


		if (!isDefined("pugenWeight"))
		{
			_rlm = _rlm.Define("pugenWeight", [this](float x, float y){
					return (x > 0 ? 1.0 : -1.0) *y;
				}, {"genWeight", "puWeight"});
		}
	}
	_jerctag = jerctag;
	_jercunctag = jercunctag;

	setupJetMETCorrection(jercfname, _jerctag);
	applyJetMETCorrections();
}

void NanoAODAnalyzerrdframe::setupObjects()
{
	// Object selection will be defined in sequence.
	// Selected objects will be stored in new vectors.
	//selectJets();
	removeOverlaps();
	selectFatJets();
}

void NanoAODAnalyzerrdframe::setupAnalysis()
{
	// Event weight for data it's always one. For MC, it depends on the sign
 	//cout<<"year===="<< _year<< "==runtype=== " <<  _runtype <<endl;
	_rlm = _rlm.Define("one", "1.0");
	if (_isData && !isDefined("evWeight"))
	{
		_rlm = _rlm.Define("evWeight", [](){
				return 1.0;
			}, {} );
	}
	/*if (_isData && !isDefined("pugenWeight"))
	{
		_rlm = _rlm.Define("pugenWeight", [](){
				return 1.0;
			}, {} );
	}*/


	//defineMoreVars();
	//defineCuts();
	//bookHists();
	setupCuts_and_Hists();
	setupTree();
}


bool NanoAODAnalyzerrdframe::readgoodjson(string goodjsonfname)
{
	auto isgoodjsonevent = [this](unsigned int runnumber, unsigned int lumisection)
		{
			auto key = std::to_string(runnumber).c_str();

			bool goodeventflag = false;


			if (jsonroot.contains(key))
			{
				for (auto &v: jsonroot[key])
				{
					if (v[0]<=lumisection && lumisection <=v[1]) goodeventflag = true;
				}
			}
			return goodeventflag;
		};

	if (goodjsonfname != "")
	{
		std::ifstream jsoninfile;
		jsoninfile.open(goodjsonfname);

		if (jsoninfile.good())
		{
			//using rapidjson
			//rapidjson::IStreamWrapper s(jsoninfile);
			//jsonroot.ParseStream(s);

			//using jsoncpp
			jsoninfile >> jsonroot;
			_rlm = _rlm.Define("goodjsonevent", isgoodjsonevent, {"run", "luminosityBlock"}).Filter("goodjsonevent");
			_jsonOK = true;
			return true;
		}
		else
		{
			cout << "Problem reading json file " << goodjsonfname << endl;
			return false;
		}
	}
	else
	{
		cout << "no JSON file given" << endl;
		return true;
	}
}

void NanoAODAnalyzerrdframe::setupJetMETCorrection(string fname, string jettag)
{
	// read from file 
	_correction_jerc = correction::CorrectionSet::from_file(fname);
	assert(_correction_jerc->validate());
	// correction type(jobconfiganalysis.py)
	_jetCorrector = _correction_jerc->compound().at(jettag);
	_jetCorrectionUnc = _correction_jerc->at(_jercunctag);
}


// Adapted from https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetRecalib.py
// and https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/JetRecalibrator.py
void NanoAODAnalyzerrdframe::applyJetMETCorrections()
{

	auto appcorrlambdaf = [this](floats jetpts, floats jetetas, floats jetAreas, floats jetrawf, float rho)->floats
	{
		floats corrfactors;
		corrfactors.reserve(jetpts.size());
		for (auto i =0; i<jetpts.size(); i++)
		{
			float rawjetpt = jetpts[i]*(1.0-jetrawf[i]);
			float corrfactor = _jetCorrector->evaluate({jetAreas[i], jetetas[i], rawjetpt, rho});
			corrfactors.emplace_back(rawjetpt * corrfactor);

		}
		return corrfactors;
	};

	auto jecuncertaintylambdaf= [this](floats jetpts, floats jetetas, floats jetAreas, floats jetrawf, float rho)->floats
		{
			floats uncertainties;
			uncertainties.reserve(jetpts.size());
			for (auto i =0; i<jetpts.size(); i++)
			{
				float rawjetpt = jetpts[i]*(1.0-jetrawf[i]);
				float corrfactor = _jetCorrector->evaluate({jetAreas[i], jetetas[i], rawjetpt, rho});
				float unc = _jetCorrectionUnc->evaluate({corrfactor*rawjetpt, jetetas[i]});
				uncertainties.emplace_back(unc);

			}
			return uncertainties;
		};

	auto metcorrlambdaf = [](float met, float metphi, floats jetptsbefore, floats jetptsafter, floats jetphis)->float
	{
		auto metx = met * cos(metphi);
		auto mety = met * sin(metphi);
		for (auto i=0; i<jetphis.size(); i++)
		{
			if (jetptsafter[i]>15.0)
			{
				metx -= (jetptsafter[i] - jetptsbefore[i])*cos(jetphis[i]);
				mety -= (jetptsafter[i] - jetptsbefore[i])*sin(jetphis[i]);
			}
		}
		return float(sqrt(metx*metx + mety*mety));
	};

	auto metphicorrlambdaf = [](float met, float metphi, floats jetptsbefore, floats jetptsafter, floats jetphis)->float
	{
		auto metx = met * cos(metphi);
		auto mety = met * sin(metphi);
		for (auto i=0; i<jetphis.size(); i++)
		{
			if (jetptsafter[i]>15.0)
			{
				metx -= (jetptsafter[i] - jetptsbefore[i])*cos(jetphis[i]);
				mety -= (jetptsafter[i] - jetptsbefore[i])*sin(jetphis[i]);
			}
		}
		return float(atan2(mety, metx));
	};

	if (_jetCorrector != 0)
	{
		_rlm = _rlm.Define("Jet_pt_corr", appcorrlambdaf, {"Jet_pt", "Jet_eta", "Jet_area", "Jet_rawFactor", "fixedGridRhoFastjetAll"});
		_rlm = _rlm.Define("Jet_pt_relerror", jecuncertaintylambdaf, {"Jet_pt", "Jet_eta", "Jet_area", "Jet_rawFactor", "fixedGridRhoFastjetAll"});
		_rlm = _rlm.Define("Jet_pt_corr_up", "Jet_pt_corr*(1.0f + Jet_pt_relerror)");
		_rlm = _rlm.Define("Jet_pt_corr_down", "Jet_pt_corr*(1.0f - Jet_pt_relerror)");
		_rlm = _rlm.Define("MET_pt_corr", metcorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr", "Jet_phi"});
		_rlm = _rlm.Define("MET_phi_corr", metphicorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr", "Jet_phi"});
		_rlm = _rlm.Define("MET_pt_corr_up", metcorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr_up", "Jet_phi"});
		_rlm = _rlm.Define("MET_phi_corr_up", metphicorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr_up", "Jet_phi"});
		_rlm = _rlm.Define("MET_pt_corr_down", metcorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr_down", "Jet_phi"});
		_rlm = _rlm.Define("MET_phi_corr_down", metphicorrlambdaf, {"MET_pt", "MET_phi", "Jet_pt", "Jet_pt_corr_down", "Jet_phi"});
	}

}

/*void NanoAODAnalyzerrdframe::selectJets()
{
	// apparently size() returns long int, which ROOT doesn't recognized for branch types
	// , so it must be cast into int if you want to save them later into a TTree
	_rlm = _rlm.Define("jetcuts", "Jet_pt>30.0 && abs(Jet_eta)<2.4 && Jet_jetId>0")
			.Define("Sel_jetpt", "Jet_pt[jetcuts]")
			.Define("Sel_jeteta", "Jet_eta[jetcuts]")
			.Define("Sel_jetphi", "Jet_phi[jetcuts]")
			.Define("Sel_jetmass", "Jet_mass[jetcuts]")
			.Define("Sel_jetbtag", "Jet_btagCSVV2[jetcuts]")
			.Define("njetspass", "int(Sel_jetpt.size())")
			.Define("jet4vecs", ::generate_4vec, {"Sel_jetpt", "Sel_jeteta", "Sel_jetphi", "Sel_jetmass"});


	_rlm = _rlm.Define("btagcuts", "Sel_jetbtag>0.8")
			.Define("Sel_bjetpt", "Sel_jetpt[btagcuts]")
			.Define("Sel_bjeteta", "Sel_jeteta[btagcuts]")
			.Define("Sel_bjetphi", "Sel_jetphi[btagcuts]")
			.Define("Sel_bjetm", "Sel_jetmass[btagcuts]")
			.Define("bnjetspass", "int(Sel_bjetpt.size())")
			.Define("bjet4vecs", ::generate_4vec, {"Sel_bjetpt", "Sel_bjeteta", "Sel_bjetphi", "Sel_bjetm"});
			;

			// calculate event weight for MC only
	if (!_isData && !isDefined("evWeight"))
	{
		_rlm = _rlm.Define("jetcutsforsf", "Jet_pt>25.0 && abs(Jet_eta)<2.4 ")
				.Define("Sel_jetforsfpt", "Jet_pt[jetcutsforsf]")
				.Define("Sel_jetforsfeta", "Jet_eta[jetcutsforsf]")
				.Define("Sel_jetforsfhad", "Jet_hadronFlavour[jetcutsforsf]")
				.Define("Sel_jetcsvv2", "Jet_btagCSVV2[jetcutsforsf]")
				.Define("Sel_jetdeepb", "Jet_btagDeepB[jetcutsforsf]");

		auto btvcentral = [this](floats &pts, floats &etas, ints &hadflav, floats &btags)->floats
		{
			return ::btvcorrection(_correction_btag1, _btvtype, "central", pts, etas, hadflav, btags); // defined in utility.cpp
		};


		_rlm = _rlm.Define("Sel_jet_deepJet_shape_central", btvcentral, {"Sel_jetforsfpt", "Sel_jetforsfeta", "Sel_jetforsfhad", "Sel_jetdeepb"});
		//_rlm = _rlm.Define("Sel_jet_deepJet_shape_central",[this](floats &pts, floats &etas, ints &hadflav, floats &btags){return ::btvcorrection(_correction_btag1, "deepJet_shape", "central", pts, etas, hadflav, btags);}
			//				, {"Sel_jetforsfpt", "Sel_jetforsfeta", "Sel_jetforsfhad", "Sel_jetdeepb"});

		// function to calculate event weight for MC events based on DeepCSV algorithm
		auto btagweightgenerator3= [this](floats &pts, floats &etas, ints &hadflav, floats &btags)->float
		{
			double bweight=1.0;

			for (auto i=0; i<pts.size(); i++)
			{
				double w = _correction_btag1->at(_btvtype)->evaluate({"central", int(hadflav[i]), fabs(float(etas[i])), float(pts[i]), float(btags[i])});
				bweight *= w;
			}
			return bweight;
		};
		_rlm = _rlm.Define("btagWeight_DeepJetrecalc", btagweightgenerator3, {"Sel_jetforsfpt", "Sel_jetforsfeta", "Sel_jetforsfhad", "Sel_jetdeepb"});
		_rlm = _rlm.Define("evWeight", "pugenWeight * btagWeight_DeepJetrecalc");
	}

}*/

void NanoAODAnalyzerrdframe::removeOverlaps()
{
	// lambda function
	// for checking overlapped jets with electrons
	auto checkoverlap = [](FourVectorVec &seljets, FourVectorVec &selele)
		{
			doubles mindrlepton;
			//cout << "selected jets size" << seljets.size() << endl;
			//cout << "selected electrons size" << selele.size() << endl;

			for (auto ajet: seljets)
			{
				std::vector<double> drelejet(selele.size());
				for (auto alepton: selele)
				{
					auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
					drelejet.emplace_back(dr);
				}
				auto mindr = selele.size()==0 ? 6.0 : *std::min_element(drelejet.begin(), drelejet.end());
				mindrlepton.emplace_back(mindr);
			}

			return mindrlepton;
		};
	//cout << "overlap removal" << endl;
	_rlm = _rlm.Define("mindrlepton", checkoverlap, {"jet4vecs","ele4vecs"});

	//cout << "redefine cleaned jets" << endl;
	_rlm = _rlm.Define("overlapcheck", "mindrlepton>0.4");

	_rlm =	_rlm.Define("Sel2_jetpt", "Sel_jetpt[overlapcheck]")
		.Define("Sel2_jeteta", "Sel_jeteta[overlapcheck]")
		.Define("Sel2_jetphi", "Sel_jetphi[overlapcheck]")
		.Define("Sel2_jetmass", "Sel_jetmass[overlapcheck]")
		.Define("Sel2_jetbtag", "Sel_jetbtag[overlapcheck]")
		.Define("ncleanjetspass", "int(Sel2_jetpt.size())")
		.Define("cleanjet4vecs", ::generate_4vec, {"Sel2_jetpt", "Sel2_jeteta", "Sel2_jetphi", "Sel2_jetmass"})
		.Define("Sel2_jetHT", "Sum(Sel2_jetpt)")
		.Define("Sel2_jetweight", "std::vector<double>(ncleanjetspass, evWeight)"); //


	_rlm = _rlm.Define("btagcuts2", "Sel2_jetbtag>0.8")
			.Define("Sel2_bjetpt", "Sel2_jetpt[btagcuts2]")
			.Define("Sel2_bjeteta", "Sel2_jeteta[btagcuts2]")
			.Define("Sel2_bjetphi", "Sel2_jetphi[btagcuts2]")
			.Define("Sel2_bjetmass", "Sel2_jetmass[btagcuts2]")
			.Define("ncleanbjetspass", "int(Sel2_bjetpt.size())")
			.Define("Sel2_bjetHT", "Sum(Sel2_bjetpt)")
			.Define("cleanbjet4vecs", ::generate_4vec, {"Sel2_bjetpt", "Sel2_bjeteta", "Sel2_bjetphi", "Sel2_bjetmass"});

}

void NanoAODAnalyzerrdframe::selectFatJets()
{
	_rlm = _rlm.Define("fatjetcuts", "FatJet_pt>400.0 && abs(FatJet_eta)<2.4 && FatJet_tau1>0.0 && FatJet_tau2>0.0 && FatJet_tau3>0.0 && FatJet_tau3/FatJet_tau2<0.5")
				.Define("Sel_fatjetpt", "FatJet_pt[fatjetcuts]")
				.Define("Sel_fatjeteta", "FatJet_eta[fatjetcuts]")
				.Define("Sel_fatjetphi", "FatJet_phi[fatjetcuts]")
				.Define("Sel_fatjetmass", "FatJet_mass[fatjetcuts]")
				.Define("nfatjetspass", "int(Sel_fatjetpt.size())")
				.Define("Sel_fatjetweight", "std::vector<double>(nfatjetspass, evWeight)")
				.Define("Sel_fatjet4vecs", ::generate_4vec, {"Sel_fatjetpt", "Sel_fatjeteta", "Sel_fatjetphi", "Sel_fatjetmass"});
}


bool NanoAODAnalyzerrdframe::helper_1DHistCreator(std::string hname, std::string title, const int nbins, const double xlow, const double xhi, std::string rdfvar, std::string evWeight, RNode *anode)
{
	//cout << "1DHistCreator " << hname  << endl;

	RDF1DHist histojets = anode->Histo1D({hname.c_str(), title.c_str(), nbins, xlow, xhi}, rdfvar, evWeight); // Fill with weight given by evWeight
	_th1dhistos[hname] = histojets;
	//histojets.GetPtr()->Print("all");
	return true;
}


// Automatically loop to create
void NanoAODAnalyzerrdframe::setupCuts_and_Hists()
{
	cout << "setting up definitions, cuts, and histograms" <<endl;

	for ( auto &c : _varinfovector)
	{
		if (c.mincutstep.length()==0) _rlm = _rlm.Define(c.varname, c.vardefinition);
	}

	for (auto &x : _hist1dinfovector)
	{
		std::string hpost = "_nocut";

		if (x.mincutstep.length()==0)
		{
			helper_1DHistCreator(std::string(x.hmodel.fName)+hpost,  std::string(x.hmodel.fTitle)+hpost, x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.varname, x.weightname, &_rlm);
		}
	}

	_rnt.setRNode(&_rlm);

	for (auto acut : _cutinfovector)
	{
		std::string cutname = "cut"+ acut.idx;
		std::string hpost = "_"+cutname;
		RNode *r = _rnt.getParent(acut.idx)->getRNode();
		auto rnext = new RNode(r->Define(cutname, acut.cutdefinition));
		*rnext = rnext->Filter(cutname);

		for ( auto &c : _varinfovector)
		{
			if (acut.idx.compare(c.mincutstep)==0) *rnext = rnext->Define(c.varname, c.vardefinition);
		}
		for (auto &x : _hist1dinfovector)
		{
			if (acut.idx.compare(0, x.mincutstep.length(), x.mincutstep)==0)
			{
				helper_1DHistCreator(std::string(x.hmodel.fName)+hpost,  std::string(x.hmodel.fTitle)+hpost, x.hmodel.fNbinsX, x.hmodel.fXLow, x.hmodel.fXUp, x.varname, x.weightname, rnext);
			}
		}
		_rnt.addDaughter(rnext, acut.idx);

	}
}

void NanoAODAnalyzerrdframe::add1DHist(TH1DModel histdef, std::string variable, std::string weight,string mincutstep)
{
	_hist1dinfovector.push_back({histdef, variable, weight, mincutstep});
}


void NanoAODAnalyzerrdframe::drawHists(RNode t)
{
	cout << "processing" <<endl;
	t.Count();
}

void NanoAODAnalyzerrdframe::addVar(varinfo v)
{
	_varinfovector.push_back(v);
}

void NanoAODAnalyzerrdframe::addVartoStore(string varname)
{
	// varname is assumed to be a regular expression.
	// e.g. if varname is "Muon_eta" then "Muon_eta" will be stored
	// if varname=="Muon_.*", then any branch name that starts with "Muon_" string will
	// be saved
	_varstostore.push_back(varname);

}

void NanoAODAnalyzerrdframe::setupTree()
{
	vector<RNodeTree *> rntends;
	_rnt.getRNodeLeafs(rntends);
	for (auto arnt: rntends)
	{
		RNode *arnode = arnt->getRNode();
		string nodename = arnt->getIndex();
		vector<string> varforthistree;
		std::map<string, int> varused;

		for (auto varname: _varstostore)
		{
			bool foundmatch = false;
			std::regex b(varname);
			for (auto a: arnode->GetColumnNames())
			{
				if (std::regex_match(a, b) && varused[a]==0)
				{
					varforthistree.push_back(a);
					varused[a]++;
					foundmatch = true;
				}
			}
			if (!foundmatch)
			{
				cout << varname << " not found at "<< nodename << endl;
			}

		}
		_varstostorepertree[nodename]  = varforthistree;
	}

}

void NanoAODAnalyzerrdframe::addCuts(string cut, string idx)
{
	_cutinfovector.push_back({cut, idx});
}


void NanoAODAnalyzerrdframe::run(bool saveAll, string outtreename)
{


	vector<RNodeTree *> rntends;
	_rnt.getRNodeLeafs(rntends);

	//_rnt.Print();
	//cout << rntends.size() << endl;
   // std::cout<< "-------------------------------------------------------------------" << std::endl;
	// on master, regex_replace doesn't work somehow
	//std::regex rootextension("\\.root");

	for (auto arnt: rntends)
	{
		string nodename = arnt->getIndex();
		//string outname = std::regex_replace(_outfilename, rootextension, "_"+nodename+".root");
		string outname = _outfilename;
		// if producing many root files due to branched selection criteria,  each root file will get a different name
		if (rntends.size()>1) outname.replace(outname.find(".root"), 5, "_"+nodename+".root");
		_outrootfilenames.push_back(outname);
		RNode *arnode = arnt->getRNode();
		std::cout<< "-------------------------------------------------------------------" << std::endl;
                cout<<"cut : " ;
                cout << arnt->getIndex();
		//cout << ROOT::RDF::SaveGraph(_rlm) << endl;
		if (saveAll) {
			arnode->Snapshot(outtreename, outname);
		}
		else {
			// use the following if you want to store only a few variables
			//arnode->Snapshot(outtreename, outname, _varstostore);
            cout << " --writing branches" << endl;
			std::cout<< "-------------------------------------------------------------------" << std::endl;
			for (auto bname: _varstostorepertree[nodename])
			{
				cout << bname << endl;
			}
			arnode->Snapshot(outtreename, outname, _varstostorepertree[nodename]);
		}
		std::cout<< "-------------------------------------------------------------------" << std::endl;
		cout << "Creating output root file :  " << endl;
		cout << outname << " ";
		cout<<endl;
		std::cout<< "-------------------------------------------------------------------" << std::endl;
		_outrootfile = new TFile(outname.c_str(), "UPDATE");
		cout << "Writing histograms...   " << endl;
		std::cout<< "-------------------------------------------------------------------" << std::endl;
		for (auto &h : _th1dhistos)
		{
			if (h.second.GetPtr() != nullptr) {
				h.second.GetPtr()->Print();
				h.second.GetPtr()->Write();
			}
		}

		_outrootfile->Write(0, TObject::kOverwrite);
		_outrootfile->Close();
	}
    std::cout<< "-------------------------------------------------------------------" << std::endl;
    std::cout << "END...  :) " << std::endl; 

}

void NanoAODAnalyzerrdframe::setParams(int year, string runtype, int datatype)
{
    /*if(debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }*/
	_year=year;
	_runtype=runtype;
	_datatype=datatype;
	

	if(_year==2016) {
        cout << "Analysing through Run 2016" << endl;
    }else if(_year==2017) {
        cout << "Analysing through Run 2017" << endl;
    }else if(_year==2018){
        cout << "Analysing through Run 2018" << endl;
    }

	if(_runtype.find("UL") != std::string::npos){
        _isUL = true;
        cout << "Ultra Legacy Selected " << endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
    }else if(_runtype.find("ReReco") != std::string::npos){
        _isReReco = true;
        cout << " ReReco  Selected!" << endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
    }
    if (!_isUL && !_isReReco){
        std::cout<< "Default run version : UL or ReReco is not selected! "<< std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
    }


	if (_datatype==0){
		_isData = false;
        std::cout << " MC input files Selected!! "<<std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;

    }else if(_datatype==1){
        _isData = true;
        std::cout << " DATA input files Selected!!" <<std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
    }if(_datatype==-1){
		std::cout<< "Default root version :checking out gen branches! "<< std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;

		if (_atree->GetBranch("genWeight") == nullptr) {
			_isData = true;
			cout << "input file is DATA" <<endl;
		}
		else
		{
			_isData = false;
			cout << "input file is MC" <<endl;
		}
	}
	TObjArray *allbranches = _atree->GetListOfBranches();
	for (int i =0; i<allbranches->GetSize(); i++)
	{
		TBranch *abranch = dynamic_cast<TBranch *>(allbranches->At(i));
		if (abranch!= nullptr){
			//cout << abranch->GetName() << endl;
			_originalvars.push_back(abranch->GetName());
		}
	}


}
//Checking HLTs in the input root file
std::string NanoAODAnalyzerrdframe::setHLT(std::string str_HLT){
    if(debug){
    std::cout<< "================================//=================================" << std::endl;
    std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
    std::cout<< "================================//=================================" << std::endl;
    }
    if(str_HLT != "" ){

        bool ctrl_HLT =isDefined(str_HLT);
        std::string output;
        if(ctrl_HLT){
            output=str_HLT;
            std::cout<<"HLT : " <<  str_HLT  << " : SUCCESSFULLY FOUND!!"<< std::endl;
        }else{
            std::cout<<"HLT : " <<  str_HLT  << " : CAN NOT BE FOUND "<< std::endl;
            std::cout<< "Check HLT branches in the input root file!!" << std::endl;
            std::cout<< "EXITING PROGRAM!!" << std::endl;

            exit(1);
        }
        return output;

    }else{ // fill the HLT names in a vector according to each year
            std::vector<string> V_output;
            if(_year==2016){
                HLTGlobalNames=HLT2016Names;
            }else if (_year==2017){
                HLTGlobalNames=HLT2017Names;
            }else if(_year==2018){
                HLTGlobalNames=HLT2018Names;
            }

            //loop on HLTs
            for (size_t i = 0; i < HLTGlobalNames.size(); i++)
            {
                /* code */
                bool ctrl_HLT = isDefined(HLTGlobalNames[i]);
                if(ctrl_HLT){
                    V_output.push_back(HLTGlobalNames[i]);
                }

            }
            std::string output_HLT;
            if(!V_output.empty()){
                for (size_t i = 0; i < V_output.size() ; i++)
                {
                    if(i!=V_output.size()-1){
                    output_HLT += V_output[i] + "==1 || " ;
                    }else{
                        output_HLT += V_output[i] + "==1 " ;
                    }
                }
            }else{
                std::cout<< " Not matched with any HLT Triggers! Please check the HLT Names in the inputfile " << std::endl;
                std::cout<< "EXITING PROGRAM!!" << std::endl;
                exit(1);
            }
            std::cout<< " HLT names =  " <<  output_HLT  << std::endl;
            return output_HLT;

    }

}
//control all branch names using in the addCuts function
std::string NanoAODAnalyzerrdframe::ctrlBranchName(std::string str_Branch){

    if(debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    bool ctrl_Branch =isDefined(str_Branch);
    std::string output;
    if(ctrl_Branch){
        output=str_Branch;
    }else{
        std::cout<<"Branch : " <<  str_Branch  << " : CAN NOT BE FOUND "<< std::endl;
        std::cout<< "Check your branches in the input root file!!" << std::endl;
        std::cout<< "EXITING PROGRAM!!" << std::endl;

        exit(1);
    }
    return output;
}



//cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
std::string NanoAODAnalyzerrdframe::ElectronID(int cutbasedID){

    if(debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
    //double Electron_eta;
    //double Electron_pt;
    if (cutbasedID==1)std::cout<< " VETO Electron ID requested    == " << cutbasedID <<std::endl;
    if (cutbasedID==2)std::cout<< " LOOSE Electron ID requested   == " << cutbasedID <<std::endl;
    if (cutbasedID==3)std::cout<< " MEDIUM Electron ID requested  == " << cutbasedID <<std::endl;
    if (cutbasedID==4)std::cout<< " TIGHT Electron ID requested   == " << cutbasedID <<std::endl;
    std::cout<< "-------------------------------------------------------------------" << std::endl;

    if (cutbasedID<0 || cutbasedID>4){
        std::cout<< "ERROR!! Wrong Electron ID requested  == " << cutbasedID << "!! Can't be applied" <<std::endl;
        std::cout<< "Please select ElectronID from 1 to 4 " <<std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
        std::cout<< "EXITING PROGRAM!!" << std::endl;
        exit(1);
    }
	/*if(_year==2018 && _isUL){
    	if(cutbasedID==2 ){
        	Electron_eta=2.5;
    		Electron_pt=10;
    	}else if (cutbasedID==3){
        	Electron_eta=2.4;
        	Electron_pt=10;
		}
	}*/

//Rdataframe look for the variables in the intput Ttree..
std::string output = Form("Electron_cutBased == %d ",cutbasedID);
//std::string output = Form("Electron_cutBased == %d &&  abs(Electron_eta)<%f && Electron_pt<%f",cutbasedID,  Electron_eta, Electron_pt);

return output;
}

std::string NanoAODAnalyzerrdframe::MuonID(int cutbasedID){
    
    if(debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
    
    //  cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
    if (cutbasedID==1)std::cout<< " Veto Muon ID requested   == " << cutbasedID <<std::endl;
    if (cutbasedID==2)std::cout<< " LOOSE Muon ID requested  == " << cutbasedID <<std::endl;
    if (cutbasedID==3)std::cout<< " MEDIUM Muon ID requested == " << cutbasedID <<std::endl;
    if (cutbasedID==4)std::cout<< " TIGHT Muon ID requested as == " << cutbasedID <<std::endl;
    std::cout<< "-------------------------------------------------------------------" << std::endl;

    if (cutbasedID<1 || cutbasedID>4){
        std::cout<< "ERROR!! Wrong Muon ID requested  == " << cutbasedID << "!! Can't be applied" <<std::endl;
        std::cout<< "Please select Muon ID from 2 to 4 " <<std::endl;
        std::cout<< "-------------------------------------------------------------------" << std::endl;
        std::cout<< "EXITING PROGRAM!!" << std::endl;

        exit(1);
    }

    string Muon_cutBased_ID;
    
    if (cutbasedID==1){
        Muon_cutBased_ID = "Muon_looseId";
        std::cout<< " VETO Muon ID requested == " << cutbasedID <<", but it doesn't exist in the nanoAOD branches. It is moved to loose MuonID. " << cutbasedID <<std::endl;
    }
    if (cutbasedID == 2){
        Muon_cutBased_ID = "Muon_looseId";
    
    }else if(cutbasedID == 3){
        Muon_cutBased_ID = "Muon_mediumId";

    }else if(cutbasedID == 4){
        Muon_cutBased_ID = "Muon_tightId";
    }
    string output;
    output = Form ("%s==true",Muon_cutBased_ID.c_str());
    return output;
}


std::string NanoAODAnalyzerrdframe::JetID(int cutbasedID){

    if(debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    if(cutbasedID<1 || cutbasedID>7){
        std::cout<< "Error Wrong JET ID requested == " << cutbasedID << "!! Can't be applied" <<std::endl;
        std::cout<< "Please select number from 1 to 7 " <<std::endl;

    }else{
	    std::cout<< " JET ID requested  == " << cutbasedID <<std::endl;
    }

    string output;
    output = Form ("Jet_jetId==%d",cutbasedID);
    return output;
}
