/*
 * BaseAnalyser.cpp
 *
 *  Created on: May 6, 2022
 *      Author: suyong
 *      Developper: cdozen
 */

#include "Math/GenVector/VectorUtil.h"
#include "BaseAnalyser.h"
#include "utility.h"

BaseAnalyser::BaseAnalyser(TTree *t, std::string outfilename)
:NanoAODAnalyzerrdframe(t, outfilename)
{
    //initiliaze the HLT names in your analyzer class
    HLT2018Names= {"HLT_PFHT380_SixJet32_DoubleBTagCSV_p075",
                    "HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0",
                    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
                    "HLT_PFJet550","HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5",
                    "HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5"};
	HLT2016Names= {"Name1","Name2"};
	HLT2017Names= {"Name1","Name2"};
}

// Define your cuts here
void BaseAnalyser::defineCuts()
{
	if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
	auto Nentry = _rlm.Count();
	// This is how you can express a range of the first 100 entries
    //_rlm = _rlm.Range(0, 100);
    //auto Nentry_100 = _rlm.Count();
	std::cout<< "-------------------------------------------------------------------" << std::endl;
    cout << "Usage of ranges:\n"
        << " - All entries: " << *Nentry << endl;
		//<< " - Entries from 0 to 100: " << *Nentry_100 << endl;
	std::cout<< "-------------------------------------------------------------------" << std::endl;

	//MinimalSelection to filter events
    addCuts("nMuon > 2 && nElectron > 0", "0");
    addCuts("nJet>0","00");
	addCuts(setHLT(),"1"); //HLT cut buy checking HLT names in the root file

}
//===============================Find Good Electrons===========================================//
//: Define Good Electrons in rdata frame
//=============================================================================================//

void BaseAnalyser::selectElectrons()
{
    cout << "select good electrons" << endl;
    if (debug){
    std::cout<< "================================//=================================" << std::endl;
    std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
    std::cout<< "================================//=================================" << std::endl;
    }
    //string veto_electron_str = Form ("Electron_pt>%f && abs(Electron_eta)<%f && %s",cutdb->Md_Electron_pT,cutdb->Md_Electron_eta, ElectronID(2).c_str());
    //_rlm = _rlm.Define("goodElectrons", veto_electron_str.c_str());

    _rlm = _rlm.Define("goodElectrons", ElectronID(2)); //without pt-eta cuts
	//_rlm = _rlm.Define("goodElectrons", "Electron_pt>30.0 && abs(Electron_eta)<2.4 && Electron_pfRelIso03_all<0.15");
    _rlm = _rlm.Define("goodElectrons_pt", "Electron_pt[goodElectrons]")
                .Define("goodElectrons_eta", "Electron_eta[goodElectrons]")
                .Define("goodElectrons_phi", "Electron_phi[goodElectrons]")
                .Define("goodElectrons_mass", "Electron_mass[goodElectrons]")
                .Define("goodElectrons_idx", ::good_idx, {"goodElectrons"})
                .Define("NgoodElectrons", "int(goodElectrons_pt.size())");

    //-------------------------------------------------------
    //generate electron 4vector from selected good Electrons
    //-------------------------------------------------------
    _rlm = _rlm.Define("goodElectron_4Vecs", ::generate_4vec, {"goodElectrons_pt", "goodElectrons_eta", "goodElectrons_phi", "goodElectrons_mass"});

}
//===============================Find Good Muons===============================================//
//: Define Good Muons in rdata frame
//=============================================================================================//
void BaseAnalyser::selectMuons()
{

    cout << "select good muons" << endl;
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    //_rlm = _rlm.Define("goodMuons", "Muon_pt>30.0 && abs(Muon_eta)<2.4 && Muon_pfRelIso03_all<0.15");
    _rlm = _rlm.Define("goodMuons", MuonID(4));
    _rlm = _rlm.Define("goodMuons_pt", "Muon_pt[goodMuons]") 
                .Define("goodMuons_eta", "Muon_eta[goodMuons]")
                .Define("goodMuons_phi", "Muon_phi[goodMuons]")
                .Define("goodMuons_mass", "Muon_mass[goodMuons]")
                .Define("goodMuons_charge", "Muon_charge[goodMuons]")
                .Define("goodMuons_idx", ::good_idx, {"goodMuons"})
                .Define("NgoodMuons", "int(goodMuons_pt.size())");

    //-------------------------------------------------------
    //generate muon 4vector from selected good Muons
    //-------------------------------------------------------
    _rlm = _rlm.Define("goodMuons_4vecs", ::generate_4vec, {"goodMuons_pt", "goodMuons_eta", "goodMuons_phi", "goodMuons_mass"});
    //_rlm = _rlm.Define("goodMuons_4vecs", ::generate_4vec, {"goodMuons_pt", "goodMuons_eta", "goodMuons_phi", "goodMuons_mass"});

}


//=================================Select Jets=================================================//
//check the twiki page :    https://twiki.cern.ch/twiki/bin/view/CMS/JetID
//to find jetId working points for the purpose of  your analysis.
    //jetId==2 means: pass tight ID, fail tightLepVeto
    //jetId==6 means: pass tight ID and tightLepVeto ID.
//=============================================================================================//
void BaseAnalyser::selectJets()
{

    cout << "select good jets" << endl;
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

    //string jetcuts_str = Form ("Jet_pt>%f && abs(Jet_eta)<%f && Jet_jetId >= %d",cutdb->Md_jet_pT,cutdb->Md_jet_eta,6 );
    //_rlm = _rlm.Define("goodJets", jetcuts_str);

   // _rlm = _rlm.Define("goodJets", "Jet_pt>30.0 && abs(Jet_eta)<2.4 && Jet_jetId >= 6");// 
    _rlm = _rlm.Define("goodJets", JetID(6)); //without pt-eta cuts
    _rlm = _rlm.Define("goodJets_pt", "Jet_pt[goodJets]")
                .Define("goodJets_eta", "Jet_eta[goodJets]")
                .Define("goodJets_phi", "Jet_phi[goodJets]")
                .Define("goodJets_mass", "Jet_mass[goodJets]")
                .Define("goodJets_idx", ::good_idx, {"goodJets"})
                .Define("goodJets_btag", "Jet_btagCSVV2[goodJets]")
                //.Define("Selected_jetbtag", "Jet_btagDeepFlavB[goodJets]") //select bjets within goodjets
                .Define("NgoodJets", "int(goodJets_pt.size())")
                .Define("goodJets_4vecs", ::generate_4vec, {"goodJets_pt", "goodJets_eta", "goodJets_phi", "goodJets_mass"});
                //.Define("goodJets_4vecs", ::generate_4vec, {"goodJets_pt", "goodJets_eta", "goodJets_phi", "goodJets_mass"});


    _rlm = _rlm.Define("btagcuts", "goodJets_btag>0.8")
			.Define("good_bjetpt", "goodJets_pt[btagcuts]")
			.Define("good_bjeteta", "goodJets_eta[btagcuts]")
			.Define("good_bjetphi", "goodJets_phi[btagcuts]")
			.Define("good_bjetmass", "goodJets_mass[btagcuts]")
			.Define("Ngood_bjets", "int(good_bjetpt.size())")
			.Define("good_bjet4vecs", ::generate_4vec, {"good_bjetpt", "good_bjeteta", "good_bjetphi", "good_bjetmass"});



}

//=============================define variables==================================================//
void BaseAnalyser::defineMoreVars()
{
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

	//addVar({"good_electron1pt", "goodElectrons_pt[0]", ""});
    addVar({"good_muon1pt", "goodMuons_pt[0]", ""});

    //selected jet candidates
    addVar({"good_jet1pt", "(goodJets_pt.size()>0) ? goodJets_pt[0] : -1", ""});
    //addVar({"Selected_jet1pt", "Selected_jetpt[0]", ""});
    addVar({"good_jet1eta", "goodJets_eta[0]", ""});
    addVar({"good_jet1mass", "goodJets_mass[0]", ""});

    //================================Store variables in tree=======================================//
    // define variables that you want to store
    //==============================================================================================//
    addVartoStore("run");
    addVartoStore("luminosityBlock");
    addVartoStore("event");
    addVartoStore("evWeight.*");

    //electron
    addVartoStore("nElectron");
    addVartoStore("Electron_charge");
    addVartoStore("NgoodElectrons");

    //muon
    addVartoStore("nMuon");
    addVartoStore("Muon_charge");
    addVartoStore("Muon_mass");
    addVartoStore("NgoodMuons");

    //jet
    addVartoStore("nJet");
    addVartoStore("Jet_pt");
    addVartoStore("NgoodJets");
    addVartoStore("goodJets_pt");

}
void BaseAnalyser::bookHists()
{
    if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
    //add x and y axis title: add1DHist( {"hnevents", "hist_title; x_axis title; y_axis title", 2, -0.5, 1.5}, "one", "evWeight", "");
    add1DHist( {"hnevents", "Number of Events", 2, -0.5, 1.5}, "one", "evWeight", "");
    add1DHist( {"hNgoodElectrons", "NumberofGoodElectrons", 5, 0.0, 5.0}, "NgoodElectrons", "evWeight", "");

    //===minimum cut step:0 : defined in defineCuts() function =====//
    //add1DHist( {"hgoodelectron1_pt", "good electron1_pt; #electron p_{T}; Entries / after ", 18, -2.7, 2.7}, "good_electron1pt", "evWeight", "0");

    add1DHist( {"hNgoodMuons", "# of good Muons", 5, 0.0, 5.0}, "NgoodMuons", "evWeight", "");
    add1DHist( {"hgood_jet1pt", "Good Jet_1 pt after cut0" , 100, 0, 1000} , "good_jet1pt", "evWeight", "0");
}
void BaseAnalyser::setTree(TTree *t, std::string outfilename)
{
	if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }

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
//================================Selected Object Definitions====================================//

void BaseAnalyser::setupObjects()
{
	// Object selection will be defined in sequence.
	// Selected objects will be stored in new vectors.
	selectElectrons();
	selectMuons();
	selectJets();
	//removeOverlaps();

}
void BaseAnalyser::setupAnalysis()
{
	if (debug){
        std::cout<< "================================//=================================" << std::endl;
        std::cout<< "Line : "<< __LINE__ << " Function : " << __FUNCTION__ << std::endl;
        std::cout<< "================================//=================================" << std::endl;
    }
	// Event weight for data it's always one. For MC, it depends on the sign
 	//cout<<"year===="<< _year<< "==runtype=== " <<  _runtype <<endl;
	_rlm = _rlm.Define("one", "1.0");
    _rlm = _rlm.Define("evWeight","one");

    //for correction define evWeights as fallows
    /*if(_isData){
        _rlm = _rlm.Define("unitGenWeight","one")
                .Define("pugenWeight","one"); // if setupcorrection in processnanoad.py then don't define here. 
        _rlm = _rlm.Define("evWeight","one");
    }*/
	//if (_isData && !isDefined("evWeight"))
	//{
	//	_rlm = _rlm.Define("evWeight", [](){
	//			return 1.0;
	//		}, {} );
	//}

	defineCuts();
	defineMoreVars();
	bookHists();
	setupCuts_and_Hists();
	setupTree();
}

