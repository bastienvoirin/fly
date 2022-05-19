/// @file TopHiggsTrileptonAnalyser.cpp
/// @brief
/// @details
/// @author
/// @version
/// @date 2022
/// @copyright

#include "TopHiggsTrileptonAnalyser.h"
#include "utility.h"
#include "Math/GenVector/VectorUtil.h"
#include <ROOT/RDataFrame.hxx>
#include <ROOT/TProcessExecutor.hxx>
#include <TStopwatch.h>

// Particle Data Group identification number of particles and antiparticles
enum PdgId {
    // Antiparticles
    ANTI_D_QUARK          =  -1,
    ANTI_U_QUARK          =  -2,
    ANTI_S_QUARK          =  -3,
    ANTI_C_QUARK          =  -4,
    ANTI_B_QUARK          =  -5,
    ANTI_T_QUARK          =  -6,
    ANTIELECTRON          = -11,
    ELECTRON_ANTINEUTRINO = -12,
    ANTIMUON              = -13,
    MUON_ANTINEUTRINO     = -14,
    ANTITAU               = -15,
    TAU_ANTINEUTRINO      = -16,
    WMINUS                = -24,

    // Particles
    D_QUARK               =   1,
    U_QUARK               =   2,
    S_QUARK               =   3,
    C_QUARK               =   4,
    B_QUARK               =   5,
    T_QUARK               =   6,
    ELECTRON              =  11,
    ELECTRON_NEUTRINO     =  12,
    MUON                  =  13,
    MUON_NEUTRINO         =  14,
    TAU                   =  15,
    TAU_NEUTRINO          =  16,
    PHOTON                =  22,
    Z0                    =  23,
    WPLUS                 =  24,
    HIGGS                 =  25
};

TopHiggsTrileptonAnalyser::TopHiggsTrileptonAnalyser(TTree *t, std::string outfilename) 
: NanoAODAnalyzerrdframe(t,outfilename)
{ // Initialise the HLT names in the analyser class
    HLT2018Names = {
        "HLT_PFHT380_SixJet32_DoubleBTagCSV_p075",
        "HLT_PFHT300PT30_QuadPFJet_75_60_45_40_TriplePFBTagCSV_3p0",
        "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5",
        "HLT_PFJet550", "HLT_PFHT400_FivePFJet_100_100_60_30_30_DoublePFBTagDeepCSV_4p5",
        "HLT_PFHT400_FivePFJet_120_120_60_30_30_DoublePFBTagDeepCSV_4p5"
    };
    HLT2017Names = {"Name1", "Name2"};
    HLT2016Names = {"Name1", "Name2"};
}

/// @brief Event selection
/// @details Cuts to be applied in order
void TopHiggsTrileptonAnalyser::defineCuts()
{
    printDebugInfo(__LINE__, __FUNCTION__);

    // Count some entries using ranges
    auto Nentry = _rlm.Count();
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    cout << "Usage of ranges:\n"
         << "  All entries: " << *Nentry << endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;


    // Basic cuts
    // 3-lepton final state (T' => t, H^0 => WW, Wb => (lν, lν), (lν, jet))
    addCuts("Selected_muon_number == 3 && abs(Selected_muon_charge_sum) == 1", "0");
    // At least 1 b-tagged jet (t => W, b => lν, jet)
    addCuts("Selected_bjet_number >= 1", "00");

    //addCuts("Selected_muon_sum_two_muons_pt > 160", "00");
    //addCuts("Selected_muon_leading_pt > 80", "000");
    //addCuts("Selected_muon_subleading_pt > 40", "0000");
    //addCuts("Selected_muon_deltaR > 1.8", "00000");
}

/// @brief Find good electrons
/// @details Define good electrons in rdata frame
void TopHiggsTrileptonAnalyser::selectElectrons()
{
    cout << "Selecting good electrons" << endl;
    printDebugInfo(__LINE__, __FUNCTION__);

    _rlm = _rlm.Define("goodElecID", ElectronID(2));
    _rlm = _rlm.Define("goodElectrons", "goodElecID && Electron_pt > 20 && abs(Electron_eta) < 2.5")
               .Define("Selected_electron_pt", "Electron_pt[goodElectrons]")
               .Define("Selected_electron_sum_all_electrons_pt", "Sum(Selected_electron_pt)")
               .Define("Selected_electron_eta", "Electron_eta[goodElectrons]")
               .Define("Selected_electron_phi", "Electron_phi[goodElectrons]")
               .Define("Selected_electron_mass", "Electron_mass[goodElectrons]")
               .Define("Selected_electron_idx", ::good_idx, {"goodElectrons"})
               .Define("Selected_electron_number", "int(Selected_electron_pt.size())");

    // Generate electron 4-vectors
    _rlm = _rlm.Define("elec4vecs", ::generate_4vec, {"Selected_electron_pt",
                                                      "Selected_electron_eta",
                                                      "Selected_electron_phi",
                                                      "Selected_electron_mass"});
}

/// @brief Find good muons
/// @details Define good muons in rdata frame
void TopHiggsTrileptonAnalyser::selectMuons()
{
    cout << "Selecting good muons" << endl;
    printDebugInfo(__LINE__, __FUNCTION__);

    // Define good muons
    _rlm = _rlm.Define("goodMuonID", MuonID(4));
    _rlm = _rlm.Define("goodMuons", "goodMuonID && Muon_pt > 20 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.05")
               .Define("Selected_muon_pt", "Muon_pt[goodMuons]")
               .Define("Selected_muon_leading_pt", "Selected_muon_pt[0]")
               .Define("Selected_muon_subleading_pt", "Selected_muon_pt[1]")
               .Define("Selected_muon_subsubleading_pt", "Selected_muon_pt[2]")
               .Define("Selected_muon_sum_two_muons_pt", "Selected_muon_pt[0] + Selected_muon_pt[1]")
               .Define("Selected_muon_sum_all_muons_pt", "Sum(Selected_muon_pt)")
               .Define("Selected_muon_eta", "Muon_eta[goodMuons]")
               .Define("Selected_muon_phi", "Muon_phi[goodMuons]")
               .Define("Selected_muon_mass", "Muon_mass[goodMuons]")
               .Define("Selected_muon_charge", "Muon_charge[goodMuons]")
               .Define("Selected_muon_charge_sum", "Selected_muon_charge[0] + Selected_muon_charge[1] + Selected_muon_charge[2]")
               .Define("Selected_muon_number", "int(Selected_muon_pt.size())")
               .Define("Selected_muon_deltaeta_01", "abs(Selected_muon_eta[0] - Selected_muon_eta[1])")
               .Define("Selected_muon_deltaeta_02", "abs(Selected_muon_eta[0] - Selected_muon_eta[2])")
               .Define("Selected_muon_deltaeta_12", "abs(Selected_muon_eta[1] - Selected_muon_eta[2])")
               .Define("Selected_muon_deltaphi_01", "abs(ROOT::VecOps::DeltaPhi(Selected_muon_phi[0], Selected_muon_phi[1]))")
               .Define("Selected_muon_deltaphi_02", "abs(ROOT::VecOps::DeltaPhi(Selected_muon_phi[0], Selected_muon_phi[2]))")
               .Define("Selected_muon_deltaphi_12", "abs(ROOT::VecOps::DeltaPhi(Selected_muon_phi[1], Selected_muon_phi[2]))")
               .Define("Selected_muon_deltaR_01", "ROOT::VecOps::DeltaR(Selected_muon_eta[0], Selected_muon_eta[1], Selected_muon_phi[0], Selected_muon_phi[1])")
               .Define("Selected_muon_deltaR_02", "ROOT::VecOps::DeltaR(Selected_muon_eta[0], Selected_muon_eta[2], Selected_muon_phi[0], Selected_muon_phi[2])")
               .Define("Selected_muon_deltaR_12", "ROOT::VecOps::DeltaR(Selected_muon_eta[1], Selected_muon_eta[2], Selected_muon_phi[1], Selected_muon_phi[2])")
               //.Define("Selected_muon_miniPFRelIso_all", "Muon_miniPFRelIso_all[goodMuons]")
               //.Define("Selected_muon_miniPFRelIso_chg", "Muon_miniPFRelIso_chg[goodMuons]")
               //.Define("Selected_muon_pfRelIso03_all", "Muon_pfRelIso03_all[goodMuons]")
               //.Define("Selected_muon_pfRelIso03_chg", "Muon_pfRelIso03_chg[goodMuons]")
               //.Define("Selected_muon_pfRelIso04_all", "Muon_pfRelIso04_all[goodMuons]")
               .Define("Selected_muon_genPartIdx", "Muon_genPartIdx[goodMuons]")
               .Define("Selected_muon_genPartFlav", "Muon_genPartFlav[goodMuons]")
               .Define("Selected_muon_pdgId", "Muon_pdgId[goodMuons]")
               .Define("Selected_muon_idx", ::good_idx, {"goodMuons"});

    // Generate muon 4-vectors
    _rlm = _rlm.Define("muon4vecs", ::generate_4vec, {"Selected_muon_pt",
                                                      "Selected_muon_eta",
                                                      "Selected_muon_phi",
                                                      "Selected_muon_mass"});

    _rlm = _rlm.Define("Vectorial_sum_three_muons", "muon4vecs[0] + muon4vecs[1] + muon4vecs[2]")
               .Define("Vectorial_sum_three_muons_mass", "Vectorial_sum_three_muons.M()");

    /// https://twiki.cern.ch/twiki/bin/view/Main/PdgId
    /// https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
    /// @brief
    /// @details
    /// @param &pdgId            Selected_muon_pdgId
    ///                          Particle Data Group identification number of
    ///                          the (3?) selected muons
    /// https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#Muon
    /// @param &genPartPdgId     GenPart_pdgId
    ///                          Particle Data Group identification number of
    ///                          the generator particles
    /// https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#GenPart
    /// @param &genPartIdx       Selected_muon_genPartIdx
    /// @param &genPartMotherIdx GenPart_genPartIdxMother
    /// https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#GenPart
    /// @return areFromW         A boolean vector specifying whether selected
    ///                          muons are indeed muons from W bosons
    auto isMuonFromW = [](ints &pdgId,
                          ints &genPartPdgId,
                          ints &genPartIdx,
                          ints &genPartMotherIdx)
    {
        bools areFromW;
        bool isFromW = false;
        unsigned int nSelectedMuons = pdgId.size();

        // Iterate over the selected muons
        for (unsigned int i = 0; i < nSelectedMuons; i++)
        {
            isFromW = false;

            // Index of the generator particle corresponding to the current selected muon
            int genIdx = genPartIdx[i];
            std::cout << "i = " << i << std::endl;

            // Match the muon to the same muon in the generator tree
            if (abs(genPartPdgId[genIdx]) == MUON)
            {
                if (pdgId[i] == -genPartPdgId[genIdx])
                {
                    std::cout << "The particle and the generator particle ";
                    std::cout << "don't have the same charge" << std::endl;
                }
            }

            // Track its mother particle
            int motherIdx = genPartMotherIdx[genIdx];
            std::cout << "motherIdx = " << motherIdx << ", ";
            std::cout << "nature = " << genPartPdgId[motherIdx] << std::endl;
            // If motherIdx is out of range, a segmentation fault will be raised
            if (motherIdx > genPartPdgId.size() || motherIdx < -1) break;

            // Check if the muon comes from a W boson
            if (abs(genPartPdgId[motherIdx]) == WPLUS)
            {
                isFromW = true;
                std::cout << "The muon comes from a W boson." << std::endl;

                // Track its grandmother particle
                int grandmotherIdx = genPartMotherIdx[motherIdx];
                std::cout << "grandmotherIdx = " << grandmotherIdx << ", ";
                std::cout << "nature = " << genPartPdgId[grandmotherIdx] << std::endl;
                // If grandmotherIdx is out of range, a segmentation fault will be raised
                if (grandmotherIdx > genPartPdgId.size() || grandmotherIdx < -1) break;

                // Make sure the muon isn't coming from a b quark
                if (abs(genPartPdgId[grandmotherIdx]) == B_QUARK)
                {
                    isFromW = false;
                    std::cout << "The muon comes from a b quark." << std::endl;
                    continue;
                }

                std::cout << "The muon comes from a W boson which ";
                std::cout << "doesn't come from a b quark!" << std::endl;
            }
            else if (abs(genPartPdgId[motherIdx]) == MUON)
            {
                // Check if the muon comes from a W boson
                isFromW = false;

                // Track its grandmother particle
                int grandmotherIdx = genPartMotherIdx[motherIdx];
                std::cout << "grandmotherIdx = " << grandmotherIdx << ", ";
                std::cout << "nature = " << genPartPdgId[grandmotherIdx] << std::endl;
                // If grandmotherIdx is out of range, a segmentation fault will be raised
                if (grandmotherIdx > genPartPdgId.size() || grandmotherIdx < -1) break;

                if (abs(genPartPdgId[grandmotherIdx]) == WPLUS)
                {
                    std::cout << "The muon indirectly comes from a W boson." << std::endl;
                    isFromW = true;
                }
            }
            else
            {
                std::cout << "NO MATCH: " << genPartPdgId[motherIdx] << std::endl;
            }

            areFromW.emplace_back(isFromW);
        }

        return areFromW;
    };

    /*
    _rlm = _rlm.Define("sel_muongenW", isMuonFromW, {"Selected_muon_pdgId",
                                                     "GenPart_pdgId",
                                                     "Selected_muon_genPartIdx",
                                                     "GenPart_genPartIdxMother"})
               .Define("MuonFromW_pt", "Selected_muon_pt[sel_muongenW]")
               .Define("MuonFromW_eta", "Selected_muon_eta[sel_muongenW]")
               .Define("MuonFromW_phi", "Selected_muon_phi[sel_muongenW]")
               .Define("MuonFromW_mass", "Selected_muon_mass[sel_muongenW]")
               .Define("MuonFromW_charge", "Selected_muon_charge[sel_muongenW]")
               .Define("MuonFromW_charge_sum", "MuonFromW_charge[0] + MuonFromW_charge[1] + MuonFromW_charge[2]")
               .Define("MuonFromW_number", "int(MuonFromW_pt.size())")
               .Define("MuonFromW_deltaeta", "abs(MuonFromW_eta[1] - MuonFromW_eta[0])")
               .Define("MuonFromW_deltaphi", "abs(ROOT::VecOps::DeltaPhi(MuonFromW_phi[0], MuonFromW_phi[1]))")
               .Define("MuonFromW_deltaR", "ROOT::VecOps::DeltaR(MuonFromW_eta[0], MuonFromW_eta[1], MuonFromW_phi[0], MuonFromW_phi[1])")
               .Define("MuonFromW_miniPFRelIso_all", "Selected_muon_miniPFRelIso_all[sel_muongenW]")
               .Define("MuonFromW_miniPFRelIso_chg", "Selected_muon_miniPFRelIso_chg[sel_muongenW]")
               .Define("MuonFromW_pfRelIso03_all", "Selected_muon_pfRelIso03_all[sel_muongenW]")
               .Define("MuonFromW_pfRelIso03_chg", "Selected_muon_pfRelIso03_chg[sel_muongenW]")
               .Define("MuonFromW_pfRelIso04_all", "Selected_muon_pfRelIso04_all[sel_muongenW]")
               .Define("MuonFromW_miniPFRelIso_all_divided_pt", "Selected_muon_miniPFRelIso_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_miniPFRelIso_chg_divided_pt", "Selected_muon_miniPFRelIso_chg[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso03_all_divided_pt", "Selected_muon_pfRelIso03_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso03_chg_divided_pt", "Selected_muon_pfRelIso03_chg[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso04_all_divided_pt", "Selected_muon_pfRelIso04_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_genPartIdx", "Selected_muon_genPartIdx[sel_muongenW]")
               .Define("MuonFromW_genPartFlav", "Selected_muon_genPartFlav[sel_muongenW]")
               .Define("MuonFromW_pdgId", "Selected_muon_pdgId[sel_muongenW]")
               .Define("p4_MuonFromW", ::generate_4vec, {"MuonFromW_pt",
                                                         "MuonFromW_eta",
                                                         "MuonFromW_phi",
                                                         "MuonFromW_mass"});
    */
}

/// @brief Select jets
/// @details Check https://twiki.cern.ch/twiki/bin/view/CMS/JetID to find
///     jetId working points for the purpose of your analysis
///     jetId == 2 means: pass tight ID and fail tightLepVeto
///     jetId == 6 means: pass tight ID and tightLepVeto
void TopHiggsTrileptonAnalyser::selectJets()
{
    cout << "Selecting good jets" << endl;
    printDebugInfo(__LINE__, __FUNCTION__);
    
    // Define good jets
    _rlm = _rlm.Define("goodJets", "Jet_pt > 30 &&  abs(Jet_eta) < 2.5 && Jet_jetId >= 6") // && Jet_btagDeepB > 0.1208")
               .Define("Selected_jet_pt", "Jet_pt[goodJets]")
               .Define("Selected_jet_sum_all_jets_pt", "Sum(Selected_jet_pt)")
               .Define("Selected_jet_eta", "Jet_eta[goodJets]")
               .Define("Selected_jet_phi", "Jet_phi[goodJets]")
               .Define("Selected_jet_mass", "Jet_mass[goodJets]")
               .Define("Selected_jet_number", "int(Selected_jet_pt.size())");

    _rlm = _rlm.Define("goodJets_btag", "goodJets && Jet_btagDeepFlavB > 0.049")
               .Define("Selected_bjet_pt", "Jet_pt[goodJets_btag]")
               .Define("Selected_bjet_sum_all_jets_pt", "Sum(Selected_bjet_pt)")
               .Define("Selected_bjet_eta", "Jet_eta[goodJets_btag]")
               .Define("Selected_bjet_phi", "Jet_phi[goodJets_btag]")
               .Define("Selected_bjet_mass", "Jet_mass[goodJets_btag]")
               .Define("Selected_bjet_number", "int(Selected_bjet_pt.size())");
}

/// @brief Check overlaps
/// @details To find the clean jets: remove jets overlapping with any lepton
///     using the lambda function
/*void TopHiggsTrileptonAnalyser::removeOverlaps()
{
    printDebugInfo(__LINE__, __FUNCTION__);

    // Lambda function for checking overlapped jets with leptons
    auto checkoverlap = [](FourVectorVec &seljets, FourVectorVec &sellep)
    {
        doubles mindrlepton;
        for (auto ajet: seljets)
        {
            auto mindr = 6.0;
            for (auto alepton : sellep)
            {
                auto dr = ROOT::Math::VectorUtil::DeltaR(ajet, alepton);
                if (dr < mindr) mindr = dr;
            }
            int out = (mindr > 0.4) ? 1 : 0;
            mindrlepton.emplace_back(out);
        }
        return mindrlepton;
    };

    // Jet-muon separation
    // Overlap removal with muon (used for b-tagging SF)
    cout << "Checking overlap with muons and jets" << endl;

    // Find b-tagged jets
    if (_year == 2016)
    {
       //_rlm = _rlm.Define("btagcuts", "Selmu_jetbtag>0.7221");
    }
    else if (_year == 2017)
    {
       //_rlm = _rlm.Define("btagcuts", "Selmu_jetbtag>0.7476");
    }
    else if (_year == 2018)
    {
       // _rlm = _rlm.Define("btagcuts", "Selmu_jetbtag>0.7100");
    }
}*/

/// @brief Event weight calculation
/// @details
/*void TopHiggsTrileptonAnalyser::calculateEvWeight()
{
    printDebugInfo(__LINE__, __FUNCTION__);
    
    // Muon scale factor
    auto muonSF = [this](floats &pt, floats &eta)->float {
        float weight = 1.0;
        if (pt.size() > 0)
        {
            for (unsigned int i=0; i<pt.size(); i++)
            {
                //float trg_SF = _muontrg->getWeight(std::abs(eta[i]), pt[i]);
                //float ID_SF = _muonid->getWeight(std::abs(eta[i]), pt[i]);
                //float Iso_SF = _muoniso->getWeight(std::abs(eta[i]), pt[i]);
                //weight *= trg_SF * ID_SF * Iso_SF;
            }
        }
        return weight;
    };
    //_rlm = _rlm.Define("evWeight_muonSF", muonSF, {"Selected_muon_pt", "Selected_muon_eta"});
}
*/

/// @brief Define variables
/// @details
void TopHiggsTrileptonAnalyser::defineMoreVars()
{
    printDebugInfo(__LINE__, __FUNCTION__);
  
    // Define St
    addVar({"St", "Sum(Selected_electron_pt) + Sum(Selected_muon_pt) + Sum(Selected_jet_pt)"}); 

    /// @brief Store variables in tree
    /// @details Define variables that you want to store in the tree
    addVartoStore("run");
    addVartoStore("luminosityBlock");
    addVartoStore("event");
    addVartoStore("evWeight.*");

    // Electrons
    addVartoStore("Selected_electron_pt");
    addVartoStore("Selected_electron_sum_all_electrons_pt");
    addVartoStore("Selected_electron_eta");
    addVartoStore("Selected_electron_phi");
    addVartoStore("Selected_electron_mass");
    addVartoStore("Selected_electron_number");
    addVartoStore("elec4vecs");
    
    // Muons
    addVartoStore("goodMuons");
    addVartoStore("Selected_muon_pt");
    addVartoStore("Selected_muon_leading_pt");
    addVartoStore("Selected_muon_subleading_pt");
    addVartoStore("Selected_muon_subsubleading_pt");
    addVartoStore("Selected_muon_sum_two_muons_pt");
    addVartoStore("Selected_muon_sum_all_muons_pt");   
    addVartoStore("Selected_muon_eta");
    addVartoStore("Selected_muon_phi");
    addVartoStore("Selected_muon_mass");
    addVartoStore("Selected_muon_charge");
    addVartoStore("Selected_muon_charge_sum");
    addVartoStore("Selected_muon_number");
    addVartoStore("Selected_muon_deltaeta_01");
    addVartoStore("Selected_muon_deltaeta_02");
    addVartoStore("Selected_muon_deltaeta_12");
    addVartoStore("Selected_muon_deltaphi_01");
    addVartoStore("Selected_muon_deltaphi_02");
    addVartoStore("Selected_muon_deltaphi_12");
    addVartoStore("Selected_muon_deltaR_01");
    addVartoStore("Selected_muon_deltaR_02");
    addVartoStore("Selected_muon_deltaR_12");
    //addVartoStore("Selected_muon_miniPFRelIso_all");
    //addVartoStore("Selected_muon_miniPFRelIso_chg");
    //addVartoStore("Selected_muon_pfRelIso03_all");
    //addVartoStore("Selected_muon_pfRelIso03_chg");
    //addVartoStore("Selected_muon_pfRelIso04_all");
    addVartoStore("Selected_muon_genPartIdx");
    addVartoStore("Selected_muon_genPartFlav");
    addVartoStore("Selected_muon_pdgId");
    addVartoStore("Selected_muon_idx");
    addVartoStore("muon4vecs");
    addVartoStore("Vectorial_sum_two_muons");
    addVartoStore("Vectorial_sum_two_muons_mass");

    // Reconstructed muons
    /*
    addVartoStore("sel_muongenW");
    addVartoStore("MuonFromW_pt");
    addVartoStore("MuonFromW_eta");
    addVartoStore("MuonFromW_phi");
    addVartoStore("MuonFromW_mass");
    addVartoStore("MuonFromW_charge");
    addVartoStore("MuonFromW_charge_sum");
    addVartoStore("MuonFromW_number");
    addVartoStore("MuonFromW_deltaeta");
    addVartoStore("MuonFromW_deltaphi");
    addVartoStore("MuonFromW_deltaR");
    //addVartoStore("MuonFromW_miniPFRelIso_all");
    //addVartoStore("MuonFromW_miniPFRelIso_chg");
    //addVartoStore("MuonFromW_pfRelIso03_all");
    //addVartoStore("MuonFromW_pfRelIso03_chg");
    //addVartoStore("MuonFromW_pfRelIso04_all");
    //addVartoStore("MuonFromW_miniPFRelIso_all_divided_pt");
    //addVartoStore("MuonFromW_miniPFRelIso_chg_divided_pt");
    //addVartoStore("MuonFromW_pfRelIso03_all_divided_pt");
    //addVartoStore("MuonFromW_pfRelIso03_chg_divided_pt");
    //addVartoStore("MuonFromW_pfRelIso04_all_divided_pt");
    addVartoStore("MuonFromW_genPartIdx");
    addVartoStore("MuonFromW_genPartFlav");
    addVartoStore("MuonFromW_pdgId");  
    addVartoStore("p4_MuonW");
    */

    // Jets
    addVartoStore("Selected_jet_pt");
    addVartoStore("Selected_jet_sum_all_jets_pt");
    addVartoStore("Selected_jet_eta");
    addVartoStore("Selected_jet_phi");
    addVartoStore("Selected_jet_mass");
    addVartoStore("Selected_jet_number");
    addVartoStore("Selected_bjet_pt");
    addVartoStore("Selected_bjet_sum_all_jets_pt");
    addVartoStore("Selected_bjet_eta");
    addVartoStore("Selected_bjet_phi");
    addVartoStore("Selected_bjet_mass");
    addVartoStore("Selected_bjet_number");

    // Others
    addVartoStore("St");
}

// TODO: clean/document
//================================Histogram Definitions===========================================//
// add1DHist(TH1DModel histdef, std::string variable, std::string weight, string mincutstep="");
//==============================================================================================//

/// add1DHist
/// @param histdef
/// @param variable
/// @param weight
/// @param mincutstep

/// Histogram definitions
void TopHiggsTrileptonAnalyser::bookHists()
{
    printDebugInfo(__LINE__, __FUNCTION__);

    add1DHist( {"hnevents", "Number of Events", 2, -0.5, 1.5}, "one", "evWeight", "");
    
    // Muons
    //add1DHist( {"hgoodelectron1_pt", "good electron1_pt", 18, -2.7, 2.7}, "good_electron1pt", "evWeight", "0");
    add1DHist({"Number_Muons", "Number of muons", 20, 0.0, 10.0}, "Selected_muon_number", "one", "0");
    add1DHist({"Pt_Muons", "Pt of muons", 40, 0.0, 200.0}, "Selected_muon_pt", "one", "");
    add1DHist({"Leading_Pt_Muons", "Pt of the leading muon", 60, 0.0, 300.0}, "Selected_muon_leading_pt", "one", "");
    add1DHist({"Subleading_Pt_Muons", "Pt of the subleading muon", 40, 0.0, 200.0}, "Selected_muon_subleading_pt", "one", "");
    add1DHist({"Subsubleading_Pt_Muons", "Pt of the subsubleading muon", 40, 0.0, 200.0}, "Selected_muon_subsubleading_pt", "one", "");
    add1DHist({"Sum_Pt_Two_Muons", "Pt of the sum of the two muons", 60, 0.0, 400.0}, "Selected_muon_sum_two_muons_pt", "one", "");
    add1DHist({"Eta_Muons", "Eta of muons", 100, -3.0, 3.0}, "Selected_muon_eta", "one", "");
    add1DHist({"Charge_Muons", "Charge of muons", 50, -2.0, 2.0}, "Selected_muon_charge", "one", "");
    add1DHist({"Deltaphi01_Muons", "Delta phi between the 2 muons", 30, 0.0, 4.0}, "Selected_muon_deltaphi_01", "one", "");
    add1DHist({"Deltaphi02_Muons", "Delta phi between the 2 muons", 30, 0.0, 4.0}, "Selected_muon_deltaphi_02", "one", "");
    add1DHist({"Deltaphi12_Muons", "Delta phi between the 2 muons", 30, 0.0, 4.0}, "Selected_muon_deltaphi_12", "one", "");
    add1DHist({"DeltaR01_Muons", "Delta R between the 2 muons", 30, 0.0, 6.0}, "Selected_muon_deltaR_01", "one", "");
    add1DHist({"DeltaR02_Muons", "Delta R between the 2 muons", 30, 0.0, 6.0}, "Selected_muon_deltaR_02", "one", "");
    add1DHist({"DeltaR12_Muons", "Delta R between the 2 muons", 30, 0.0, 6.0}, "Selected_muon_deltaR_12", "one", "");
    //add1DHist({"Sum_Mass_Two_Muons", "Invariant mass of the two muons", 20, 0.0, 1.0}, "Vectorial_sum_two_muons_mass", "one", "");
    //add1DHist({"miniPFRelIso_all_Muons", "miniPFRel_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "Selected_muon_miniPFRelIso_all", "one", "");
    //add1DHist({"miniPFRelIso_chg_Muons", "miniPFRel_chg isolation variable between the 2 muons", 50, 0.0, 10.0}, "Selected_muon_miniPFRelIso_chg", "one", "");
    //add1DHist({"pfRelIso03_all_Muons", "pfRelIso03_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "Selected_muon_pfRelIso03_all", "one", "");
    //add1DHist({"pfRelIso03_chg_Muons", "pfRelIso03_chg isolation variable between the 2 muons", 50, 0.0, 10.0}, "Selected_muon_pfRelIso03_chg", "one", "");
    //add1DHist({"pfRelIso04_all_Muons", "pfRelIso04_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "Selected_muon_pfRelIso04_all", "one", "");

    // Reconstructed muons
    /*
    add1DHist({"Number_MuonsFromW", "Number of reconstructed muons", 20, 0.0, 10.0}, "MuonFromW_number", "one", "0");
    add1DHist({"Pt_MuonsFromW", "Pt of reconstructed muons", 20, 0.0, 100.0}, "MuonFromW_pt", "one", "0");
    add1DHist({"Eta_MuonsFromW", "Eta of reconstructed muons", 100, -3.0, 3.0}, "MuonFromW_eta", "one", "0");
    add1DHist({"Charge_MuonsFromW","Charge of reconstructed muons", 50, -2.0, 2.0}, "MuonFromW_charge", "one", "0");
    add1DHist({"Deltaphi_MuonsFromW", "Delta phi between the 2 reconstructed muons", 30, 0.0, 4.0}, "MuonFromW_deltaphi", "one", "0");
    add1DHist({"DeltaR_MuonsFromW", "Delta R between the 2 reconstructed muons", 30, 0.0, 6.0}, "MuonFromW_deltaR", "one", "0");
    */

    // Others
    //add1DHist({"St", "St", 50, 0.0, 1500.0}, "St", "one", "");
}

/// @brief
/// @details
/// @param t
/// @param outfilename

void TopHiggsTrileptonAnalyser::setTree(TTree *t, std::string outfilename)
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


/// @brief Selected object definitions
/// @details
void TopHiggsTrileptonAnalyser::setupObjects()
{
	// Object selection will be defined in sequence.
	// Selected objects will be stored in new vectors.
	selectElectrons();
	selectMuons();
	selectJets();
	//removeOverlaps();

}

/// @brief
/// @details
void TopHiggsTrileptonAnalyser::setupAnalysis()
{
	printDebugInfo(__LINE__, __FUNCTION__);
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

/// @brief Prints debug info
/// @details Prints a line and function for debugging purposes
/// @param line
/// @param func
void TopHiggsTrileptonAnalyser::printDebugInfo(int line, std::string func)
{
    if (debug)
    {
        std::cout << "=============================//=============================" << std::endl;
        std::cout << "Line: " << line << std::endl;
        std::cout << "Func: " << func << std::endl;
        std::cout << "=============================//=============================" << std::endl;
    }
}
