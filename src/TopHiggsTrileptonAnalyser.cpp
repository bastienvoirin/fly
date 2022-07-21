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
: NanoAODAnalyzerrdframe(t, outfilename)
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
    //_rlm = _rlm.Range(0, 1000);
    Nentry = _rlm.Count();
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    cout << "Usage of ranges:\n" << "    All entries: " << *Nentry << endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    // Preliminary cut: generator-level signal selection (3 leptons from W bosons)
    addCuts("IsSignal == 0 || (IsSignal == 1 && Sum(IsFromW) >= 3)", "0");

    // Basic selection: 3 leptons, at least 2 oppositely charged leptons, at least 2 oppositely charged same flavor leptons, remove reconstruction issue
    //(jul21a_benj_id)//addCuts("(sel_lepton_number == 3) && (abs(sel_lepton_charge_sum) == 1) && (ll_flav_orderby_dr_delta_index > 0)", "00");
    addCuts("(DeltaR_01 > 0.05) && (DeltaR_02 > 0.05) && (DeltaR_12 > 0.05)", "000");

    addCuts("sel_lepton_sum_3_lepton_pt >= 160.0", "0000");
    addCuts("ll_0_orderby_dr_DeltaR <= 0.8", "00000"); // Should eliminate many ttbar events
    addCuts("sel_bjet_number >= 1", "000000"); // At least 1 b-tagged jet: t => (W, b) => (lν, b-jet); relevant for WZ background, not for ttbar
    addCuts("abs(ll_flav_0_orderby_dr_mass - 91.2) >= 10", "0000000"); // Should eliminate the Z peak of WZ events
}

/// @brief
/// @details
void TopHiggsTrileptonAnalyser::selectSignal(bool isSignal)
{
    std::cout << "################ selectSignal(isSignal = " << isSignal << ") ################" << std::endl;

    _rlm = _rlm.Define("IsSignal",        isSignal ? "1" : "0")
               .Define("GenMu",           "abs(GenPart_pdgId) == 13") // generator-level muons
               .Define("GenEl",           "abs(GenPart_pdgId) == 11") // generator-level electrons
               .Define("MotherGenMu",     "GenPart_genPartIdxMother[GenMu]") // mother particles of generator-level muons
               .Define("MotherGenEl",     "GenPart_genPartIdxMother[GenEl]") // mother particles of generator-level electrons
               .Define("MotherPartPdgId", "ROOT::VecOps::Take(GenPart_pdgId, ROOT::VecOps::Concatenate(MotherGenMu, MotherGenEl))") // PDG Id of the mother particles of generator-level muons
               .Define("IsFromW",         "abs(MotherPartPdgId) == 24"); // is the generator-level muon from a W boson?
}

/// @brief Find good electrons
/// @details Define good electrons in rdata frame
void TopHiggsTrileptonAnalyser::selectElectrons()
{
    cout << "Selecting good electrons" << endl;
    printDebugInfo(__LINE__, __FUNCTION__);

    _rlm = _rlm.Define("goodElecID",              ElectronID(4)); // 2: loose, 3: medium, 4: tight
  //_rlm = _rlm.Define("goodElecBenj",            "goodElecID && Electron_cutBased_HEEP == 1 && Electron_pt > 25 && abs(Electron_eta) < 2.5 && Electron_sip3d < 3 && Electron_miniPFRelIso_all < 0.05")
    _rlm = _rlm.Define("goodElectrons",           "goodElecID && Electron_pt > 25 && abs(Electron_eta) < 2.5 && Electron_sip3d < 3 && Electron_miniPFRelIso_all < 0.05")//"goodElecBenj && Electron_pt > 20 && abs(Electron_eta) < 2.4 && abs(Electron_dxy) < 0.05 && abs(Electron_dz) < 0.1 && Electron_sip3d < 8 && Electron_miniPFRelIso_all < 0.4 && Electron_sieie < 0.03 && Electron_hoe < 0.1 && Electron_eInvMinusPInv > -0.04 && Electron_lostHits == 0 && Electron_convVeto == 1 && Electron_mvaFall17V2noIso_WPL == 1")
               .Define("sel_el_pt",               "Electron_pt[goodElectrons]")
               .Define("sel_el_leading_pt",       "sel_el_pt[0]")
               .Define("sel_el_subleading_pt",    "sel_el_pt[1]")
               .Define("sel_el_subsubleading_pt", "sel_el_pt[2]")
               .Define("sel_el_sum_two_el_pt",    "sel_el_pt[0] + sel_el_pt[1]")
               .Define("sel_el_sum_all_el_pt",    "Sum(sel_el_pt)")
               .Define("sel_el_eta",              "Electron_eta[goodElectrons]")
               .Define("sel_el_phi",              "Electron_phi[goodElectrons]")
               .Define("sel_el_mass",             "Electron_mass[goodElectrons]")
               .Define("sel_el_charge",           "Electron_charge[goodElectrons]")
               .Define("sel_el_charge_sum",       "sel_el_charge[0] + sel_el_charge[1] + sel_el_charge[2]")
               .Define("sel_el_number",           "int(sel_el_pt.size())")
               .Define("sel_el_genPartIdx",       "Electron_genPartIdx[goodElectrons]")
               .Define("sel_el_genPartFlav",      "Electron_genPartFlav[goodElectrons]")
               .Define("sel_el_pdgId",            "Electron_pdgId[goodElectrons]");

    // Generate electron 4-vectors
    _rlm = _rlm.Define("elec4vecs", ::generate_4vec, {"sel_el_pt",
                                                      "sel_el_eta",
                                                      "sel_el_phi",
                                                      "sel_el_mass"});
}

/// @brief Find good muons
/// @details Define good muons in rdata frame
void TopHiggsTrileptonAnalyser::selectMuons()
{
    cout << "Selecting good muons" << endl;
    printDebugInfo(__LINE__, __FUNCTION__);

    // TODO: https://root.cern/doc/master/group__vecops.html

    // Define good muons
    _rlm = _rlm.Define("goodMuonIDLoose",         MuonID(2)) // 2: loose, 3: medium, 4: tight
               .Define("goodMuonID",              MuonID(4))
               .Define("goodMuonsLoose",          "goodMuonIDLoose && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.40")
               .Define("goodMuons",               "goodMuonsLoose && goodMuonID && Muon_pt > 20 && abs(Muon_eta) < 2.4 && Muon_miniPFRelIso_all < 0.05 && Muon_sip3d < 3")
  //_rlm = _rlm.Define("goodMuonsLoose",          "goodMuonID && Muon_pt > 10 && abs(Muon_eta) < 2.4 && Muon_looseId && Muon_miniPFRelIso_all < 0.40")
             //.Define("goodMuons",               "goodMuonsBenj && goodMuonID && Muon_pt > 20 && abs(Muon_eta) < 2.4 && abs(Muon_dxy) < 0.05 && abs(Muon_dz) < 0.1 && Muon_sip3d < 3 && Muon_miniPFRelIso_all < 0.05 && Muon_mvaTTH > 0.85")
               .Define("sel_mu_pt",               "Muon_pt[goodMuons]")
               .Define("sel_mu_leading_pt",       "sel_mu_pt[0]")
               .Define("sel_mu_subleading_pt",    "sel_mu_pt[1]")
               .Define("sel_mu_subsubleading_pt", "sel_mu_pt[2]")
               .Define("sel_mu_sum_two_muons_pt", "sel_mu_pt[0] + sel_mu_pt[1]")
               .Define("sel_mu_sum_all_muons_pt", "Sum(sel_mu_pt)")
               .Define("sel_mu_eta",              "Muon_eta[goodMuons]")
               .Define("sel_mu_phi",              "Muon_phi[goodMuons]")
               .Define("sel_mu_mass",             "Muon_mass[goodMuons]")
               .Define("sel_mu_charge",           "Muon_charge[goodMuons]")
               .Define("sel_mu_charge_sum",       "sel_mu_charge[0] + sel_mu_charge[1] + sel_mu_charge[2]")
               .Define("sel_mu_number",           "int(sel_mu_pt.size())")
               .Define("sel_mu_genPartIdx",       "Muon_genPartIdx[goodMuons]")
               .Define("sel_mu_genPartFlav",      "Muon_genPartFlav[goodMuons]")
               .Define("sel_mu_pdgId",            "Muon_pdgId[goodMuons]");
             //.Define("sel_mu_deltaeta_01", "abs(sel_mu_eta[0] - sel_mu_eta[1])")
             //.Define("sel_mu_deltaeta_02", "abs(sel_mu_eta[0] - sel_mu_eta[2])")
             //.Define("sel_mu_deltaeta_12", "abs(sel_mu_eta[1] - sel_mu_eta[2])")
             //.Define("sel_mu_deltaphi_01", "abs(ROOT::VecOps::DeltaPhi(sel_mu_phi[0], sel_mu_phi[1]))")
             //.Define("sel_mu_deltaphi_02", "abs(ROOT::VecOps::DeltaPhi(sel_mu_phi[0], sel_mu_phi[2]))")
             //.Define("sel_mu_deltaphi_12", "abs(ROOT::VecOps::DeltaPhi(sel_mu_phi[1], sel_mu_phi[2]))")
             //.Define("sel_mu_dr_01", "ROOT::VecOps::DeltaR(sel_mu_eta[0], sel_mu_eta[1], sel_mu_phi[0], sel_mu_phi[1])")
             //.Define("sel_mu_dr_02", "ROOT::VecOps::DeltaR(sel_mu_eta[0], sel_mu_eta[2], sel_mu_phi[0], sel_mu_phi[2])")
             //.Define("sel_mu_dr_12", "ROOT::VecOps::DeltaR(sel_mu_eta[1], sel_mu_eta[2], sel_mu_phi[1], sel_mu_phi[2])")
             //.Define("sel_mu_dr_min", "min(min(sel_mu_dr_01, sel_mu_dr_02), sel_mu_dr_12)")
             //.Define("mu_01", "(sel_mu_charge[0]*sel_mu_charge[1] + 1) / 2")
             //.Define("mu_02", "(sel_mu_charge[0]*sel_mu_charge[2] + 1) / 2")
             //.Define("mu_12", "(sel_mu_charge[1]*sel_mu_charge[2] + 1) / 2")
             //.Define("sel_mu_dr_min_neutral", "min(sel_mu_dr_01*mu_12 + sel_mu_dr_02*mu_01 + sel_mu_dr_12*mu_02, sel_mu_dr_02*mu_12 + sel_mu_dr_12*mu_01 + sel_mu_dr_01*mu_02)")
             //.Define("sel_mu_idx", ::good_idx, {"goodMuons"})
             //.Define("sel_mu_idx", "ROOT::VecOps::Nonzero(goodMuons)");

    // Generate muon 4-vectors
    /*_rlm = _rlm.Define("muon4vecs", ::generate_4vec, {"sel_mu_pt",
                                                      "sel_mu_eta",
                                                      "sel_mu_phi",
                                                      "sel_mu_mass"})
               .Define("dimu_0_orderby_pt_mu_0", "0")
               .Define("dimu_0_orderby_pt_mu_1", "2-1*abs(Muon_charge[0]-Muon_charge[1])/2")
               .Define("dimu_1_orderby_pt_mu_0", "1*abs(Muon_charge[1]-Muon_charge[2])/2")
               .Define("dimu_1_orderby_pt_mu_1", "2")
               .Define("dimu_0_orderby_pt_vectsum", "muon4vecs[dimu_0_orderby_pt_mu_0] + muon4vecs[dimu_0_orderby_pt_mu_1]")
               .Define("dimu_1_orderby_pt_vectsum", "muon4vecs[dimu_1_orderby_pt_mu_0] + muon4vecs[dimu_1_orderby_pt_mu_1]")
               .Define("dimu_0_orderby_pt_mass", "dimu_0_orderby_pt_vectsum.M()")
               .Define("dimu_1_orderby_pt_mass", "dimu_1_orderby_pt_vectsum.M()");*/

    //_rlm = _rlm.Define("Vectorial_sum_three_muons", "muon4vecs[0] + muon4vecs[1] + muon4vecs[2]")
    //           .Define("Vectorial_sum_three_muons_mass", "Vectorial_sum_three_muons.M()");
               
    //_rlm = _rlm.Define("dimuon_candidates", ::dimuon, {"muon4vecs", "sel_mu_charge"});
    /*_rlm = _rlm.Define("dimu_orderby_dr", ::nearest, {"sel_mu_eta", "sel_mu_eta", "sel_mu_phi", "sel_mu_phi", "sel_mu_charge", "sel_mu_charge"}) // TODO: vectorize
               .Define("dimu_orderby_dr_DeltaR", "ROOT::VecOps::DeltaR(sel_mu_eta[dimu_orderby_dr[0]], sel_mu_eta[dimu_orderby_dr[1]], sel_mu_phi[dimu_orderby_dr[0]], sel_mu_phi[dimu_orderby_dr[1]])")
               .Define("dimu_orderby_dr_DeltaPhi", "abs(ROOT::VecOps::DeltaPhi(sel_mu_phi[dimu_orderby_dr[0]], sel_mu_phi[dimu_orderby_dr[1]]))")
               .Define("dimu_0_orderby_dr_vectsum", "muon4vecs[dimu_orderby_dr[0]] + muon4vecs[dimu_orderby_dr[1]]")
               .Define("dimu_1_orderby_dr_vectsum", "muon4vecs[dimu_orderby_dr[2]] + muon4vecs[dimu_orderby_dr[3]]")
               .Define("dimu_0_orderby_dr_mass", "dimu_0_orderby_dr_vectsum.M()")
               .Define("dimu_1_orderby_dr_mass", "dimu_1_orderby_dr_vectsum.M()");*/

    /*_rlm = _rlm.Define("dimu", "sel_mu_idx.size() >= 2 && ROOT::VecOps::Combinations(sel_mu_idx, 2) || [[], []]")
               .Define("dimu_opp_charge_0", "dimu[0][Muon_charge[dimu[0]] != Muon_charge[dimu[1]]]")
               .Define("dimu_opp_charge_1", "dimu[1][Muon_charge[dimu[0]] != Muon_charge[dimu[1]]]");*/
               //.Define("dimu_DeltaR", "ROOT::VecOps::DeltaR(sel_mu_eta[dimu_opp_charge_0], sel_mu_eta[dimu_opp_charge_1], sel_mu_phi[dimu_opp_charge_0], sel_mu_phi[dimu_opp_charge_1])")
               //.Define("dimu_opp_charge_dr_0", "ROOT::VecOps::Take(dimu_opp_charge_0, ROOT::VecOps::Argsort(dimu_DeltaR))")
               //.Define("dimu_opp_charge_dr_1", "ROOT::VecOps::Take(dimu_opp_charge_1, ROOT::VecOps::Argsort(dimu_DeltaR))");

    /// https://twiki.cern.ch/twiki/bin/view/Main/PdgId
    /// https://pdg.lbl.gov/2019/reviews/rpp2018-rev-monte-carlo-numbering.pdf
    /// @brief
    /// @details
    /// @param &pdgId            sel_mu_pdgId
    ///                          Particle Data Group identification number of
    ///                          the (3?) selected muons
    /// https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#muon
    /// @param &genPartPdgId     GenPart_pdgId
    ///                          Particle Data Group identification number of
    ///                          the generator particles
    /// https://cms-nanoaod-integration.web.cern.ch/integration/master-106X/mc102X_doc.html#GenPart
    /// @param &genPartIdx       sel_mu_genPartIdx
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

    
    /*_rlm = _rlm.Define("sel_muongenW", isMuonFromW, {"sel_mu_pdgId",
                                                       "GenPart_pdgId",
                                                       "sel_mu_genPartIdx",
                                                       "GenPart_genPartIdxMother"})
               .Define("MuonFromW_pt", "sel_mu_pt[sel_muongenW]")
               .Define("MuonFromW_eta", "sel_mu_eta[sel_muongenW]")
               .Define("MuonFromW_phi", "sel_mu_phi[sel_muongenW]")
               .Define("MuonFromW_mass", "sel_mu_mass[sel_muongenW]")
               .Define("MuonFromW_charge", "sel_mu_charge[sel_muongenW]")
               .Define("MuonFromW_charge_sum", "MuonFromW_charge[0] + MuonFromW_charge[1] + MuonFromW_charge[2]")
               .Define("MuonFromW_number", "int(MuonFromW_pt.size())")
               .Define("MuonFromW_deltaeta", "abs(MuonFromW_eta[1] - MuonFromW_eta[0])")
               .Define("MuonFromW_deltaphi", "abs(ROOT::VecOps::DeltaPhi(MuonFromW_phi[0], MuonFromW_phi[1]))")
               .Define("MuonFromW_deltaR", "ROOT::VecOps::DeltaR(MuonFromW_eta[0], MuonFromW_eta[1], MuonFromW_phi[0], MuonFromW_phi[1])")
               .Define("MuonFromW_miniPFRelIso_all", "sel_mu_miniPFRelIso_all[sel_muongenW]")
               .Define("MuonFromW_miniPFRelIso_chg", "sel_mu_miniPFRelIso_chg[sel_muongenW]")
               .Define("MuonFromW_pfRelIso03_all", "sel_mu_pfRelIso03_all[sel_muongenW]")
               .Define("MuonFromW_pfRelIso03_chg", "sel_mu_pfRelIso03_chg[sel_muongenW]")
               .Define("MuonFromW_pfRelIso04_all", "sel_mu_pfRelIso04_all[sel_muongenW]")
               .Define("MuonFromW_miniPFRelIso_all_divided_pt", "sel_mu_miniPFRelIso_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_miniPFRelIso_chg_divided_pt", "sel_mu_miniPFRelIso_chg[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso03_all_divided_pt", "sel_mu_pfRelIso03_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso03_chg_divided_pt", "sel_mu_pfRelIso03_chg[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_pfRelIso04_all_divided_pt", "sel_mu_pfRelIso04_all[sel_muongenW]/MuonFromW_pt")
               .Define("MuonFromW_genPartIdx", "sel_mu_genPartIdx[sel_muongenW]")
               .Define("MuonFromW_genPartFlav", "sel_mu_genPartFlav[sel_muongenW]")
               .Define("MuonFromW_pdgId", "sel_mu_pdgId[sel_muongenW]")
               .Define("p4_MuonFromW", ::generate_4vec, {"MuonFromW_pt",
                                                         "MuonFromW_eta",
                                                         "MuonFromW_phi",
                                                         "MuonFromW_mass"});*/
}

void TopHiggsTrileptonAnalyser::selectLeptons()
{
    printDebugInfo(__LINE__, __FUNCTION__);

    _rlm = _rlm.Define("sel_el_mu_pt",                "ROOT::VecOps::Concatenate(sel_el_pt, sel_mu_pt)")
               .Define("sel_lepton_argsort_pt",       "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(sel_el_mu_pt))")
               .Define("sel_lepton_pt",               "ROOT::VecOps::Take(sel_el_mu_pt, sel_lepton_argsort_pt)")
               .Define("sel_lepton_leading_pt",       "sel_lepton_pt[0]")
               .Define("sel_lepton_subleading_pt",    "sel_lepton_pt[1]")
               .Define("sel_lepton_subsubleading_pt", "sel_lepton_pt[2]")
               .Define("sel_lepton_sum_2_lepton_pt",  "sel_lepton_pt[0] + sel_lepton_pt[1]")
               .Define("sel_lepton_sum_3_lepton_pt",  "sel_lepton_pt[0] + sel_lepton_pt[1] + sel_lepton_pt[2]")
               .Define("sel_lepton_eta",              "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_eta,         sel_mu_eta),         sel_lepton_argsort_pt)")
               .Define("sel_lepton_phi",              "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_phi,         sel_mu_phi),         sel_lepton_argsort_pt)")
               .Define("sel_lepton_mass",             "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_mass,        sel_mu_mass),        sel_lepton_argsort_pt)")
               .Define("sel_lepton_charge",           "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_charge,      sel_mu_charge),      sel_lepton_argsort_pt)")
               .Define("sel_lepton_genPartIdx",       "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_genPartIdx,  sel_mu_genPartIdx),  sel_lepton_argsort_pt)")
               .Define("sel_lepton_genPartFlav",      "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_genPartFlav, sel_mu_genPartFlav), sel_lepton_argsort_pt)")
               .Define("sel_lepton_pdgId",            "ROOT::VecOps::Take(ROOT::VecOps::Concatenate(sel_el_pdgId,       sel_mu_pdgId),       sel_lepton_argsort_pt)")
               .Define("sel_lepton_charge_sum",       "Sum(sel_lepton_charge)")
               .Define("sel_lepton_number",           "int(sel_lepton_pt.size())")
               .Define("sel_lepton_0_abs_pdgId",      "abs(sel_lepton_pdgId[0])")
               .Define("sel_lepton_1_abs_pdgId",      "abs(sel_lepton_pdgId[1])")
               .Define("sel_lepton_2_abs_pdgId",      "abs(sel_lepton_pdgId[2])");

    _rlm = _rlm.Define("lepton4vecs", ::generate_4vec, {"sel_lepton_pt",
                                                        "sel_lepton_eta",
                                                        "sel_lepton_phi",
                                                        "sel_lepton_mass"})
               .Define("ll_0_orderby_pt_lepton_0", "0")
               .Define("ll_0_orderby_pt_lepton_1", "2-1*abs(sel_lepton_charge[0]-sel_lepton_charge[1])/2")
               .Define("ll_1_orderby_pt_lepton_0", "1*abs(sel_lepton_charge[1]-sel_lepton_charge[2])/2")
               .Define("ll_1_orderby_pt_lepton_1", "2")
               .Define("ll_0_orderby_pt_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_0_orderby_pt_lepton_0], sel_lepton_eta[ll_0_orderby_pt_lepton_1], sel_lepton_phi[ll_0_orderby_pt_lepton_0], sel_lepton_phi[ll_0_orderby_pt_lepton_1])")
               .Define("ll_1_orderby_pt_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_1_orderby_pt_lepton_0], sel_lepton_eta[ll_1_orderby_pt_lepton_1], sel_lepton_phi[ll_1_orderby_pt_lepton_0], sel_lepton_phi[ll_1_orderby_pt_lepton_1])")
               .Define("ll_0_orderby_pt_vectsum", "lepton4vecs[ll_0_orderby_pt_lepton_0] + lepton4vecs[ll_0_orderby_pt_lepton_1]")
               .Define("ll_1_orderby_pt_vectsum", "lepton4vecs[ll_1_orderby_pt_lepton_0] + lepton4vecs[ll_1_orderby_pt_lepton_1]")
               .Define("ll_0_orderby_pt_mass", "ll_0_orderby_pt_vectsum.M()")
               .Define("ll_1_orderby_pt_mass", "ll_1_orderby_pt_vectsum.M()");

    _rlm = _rlm.Define("ll_orderby_dr", ::nearest, {"sel_lepton_eta",
                                                    "sel_lepton_phi",
                                                    "sel_lepton_charge"})
               .Define("ll_0_orderby_dr_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_orderby_dr[0]], sel_lepton_eta[ll_orderby_dr[1]], sel_lepton_phi[ll_orderby_dr[0]], sel_lepton_phi[ll_orderby_dr[1]])")
               .Define("ll_1_orderby_dr_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_orderby_dr[2]], sel_lepton_eta[ll_orderby_dr[3]], sel_lepton_phi[ll_orderby_dr[2]], sel_lepton_phi[ll_orderby_dr[3]])")
               .Define("ll_0_orderby_dr_DeltaR_z", "ll_0_orderby_dr_DeltaR < 0.025")
               .Define("ll_0_orderby_dr_DeltaPhi", "abs(ROOT::VecOps::DeltaPhi(sel_lepton_phi[ll_orderby_dr[0]], sel_lepton_phi[ll_orderby_dr[1]]))")
               .Define("ll_1_orderby_dr_DeltaPhi", "abs(ROOT::VecOps::DeltaPhi(sel_lepton_phi[ll_orderby_dr[2]], sel_lepton_phi[ll_orderby_dr[3]]))")
               .Define("ll_0_orderby_dr_DeltaEta", "abs(sel_lepton_eta[ll_orderby_dr[0]] - sel_lepton_eta[ll_orderby_dr[1]])")
               .Define("ll_1_orderby_dr_DeltaEta", "abs(sel_lepton_eta[ll_orderby_dr[2]] - sel_lepton_eta[ll_orderby_dr[3]])")
               .Define("ll_0_orderby_dr_vectsum", "lepton4vecs[ll_orderby_dr[0]] + lepton4vecs[ll_orderby_dr[1]]")
               .Define("ll_1_orderby_dr_vectsum", "lepton4vecs[ll_orderby_dr[2]] + lepton4vecs[ll_orderby_dr[3]]")
               .Define("ll_0_orderby_dr_mass", "ll_0_orderby_dr_vectsum.M()")
               .Define("ll_1_orderby_dr_mass", "ll_1_orderby_dr_vectsum.M()")
               .Define("ll_orderby_dr_0", "ll_orderby_dr[0]")
               .Define("ll_orderby_dr_1", "ll_orderby_dr[1]")
               .Define("ll_orderby_dr_delta_index", "ll_orderby_dr[1] - ll_orderby_dr[0]");

    _rlm = _rlm.Define("ll_flav_orderby_dr", ::nearestSameFlavor, {"sel_lepton_eta",
                                                                   "sel_lepton_phi",
                                                                   "sel_lepton_charge",
                                                                   "sel_lepton_pdgId"})
               .Define("ll_flav_0_orderby_dr_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_flav_orderby_dr[0]], sel_lepton_eta[ll_flav_orderby_dr[1]], sel_lepton_phi[ll_flav_orderby_dr[0]], sel_lepton_phi[ll_flav_orderby_dr[1]])")
             //.Define("ll_flav_1_orderby_dr_DeltaR", "ROOT::VecOps::DeltaR(sel_lepton_eta[ll_flav_orderby_dr[2]], sel_lepton_eta[ll_flav_orderby_dr[3]], sel_lepton_phi[ll_flav_orderby_dr[2]], sel_lepton_phi[ll_flav_orderby_dr[3]])")
               .Define("ll_flav_0_orderby_dr_DeltaPhi", "abs(ROOT::VecOps::DeltaPhi(sel_lepton_phi[ll_flav_orderby_dr[0]], sel_lepton_phi[ll_flav_orderby_dr[1]]))")
             //.Define("ll_flav_1_orderby_dr_DeltaPhi", "abs(ROOT::VecOps::DeltaPhi(sel_lepton_phi[ll_flav_orderby_dr[2]], sel_lepton_phi[ll_flav_orderby_dr[3]]))")
               .Define("ll_flav_0_orderby_dr_DeltaEta", "abs(sel_lepton_eta[ll_flav_orderby_dr[0]] - sel_lepton_eta[ll_flav_orderby_dr[1]])")
             //.Define("ll_flav_1_orderby_dr_DeltaEta", "abs(ROOT::VecOps::DeltaEta(sel_lepton_eta[ll_flav_orderby_dr[2]], sel_lepton_eta[ll_flav_orderby_dr[3]]))")
               .Define("ll_flav_0_orderby_dr_vectsum", "lepton4vecs[ll_flav_orderby_dr[0]] + lepton4vecs[ll_flav_orderby_dr[1]]")
             //.Define("ll_flav_1_orderby_dr_vectsum", "lepton4vecs[ll_flav_orderby_dr[2]] + lepton4vecs[ll_flav_orderby_dr[3]]")
               .Define("ll_flav_0_orderby_dr_mass", "ll_flav_0_orderby_dr_vectsum.M()")
               .Define("ll_flav_0_orderby_dr_DeltaR_z", "ll_flav_0_orderby_dr_DeltaR < 0.025")
               .Define("ll_flav_orderby_dr_0", "ll_flav_orderby_dr[0]")
               .Define("ll_flav_orderby_dr_1", "ll_flav_orderby_dr[1]")
               .Define("ll_flav_orderby_dr_delta_index", "ll_flav_orderby_dr[1] - ll_flav_orderby_dr[0]")
             //.Define("ll_flav_1_orderby_dr_mass", "ll_flav_1_orderby_dr_vectsum.M()");
               .Define("DeltaR_01", "ROOT::VecOps::DeltaR(sel_lepton_eta[0], sel_lepton_eta[1], sel_lepton_phi[0], sel_lepton_phi[1])")
               .Define("DeltaR_02", "ROOT::VecOps::DeltaR(sel_lepton_eta[0], sel_lepton_eta[2], sel_lepton_phi[0], sel_lepton_phi[2])")
               .Define("DeltaR_12", "ROOT::VecOps::DeltaR(sel_lepton_eta[1], sel_lepton_eta[2], sel_lepton_phi[1], sel_lepton_phi[2])");
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
    _rlm = _rlm.Define("goodJets",                "Jet_pt > 30 && abs(Jet_eta) < 2.5 && Jet_jetId >= 6") // && Jet_btagDeepB > 0.1208")
               .Define("sel_jet_pt",              "Jet_pt[goodJets]")
               .Define("sel_jet_sum_all_jets_pt", "Sum(sel_jet_pt)")
               .Define("sel_jet_eta",             "Jet_eta[goodJets]")
               .Define("sel_jet_phi",             "Jet_phi[goodJets]")
               .Define("sel_jet_mass",            "Jet_mass[goodJets]")
               .Define("sel_jet_number",           "int(sel_jet_pt.size())");

    _rlm = _rlm.Define("goodJets_btag",            "goodJets && abs(Jet_eta) < 2.5 && Jet_btagDeepFlavB > 0.71") // DeepJet loose = 0.0490, DeepJet medium = 0.2783 DeepJet tight WP = 0.7100
               .Define("sel_bjet_pt",              "Jet_pt[goodJets_btag]")
               .Define("sel_bjet_sum_all_jets_pt", "Sum(sel_bjet_pt)")
               .Define("sel_bjet_eta",             "Jet_eta[goodJets_btag]")
               .Define("sel_bjet_phi",             "Jet_phi[goodJets_btag]")
               .Define("sel_bjet_mass",            "Jet_mass[goodJets_btag]")
               .Define("sel_bjet_pt_leading",      "(Jet_pt[goodJets_btag])[0]")
               .Define("sel_bjet_eta_leading",     "(Jet_eta[goodJets_btag])[0]")
               .Define("sel_bjet_mass_leading",    "(Jet_mass[goodJets_btag])[0]")
               .Define("sel_bjet_number",          "int(sel_bjet_pt.size())");

    _rlm = _rlm.Define("jet4vecs", ::generate_4vec, {"sel_bjet_pt",
                                                     "sel_bjet_eta",
                                                     "sel_bjet_phi",
                                                     "sel_bjet_mass"});
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
    //_rlm = _rlm.Define("evWeight_muonSF", muonSF, {"sel_mu_pt", "sel_mu_eta"});
}
*/

/// @brief Define variables
/// @details
void TopHiggsTrileptonAnalyser::defineMoreVars()
{
    printDebugInfo(__LINE__, __FUNCTION__);
  
    // Define St
    addVar({"St", "Sum(sel_el_pt) + Sum(sel_mu_pt) + Sum(sel_jet_pt)"});
    
    //addVar({"dimuon_deltaphi", "abs(dimuon_candidates[0])", ""});
    //addVar({"dimuon_deltaR", "dimuon_candidates[1]", ""});
    //addVartoStore("dimuon_deltaphi");
    //addVartoStore("dimuon_deltaR");

    /// @brief Store variables in tree
    /// @details Define variables that you want to store in the tree
    addVartoStore("run");
    addVartoStore("luminosityBlock");
    addVartoStore("event");
    addVartoStore("evWeight.*");

    // Electrons
    addVartoStore("sel_el_number");
    
    // Muons
    addVartoStore("sel_mu_number");

    addVartoStore("sel_lepton_pt");
    addVartoStore("sel_lepton_leading_pt");
    addVartoStore("sel_lepton_subleading_pt");
    addVartoStore("sel_lepton_subsubleading_pt");
    addVartoStore("sel_lepton_sum_2_lepton_pt");
    addVartoStore("sel_lepton_sum_3_lepton_pt");
    addVartoStore("sel_lepton_eta");
    addVartoStore("sel_lepton_phi");
    addVartoStore("sel_lepton_mass");
    addVartoStore("sel_lepton_charge");
    addVartoStore("sel_lepton_genPartIdx");
    addVartoStore("sel_lepton_genPartFlav");
    addVartoStore("sel_lepton_pdgId");
    addVartoStore("sel_lepton_charge_sum");
    addVartoStore("sel_lepton_number");

    addVartoStore("ll_0_orderby_pt_mu_0");
    addVartoStore("ll_0_orderby_pt_mu_1");
    addVartoStore("ll_1_orderby_pt_mu_0");
    addVartoStore("ll_1_orderby_pt_mu_1");
    addVartoStore("ll_0_orderby_pt_vectsum");
    addVartoStore("ll_1_orderby_pt_vectsum");
    addVartoStore("ll_0_orderby_pt_mass");
    addVartoStore("ll_1_orderby_pt_mass");

    addVartoStore("ll_0_orderby_dr_DeltaR");
    addVartoStore("ll_1_orderby_dr_DeltaR");
    addVartoStore("ll_0_orderby_dr_DeltaPhi");
    addVartoStore("ll_1_orderby_dr_DeltaPhi");
    addVartoStore("ll_0_orderby_dr_vectsum");
    addVartoStore("ll_1_orderby_dr_vectsum");
    addVartoStore("ll_0_orderby_dr_mass");
    addVartoStore("ll_1_orderby_dr_mass");

    addVartoStore("ll_flav_0_orderby_dr_DeltaR");
    addVartoStore("ll_flav_1_orderby_dr_DeltaR");
    addVartoStore("ll_flav_0_orderby_dr_DeltaPhi");
    addVartoStore("ll_flav_1_orderby_dr_DeltaPhi");
    addVartoStore("ll_flav_0_orderby_dr_vectsum");
    addVartoStore("ll_flav_1_orderby_dr_vectsum");
    addVartoStore("ll_flav_0_orderby_dr_mass");
    addVartoStore("ll_flav_1_orderby_dr_mass");

    // Reconstructed muons
    //addVartoStore("sel_muongenW");
    //addVartoStore("MuonFromW_pt");
    //addVartoStore("MuonFromW_eta");
    //addVartoStore("MuonFromW_phi");
    //addVartoStore("MuonFromW_mass");
    //addVartoStore("MuonFromW_charge");
    //addVartoStore("MuonFromW_charge_sum");
    //addVartoStore("MuonFromW_number");
    //addVartoStore("MuonFromW_deltaeta");
    //addVartoStore("MuonFromW_deltaphi");
    //addVartoStore("MuonFromW_deltaR");
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
    //addVartoStore("MuonFromW_genPartIdx");
    //addVartoStore("MuonFromW_genPartFlav");
    //addVartoStore("MuonFromW_pdgId");  
    //saddVartoStore("p4_MuonW");

    // Jets
    addVartoStore("sel_jet_pt");
    addVartoStore("sel_jet_sum_all_jets_pt");
    addVartoStore("sel_jet_eta");
    addVartoStore("sel_jet_phi");
    addVartoStore("sel_jet_mass");
    addVartoStore("sel_jet_number");
    addVartoStore("sel_bjet_pt");
    addVartoStore("sel_bjet_sum_all_jets_pt");
    addVartoStore("sel_bjet_eta");
    addVartoStore("sel_bjet_phi");
    addVartoStore("sel_bjet_mass");
    addVartoStore("sel_bjet_number");
    addVartoStore("jet4vecs");

    _rlm = _rlm.Define("Vectorial_sum_bl",      "lepton4vecs[3 - ll_orderby_dr[0] - ll_orderby_dr[1]] + jet4vecs[0]")
               .Define("Vectorial_sum_bl_mass", "Vectorial_sum_bl.M()");

    _rlm = _rlm.Define("Vectorial_sum_bl_flav",      "lepton4vecs[3 - ll_flav_orderby_dr[0] - ll_flav_orderby_dr[1]] + jet4vecs[0]")
               .Define("Vectorial_sum_bl_flav_mass", "Vectorial_sum_bl_flav.M()");

    addVartoStore("Vectorial_sum_bl");
    addVartoStore("Vectorial_sum_bl_mass");

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

    
    add1DHist({"hnevents", "Number of events;Number of events;Events", 2, -0.5, 1.5}, "one", "evWeight", "");
    
    add1DHist({"Number_Muons",                    "Number of muons;Number of muons;Events",                                                                     6, -0.5,    5.5},  "sel_mu_number",                  "one", "");
    add1DHist({"Number_Electrons",                "Number of electrons;Number of electrons;Events",                                                             6, -0.5,    5.5},  "sel_el_number",                  "one", "");
    add1DHist({"Number_Leptons",                  "Number of leptons;Number of leptons;Events",                                                                 6, -0.5,    5.5},  "sel_lepton_number",              "one", "");
    add1DHist({"Pt_Leptons",                      "p_{T} of leptons;Leptons p_{T} (GeV);Events",                                                               50,  0.0,  250.0},  "sel_lepton_pt",                  "one", "00");
    add1DHist({"Leading_Pt_Leptons",              "p_{T} of leading leptons;Leading lepton p_{T} (GeV);Events",                                                70,  0.0,  350.0},  "sel_lepton_leading_pt",          "one", "00");
    add1DHist({"Subleading_Pt_Leptons",           "p_{T} of subleading leptons;Subleading lepton p_{T} (GeV);Events",                                          50,  0.0,  250.0},  "sel_lepton_subleading_pt",       "one", "00");
    add1DHist({"Subsubleading_Pt_Leptons",        "p_{T} of subsubleading leptons;Subsubleading lepton p_{T} (GeV);Events",                                    40,  0.0,  200.0},  "sel_lepton_subsubleading_pt",    "one", "00");
    add1DHist({"Sum_Pt_3_Leptons",                "Sum of the p_{T} of the 3 leptons;Total p_{T} (GeV);Events",                                                64,  0.0,  640.0},  "sel_lepton_sum_3_lepton_pt",     "one", "00");
    add1DHist({"Charge_3_Leptons",                "Total charge of the 3 leptons;Total charge of the 3 leptons;Events",                                         7, -3.5,    3.5},  "sel_lepton_charge_sum",          "one", "");
    add1DHist({"St",                              "s_{T};s_{T} (GeV);Events",                                                                                  70,  0.0, 1750.0},  "St",                             "one", "00");
    add1DHist({"Mass_bl",                         "Mass of the bl (dilepton = #DeltaR_{min});Mass of the bl (GeV);Events",                                     60,  0.0,  300.0},  "Vectorial_sum_bl_mass",          "one", "000000");
    add1DHist({"Mass_bl_flav",                    "Mass of the bl (dilepton = #DeltaR_{min}, same flav.);Mass of the bl same flav. (GeV);Events",              60,  0.0,  300.0},  "Vectorial_sum_bl_flav_mass",     "one", "000000");
    
    add1DHist({"Sel_bjet_sum_all_jets_pt",        "Sum of the p_{T} of the b-tagged jets;Sum of the p_{T} of the b-tagged jets;Events",                        64, -0.5,  640.0},  "sel_bjet_sum_all_jets_pt",       "one", "");
    add1DHist({"Sel_bjet_pt",                     "p_{T} of the b-tagged jets;p_{T} of the b-tagged jets;Events",                                              70, -0.5,  350.0},  "sel_bjet_pt",                    "one", "000000");
    add1DHist({"Sel_bjet_eta",                    "#eta of the b-tagged jets;#eta of the b-tagged jets;Events",                                                60,  0.0,    3.0},  "sel_bjet_eta",                   "one", "000000");
    add1DHist({"Sel_bjet_mass",                   "Mass of the b-tagged jets;Mass of the b-tagged jets;Events",                                               100, -0.5,  100.0},  "sel_bjet_mass",                  "one", "000000");
    add1DHist({"Sel_bjet_pt_leading",             "p_{T} of the leading b-tagged jet;p_{T} of the leading b-tagged jet;Events",                                70, -0.5,  350.0},  "sel_bjet_pt_leading",            "one", "000000");
    add1DHist({"Sel_bjet_eta_leading",            "#eta of the leading b-tagged jet;#eta of the leading b-tagged jet;Events",                                  60,  0.0,    3.0},  "sel_bjet_eta_leading",           "one", "000000");
    add1DHist({"Sel_bjet_mass_leading",           "Mass of the leading b-tagged jet;Mass of the leading b-tagged jet;Events",                                 100,  0.0,  100.0},  "sel_bjet_mass_leading",          "one", "000000");
    add1DHist({"Sel_bjet_number",                 "Number of b-tagged jets;Number of b-tagged jets;Events",                                                     6, -0.5,    5.5},  "sel_bjet_number",                "one", "");

    add1DHist({"DeltaIndex_lepton_dilepton",      "Index(lepton_{1}) - Index(lepton_{0});Index(lepton_{1}) - Index(lepton_{0});Events",                         5, -0.25,   2.25}, "ll_orderby_dr_delta_index",      "one", "00");
    add1DHist({"DeltaIndex_lepton_dilepton_flav", "Index(lepton_{1}) - Index(lepton_{0}) same flav.;Index(lepton_{1}) - Index(lepton_{0}) same flav.;Events",   5, -0.25,   2.25}, "ll_flav_orderby_dr_delta_index", "one", "00");

  //add1DHist({"Index_lepton_0_dilepton",         "Index of lepton_{0} (dilepton = #DeltaR_{min});Index of lepton_{0};Events",                                  7, -1.25,   2.25}, "ll_orderby_dr_0",                "one", "00");
  //add1DHist({"Index_lepton_1_dilepton",         "Index of lepton_{1} (dilepton = #DeltaR_{min});Index of lepton_{1};Events",                                  7, -1.25,   2.25}, "ll_orderby_dr_1",                "one", "00");
  //add1DHist({"Index_lepton_0_dilepton_flav",    "Index of lepton_{0} (dilepton = #DeltaR_{min}, same flav.);Index of lepton_{0} same flav.;Events",           7, -1.25,   2.25}, "ll_flav_orderby_dr_0",           "one", "00");
  //add1DHist({"Index_lepton_1_dilepton_flav",    "Index of lepton_{1} (dilepton = #DeltaR_{min}, same flav.);Index of lepton_{1} same flav.;Events",           7, -1.25,   2.25}, "ll_flav_orderby_dr_1",           "one", "00");

  //add1DHist({"DeltaR_min_dilepton_zero",        "#DeltaR_{min} < 0.025?;#DeltaR_{min} < 0.025?;Events",                                                       3, -0.25,   1.25}, "ll_0_orderby_dr_DeltaR_z",       "one", "00");
  //add1DHist({"DeltaR_min_dilepton_flav_zero",   "#DeltaR_{min, same flav.} < 0.025?;#DeltaR_{min, same flav.} < 0.025?;Events",                               3, -0.25,   1.25}, "ll_flav_0_orderby_dr_DeltaR_z",  "one", "00");
    add1DHist({"DeltaR_min_dilepton",             "#DeltaR_{min};#DeltaR_{min};Events",                                                                        80,  0.0,    4.0},  "ll_0_orderby_dr_DeltaR",         "one", "00");
    add1DHist({"DeltaR_min_dilepton_flav",        "#DeltaR_{min, same flav.};#DeltaR_{min, same flav.};Events",                                                80,  0.0,    4.0},  "ll_flav_0_orderby_dr_DeltaR",    "one", "00");
    add1DHist({"DeltaR_min_dilepton_zoom",        "#DeltaR_{min};#DeltaR_{min};Events",                                                                       100,  0.0,    0.1},  "ll_0_orderby_dr_DeltaR",         "one", "00");
    add1DHist({"DeltaR_min_dilepton_flav_zoom",   "#DeltaR_{min, same flav.};#DeltaR_{min, same flav.};Events",                                               100,  0.0,    0.1},  "ll_flav_0_orderby_dr_DeltaR",    "one", "00");
  //add1DHist({"DeltaR_dilepton_with_leading",    "#Delta R(dilepton including leading lepton);#DeltaR;Events",                                                55,  0.0,    5.5},  "ll_0_orderby_pt_DeltaR",         "one", "00");
  //add1DHist({"DeltaR_dilepton_without_leading", "#Delta R(dilepton excluding leading lepton);#DeltaR;Events",                                                55,  0.0,    5.5},  "ll_1_orderby_pt_DeltaR",         "one", "00");
    add1DHist({"DeltaPhi_dilepton",               "#Delta#phi (dilepton = #DeltaR_{min});#Delta#phi;Events",                                                   35,  0.0,    3.5},  "ll_0_orderby_dr_DeltaPhi",       "one", "00");
    add1DHist({"DeltaPhi_dilepton_flav",          "#Delta#phi (dilepton = #DeltaR_{min}, same flav.);#Delta#phi same flav.;Events",                            35,  0.0,    3.5},  "ll_flav_0_orderby_dr_DeltaPhi",  "one", "00");
    add1DHist({"DeltaEta_dilepton",               "#Delta#eta (dilepton = #DeltaR_{min});#Delta#eta;Events",                                                  100,  0.0,    2.5},  "ll_0_orderby_dr_DeltaEta",       "one", "00");
    add1DHist({"DeltaEta_dilepton_flav",          "#Delta#eta (dilepton = #DeltaR_{min}, same flav.);#Delta#eta same flav.;Events",                           100,  0.0,    2.5},  "ll_flav_0_orderby_dr_DeltaEta",  "one", "00");
    add1DHist({"Mass_dilepton",                   "Mass of the dilepton (#DeltaR_{min});Mass of the dilepton (GeV);Events",                                   100,  0.0,  200.0},  "ll_0_orderby_dr_mass",           "one", "00");
    add1DHist({"Mass_dilepton_flav",              "Mass of the dilepton (#DeltaR_{min}, same flav.);Mass of the dilepton same flav. (GeV);Events",            100,  0.0,  200.0},  "ll_flav_0_orderby_dr_mass",      "one", "00"); 

    //add1DHist({"DeltaR_min_Dimuon", "DeltaR min;Minimum DeltaR between 2 muons of opposite charges;Events", 30, 0.0, 6.0}, "dimuon_deltaR", "one", "");
    //add1DHist({"Sum_Mass_Two_Muons", "Invariant mass of the two muons", 20, 0.0, 1.0}, "Vectorial_sum_two_muons_mass", "one", "");
    //add1DHist({"miniPFRelIso_all_Muons", "miniPFRel_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "sel_mu_miniPFRelIso_all", "one", "");
    //add1DHist({"miniPFRelIso_chg_Muons", "miniPFRel_chg isolation variable between the 2 muons", 50, 0.0, 10.0}, "sel_mu_miniPFRelIso_chg", "one", "");
    //add1DHist({"pfRelIso03_all_Muons", "pfRelIso03_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "sel_mu_pfRelIso03_all", "one", "");
    //add1DHist({"pfRelIso03_chg_Muons", "pfRelIso03_chg isolation variable between the 2 muons", 50, 0.0, 10.0}, "sel_mu_pfRelIso03_chg", "one", "");
    //add1DHist({"pfRelIso04_all_Muons", "pfRelIso04_all isolation variable between the 2 muons", 50, 0.0, 10.0}, "sel_mu_pfRelIso04_all", "one", "");

    // Reconstructed muons
    /*
    add1DHist({"Number_MuonsFromW", "Number of reconstructed muons", 20, 0.0, 10.0}, "MuonFromW_number", "one", "0");
    add1DHist({"Pt_MuonsFromW", "Pt of reconstructed muons", 20, 0.0, 100.0}, "MuonFromW_pt", "one", "0");
    add1DHist({"Eta_MuonsFromW", "Eta of reconstructed muons", 100, -3.0, 3.0}, "MuonFromW_eta", "one", "0");
    add1DHist({"Charge_MuonsFromW","Charge of reconstructed muons", 50, -2.0, 2.0}, "MuonFromW_charge", "one", "0");
    add1DHist({"Deltaphi_MuonsFromW", "Delta phi between the 2 reconstructed muons", 30, 0.0, 4.0}, "MuonFromW_deltaphi", "one", "0");
    add1DHist({"DeltaR_MuonsFromW", "Delta R between the 2 reconstructed muons", 30, 0.0, 6.0}, "MuonFromW_deltaR", "one", "0");
    */
}

/// @brief
/// @details
/// @param t
/// @param outfilename

void TopHiggsTrileptonAnalyser::setTree(TTree *t, std::string outfilename)
{
	if (debug){
        std::cout << "=============================//=============================" << std::endl;
        std::cout << "Line: " << __LINE__ << std::endl;
        std::cout << "Func: " << __FUNCTION__ << std::endl;
        std::cout << "=============================//=============================" << std::endl;
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
void TopHiggsTrileptonAnalyser::setupObjects(bool isSignal)
{
    std::cout << "################ setupObjects(isSignal = " << isSignal << ") ################" << std::endl;
	// Object selection will be defined in sequence.
	// Selected objects will be stored in new vectors.
	selectElectrons();
	selectMuons();
    selectLeptons();
	selectJets();
    selectSignal(isSignal);
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
    _rlm = _rlm.Define("evWeight", "one");

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
    return;
    if (debug)
    {
        std::cout << "=============================//=============================" << std::endl;
        std::cout << "Line: " << line << std::endl;
        std::cout << "Func: " << func << std::endl;
        std::cout << "=============================//=============================" << std::endl;
    }
}
