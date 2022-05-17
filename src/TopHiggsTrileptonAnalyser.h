#ifndef TOPHIGGSTRILEPTONANALYSER_H_
#define TOPHIGGSTRILEPTONANALYSER_H_

#include "NanoAODAnalyzerrdframe.h"

class TopHiggsTrileptonAnalyser: public NanoAODAnalyzerrdframe
{
    public:
        TopHiggsTrileptonAnalyser(
            TTree *t,
            std::string outfilename
        );
        void defineCuts();
        void defineMoreVars();
        void bookHists();
        void setTree(
            TTree *t,
            std::string outfilename
        );
        void setupAnalysis();
        void setupObjects(); 
        void selectElectrons();
        void selectMuons();
        void selectJets();
        //void removeOverlaps();
        //void selectFatJets();
        //void calculateEvWeight();
        void printDebugInfo(
            int line,
            std::string func
        );

        bool debug = true;

        //private:
        //    std::string year;
        //    std::string runtype;
        //    std::string syst;
};

#endif // TOPHIGGSTRILEPTONANALYSER_H_
