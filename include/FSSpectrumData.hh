#ifndef FSSPECTRUMDATA_HH
#define FSSPECTRUMDATA_HH

#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#define numBins 45

class FSSpectrumData
{
    public:
        FSSpectrumData();
        virtual ~FSSpectrumData();
        void AddData(G4HadFinalState* result);
        void AddSecondary(G4DynamicParticle &nSec);
        void AddDelayed(G4DynamicParticle &nDel);
        void AddPhoton(G4DynamicParticle &pSec);
        void AddPrimary(G4HadFinalState &nPrim);

        std::vector<double>& GetNSecMomDirPhi() {return nSecMomDirPhi;}
        std::vector<double>& GetNDelMomDirPhi() {return nDelMomDirPhi;}
        std::vector<double>& GetPSecMomDirPhi() {return pSecMomDirPhi;}
        std::vector<double>& GetNPrimMomDirPhi() {return nPrimMomDirPhi;}
        std::vector<double>& GetNSecMomDirTheta() {return nSecMomDirTheta;}
        std::vector<double>& GetNDelMomDirTheta() {return nDelMomDirTheta;}
        std::vector<double>& GetPSecMomDirTheta() {return pSecMomDirTheta;}
        std::vector<double>& GetNPrimMomDirTheta() {return nPrimMomDirTheta;}
        std::vector<double>& GetNSecKEn() {return nSecKEn;}
        std::vector<double>& GetNDelKEn() {return nDelKEn;}
        std::vector<double>& GetPSecKEn() {return pSecKEn;}
        std::vector<double>& GetNPrimKEn() {return nPrimKEn;}
        std::vector<double>& GetNSecYield() {return nSecYield;}
        std::vector<double>& GetNDelYield() {return nDelYield;}
        std::vector<double>& GetPSecYield() {return pSecYield;}

        void SetSecYield(int yield)
        {
            nSecYield.push_back(yield);
        }
        void SetDelYield(int yield)
        {
            nDelYield.push_back(yield);
        }
        void SetPhYield(int yield)
        {
            pSecYield.push_back(yield);
        }
        double CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, bool *relevantData);
        double CompareHist(std::stringstream &stream, std::vector<double> &g4ndlData, std::vector<double> &mcnpData);
        void SetDataStream( std::string filename , std::stringstream& ss, bool overWrite );
        void Clear()
        {
            nSecMomDirPhi.clear();
            nDelMomDirPhi.clear();
            pSecMomDirPhi.clear();
            nPrimMomDirPhi.clear();
            nSecMomDirTheta.clear();
            nDelMomDirTheta.clear();
            pSecMomDirTheta.clear();
            nPrimMomDirTheta.clear();
            nSecKEn.clear();
            nDelKEn.clear();
            pSecKEn.clear();
            nPrimKEn.clear();
            nSecYield.clear();
            nDelYield.clear();
            pSecYield.clear();
        }
    protected:

    private:
    std::vector<double> nSecMomDirPhi,  nSecMomDirTheta, nDelMomDirPhi, nDelMomDirTheta, pSecMomDirPhi, pSecMomDirTheta, nPrimMomDirPhi, nPrimMomDirTheta;
    std::vector<double> nSecKEn, nDelKEn, pSecKEn, nPrimKEn;
    std::vector<double> nSecYield, nDelYield, pSecYield;
};

#endif // FSSPECTRUMDATA_HH
