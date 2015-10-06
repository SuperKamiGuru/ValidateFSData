#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <iomanip>

//Models
#include "include/VFDNeutronHPElasticFS.hh"
#include "include/VFDNeutronHPInelastic.hh"
#include "include/VFDNeutronHPCaptureFS.hh"
#include "include/VFDNeutronHPFSFissionFS.hh"
#include "include/VFDNeutronHPFissionBaseFS.hh"
#include "include/VFDNeutronHPFinalState.hh"

//Data Containers
#include "include/IsotopeMass.hh"
#include "include/FSSpectrumData.hh"
#include "include/ElementNames.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"

#define numPointPerOrder 20
#define numFS 43
#define numSamples 300
#define numBins 45

// need to repeat FS application multiple times for the same incoming energy, for processes that use sampling and random techniques (ie everything except elastic) and then compare
// distributions, this can easily be done by setting the yield to 50 for these reactions, while seperately extracting the real yield using new function
// or just reapeat the process multiple timesS
// could get rid of unnessary code in FS classes to speed up calculation,  no need to calc dopp broad nuc or use so many data containers, could make our own data container that stores all of the data we care about
// fix fission fs
// added in fission ff/

// we changed the max tries from 1000 to 50 in VFDNeutronHPPhotonDist.cc:372 to improve efficiency since if it doesn't work after 50 tries it probably won't work for a 1000
// *we took out the function adjust_final_state on line 709 of VFDNeutronHPInelasticCompFS since we are having trouble initiating the ion table, fix this before comparing results to MCNP

// 080808 Something unexpected is happen in G4NeutronHPLabAngularEnergy is caused by the energy regime of the angular-energy not extending down to the minimum energy supplied in the cross-section data, this is not our fault
// this is the fault of the G4NDL data, and seems to mainly occur for the isotopes of berylium
using namespace std;

void ExtractDir(G4String dirName, VFDNeutronHPFinalState ***isoFSData, string **isoNameList, int &vecSize, VFDNeutronHPFinalState *fsType, std::pair <double,double> **enerBound);
double CompareData(string outDir, string *isoNameListG4NDL, VFDNeutronHPFinalState **isoFSDataG4NDL, int fsSizeG4NDL, std::pair <double,double> *enerBoundG4NDL, string *isoNameListMCNP, VFDNeutronHPFinalState **isoFSDataMCNP, int fsSizeMCNP, std::pair <double,double> *enerBoundMCNP, bool sampling, bool *relevantData);
double CalcDiff(string outFileName, VFDNeutronHPFinalState *isoFSDataG4NDL, std::pair <double,double> enerBoundG4NDL, VFDNeutronHPFinalState *isoFSDataMCNP, std::pair <double,double> enerBoundMCNP, bool sampling, bool *relevantData);
void ExtractZA(string fileName, int &Z, int &A);
bool DirectoryExists( const char* pzPath );
void GetDataStream( string filename, std::stringstream& ss);
void SetDataStream( string filename , std::stringstream& ss, bool overWrite );

double *incEnerVec;
int incEnerSize;
//G4Track *nTrack;
G4DynamicParticle *nParticle;

int main (int argc, char** argv)
{
    G4Neutron::NeutronDefinition();
    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *neutron = pTable->FindParticle("neutron");
//    G4Geantino::GeantinoDefinition();
//    G4GenericIon::GenericIonDefinition();
//    G4Neutron::NeutronDefinition();
//    G4Deuteron::DeuteronDefinition();
//    G4Triton::TritonDefinition();
//    G4Gamma::GammaDefinition();
//    G4Proton::ProtonDefinition();
//    G4AntiProton::AntiProtonDefinition();
//    G4PionPlus::PionPlusDefinition();
//    G4PionZero::PionZeroDefinition();
//    G4PionMinus::PionMinusDefinition();
//    G4KaonPlus::KaonPlusDefinition();
//    G4KaonZero::KaonZeroDefinition();
//    G4KaonMinus::KaonMinusDefinition();
//    G4KaonZeroShort::KaonZeroShortDefinition();
//    G4KaonZeroLong::KaonZeroLongDefinition();
//    G4Alpha::AlphaDefinition();
//    G4Electron::ElectronDefinition();
//    G4Positron::PositronDefinition();
//    G4ParticleDefinition *particle = G4ParticleTable::GetParticleTable()->GetGenericIon();
//    G4ProcessManager *pmanager = new G4ProcessManager(particle);
//    particle->SetProcessManager(pmanager);
//    pTable->GetIonTable()->CreateAllIsomer();
//    pTable->SetReadiness(true);

    nParticle = new G4DynamicParticle(neutron,G4ThreeVector(0.,0.,1.), 10.);

    //nTrack = new G4Track(apValueDynamicParticle, 0., G4ThreeVector(0.,0.,0.));

    int fsSizeG4NDL=0, fsSizeMCNP=0;
    double tempDiff;
    vector<double> sumDiffVec;
    string G4NDLDir, MCNPDir, outDir, compareCS="false", dirName, temp;
    string *isoNameListG4NDL, *isoNameListMCNP;
    ElementNames* elementNames;
    IsotopeMass* isoMasses;
    elementNames->SetElementNames();
    isoMasses->SetIsotopeMass();
    string fsDirNameList[numFS] = {"Elastic/FS/", "Capture/FS/", "Fission/FS/", "Fission/FC/", "Fission/SC/", "Fission/TC/", "Fission/LC/", "Inelastic/F01/",
                                    "Inelastic/F02/", "Inelastic/F03/", "Inelastic/F04/", "Inelastic/F05/", "Inelastic/F06/", "Inelastic/F07/", "Inelastic/F08/", "Inelastic/F09/",
                                    "Inelastic/F10/", "Inelastic/F11/", "Inelastic/F12/", "Inelastic/F13/", "Inelastic/F14/", "Inelastic/F15/", "Inelastic/F16/", "Inelastic/F17/",
                                    "Inelastic/F18/", "Inelastic/F19/", "Inelastic/F20/", "Inelastic/F21/", "Inelastic/F22/", "Inelastic/F23/", "Inelastic/F24/", "Inelastic/F25/",
                                    "Inelastic/F26/", "Inelastic/F27/", "Inelastic/F28/", "Inelastic/F29/", "Inelastic/F30/", "Inelastic/F31/", "Inelastic/F32/", "Inelastic/F33/",
                                    "Inelastic/F34/", "Inelastic/F35/","Inelastic/F36/"};
    std::pair<double,double> *enerBoundG4NDL, *enerBoundMCNP;
    bool *relevantData[numFS];
    bool relDataTemp[5][14] = {{0,0,0,0,0,0,1,1,0,0,0,0,0,0}, {0,0,0,0,1,1,0,0,0,0,1,0,0,1}, {1,1,1,1,1,1,1,1,1,1,1,1,1,1}, {1,1,0,0,0,0,1,1,1,0,0,1,0,0},{1,1,0,0,1,1,1,1,1,0,1,1,0,1}};
    relevantData[0] = relDataTemp[0];
    relevantData[1] = relDataTemp[1];
    relevantData[2] = relDataTemp[2];
    for(int i=3; i<7; i++)
    {
        relevantData[i] = relDataTemp[3];
    }
    for(int i=7; i<numFS; i++)
    {
        relevantData[i] = relDataTemp[4];
    }

    VFDNeutronHPFinalState* fsTypeList[numFS] = {new VFDNeutronHPElasticFS, new VFDNeutronHPCaptureFS, new VFDNeutronHPFSFissionFS, new VFDNeutronHPFissionBaseFS, new VFDNeutronHPFissionBaseFS,
    new VFDNeutronHPFissionBaseFS, new VFDNeutronHPFissionBaseFS, new VFDNeutronHPNInelasticFS, new VFDNeutronHPNXInelasticFS
    , new VFDNeutronHP2NDInelasticFS, new VFDNeutronHP2NInelasticFS, new VFDNeutronHP3NInelasticFS, new VFDNeutronHPNAInelasticFS, new VFDNeutronHPN3AInelasticFS
    , new VFDNeutronHP2NAInelasticFS, new VFDNeutronHP3NAInelasticFS, new VFDNeutronHPNPInelasticFS, new VFDNeutronHPN2AInelasticFS, new VFDNeutronHP2N2AInelasticFS
    , new VFDNeutronHPNDInelasticFS, new VFDNeutronHPNTInelasticFS, new VFDNeutronHPNHe3InelasticFS, new VFDNeutronHPND2AInelasticFS, new VFDNeutronHPNT2AInelasticFS
    , new VFDNeutronHP4NInelasticFS, new VFDNeutronHP2NPInelasticFS, new VFDNeutronHP3NPInelasticFS, new VFDNeutronHPN2PInelasticFS, new VFDNeutronHPNPAInelasticFS
    , new VFDNeutronHPPInelasticFS, new VFDNeutronHPDInelasticFS, new VFDNeutronHPTInelasticFS, new VFDNeutronHPHe3InelasticFS, new VFDNeutronHPAInelasticFS
    , new VFDNeutronHP2AInelasticFS, new VFDNeutronHP3AInelasticFS, new VFDNeutronHP2PInelasticFS, new VFDNeutronHPPAInelasticFS, new VFDNeutronHPD2AInelasticFS
    , new VFDNeutronHPT2AInelasticFS, new VFDNeutronHPPDInelasticFS, new VFDNeutronHPPTInelasticFS, new VFDNeutronHPDAInelasticFS };

    VFDNeutronHPFinalState **isoFSDataG4NDL, **isoFSDataMCNP;
    //bool compCS=0;
    stringstream stream;
    if(argc == 5)
    {
        stream << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
        stream >> G4NDLDir >> MCNPDir >> outDir >> compareCS;
    }
    if(argc == 4)
    {
        stream << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3];
        stream >> G4NDLDir >> MCNPDir >> outDir;
    }
    else
    {
        cout << "Please input the G4NDL data directory, the directory contianing the converted MCNP data, the output directory, and the compare CS data flag" << endl;
        for(int i=0; i<numFS; i++)
        {
            delete fsTypeList[i];
        }
        return 0;
    }

//    if(compareCS=="true"||compareCS=="1"||compareCS=="True")
//        compCS=1;

    // creates an incoming neutron energy vector with an equal number of points of each order of magnitude
//    int numStep=(numPoints-1)/13;
//    double enStep=9.0e-11/(numStep);
//    incEnerVec[0]=1.0e-11;
//
//
//    for(int i=0; i<13; i++)
//    {
//        for(int step=0; step<numStep; step++)
//        {
//            incEnerVec[i*numStep+step+1] = incEnerVec[i*numStep+step] + enStep;
//        }
//        if(i==11)
//        {
//            enStep=1.0e+1/(numStep);
//        }
//        else
//        {
//            enStep*=10;
//        }
//    }

    if(G4NDLDir[G4NDLDir.size()-1]!='/')
        G4NDLDir.push_back('/');

    if(MCNPDir[MCNPDir.size()-1]!='/')
        MCNPDir.push_back('/');

    if(outDir[outDir.size()-1]!='/')
        outDir.push_back('/');

    for(int i=0; i<numFS; i++)
    {
        cout << "\n\nExtracting Directory " << G4NDLDir+fsDirNameList[i] << endl;
        ExtractDir(G4NDLDir+fsDirNameList[i], &isoFSDataG4NDL, &isoNameListG4NDL, fsSizeG4NDL, fsTypeList[i], &enerBoundG4NDL);
        cout << "\nExtracting Directory " << MCNPDir+fsDirNameList[i] << endl;
        ExtractDir(MCNPDir+fsDirNameList[i], &isoFSDataMCNP, &isoNameListMCNP, fsSizeMCNP, fsTypeList[i], &enerBoundMCNP);
        cout << "\nComparing FS Data " << G4NDLDir+fsDirNameList[i] << endl;
        sumDiffVec.push_back(CompareData(outDir+fsDirNameList[i], isoNameListG4NDL, isoFSDataG4NDL, fsSizeG4NDL, enerBoundG4NDL, isoNameListMCNP, isoFSDataMCNP, fsSizeMCNP, enerBoundMCNP, i!=0, relevantData[i]));
        for(int j=0; j<fsSizeG4NDL; j++)
        {
            delete isoFSDataG4NDL[j];
        }
        for(int j=0; j<fsSizeMCNP; j++)
        {
            delete isoFSDataMCNP[j];
        }
        if(fsSizeG4NDL>0)
        {
            delete [] isoFSDataG4NDL;
            delete [] isoNameListG4NDL;
            delete [] enerBoundG4NDL;
        }
        if(fsSizeMCNP>0)
        {
            delete [] isoFSDataMCNP;
            delete [] isoNameListMCNP;
            delete [] enerBoundMCNP;
        }
    }
    for(int i=0; i<numFS; i++)
    {
        delete fsTypeList[i];
    }
    //delete nTrack;
    delete nParticle;

    for(int i=0; i<int(sumDiffVec.size()); i++)
    {
        for(int j=i; j<int(sumDiffVec.size()); j++)
        {
            if(sumDiffVec[i]<sumDiffVec[j])
            {
                tempDiff=sumDiffVec[i];
                temp=fsDirNameList[i];
                sumDiffVec[i]=sumDiffVec[j];
                fsDirNameList[i]=fsDirNameList[j];
                sumDiffVec[j]=tempDiff;
                fsDirNameList[j]=temp;
            }
        }
    }
    stream.str("");
    stream.clear();
    stream << "The FS directories, ordered from greatest to lowest discrepancy\n" << std::endl;

    for(int i=0; i<int(sumDiffVec.size()); i++)
    {
        stream << fsDirNameList[i] << "has summed difference squared of " << sumDiffVec[i] << std::endl;
    }
    SetDataStream( outDir+"FileDiffList.txt", stream, true);

    elementNames->ClearStore();
    isoMasses->ClearStore();
    return 0;
}

void ExtractDir(G4String dirName, VFDNeutronHPFinalState ***isoFSData, string **isoNameList, int &vecSize, VFDNeutronHPFinalState *fsType, std::pair <double,double> **enerBound)
{
    int i, j;
    int Z, A, pos1, pos2;
    bool ahead;
    string isoName;
    vector<string> isoNameVec;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dirName.c_str())) != NULL)
    {
        while ((ent = readdir (dir)) != NULL)
        {
            if((string(ent->d_name)!="..")&&(string(ent->d_name)!="."))
            {
                isoName = ent->d_name;
                ExtractZA(isoName, Z, A);
                if((Z!=-1)&&(A!=-1))
                {
                    ahead=false;
                    for(i=0; i<int(isoNameVec.size()); i++)
                    {
                        for(j=0; j<int(min(isoName.length(), isoNameVec[i].length())); j++)
                        {
                            if(isoName[j]>isoNameVec[i][j])
                            {
                                break;
                            }
                            if(isoName[j]<isoNameVec[i][j])
                            {
                                ahead=true;
                                break;
                            }
                            if((isoName.length()<isoNameVec[i].length())&&(j==int(isoName.length()-1)))
                            {
                                ahead=true;
                            }
                        }
                        if(ahead)
                        {
                            break;
                        }
                    }
                    isoNameVec.insert(isoNameVec.begin()+i,isoName);
                }

            }
        }
    }

    pos1 = dirName.find_last_of('/');
    pos2 = dirName.find_last_of('/', pos1-1);
    G4String theFSType = dirName.substr(pos2+1,pos1-pos2);
    dirName = dirName.substr(0,pos2+1);

    vecSize=isoNameVec.size();
    isoFSData[0] = new VFDNeutronHPFinalState *[vecSize];
    isoNameList[0] = new string [vecSize];
    enerBound[0] = new std::pair<double,double> [vecSize];

    IsotopeMass* isoMasses;
    int total;
    double minE, maxE, dummy;
    stringstream stream;

    for(i=0; i<int(isoNameVec.size()); i++)
    {
        isoFSData[0][i] = fsType->New();
        isoNameList[0][i] = isoNameVec[i];
        ExtractZA(isoNameList[0][i], Z, A);
        isoFSData[0][i]->Init(A, Z, 0, dirName, theFSType);

        if(isoFSData[0][i]->HasXsec())
        {
            enerBound[0][i].first = (1+1.0e-5)*isoFSData[0][i]->GetMinEnergy();
            enerBound[0][i].second = (1-1.0e-5)*isoFSData[0][i]->GetMaxEnergy();
        }
        else
        {

            GetDataStream(dirName+"CrossSection/"+isoNameList[0][i], stream);
            stream >> dummy >> dummy >> total;
            stream >> minE >> dummy;
            for(j=0; j<(total-2); j++)
            {
                stream >> dummy >> dummy;
            }
            stream >> maxE >> dummy;

            enerBound[0][i].first = (1+1.0e-5)*minE/1.0e+06;
            enerBound[0][i].second = (1-1.0e-5)*maxE/1.0e+06;

            stream.str("");
            stream.clear();
        }

        VFDNeutronHPFissionBaseFS *fisBaseFS = dynamic_cast<VFDNeutronHPFissionBaseFS*>(isoFSData[0][i]);
        if(fisBaseFS)
        {
            fisBaseFS->SetMass(isoMasses->GetIsotopeMass(Z,A));
        }
    }
    isoNameVec.clear();
}

double CompareData(string outDir, string *isoNameListG4NDL, VFDNeutronHPFinalState **isoFSDataG4NDL, int fsSizeG4NDL, std::pair <double,double> *enerBoundG4NDL, string *isoNameListMCNP, VFDNeutronHPFinalState **isoFSDataMCNP, int fsSizeMCNP, std::pair <double,double> *enerBoundMCNP, bool sampling, bool *relevantData)
{
    bool match;
    string temp;
    vector<double> sumErrorVec;
    double sumDiff=0., tempDiff;
    stringstream stream;
    int zG4NDL, zMCNP, aG4NDL, aMCNP;
    stream << "The FS files that have been compared in this directory, ordered from greatest to lowest discrepancy\n" << std::endl;

    if(!(DirectoryExists((outDir).c_str())))
    {
        system( ("mkdir -p -m=666 "+outDir).c_str());
        if(DirectoryExists((outDir).c_str()))
        {
            cout << "created directory " << outDir << "\n" << endl;
        }
        else
        {
            cout << "\nError: could not create directory " << outDir << "\n" << endl;
            return -1;
        }
    }

    cout << "[";
    cout.flush();
    int nBar=0;
    for(int j=0; j<fsSizeG4NDL; j++)
    {
        match=false;
        ExtractZA(isoNameListG4NDL[j], zG4NDL, aG4NDL);
        if(j<fsSizeMCNP)
        {
            for(int k=j; k<fsSizeMCNP; k++)
            {
                ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                {
                    sumErrorVec.push_back(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                    match=true;
                    break;
                }
            }
            if(!match)
            {
                for(int k=0; k<j; k++)
                {
                    ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                    if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                    {
                        sumErrorVec.push_back(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                        match=true;
                        break;
                    }
                }
            }
        }
        else
        {
            for(int k=0; k<fsSizeMCNP; k++)
            {
                ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                {
                    sumErrorVec.push_back(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                    match=true;
                    break;
                }
            }
        }
        while(nBar<int(j*80/fsSizeG4NDL))
        {
            nBar++;
            cout << "#";
            cout.flush();
        }
    }
    cout << "#]" << endl;

    for(int i=0; i<int(sumErrorVec.size()); i++)
    {
        for(int j=i; j<int(sumErrorVec.size()); j++)
        {
            if(sumErrorVec[i]<sumErrorVec[j])
            {
                tempDiff=sumErrorVec[i];
                temp=isoNameListG4NDL[i];
                sumErrorVec[i]=sumErrorVec[j];
                isoNameListG4NDL[i]=isoNameListG4NDL[j];
                sumErrorVec[j]=tempDiff;
                isoNameListG4NDL[j]=temp;
            }
        }
    }

    for(int i=0; i<int(sumErrorVec.size()); i++)
    {
        stream << isoNameListG4NDL[i] << "has summed difference squared of " << sumErrorVec[i] << std::endl;
        sumDiff+=sumErrorVec[i];
    }

    SetDataStream( outDir+"FileDiffList.txt", stream, true);

    return sumDiff;
}

double CalcDiff(string outFileName, VFDNeutronHPFinalState *isoFSDataG4NDL, std::pair <double,double> enerBoundG4NDL, VFDNeutronHPFinalState *isoFSDataMCNP, std::pair <double,double> enerBoundMCNP, bool sampling, bool *relevantData)
{
    stringstream stream;
    G4HadProjectile *projectile;
    int points=1;
    double totalDiff=0.;
    FSSpectrumData fsDataG4NDL, fsDataMCNP;
    G4HadFinalState *resultG4NDLVec, *resultMCNPVec;

    double emax, emin, enStep;
    if(enerBoundG4NDL.first<enerBoundMCNP.first)
    {
        emin=enerBoundMCNP.first;
    }
    else
    {
        emin=enerBoundG4NDL.first;
    }

    if(enerBoundG4NDL.second>enerBoundMCNP.second)
    {
        emax=enerBoundMCNP.second;
    }
    else
    {
        emax=enerBoundG4NDL.second;
    }

    int numLoops=(std::ceil(std::log10(emax/emin)));
    incEnerSize = numPointPerOrder*numLoops+1;
    incEnerVec = new double[incEnerSize];
    incEnerVec[0]=emin;

    for(int i=0; i<numLoops; i++)
    {
        (incEnerVec[i*numPointPerOrder]*10<emax)? (enStep=9*incEnerVec[i*numPointPerOrder]/numPointPerOrder) : (enStep=(emax-incEnerVec[i*numPointPerOrder])/numPointPerOrder);
        for(int step=0; step<numPointPerOrder; step++)
        {
            incEnerVec[i*numPointPerOrder+step+1] = incEnerVec[i*numPointPerOrder+step] + enStep;
        }
    }

    stream << "Comparison of G4NDL Data to the Converted MCNP data using\n" << incEnerSize << " incoming energy points\n" << numBins << " bins per histogram\n" << std::endl;
    for(int i=0; i<14; i++)
    {
        stream << relevantData[i] << '\t';
    }
    stream << '\n' << endl;

    SetDataStream( outFileName, stream, true);


    for(int i=0; i<incEnerSize; i++)
    {
        nParticle->SetKineticEnergy(incEnerVec[i]);
        projectile = new G4HadProjectile(*nParticle);

        fsDataG4NDL.Clear();
        fsDataMCNP.Clear();

        if(sampling)
        {
            points=numSamples;
        }
        for(int l=0; l<points; l++)
        {
            resultG4NDLVec = isoFSDataG4NDL->ApplyYourself(*projectile);
            resultMCNPVec = isoFSDataMCNP->ApplyYourself(*projectile);
            if((resultG4NDLVec==0)||(resultMCNPVec==0))
            {
                stream.clear();
                stream.str("");
                stream << "No final state data to compare!" << std::endl;
                SetDataStream( outFileName, stream, false);
                delete projectile;
                return 0.;
            }

            fsDataG4NDL.AddData(resultG4NDLVec);
            fsDataMCNP.AddData(resultMCNPVec);
        }
        totalDiff+=fsDataG4NDL.CompareFSData(outFileName, fsDataMCNP, relevantData);

        delete projectile;
    }
    delete [] incEnerVec;
    return totalDiff/points;
}



//ExtractZA
//extracts the Z and the A isotope numbers from the file name
void ExtractZA(string fileName, int &Z, int &A)
{
    std::size_t startPos=0;
    stringstream ss;
    ElementNames* elementNames;
    IsotopeMass* isoMass;
    while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
        startPos++;

    if(startPos==fileName.length())
    {
        //cout << "### File Name Does Not Contian a Z or an A Value " << fileName << " is Invalid for Broadening ###" << std::endl;
        Z=A=-1;
    }
    else
    {
    ////
        std::size_t found1 = fileName.find_first_of('_', startPos);
        if (found1==std::string::npos)
        {
            /*cout << "### File Name Does Not Contian a '_', two are needed, one to seperate the Z and A value, \
            and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
            Z=A=-1;
        }
        else
        {
            std::size_t found2 = fileName.find_first_of('_', found1+1);
            if (found2==std::string::npos)
            {
                /*cout << "### File Name Does Not Contian a second '_', two are needed, one to seperate the Z and A value, \
                and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
                Z=A=-1;
            }
            else
            {

                ss.str(fileName.substr(startPos, found1));
                ss >> Z;
                ss.str("");
                ss.clear();
                if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                    ss.str(fileName.substr(found1+1, found2-found1-3));
                else
                    ss.str(fileName.substr(found1+1, found2-found1-1));
                ss >> A;
                ss.str("");
                ss.clear();
                ss.str(fileName.substr(found2+1));
                if (!(elementNames->CheckName(ss.str(), Z)))
                {
                    //cout << "### " << fileName << " does not include the correct element name at the end ###" << std::endl;
                    Z=A=-1;
                }
                ss.str("");
                ss.clear();
                if((Z>0)&&(A==0))
                {
                    isoMass->GetMostNatAbunIso(Z, A);
                }
            }

        }

    }
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        closedir (pDir);
    }

    return bExists;
}

void GetDataStream( string filename, std::stringstream& ss)
{
   string* data=NULL;
   std::ifstream* in=NULL;
   //string compfilename(filename);

// Use regular text file
  std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
  if ( thefData.good() )
  {
     int file_size = thefData.tellg();
     thefData.seekg( 0 , std::ios::beg );
     char* filedata = new char[ file_size ];
     while ( thefData )
     {
        thefData.read( filedata , file_size );
     }
     thefData.close();
     data = new string ( filedata , file_size );
     delete [] filedata;
  }
  else
  {
// found no data file
//                 set error bit to the stream
     ss.setstate( std::ios::badbit );
     cout << std::endl << "### failed to open ascii file " << filename << " ###" << std::endl;
  }

   if (data != NULL)
   {
        ss.str(*data);
        if(data[0][data->size()-1]!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

   if(in!=NULL)
   {
        in->close();
        delete in;
   }

    if(data!=NULL)
        delete data;
}


void SetDataStream( string filename , std::stringstream& ss, bool overWrite )
{
    // Use regular text file
    string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.erase((compfilename.length()-2),2);
    }

    std::ofstream *out;
    if(overWrite)
        out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::trunc );
    else
        out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::app );
    if ( ss.good() )
    {
         ss.seekg( 0 , std::ios::end );
         int file_size = ss.tellg();
         ss.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( ss ) {
            ss.read( filedata , file_size );
            if(!file_size)
            {
                cout << "\n #### Error the size of the stringstream is invalid ###" << std::endl;
                break;
            }
         }
         out->write(filedata, file_size);
        if (out->fail())
        {
            cout << std::endl << "writing the ascii data to the output file " << compfilename << " failed" << std::endl
                 << " may not have permission to delete an older version of the file" << std::endl;
        }

         delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
     ss.setstate( std::ios::badbit );

     cout << std::endl << "### failed to write to ascii file " << compfilename << " ###" << std::endl;
    }
    out->close();
    delete out;
   ss.str("");
}
