#include "../include/FSSpectrumData.hh"


FSSpectrumData::FSSpectrumData()
{
    //ctor
}

FSSpectrumData::~FSSpectrumData()
{
    //dtor
}

void FSSpectrumData::AddData(G4HadFinalState* result)
{
    G4DynamicParticle *particle;
    int numSec, pSecY, nSecY, nDelY;
    pSecY = nSecY = nDelY = 0;

    numSec=result->GetNumberOfSecondaries();
    AddPrimary(result[0]);

    for(int s=0; s<numSec; s++)
    {
        particle=result->GetSecondary(s)->GetParticle();
        if(particle->GetParticleDefinition()->GetParticleName()=="gamma")
        {
            pSecY++;
            AddPhoton(particle[0]);
        }
        else if(particle->GetParticleDefinition()->GetParticleName()=="neutron")
        {
            if(particle->GetProperTime()>0.)
            {
                nDelY++;
                AddDelayed(particle[0]);
            }
            else
            {
                nSecY++;
                AddSecondary(particle[0]);
            }
        }
        delete particle;
    }
    SetSecYield(nSecY);
    SetDelYield(pSecY);
    SetPhYield(nDelY);

    result->Clear();
}

void FSSpectrumData::AddSecondary(G4DynamicParticle &nSec)
{
    nSecMomDirPhi.push_back(nSec.GetMomentumDirection().getPhi());
    nSecMomDirTheta.push_back(nSec.GetMomentumDirection().getTheta());
    nSecKEn.push_back(nSec.GetKineticEnergy());
}

void FSSpectrumData::AddDelayed(G4DynamicParticle &nDel)
{
    nDelMomDirPhi.push_back(nDel.GetMomentumDirection().getPhi());
    nDelMomDirTheta.push_back(nDel.GetMomentumDirection().getTheta());
    nDelKEn.push_back(nDel.GetKineticEnergy());
}

void FSSpectrumData::AddPhoton(G4DynamicParticle &pSec)
{
    pSecMomDirPhi.push_back(pSec.GetMomentumDirection().getPhi());
    pSecMomDirTheta.push_back(pSec.GetMomentumDirection().getTheta());
    pSecKEn.push_back(pSec.GetKineticEnergy());
}

void FSSpectrumData::AddPrimary(G4HadFinalState &nPrim)
{
    nPrimMomDirPhi.push_back(nPrim.GetMomentumChange().getPhi());
    nPrimMomDirTheta.push_back(nPrim.GetMomentumChange().getTheta());
    nPrimKEn.push_back(nPrim.GetEnergyChange());
}

double FSSpectrumData::CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, bool *relevantData)
{
    std::stringstream stream;
    double totalDiff=0.;

    if(relevantData[0])
    {
        stream << "Comparing the secondary neutron x-y angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nSecMomDirPhi, mcnpData.GetNSecMomDirPhi());
    }

    if(relevantData[1])
    {
        stream << "Comparing the secondary neutron z-xy angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nSecMomDirTheta, mcnpData.GetNSecMomDirTheta());
    }

    if(relevantData[2])
    {
        stream << "Comparing the delayed neutron x-y angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nDelMomDirPhi, mcnpData.GetNDelMomDirPhi());
    }

    if(relevantData[3])
    {
        stream << "Comparing the delayed neutron z-xy angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nDelMomDirTheta, mcnpData.GetNDelMomDirTheta());
    }

    if(relevantData[4])
    {
        stream << "Comparing the secondary photon x-y angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, pSecMomDirPhi, mcnpData.GetPSecMomDirPhi());
    }

    if(relevantData[5])
    {
        stream << "Comparing the secondary photon z-xy angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, pSecMomDirTheta, mcnpData.GetPSecMomDirTheta());
    }

    if(relevantData[6])
    {
        stream << "Comparing the primary neutron x-y angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nPrimMomDirPhi, mcnpData.GetNPrimMomDirPhi());
    }

    if(relevantData[7])
    {
        stream << "Comparing the primary neutron z-xy angular direction\n" << std::endl;
        totalDiff+=CompareHist(stream, nPrimMomDirTheta, mcnpData.GetNPrimMomDirTheta());
    }

    if(relevantData[8])
    {
        stream << "Comparing the secondary neutron kinetic energy\n" << std::endl;
        totalDiff+=CompareHist(stream, nSecKEn, mcnpData.GetNSecKEn());
    }

    if(relevantData[9])
    {
        stream << "Comparing the delayed neutron kinetic energy\n" << std::endl;
        totalDiff+=CompareHist(stream, nDelKEn, mcnpData.GetNDelKEn());
    }

    if(relevantData[10])
    {
        stream << "Comparing the secondary photon kinetic energy\n" << std::endl;
        totalDiff+=CompareHist(stream, pSecKEn, mcnpData.GetPSecKEn());
    }

    if(relevantData[11])
    {
        stream << "Comparing the secondary neutron yield\n" << std::endl;
        totalDiff+=CompareHist(stream, nSecYield, mcnpData.GetNSecYield());
    }

    if(relevantData[12])
    {
        stream << "Comparing the delayed neutron yield\n" << std::endl;
        totalDiff+=CompareHist(stream, nDelYield, mcnpData.GetNDelYield());
    }

    if(relevantData[13])
    {
        stream << "Comparing the secondary photon yield\n" << std::endl;
        totalDiff+=CompareHist(stream, pSecYield, mcnpData.GetPSecYield());
    }

    SetDataStream( outFileName, stream, false );
    return totalDiff;
}

double FSSpectrumData::CompareHist(std::stringstream &stream, std::vector<double> &g4ndlData, std::vector<double> &mcnpData)
{
    double g4ndlHist[numBins]={0.};
    double mcnpHist[numBins]={0.};
    double diffHist[numBins]={0.};
    double binBounds[numBins+1];
    double maxNum, minNum, sumDiff=0.;

    stream.fill(' ');
    stream.precision(6);

    if(g4ndlData.size()>0)
        minNum=maxNum=g4ndlData[0];
    else if(mcnpData.size()>0)
        minNum=maxNum=mcnpData[0];
    else
    {
        return 0.;
    }

    for(int i=0; i<int(g4ndlData.size()); i++)
    {
        if(maxNum<g4ndlData[i])
        {
            maxNum=g4ndlData[i];
        }
        if(minNum>g4ndlData[i])
        {
            minNum=g4ndlData[i];
        }
    }

    for(int i=0; i<int(mcnpData.size()); i++)
    {
        if(maxNum<mcnpData[i])
        {
            maxNum=mcnpData[i];
        }
        if(minNum>mcnpData[i])
        {
            minNum=mcnpData[i];
        }
    }

    for(int i=0; i<numBins+1; i++)
    {
        binBounds[i]=(maxNum-minNum)*i/numBins+minNum;
    }

    for(int i=0; i<int(g4ndlData.size()); i++)
    {
        for(int j=1; j<numBins+1; j++)
        {
            if(g4ndlData[i]<=binBounds[j])
            {
                g4ndlHist[j-1]++;
                break;
            }
        }
    }

    if(g4ndlData.size()>0)
    {
        for(int j=0; j<numBins; j++)
        {
            g4ndlHist[j]/=g4ndlData.size();
        }
    }

    for(int i=0; i<int(mcnpData.size()); i++)
    {
        for(int j=1; j<numBins+1; j++)
        {
            if(mcnpData[i]<=binBounds[j])
            {
                mcnpHist[j-1]++;
                break;
            }
        }
    }

    if(mcnpData.size()>0)
    {
        for(int j=0; j<numBins; j++)
        {
            mcnpHist[j]/=mcnpData.size();
        }
    }

    stream << "G4NDL Hist" << std::endl;
    for(int j=0; j<numBins; j++)
    {
        stream << std::setw(14) << std::right << binBounds[j] << std::setw(14) << std::right << g4ndlHist[j];
        if(j==numBins-1)
        {
            stream << std::setw(14) << std::right << binBounds[numBins];
        }
        if(((j+1)%3==0)||(j==numBins-1))
            stream << '\n';
    }

    stream << "MCNP Hist" << std::endl;
    for(int j=0; j<numBins; j++)
    {
        stream << std::setw(14) << std::right << binBounds[j] << std::setw(14) << std::right << mcnpHist[j];
        if(j==numBins-1)
        {
            stream << std::setw(14) << std::right << binBounds[numBins];
        }
        if(((j+1)%3==0)||(j==numBins-1))
            stream << '\n';
    }

    stream << "Diff Sq Hist" << std::endl;
    for(int j=0; j<numBins; j++)
    {
        diffHist[j] = std::pow(mcnpHist[j]-g4ndlHist[j],2);
        sumDiff+=diffHist[j];
        stream << std::setw(14) << std::right << binBounds[j] << std::setw(14) << std::right << diffHist[j];
        if(j==numBins-1)
        {
            stream << std::setw(14) << std::right << binBounds[numBins];
        }
        if(((j+1)%3==0)||(j==numBins-1))
            stream << '\n';
    }

    stream << '\n';
    return sumDiff;
}

void FSSpectrumData::SetDataStream( std::string filename , std::stringstream& ss, bool overWrite )
{
    // Use regular text file
    std::string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.erase(compfilename.size()-2, 2);
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
                std::cout << "\n #### Error the size of the stringstream is invalid ###" << std::endl;
                break;
            }
         }
         out->write(filedata, file_size);
        if (out->fail())
        {
            std::cout << std::endl << "writing the ascii data to the output file " << compfilename << " failed" << std::endl
                 << " may not have permission to delete an older version of the file" << std::endl;
        }

         delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
     ss.setstate( std::ios::badbit );

     std::cout << std::endl << "### failed to write to ascii file " << compfilename << " ###" << std::endl;
    }
    out->close();
    delete out;
   ss.str("");
}
