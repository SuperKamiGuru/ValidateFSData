//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
#ifndef G4NeutronHPFissionSpectrum_h
#define G4NeutronHPFissionSpectrum_h 1

#include <fstream>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4NeutronHPVector.hh"
#include "G4VNeutronHPEDis.hh"

// we will need a List of these .... one per term.

class VFDNeutronHPFissionSpectrum : public G4VNeutronHPEDis
{
  public:
  VFDNeutronHPFissionSpectrum()
  {
    expm1 = std::exp(-1.);
  }
  ~VFDNeutronHPFissionSpectrum()
  {
  }

  inline void Init(std::istream & aDataFile)
  {
    theFractionalProb.Init(aDataFile, CLHEP::eV);
    theThetaDist.Init(aDataFile, CLHEP::eV);
  }

  inline G4double GetFractionalProbability(G4double anEnergy)
  {
    return theFractionalProb.GetY(anEnergy);
  }

  inline G4double Sample(G4double anEnergy)
  {
    G4double theta = theThetaDist.GetY(anEnergy);
    // here we need to sample Maxwells distribution, if
    // need be.
    G4double result, cut;
    G4double range =5/CLHEP::eV;
    if(theta>range)
    {
        G4cout << "Error theta exceeds our set maximum energy value of range=5MeV" << G4endl;
    }
    G4double max = Maxwell((theta)/2., theta);
    G4double value;
    do
    {
      result = range*G4UniformRand();
      value = Maxwell(result, theta);
      cut = G4UniformRand();
    }
    while(cut > value/max);
    return result;
  }

  private:

  // this is the function to sample from.
  inline G4double Maxwell(G4double anEnergy, G4double theta)
  {
    G4double result = std::sqrt(anEnergy)*std::exp(-anEnergy/theta);
    return result;
  }

  private:

  G4double expm1;

  G4NeutronHPVector theFractionalProb;

  G4NeutronHPVector theThetaDist;

};

#endif
