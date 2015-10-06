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
 // Hadronic Process: High Precision low E neutron tracking
 // original by H.P. Wellisch, TRIUMF, 14-Feb-97
 // Builds and has the Cross-section data for one material.

#ifndef VFDNeutronHPInelastic_h
#define VFDNeutronHPInelastic_h 1

// Class Description
// Final state production model for a high precision (based on evaluated data
// libraries) description of neutron inelastic scattering below 20 MeV;
// 36 exclusive final states are consideded.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with
// the corresponding process.
// Class Description - End

#include "globals.hh"
#include "VFDNeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "VFDNeutronHPChannelList.hh"

#include "VFDNeutronHP2AInelasticFS.hh"
#include "VFDNeutronHP2N2AInelasticFS.hh"
#include "VFDNeutronHP2NAInelasticFS.hh"
#include "VFDNeutronHP2NDInelasticFS.hh"
#include "VFDNeutronHP2NInelasticFS.hh"
#include "VFDNeutronHP2NPInelasticFS.hh"
#include "VFDNeutronHP2PInelasticFS.hh"
#include "VFDNeutronHP3AInelasticFS.hh"
#include "VFDNeutronHP3NAInelasticFS.hh"
#include "VFDNeutronHP3NInelasticFS.hh"
#include "VFDNeutronHP3NPInelasticFS.hh"
#include "VFDNeutronHP4NInelasticFS.hh"
#include "VFDNeutronHPAInelasticFS.hh"
#include "VFDNeutronHPD2AInelasticFS.hh"
#include "VFDNeutronHPDAInelasticFS.hh"
#include "VFDNeutronHPDInelasticFS.hh"
#include "VFDNeutronHPHe3InelasticFS.hh"
#include "VFDNeutronHPN2AInelasticFS.hh"
#include "VFDNeutronHPN2PInelasticFS.hh"
#include "VFDNeutronHPN3AInelasticFS.hh"
#include "VFDNeutronHPNAInelasticFS.hh"
#include "VFDNeutronHPND2AInelasticFS.hh"
#include "VFDNeutronHPNDInelasticFS.hh"
#include "VFDNeutronHPNHe3InelasticFS.hh"
#include "VFDNeutronHPNInelasticFS.hh"
#include "VFDNeutronHPNPAInelasticFS.hh"
#include "VFDNeutronHPNPInelasticFS.hh"
#include "VFDNeutronHPNT2AInelasticFS.hh"
#include "VFDNeutronHPNTInelasticFS.hh"
#include "VFDNeutronHPNXInelasticFS.hh"
#include "VFDNeutronHPPAInelasticFS.hh"
#include "VFDNeutronHPPDInelasticFS.hh"
#include "VFDNeutronHPPInelasticFS.hh"
#include "VFDNeutronHPPTInelasticFS.hh"
#include "VFDNeutronHPT2AInelasticFS.hh"
#include "VFDNeutronHPTInelasticFS.hh"

class VFDNeutronHPInelastic : public G4HadronicInteraction
{
  public:

  VFDNeutronHPInelastic();

  ~VFDNeutronHPInelastic();

  G4HadFinalState * ApplyYourself(const G4HadProjectile & aTrack, G4Nucleus & aTargetNucleus);

  virtual const std::pair<G4double, G4double> GetFatalEnergyCheckLevels() const;

   public:
      G4int GetVerboseLevel() const;
      void SetVerboseLevel( G4int );

  private:

  G4double * xSec;
  //VFDNeutronHPChannelList * theInelastic; // one List per element
      std::vector<VFDNeutronHPChannelList*> theInelastic; // one List per element
  G4String dirName;
  G4int numEle;
      void addChannelForNewElement();

  private:

   VFDNeutronHP2AInelasticFS the2AFS;
   VFDNeutronHP2N2AInelasticFS the2N2AFS;
   VFDNeutronHP2NAInelasticFS the2NAFS;
   VFDNeutronHP2NDInelasticFS the2NDFS;
   VFDNeutronHP2NInelasticFS the2NFS;
   VFDNeutronHP2NPInelasticFS the2NPFS;
   VFDNeutronHP2PInelasticFS the2PFS;
   VFDNeutronHP3AInelasticFS the3AFS;
   VFDNeutronHP3NAInelasticFS the3NAFS;
   VFDNeutronHP3NInelasticFS the3NFS;
   VFDNeutronHP3NPInelasticFS the3NPFS;
   VFDNeutronHP4NInelasticFS the4NFS;
   VFDNeutronHPAInelasticFS theAFS;
   VFDNeutronHPD2AInelasticFS theD2AFS;
   VFDNeutronHPDAInelasticFS theDAFS;
   VFDNeutronHPDInelasticFS theDFS;
   VFDNeutronHPHe3InelasticFS theHe3FS;
   VFDNeutronHPN2AInelasticFS theN2AFS;
   VFDNeutronHPN2PInelasticFS theN2PFS;
   VFDNeutronHPN3AInelasticFS theN3AFS;
   VFDNeutronHPNAInelasticFS theNAFS;
   VFDNeutronHPND2AInelasticFS theND2AFS;
   VFDNeutronHPNDInelasticFS theNDFS;
   VFDNeutronHPNHe3InelasticFS theNHe3FS;
   VFDNeutronHPNInelasticFS theNFS;
   VFDNeutronHPNPAInelasticFS theNPAFS;
   VFDNeutronHPNPInelasticFS theNPFS;
   VFDNeutronHPNT2AInelasticFS theNT2AFS;
   VFDNeutronHPNTInelasticFS theNTFS;
   VFDNeutronHPNXInelasticFS theNXFS;
   VFDNeutronHPPAInelasticFS thePAFS;
   VFDNeutronHPPDInelasticFS thePDFS;
   VFDNeutronHPPInelasticFS thePFS;
   VFDNeutronHPPTInelasticFS thePTFS;
   VFDNeutronHPT2AInelasticFS theT2AFS;
   VFDNeutronHPTInelasticFS theTFS;
};

#endif
