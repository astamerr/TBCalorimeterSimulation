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
/// \file EventAction.hh
/// \brief Definition of the B4c::EventAction class

#ifndef B4cEventAction_h
#define B4cEventAction_h 1

#include "G4UserEventAction.hh"

#include "CalorHit.hh"

#include "globals.hh"
#include <vector>

namespace B4c
{

/// Event action class
///
/// In EndOfEventAction(), it prints the accumulated quantities of the energy
/// deposit and track lengths of charged particles in Absober and Gap layers
/// stored in the hits collections.

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override;

  void  BeginOfEventAction(const G4Event* event) override;
  void    EndOfEventAction(const G4Event* event) override;
  
  G4double firstHadInteractionDepth;
  G4double firstHadInteractionTime;
  G4double leakageEnergy;
  
  //vectors for gap hit
  std::vector<G4double> GapHitEdep;
  std::vector<G4int> GapHitParticle;
  std::vector<G4double> GapHitX;
  std::vector<G4double> GapHitY;
  std::vector<G4double> GapHitZ;
  std::vector<G4int> GapHitLayer;
  std::vector<G4int> GapHitRow;
  std::vector<G4int> GapHitColumn;
  std::vector<G4double> GapHitTime;
  
  //vectors for absorber hit
  std::vector<G4double> AbsoHitEdep;
  std::vector<G4int> AbsHitParticle;
  std::vector<G4double> AbsHitX;
  std::vector<G4double> AbsHitY;
  std::vector<G4double> AbsHitZ;
  std::vector<G4int> AbsoHitLayer;
  std::vector<G4int> AbsHitRow;
  std::vector<G4int> AbsHitColumn;
  std::vector<G4double> AbsHitTime;
  
  //std::vector<G4int>& GetHitLayer() {return fLayerID; }

private:
  // methods
  CalorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double gapTrackLength) const;

  // data members
  G4int fAbsHCID = -1;
  G4int fGapHCID = -1;
  
};

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


