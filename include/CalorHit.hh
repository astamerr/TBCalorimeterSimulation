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
/// \file CalorHit.hh
/// \brief Definition of the B4c::CalorHit class

#ifndef B4cCalorHit_h
#define B4cCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"


namespace B4c
{

/// Calorimeter hit class
///
/// It defines data members to store the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class CalorHit : public G4VHit
{
  public:
    CalorHit();
    CalorHit(G4int iCol,G4int iRow, G4int iLayer);
    CalorHit(const CalorHit&) = default;
    ~CalorHit() override;

    // operators
    CalorHit& operator=(const CalorHit&) = default;
    G4bool operator==(const CalorHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    void Draw()  override{}
    void Print() override;
    
    void SetParticleID(G4int i) {fParticleID = i;}
    G4int GetParticleID() const { return fParticleID;}
    
    void SetColumnID(G4int z) { fColumnID = z; }
    G4int GetColumnID() const { return fColumnID; }

    void SetRowID(G4int z) { fRowID = z; }
    G4int GetRowID() const { return fRowID; }
    
    void SetLayerID(G4int z) { fLayerID = z; }
    G4int GetLayerID() const { return fLayerID; }
    
    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }
    
    void SetTime(G4double t) { fTime = t; }
    G4double GetTime() const { return fTime; }
    
    void SetTimeLoc(G4double t) { fTimeLoc = t; }
    G4double GetTimeLoc() const { return fTimeLoc; }

    //void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    //G4double GetEdep() const { return fEdep; }

    // methods to handle data
    void Add(G4double de, G4double dl);

    // get methods
    G4double GetEdep() const;
    G4double GetTrackLength() const;

  private:
    G4int fParticleID = 9999;
    G4int fColumnID = -1;
    G4int fRowID = -1;
    G4int fLayerID = -1;
    G4double fEdep = 0.;        ///< Energy deposit in the sensitive volume
    G4double fTrackLength = 0.; ///< Track length in the  sensitive volume
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    G4double fTime = 0.;
    G4double fTimeLoc = 0.;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using CalorHitsCollection = G4THitsCollection<CalorHit>;

extern G4ThreadLocal G4Allocator<CalorHit>* CalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* CalorHit::operator new(size_t)
{
  if (!CalorHitAllocator) {
    CalorHitAllocator = new G4Allocator<CalorHit>;
  }
  void *hit;
  hit = (void *) CalorHitAllocator->MallocSingle();
  return hit;
}

inline void CalorHit::operator delete(void *hit)
{
  if (!CalorHitAllocator) {
    CalorHitAllocator = new G4Allocator<CalorHit>;
  }
  CalorHitAllocator->FreeSingle((CalorHit*) hit);
}

inline void CalorHit::Add(G4double de, G4double dl) {
  fEdep += de;
  fTrackLength += dl;
}

inline G4double CalorHit::GetEdep() const {
  return fEdep;
}

inline G4double CalorHit::GetTrackLength() const {
  return fTrackLength;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
