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
// $Id: B1DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps

  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 2*m;
  G4double world_sizeZ  = 2*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  // Steel Housing

  G4Material* steelHousing_mat = nist->FindOrBuildMaterial("G4_A-150_TISSUE");
  G4ThreeVector pos1 = G4ThreeVector(0, 2*cm, 2*cm);

  // Conical section shape

  G4double steelHousing_rmina =  0.*m, steelHousing_rmaxa = 1.5*m;
  G4double steelHousing_rminb =  0.*m, steelHousing_rmaxb = 1.5*m;
  G4double steelHousing_hz = 1.5*m;
  G4double steelHousing_phimin = 0.*deg, steelHousing_phimax = 360.*deg;
  G4Cons* solidsteelHousing =
    new G4Cons("steelHousing",
    steelHousing_rmina, steelHousing_rmaxa, steelHousing_rminb, steelHousing_rmaxb, steelHousing_hz,
    steelHousing_phimin, steelHousing_phimax);

  G4LogicalVolume* logicsteelHousing =
    new G4LogicalVolume(solidsteelHousing,         //its solid
                        steelHousing_mat,          //its material
                        "steelHousing");           //its name

G4VisAttributes visibility;
visibility.SetForceWireframe(true);
logicsteelHousing->SetVisAttributes(visibility);


  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicsteelHousing,         //its logical volume
                    "steelHousing",            //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  // Liquid Argom

  G4double density = 1.390*g/cm3;
  G4double z = 18., a = 39.95*g/mole;
  G4Material* lAr_mat = new G4Material ("liquidArgon", z, a, density);
  G4ThreeVector pos3 = G4ThreeVector(0, 0, 0);

  // Conical section shape

                    G4double lAr_rmina =  0.*cm, lAr_rmaxa = 1.*m;
                    G4double lAr_rminb =  0.*cm, lAr_rmaxb = 1.*m;
                    G4double lAr_hz = 1.*m;
                    G4double lAr_phimin = 0.*deg, lAr_phimax = 360.*deg;
                    G4Cons* solidlAr =
                      new G4Cons("lAr",
                      lAr_rmina, lAr_rmaxa, lAr_rminb, lAr_rmaxb, lAr_hz,
                      lAr_phimin, lAr_phimax);

                    G4LogicalVolume* logiclAr =
                      new G4LogicalVolume(solidlAr,         //its solid
                                          lAr_mat,          //its material
                                          "lAr");           //its name

                    new G4PVPlacement(0,                    //no rotation
                                      pos3,                 //at position
                                      logiclAr,             //its logical volume
                                      "lAr",                //its name
                                      logicsteelHousing,                //its mother  volume
                                      false,                   //no boolean operation
                                      0,                       //copy number
                                      checkOverlaps);          //overlaps checking




  return physWorld;
}
