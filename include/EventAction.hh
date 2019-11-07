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
/// \file electromagnetic/TestEm18/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 82401 2014-06-18 14:43:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"
#include "G4Step.hh"

#include "Randomize.hh"
#include <iomanip>

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction(DetectorConstruction*, RunAction*);
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void   EndOfEventAction(const G4Event*);

    
    void AddEnergyDeposit1(G4double edep)                     {fEnergyDeposit1  += edep;};
    void AddEnergyDeposit2(G4double edep)                     {fEnergyDeposit2  += edep;};
    void AddPrimaryTrackLength(G4double track)                {fPrimaryTrackLength  += track;};
    void AddSecondary1(G4double ekin)                         {fEnergySecondary1  += ekin;};
    void AddSecondary2(G4double ekin)                         {fEnergySecondary2  += ekin;};
    void AddSecondaryxPolarization(G4double spolarization)    {fSecondaryxPolarization  += spolarization;};
    void AddSecondaryyPolarization(G4double spolarization)    {fSecondaryyPolarization  += spolarization;};
    void AddSecondaryzPolarization(G4double spolarization)    {fSecondaryzPolarization  += spolarization;};
    void AddSecondaryTrackLength(G4double track)              {fSecondaryTrackLength  += track;};
    void AddSecondaryDetTrackLength(G4double track)           {fSecondaryDetTrackLength  += track;};
    void AddTertiary1(G4double ekin)                          {fEnergyTertiary1  += ekin;};
    void AddTertiary2(G4double ekin)                          {fEnergyTertiary2  += ekin;};
    void AddTertiaryxPolarization(G4double spolarization)     {fTertiaryxPolarization  += spolarization;};
    void AddTertiaryyPolarization(G4double spolarization)     {fTertiaryyPolarization  += spolarization;};
    void AddTertiaryzPolarization(G4double spolarization)     {fTertiaryzPolarization  += spolarization;};
    void AddTertiaryTrackLength(G4double track)               {fTertiaryTrackLength  += track;};
    void TrackCheck(G4double trackcheck)                      {fTrack1 = fTrack2; fTrack2 = trackcheck; 
                                                               if (fTrack1 > fTrack2) {
                                                                   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
                                                                   analysisManager->FillH1(31, fTrack1);
                                                                   //G4cout 
     									//<< "\n Last track length of secondary created in detector calculated = " 
     									//<< G4BestUnit(fTrack1, "Length")
     									//<< G4endl;
                                                               }
							      };
    void AddEnergyStrip1a(G4double edep, G4int i)  {fEnergyStrip1a[i] += edep;};
    void AddEnergyStrip1b(G4double edep, G4int i)  {fEnergyStrip1b[i] += edep;};
    void AddEnergyStrip2a(G4double edep, G4int i)  {fEnergyStrip2a[i] += edep;};
    void AddEnergyStrip2b(G4double edep, G4int i)  {fEnergyStrip2b[i] += edep;};

    void AddEnergyPixel(G4double edep, G4int i, G4int j, G4int k, G4int l)  {fEnergyPixel[i][j][k][l] += edep;};
    void AddEnergyPixel1(G4double edep, G4int k, G4int l)  {fEnergyPixel1[k][l] += edep;};

    void AddMomentumDirection1(G4double dirx, G4double diry, G4double dirz)  {fMomDir1x += dirx; fMomDir1y += diry; fMomDir1z += dirz;};
    void AddMomentumDirection2(G4double dirx, G4double diry, G4double dirz)  {fMomDir2x += dirx; fMomDir2y += diry; fMomDir2z += dirz;};

    void AddPointDet1entReal(G4ThreeVector vec)  {fPointDet1entReal += vec;};
    void AddPointDet2entReal(G4ThreeVector vec)  {fPointDet2entReal += vec;};

    void AddPointPix1entReal(G4ThreeVector vec)  {fPointPix1entReal += vec;};
    void AddPointPix2entReal(G4ThreeVector vec)  {fPointPix2entReal += vec;};

    void AddPointPix3rentReal(G4ThreeVector vec)  {fPointPix3rentReal += vec;};
    void AddPointPix4rentReal(G4ThreeVector vec)  {fPointPix4rentReal += vec;};

    void AddA1(G4ThreeVector vec)  {A1 += vec;};
    void AddA(G4ThreeVector vec)   {A += vec;};
    void AddA2(G4ThreeVector vec)  {A2 += vec;};
    void AddB1(G4ThreeVector vec)  {B1 += vec;};
    void AddB(G4ThreeVector vec)   {B += vec;};
    void AddB2(G4ThreeVector vec)  {B2 += vec;};
    void AddC1(G4ThreeVector vec)  {C1 += vec;};
    void AddC(G4ThreeVector vec)   {C += vec;};
    void AddC2(G4ThreeVector vec)  {C2 += vec;};
    void AddD1(G4ThreeVector vec)  {D1 += vec;};
    void AddD(G4ThreeVector vec)   {D += vec;};
    void AddD2(G4ThreeVector vec)  {D2 += vec;};
    void AddE1(G4ThreeVector vec)  {E1 += vec;};
    void AddE(G4ThreeVector vec)   {E += vec;};
    void AddE2(G4ThreeVector vec)  {E2 += vec;};
    void AddF1(G4ThreeVector vec)  {F1 += vec;};
    void AddF(G4ThreeVector vec)   {F += vec;};
    void AddF2(G4ThreeVector vec)  {F2 += vec;};
    void AddG1(G4ThreeVector vec)  {G1 += vec;};
    void AddG(G4ThreeVector vec)   {G += vec;};
    void AddG2(G4ThreeVector vec)  {G2 += vec;};
    void AddH1(G4ThreeVector vec)  {H1 += vec;};
    void AddH(G4ThreeVector vec)   {H += vec;};
    void AddH2(G4ThreeVector vec)  {H2 += vec;};
    void AddI1(G4ThreeVector vec)  {I1 += vec;};
    void AddI(G4ThreeVector vec)   {I += vec;};
    void AddI2(G4ThreeVector vec)  {I2 += vec;};
    void AddJ1(G4ThreeVector vec)  {J1 += vec;};
    void AddJ(G4ThreeVector vec)   {J += vec;};
    void AddJ2(G4ThreeVector vec)  {J2 += vec;};
    void AddBr1(G4ThreeVector vec)  {Br1 += vec;};
    void AddBr2(G4ThreeVector vec)  {Br2 += vec;};
    void AddCr1(G4ThreeVector vec)  {Cr1 += vec;};
    void AddCr2(G4ThreeVector vec)  {Cr2 += vec;};

    void AddAReal(G4ThreeVector vec)   {AReal += vec;};
    void AddBReal(G4ThreeVector vec)   {BReal += vec;};
    void AddCReal(G4ThreeVector vec)   {CReal += vec;};
    void AddDReal(G4ThreeVector vec)   {DReal += vec;};
    void AddEReal(G4ThreeVector vec)   {EReal += vec;};
    void AddFReal(G4ThreeVector vec)   {FReal += vec;};
    void AddGReal(G4ThreeVector vec)   {GReal += vec;};
    void AddHReal(G4ThreeVector vec)   {HReal += vec;};
    void AddIReal(G4ThreeVector vec)   {IReal += vec;};
    void AddJReal(G4ThreeVector vec)   {JReal += vec;};
    void AddBr1Real(G4ThreeVector vec)  {Br1Real += vec;};
    void AddBr2Real(G4ThreeVector vec)  {Br2Real += vec;};
    void AddCr1Real(G4ThreeVector vec)  {Cr1Real += vec;};
    void AddCr2Real(G4ThreeVector vec)  {Cr2Real += vec;};

        
  private:
    DetectorConstruction* fDetectorconstruction;
    RunAction*    fRunAction;
    
    G4double      fEnergyDeposit1;
    G4double      fEnergyDeposit2;
    G4double      fPrimaryTrackLength;
    G4double      fEnergySecondary1;
    G4double      fEnergySecondary2;       
    G4double      fSecondaryxPolarization;  
    G4double      fSecondaryyPolarization;  
    G4double      fSecondaryzPolarization; 
    G4double      fSecondaryTrackLength; 
    G4double      fSecondaryDetTrackLength;
    G4double      fEnergyTertiary1;
    G4double      fEnergyTertiary2;   
    G4double      fTertiaryxPolarization;   
    G4double      fTertiaryyPolarization; 
    G4double      fTertiaryzPolarization; 
    G4double      fTertiaryTrackLength;
    G4double      fEnergyStrip1a[1017];
    G4double      fEnergyStrip1b[1017];
    G4double      fEnergyStrip2a[1017];
    G4double      fEnergyStrip2b[1017];
    G4double      fWeightStrip1a[1017];
    G4double      fWeightStrip1b[1017];
    G4double      fWeightStrip2a[1017];
    G4double      fWeightStrip2b[1017];
    G4double      fMomDir1x, fMomDir1y, fMomDir1z;
    G4double      fMomDir2x, fMomDir2y, fMomDir2z;
    G4double      fEnergyPixel[3][17][81][53];
    G4double      fEnergyPixel1[217][865];
    G4double      fWeightPixel1[217][865];
    //G4double      fWeightPixel1[17][81][53];
    G4double      fWeightPixel2[17][81][53];

    G4double 	  fStripCenterNRX, fStripCenterNRY, fStripCenterNRZ; //Strip center position if the DUT wasn't rotated
    G4double 	  fStripCenterX, fStripCenterY, fStripCenterZ;

    G4double 	  fPixelCenterNRX, fPixelCenterNRY, fPixelCenterNRZ, fPixelCenterNRX1, fPixelCenterNRY1, fPixelCenterNRZ1; //Pixel center position if the DUT wasn't rotated
    G4double 	  fPixelCenterX, fPixelCenterY, fPixelCenterZ;

    G4double	  fDiff1x, fDiff1y, fDiff2x, fDiff2y;

    G4double      Det1SizeZ, Det2SizeZ, Dist, Strip1Depth, Strip1Length, Strip2Depth, Strip2Length, StripDist, StripWidth, StripPitch, posEndArm1, posBeginningArm2, BPIXSizeZ, pixelDepth, pixelX, pixelY, pixelHalf, ROChor, ROCvert, XangleDUT, XangleBPIX, YangleBPIX, ElField1, ElField2;

    G4int 	  totalNbStrips, stripNo, totalNbHitStrips;
    G4int         fCharge1, fCharge2, fCS1a, fCS1b, fCS2a, fCS2b, fChargePix, fChargeStrip1, fChargeStrip2, fChargePix1, fChargePix2, fChargePixel1, fChargePixel2;
    G4int         fHitSensor1, fHitSensor2;
    G4int	  fHitPixelDet1, fHitPixelDet2;
    G4int	  fPixThreshold;

    G4int         fNbHitsStrip1a[1017];
    G4int         fNbHitsStrip1b[1017];
    G4int         fNbHitsStrip2a[1017];
    G4int         fNbHitsStrip2b[1017];
    G4int         fChargeStrip1a[1017];
    G4int         fChargeStrip1b[1017];
    G4int         fChargeStrip2a[1017];
    G4int         fChargeStrip2b[1017];
    G4int         fChargePixel[3][17][81][53];
    G4int         fChargePixel1Matrix[217][865];
    G4int         fNbHitsPixel[3][17][81][53];

    G4int         fHitsMultiplicity1b;
    G4int         fHitsMultiplicity2b;
    G4int         fHitsMultiplicityPix[3];

    G4ThreeVector fPointPix1ent, fPointPix1entReal, fPointPix1mid, fPointPix1ex, fPointDet1ent, fPointDet1entReal, fPointDet1mid, fPointDet1ex, fPointDet2ent, fPointDet2entReal, fPointDet2mid, fPointDet2ex, fPointPix2ent, fPointPix2entReal,  fPointPix2mid, fPointPix2ex, fPointPix3rentReal, fPointPix4rentReal;



    G4ThreeVector A1, A, A2, At, B1, B, B2, Bt, C1, C, C2, Ct, D1, D, D2, Dt, E1, E, E2, Et, F1, F, F2, Ft, G1, G, G2, Gt, H1, H, H2, Ht, I1, I, I2, It, J1, J, J2, Jt, Rav, Br1, Br1t, Br2, Br2t, Cr1, Cr1t, Cr2, Cr2t;
    G4ThreeVector AReal, BReal, CReal, DReal, EReal, FReal, GReal, HReal, IReal, JReal, Br1Real, Br2Real, Cr1Real, Cr2Real;

    G4double M[3][3], MA[3][3], MB[3][3], MC[3][3], MD[3][3], MG[3][3], MH[3][3], MI[3][3], MJ[3][3];
    G4double p1, eig1, eig2, eig3, traceM, q, p2, p, detBa, r, phi;
    G4double Iu[3][3], Ba[3][3], Ga[3][3];

    G4double tB, tC;
    G4double Ga11, Ga00, lambdaVectorX, lambdaVectorY, lambdaVectorZ;

    //primary track's deflection angle in radians
    G4double ftheta;

    G4double lambda;
    G4double lambdaVector[3];

    //Let A = fPointPix1mid, B = fPointDet1mid, C = fPointDet2mid, D = fPointPix2mid. Let dB = distance between B and AD and dC = distance between C and AD.
    //Let B1 be the AD point with same z as B. Let C1 be the AD point with the same z as C.
    G4double dB, dC;
    G4double dBAD, dCAD, dAD;
    G4ThreeVector xBA, xBD, xCA, xCD, xAD;
    G4ThreeVector xBAD, xCAD;

    G4double B1Bx, B1By, C1Cx, C1Cy;

    G4double      fTrack1;
    G4double      fTrack2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
