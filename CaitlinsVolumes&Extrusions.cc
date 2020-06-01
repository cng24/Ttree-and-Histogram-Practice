#include "G4ios.hh"
#include "globals.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4EllipticalTube.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4RunManager.hh"

#include "G4SDManager.hh"
#include "MMG4DetectorConstruction.hh"
#include "MMG4DetectorMessenger.hh"
#include "MMG4Materials.hh"
#include "MMG4PhotonDetSD.hh"

#include "G4UserLimits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//TESTCOMMENT
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MMG4DetectorConstruction::MMG4DetectorConstruction()
: fMaterials(NULL), fLogicHole(NULL), fLogicWorld(NULL),
fPhysiWorld(NULL), fPhysiHole(NULL)
{
  //fDetectorMessenger helps in controlling the parameters of the detector geometry
  fDetectorMessenger = new MMG4DetectorMessenger(this);
  
  // The number of layers of cladding to be used over the WLS fiber
  fNumOfCladLayers = 0;
  
  fSurfaceRoughness = 1;
  
  fMirrorToggle = true;
  fMirrorPolish = 1.;
  fMirrorReflectivity = 1.;
  
  fSiPMPolish = 1.;
  fSiPMReflectivity = 0.;
  
  fExtrusionPolish = 1.;
  fExtrusionReflectivity = 1.;
  
  fXYRatio = 1.0;
  
  fMMG4fiberLength = 250*mm
  fFiberLength1 = 100*mm; //will be 514.15mm with semicircle
  fMMG4fiberY = 20*cm;
  fMMG4fiberRY  = 0.5*mm;
  fMMG4fiberOrigin = 0*cm;
  
  fSiPMShape = "Square";
  fSiPMHalfL = fMMG4fiberRY;
  fSiPMDist  = 0.00*mm;
  fSiPMTheta = 0.0*deg;
  fSiPMThickness     = 1*mm;
  
  fClrfiberLength  = fSiPMThickness + 10.*nm;
  fMirrorZ    = 0.1*mm;
  
  fPanelWidth        = 0.254*cm;
  fPanelBase          = 25*cm;
  
  fCoatingThickness = 0.21*mm; //thickness of reflector shield
  
  
  
fGrooveDepth = 1.1*mm; //I make it slightly larger than Fiber so we can have the sealant
fFiberLength1 = 100*mm;  
fpRming = 100 * mm
fpRmaxg = 103 * mm  
fpDzg = fGrooveDepth/2 * mm 
fpDPhig = 0 * Degree
fpSPhig = 180 * Degree 
fpRminf = 0 * mm
fpRmaxf = 1.00 * mm 
fpRtorf = 100 * mm 
fpDPhif = 0 * Degree
fpSPhif = 180 * Degree
fSiPMbase = 3*mm
fSiPMy = 12.5*mm //random location in the 25 mm zone between edge and fiber
  
  
//  fCoatingRadius    = 0.*mm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MMG4DetectorConstruction::~MMG4DetectorConstruction()
{
  if (fDetectorMessenger) delete fDetectorMessenger;
  if (fMaterials)         delete fMaterials;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MMG4DetectorConstruction::Construct()
{
  if (fPhysiWorld) {
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    G4LogicalSkinSurface::CleanSurfaceTable();
    G4LogicalBorderSurface::CleanSurfaceTable();
  }
  
  fMaterials = MMG4Materials::GetInstance();
  
  UpdateGeometryParameters();
  
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* MMG4DetectorConstruction::ConstructDetector()
{
  //--------------------------------------------------
  // Construct the world volume
  //--------------------------------------------------
  
  G4VSolid* solidWorld =
  new G4Box("World", fWorldSizeX, fWorldSizeY, fWorldSizeZ);
  
  // Build a logic world with air as the material
  fLogicWorld = new G4LogicalVolume(solidWorld,
                                    FindMaterial("G4_AIR"),
                                    "World");
  
  // Build a phsics world
  fPhysiWorld = new G4PVPlacement(0,
                                  G4ThreeVector(),
                                  fLogicWorld,
                                  "World",
                                  0,
                                  false,
                                  0);
  
  //--------------------------------------------------
  // Extrusion: Extrude a volume the size fo scintillator to place scintillator inside
  //--------------------------------------------------
  
  G4VSolid* solidExtrusion =
  new G4Box("Extrusion",GetPanelBase(),GetPanelBase(),GetPanelWidth());
  
  // Make a logical volume using the coating material defined in MMG4Materials
  G4LogicalVolume* logicExtrusion =
  new G4LogicalVolume(solidExtrusion,
                      FindMaterial("Coating"),
                      "Extrusion");
  
  //Build optical surface using predefined opticalSurfaceModel, surface finish, and the material
  G4OpticalSurface* TiO2Surface = new G4OpticalSurface("TiO2Surface",
                                                       glisur,
                                                       ground,
                                                       dielectric_metal,
                                                       fExtrusionPolish);
  
  G4MaterialPropertiesTable* TiO2SurfaceProperty =
  new G4MaterialPropertiesTable();
  
  G4double p_TiO2[] = {2.00*eV, 3.47*eV};
  const G4int nbins = sizeof(p_TiO2)/sizeof(G4double);
  
  G4double refl_TiO2[] = {fExtrusionReflectivity,fExtrusionReflectivity};
  assert(sizeof(refl_TiO2) == sizeof(p_TiO2));
  G4double effi_TiO2[] = {0, 0};
  assert(sizeof(effi_TiO2) == sizeof(p_TiO2));
  
  TiO2SurfaceProperty -> AddProperty("REFLECTIVITY",p_TiO2,refl_TiO2,nbins);
  TiO2SurfaceProperty -> AddProperty("EFFICIENCY",p_TiO2,effi_TiO2,nbins);
  
  TiO2Surface -> SetMaterialPropertiesTable(TiO2SurfaceProperty);
  
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicExtrusion,
                    "Extrusion",
                    fLogicWorld,
                    false,
                    0);
  
  new G4LogicalSkinSurface("TiO2Surface",logicExtrusion,TiO2Surface);
  
  //--------------------------------------------------
  // Scintillator
  //--------------------------------------------------
  
  G4VSolid* solidScintillator = new G4Box("Scintillator",
                                          GetPanelBase()-GetCoatingThickness(),
                                          GetPanelBase()-GetCoatingThickness(),
                                          GetPanelWidth()-GetCoatingThickness());
  
  G4LogicalVolume* logicScintillator =
  new G4LogicalVolume(solidScintillator,
                      FindMaterial("Polystyrene"),
                      "Scintillator");
  
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logicScintillator,
                    "Scintillator",
                    logicExtrusion,
                    false,
                    0);
  
  if (GetCoatingThickness() > 0.*mm) {
    G4Box* outerSolidScintCoating = new G4Box("OuterCoatingOfPanel",
                                         GetPanelBase(),
                                         GetPanelBase(),
                                         GetPanelWidth());
    
    G4Box* innerSolidScintCoating = new G4Box("InnerCoatingOfPanel",
                                         GetPanelBase()-GetCoatingThickness(),
                                         GetPanelBase()-GetCoatingThickness(),
                                         GetPanelWidth()-GetCoatingThickness());
    
    G4SubtractionSolid* solidScintCoating = new G4SubtractionSolid("ScintCoating",
                                                 outerSolidScintCoating,
                                                 innerSolidScintCoating);
    G4LogicalVolume* logicScintCoating =
    new G4LogicalVolume(solidScintCoating,
                        FindMaterial("Polystyrene"),
                        "ScintCoating");
  
    new G4PVPlacement(0,
                      G4ThreeVector(0,0,0),
                      logicScintCoating,
                      "ScintCoating",
                      logicExtrusion,
                      false,
                      0);
  }
  
 //********************GROOVE EXTRUSION AND CONSTRUCTION********************************************************
  
//Groove Extrusion
//put into corresponding scintillator logical volume defined previously
//Groove1 Extrusion
G4VSolid* Groove1Extrusion = new G4Box("Extrusion",fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove1logicExtrusion = new G4LogicalVolume(Groove1Extrusion, FindMaterial("G4_AIR"), "Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), -1*fFiberLength1/2, (fPanelLength/2 - 25*mm)),
                    Groove1logicExtrusion, "Extrusion", logicScintillator, false, 0);

//Groove2 SemiCircle Extrusion
G4VSolid* SemicircGrooveExtrusion = new G4Tubs("Semicircle Groove Extrusion", fpRming, fpRmaxg, fpDzg, fpSPhig, fpDPhig);
G4VLogicalVolume* SemicricGrooveLogExt = new G4LogicalVolume(SemicircGrooveExtrusion, FindMaterial("G4_AIR"), 
      "Semicircle Groove Extrusion Logic");
G4VPhysicalVolume* SemicircGroovePhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2),0,0), 
      SemicircGrooveLogExt, "Semicircle Groove Extrusion Physics", logicScintillator, false, 0);

//Groove3 Extrusion
G4VSolid* Groove3Extrusion = new G4Box("Extrusion",fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove3logicExtrusion =  new G4LogicalVolume(Groove3Extrusion,FindMaterial("G4_AIR"),"Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), (-1*fFiberLength1/2), -1*(fPanelLength/2 - 25*mm)),
                    Groove3logicExtrusion, "Extrusion",logicScintillator,false,0);
    
//Groove4/Connectorized Groove Extrusion
G4VSolid* Groove4Extrusion = new G4Box("Extrusion", 1.5*mm, 12.5*mm, 1.5*mm); //3mm groove depth for connector
G4LogicalVolume* Groove4logicExtrusion =  new G4LogicalVolume(Groove4Extrusion,FindMaterial("G4_AIR"),"Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0, 
                    G4ThreeVector((fPanelWidth/2 - 1.5*mm), -1*(fFiberLength1 + 12.5*mm), -1*(fPanelLength/2 - 25*mm)),
                    Groove4logicExtrusion, "Extrusion",logicScintillator,false,0);


//Groove Construction 
//put into corresponding groove extrusion logical volume for G4PVPlacement
//Groove1 Construction 
G4VSolid* Groove1 = new G4Box("Groove for Fiber 1", fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove1Log = new G4LogicalVolume(Groove1, FindMaterial("PMMA"), "Groove for Fiber 1", 0,0,0);
G4VPhysicalVolume* Groove1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), 
      (-1*fFiberLength1/2), (fPanelLength/2 - 25*mm)), Groove1Log, "Groove1Phys", Groove1logicExtrusion, false, 0);

//Groove2 SemiCircle Construction
G4Vsolid* SemicircGroove = new G4Tubs("Semicircle Groove", fpRming, fpRmaxg, fpDzg, fpSPhig, fpDPhig);
G4VLogicalVolume* SemicricGrooveLog = new G4LogicalVolume(SemicircGroove, FindMaterial("PMMA"), 
      "Semicircle Groove", 0,0,0);
G4VPhysicalVolume* SemicircGroovePhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2),0,0), 
      SemicircGrooveLog, "SemicircGroove Physics", SemicricGrooveLogExt, false, 0);

//Groove3 Construction
G4VSolid* Groove3 = new G4Box("Groove for Fiber 3", fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove3Log = new G4LogicalVolume(Groove3, FindMaterial("PMMA"), "Groove for Fiber 3", 0,0,0);
G4VPhysicalVolume* Groove3Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), 
       -1*fFiberLength1/2, -1*(fPanelLength/2 - 25*mm)), Groove3Log, "Groove3Phys", Groove3logicExtrusion, false, 0);

//Groove4/Connectorized Groove Construction
G4VSolid* Groove4 = new G4Box("Groove for Connector", 1.5*mm, 12.5*mm, 1.5*mm); //3mm groove depth for connector
G4LogicalVolume* Groove4Log =  new G4LogicalVolume(Groove4,FindMaterial("PMMA"),"Groove for Fiber 4");
G4VPhysicalVolume* Groove4Phys = new G4PVPlacement(0, 
                    G4ThreeVector((fPanelWidth/2 - 1.5*mm), -1*(fFiberLength1 + 12.5*mm), -1*(fPanelLength/2 - 25*mm)),
                    Groove4Log, "Connectorized Groove or Fiber 4 Phys",Groove4logicExtrusion,false,0);


  
 

  
  //--------------------------------------------------
  // Fiber
  //--------------------------------------------------
  
  ConstructFiber();
  
  //--------------------------------------------------
  // End of Construction
  //--------------------------------------------------
  
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::ConstructFiber()
{
  if (!(fLogicHole) || !(fPhysiHole) ) {
    std::ostringstream o;
    o << "The Fiber Hole has not been constructed";
    G4Exception("MMG4DetectorConstruction::ConstructFiber","",
                FatalException,o.str().c_str());
  }
  

  // Pointers to the most recently constructed volume
  G4LogicalVolume* logicPlacement = fLogicHole;
  G4VPhysicalVolume* physiPlacement = fPhysiHole;
  
  //--------------------------------------------------
  // Fiber Construction
  //--------------------------------------------------
  /* Commented out cladding stuff because I'm not sure what it is/why we are making volumes for it.
  
  
  // Boundary Surface Properties
  G4OpticalSurface* opSurface = NULL;
  
  if (fSurfaceRoughness < 1.)
    opSurface = new G4OpticalSurface("RoughSurface",          // Surface Name
                                     glisur,                  // SetModel
                                     ground,                  // SetFinish
                                     dielectric_dielectric,   // SetType
                                     fSurfaceRoughness);      // SetPolish
  
  G4LogicalVolume   *logicClad1, *logicClad2;
  G4VPhysicalVolume *physiClad1, *physiClad2;
  
  // Determine the number of cladding layers to be built
  switch ( fNumOfCladLayers ) {
      
    case 2:
      
      //--------------------------------------------------
      // Cladding 2
      //--------------------------------------------------
      
      G4VSolid* solidClad2;
      
      if (fXYRatio == 1.)
        solidClad2 = new G4Tubs("Clad2",0.,fClad2RX,fClad2Length,0.0*rad,twopi*rad);
      else
        solidClad2 = new G4EllipticalTube("Clad2",fClad2RX,fClad2RY,fClad2Length);
      
      logicClad2  = new G4LogicalVolume(solidClad2,
                                        FindMaterial("FPethylene"),
                                        "Clad2");
      
      physiClad2 = new G4PVPlacement(0,
                                     G4ThreeVector(),
                                     logicClad2,
                                     "Clad2",
                                     logicPlacement,
                                     false,
                                     0);
      
      // Place the rough surface only if needed
      if (opSurface) {
        new G4LogicalBorderSurface("surfaceClad2Out",
                                   physiClad2,
                                   physiPlacement,
                                   opSurface);
        new G4LogicalBorderSurface("surfaceClad2In",
                                   physiPlacement,
                                   physiClad2,
                                   opSurface);
      }
      
      logicPlacement = logicClad2;
      physiPlacement = physiClad2;
      break;
      
    case 1:
      
      //--------------------------------------------------
      // Cladding 1
      //--------------------------------------------------
      
      G4VSolid* solidClad1;
      
      if (fXYRatio == 1.)
        solidClad1 = new G4Tubs("Clad1",0.,fClad1RX,fClad1Length,0.0*rad,twopi*rad);
      else
        solidClad1 = new G4EllipticalTube("Clad1",fClad1RX,fClad1RY,fClad1Length);
      
      logicClad1 = new G4LogicalVolume(solidClad1,
                                       FindMaterial("Pethylene"),
                                       "Clad1");
      
      physiClad1 = new G4PVPlacement(0,
                                     G4ThreeVector(0.0,0.0,0.0),
                                     logicClad1,
                                     "Clad1",
                                     logicPlacement,
                                     false,
                                     0);
      
      // Place the rough surface only if needed
      if (opSurface) {
        new G4LogicalBorderSurface("surfaceClad1Out",
                                   physiClad1,
                                   physiPlacement,
                                   opSurface);
        new G4LogicalBorderSurface("surfaceClad1In",
                                   physiPlacement,
                                   physiClad1,
                                   opSurface);
      }
      
      logicPlacement = logicClad1;
      physiPlacement = physiClad1;
      break;
      
    default:
    
    */ 
      
      //******************************************************BECAUSE FNUMCLADLAYERS=0, START HERE******************

      //--------------------------------------------------
      // MMG4 Fiber
      //--------------------------------------------------
      
      //Fiber Extrusion
//put into corresponding grooves logical volume
//Fiber1 Holes inside Groove Construction logic volumes (?)...change fLogicHole
G4VSolid* Fiber1Hole = new G4Tubs("Fiber1 Hole",0.0*cm,fFiberRadius, fFiberLength1,0.*deg, 360.*deg);
G4LogicalVolume* Fiber1HoleLog = new G4LogicalVolume(Fiber1Hole, FindMaterial("G4_AIR"), "Fiber1 Logical Hole"); 
    fPhysiHole = new G4PVPlacement(g4rot, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
                                                 (-1*fFiberLength1/2),(fPanelLength/2 - 25*mm)),Fiber1HoleLog,"Fiber1 Hole",
                                   Groove1Log, //want to put the hole inside the groove1 constructed logical volume
                                   false, 0);

//Fiber2 Semicircle Hole
G4VSolid* SemicircSolidHole = new G4Torus("Semicircle Fiber Region Hole", fpRmin, fpRmaz, fpRtor, fDPhif, fpSPhif);
G4VLogicalVolume* SemicircLogHole = new G4LogicalVolume(SemicircSolidHole, FindMaterial("G4_AIR"), 
                                                    "Semicircle Fiber Region LogicHole");
G4VPhysicalVolume* SemicircPhysHole = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 0,0), 
                                                        SemicircLogHole, "SemicircPhys", SemicricGrooveLog, false, 0);

//Fiber3 Hole
G4VSolid* Fiber3Hole = new G4Tubs("Fiber3 Hole", 0.0*cm,fFiberRadius, fFiberLength1,0.*deg,360.*deg);
G4LogicalVolume* Fiber3HoleLog = new G4LogicalVolume(Fiber3Hole,FindMaterial("G4_AIR"),"Fiber3 Logical Hole");
    fPhysiHole = new G4PVPlacement(g4rot, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), (-1*fFiberLength1/2),
                                   -1*(fPanelLength/2 - 25*mm)),Fiber3HoleLog, "Fiber3 Hole",
                                   Groove3Log, false, 0);
      
      
      
      
//Fiber Construction
//put into corresponding Fiber hole logical volume
//Fiber1
G4VSolid* Fiber1Solid = new G4Tubs("Fiber region 1", fFiberRx, fFiberLength1, 0.0*rad, twopi*rad);
G4LogicalVolume* Fiber1Log = new G4LogicalVolume(Fiber1Solid, FindMaterial("PMMA"), "Fiber1Log");
G4vPhysicalVolume* Fiber1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
      -1*fFiberLength1/2, (fPanelLength/2 - 25*mm)), Fiber1Log, "Fiber1 Physical Volume", Fiber1HoleLog, false, 0);

//Fiber2, Fiber Semicircle
G4VSolid* SemicircSolid = new G4Torus("Semicircle Fiber Region", fpRmin, fpRmaz, fpRtor, fDPhif, fpSPhif);
G4VLogicalVolume* SemicircLog = new G4LogicalVolume(SemicircSolid, FindMaterial("PMMA"), "Semicircle Fiber Region", 
       0,0,0);
G4VPhysicalVolume* SemicircPhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius),
       0,0), SemicircLog, "SemicircPhys", SemicircLogHole, false, 0);

//Fiber3
G4VSolid* Fiber1Solid = new G4Tubs("Fiber region 1", fFiberRx, fFiberLength1, 0.0*rad, twopi*rad);
G4LogicalVolume* Fiber1Log = new G4LogicalVolume(Fiber1Solid, FindMaterial("PMMA"), "Fiber1Log");
G4vPhysicalVolume* Fiber1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
      -1*fFiberLength1/2, -1*(fPanelLength/2 - 25*mm)), Fiber1Log, "Fiber1 Physical Volume", Fiber3HoleLog, false, 0);
      
      

      
      // Place the rough surface only if needed
   /*
      if (opSurface) {
        new G4LogicalBorderSurface("surfaceMMG4Out",
                                   physiMMG4fiber,
                                   physiPlacement,
                                   opSurface);
        new G4LogicalBorderSurface("surfaceMMG4In",
                                   physiPlacement,
                                   physiMMG4fiber,
                                   opSurface);
      }
  }
  
  */
  
  //--------------------------------------------------
  // Mirror for reflection at one of the end
  //--------------------------------------------------
  
  // Place the mirror only if the user wants the mirror
  /*
  if (fMirrorToggle) {
    
    G4VSolid* solidMirror = new G4Box("Mirror",
                                      fMirrorRmax,
                                      fMirrorRmax,
                                      fMirrorZ);
    
    G4LogicalVolume* logicMirror = new G4LogicalVolume(solidMirror,
                                                       FindMaterial("G4_Al"),
                                                       "Mirror");
    
    G4OpticalSurface* mirrorSurface = new G4OpticalSurface("MirrorSurface",
                                                           glisur,
                                                           ground,
                                                           dielectric_metal,
                                                           fMirrorPolish);
    
    G4MaterialPropertiesTable* mirrorSurfaceProperty =
    new G4MaterialPropertiesTable();
    
    G4double p_mirror[] = {2.00*eV, 3.47*eV};
    const G4int nbins = sizeof(p_mirror)/sizeof(G4double);
    G4double refl_mirror[] = {fMirrorReflectivity,fMirrorReflectivity};
    assert(sizeof(refl_mirror) == sizeof(p_mirror));
    G4double effi_mirror[] = {0, 0};
    assert(sizeof(effi_mirror) == sizeof(effi_mirror));
    
    mirrorSurfaceProperty->
    AddProperty("REFLECTIVITY",p_mirror,refl_mirror,nbins);
    mirrorSurfaceProperty->
    AddProperty("EFFICIENCY",p_mirror,effi_mirror,nbins);
    
    mirrorSurface -> SetMaterialPropertiesTable(mirrorSurfaceProperty);
    
    G4RotationMatrix* g4rot = new G4RotationMatrix();
    *g4rot = StringToRotationMatrix("Y90");
    *g4rot = g4rot->inverse();
    new G4PVPlacement(g4rot,
                      G4ThreeVector(0.0,0.0,fMirrorOrigin),
                      logicMirror,
                      "Mirror",
                      fLogicWorld,
                      false,
                      0);
    
    new G4LogicalSkinSurface("MirrorSurface",logicMirror,mirrorSurface);
  }
  */
  
  //--------------------------------------------------
  // Coupling at the read-out end
  //--------------------------------------------------
  
  // Clear Fiber (Coupling Layer)
  /*
  G4VSolid* solidCouple = new G4Box("Couple",fCoupleRX,fCoupleRY,fCoupleLength);
  
  G4LogicalVolume*   logicCouple = new G4LogicalVolume(solidCouple,
                                                       FindMaterial("G4_AIR"),
                                                       "Couple");
  
  G4RotationMatrix* g4rot = new G4RotationMatrix();
  *g4rot = StringToRotationMatrix("Y90");
  *g4rot = g4rot->inverse();
  if (*g4rot == G4RotationMatrix()) g4rot = NULL;
  
  new G4PVPlacement(g4rot,
                    G4ThreeVector(fCoupleOrigin,fMMG4fiberY,0.0),
                    logicCouple,
                    "Couple",
                    fLogicWorld,
                    false,
                    0);  */
  
  //--------------------------------------------------
  // A logical layer in front of PhotonDet
  //--------------------------------------------------
  
  // Purpose: Preventing direct dielectric to metal contact
      
  //********************************************************SiPM PLACEMENT************************************
  
      
  //Pranava's Original SiPM hole placement code
      /*
  // Check for valid placement of PhotonDet
  if (fSiPMTheta > std::atan(fSiPMDist / fSiPMHalfL)) {
    
    fSiPMTheta = 0;
    fSiPMOriginX  = std::sin(fSiPMTheta) * (fSiPMDist + fClrfiberLength);
    fSiPMOriginZ  = -fCoupleLength+std::cos(fSiPMTheta)*(fSiPMDist+fClrfiberLength);
    G4cerr << "Invalid alignment.  Alignment Reset to 0" << G4endl;
  }
  
  // Clear Fiber (Coupling Layer)
  G4VSolid* solidClrfiber;
  
  if ( fSiPMShape == "Square" )
    solidClrfiber =
    new G4Box("ClearFiber",fClrfiberHalfL,fClrfiberHalfL,fClrfiberLength);
  else
    solidClrfiber =
    new G4Tubs("ClearFiber",0.,fClrfiberHalfL,fClrfiberLength,0.0*rad,twopi*rad);
  
  G4LogicalVolume*   logicClrfiber =
  new G4LogicalVolume(solidClrfiber,
                      FindMaterial("G4_AIR"),
                      "ClearFiber");
  
  new G4PVPlacement(new G4RotationMatrix(CLHEP::HepRotationY(-fSiPMTheta)),
                    G4ThreeVector(fSiPMOriginX,0.0,fSiPMOriginZ),
                    logicClrfiber,
                    "ClearFiber",
                    logicCouple,
                    false,
                    0);
                    
                    */
//SiPM Extrusion
//put into World Logical Volume I guess? Because want an upright square?
  G4VSolid* SiPMExtrusion = new G4Box("SiPM Solid Extrusion", fSiPMbase/2, fSiPMthickness/2, fSiPMbase/2);
  G4LogicalVolume* SiPMLogExtrusion = new G4LogicalVolume(SiPMExtrusion, FindMaterial("G4_AIR"), "Extrusion");
G4vPhysicalVolume* SiPMPhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 + fSiPMbase/2), -1*(fFiberLength1 + fSiPMy),(fPanelLength/2 - 25*mm)),
                    SiPMLogExtrusion, "Extrusion", fLogicWorld, false, 0);

//SiPM construction
//put into SiPM Extrusion Logical Volume
G4VSolid* SiPMsolid = new G4Box("SiPM Solid", fSiPMbase/2, fSiPMthickness/2, fSiPMbase/2);
G4VLogicalVolume* SiPMlog = new G4LogicalVolume(SiPMsolid, FindMaterial("G4_Al"), "SiPM Solid", 0,0,0);
G4VPhysicalVolume* SiPMPhys = new G4PVPlacement(0, 
                     G4ThreeVector((fPanelWidth/2 + fSiPMbase/2), -1*(fFiberLength1 + fSiPMy),(fPanelLength/2 - 25*mm)), 
                     SiPMlog, "SiPMPhys", SiPMLogExtrusion, false, 0); 
      
      
      
  //--------------------------------------------------
  // PhotonDet (Sensitive Detector)
  //--------------------------------------------------
  
  // Physical Construction
  // PhotonDet Surface Properties
      /*GET THESE FROM GABE FOR THE SIPM
  G4OpticalSurface* photonDetSurface = new G4OpticalSurface("PhotonDetSurface",
                                                            glisur,
                                                            ground,
                                                            dielectric_metal,
                                                            fSiPMPolish);
  
  G4MaterialPropertiesTable* photonDetSurfaceProperty =
  new G4MaterialPropertiesTable();
  
  G4double p_SiPM[] = {2.00*eV, 3.47*eV};
  const G4int nbins = sizeof(p_SiPM)/sizeof(G4double);
  G4double refl_SiPM[] = {fSiPMReflectivity,fSiPMReflectivity};
  assert(sizeof(refl_SiPM) == sizeof(p_SiPM));
  G4double effi_SiPM[] = {1, 1};
  assert(sizeof(effi_SiPM) == sizeof(p_SiPM));
  
  photonDetSurfaceProperty->AddProperty("REFLECTIVITY",p_SiPM,refl_SiPM,nbins);
  photonDetSurfaceProperty->AddProperty("EFFICIENCY",p_SiPM,effi_SiPM,nbins);
  
  photonDetSurface->SetMaterialPropertiesTable(photonDetSurfaceProperty);
  
  new G4LogicalSkinSurface("PhotonDetSurface",logicPhotonDet,photonDetSurface);
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::ConstructSDandField()
{
  if (!fSiPMSD.Get()) {
    G4String SiPMSDName = "MMG4/PhotonDet";
    MMG4PhotonDetSD* SiPMSD = new MMG4PhotonDetSD(SiPMSDName);
    G4SDManager::GetSDMpointer()->AddNewDetector(SiPMSD);
    fSiPMSD.Put(SiPMSD);
  }
  SetSensitiveDetector("PhotonDet_LV", fSiPMSD.Get(), true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::UpdateGeometryParameters()
{
  fMMG4fiberRX  = fXYRatio * fMMG4fiberRY;
  
  //PTS: should be made configurable
  fClad1RX = fMMG4fiberRX + 0.03*fMMG4fiberRX;
  fClad1RY = fMMG4fiberRY + 0.03*fMMG4fiberRY;
  fClad1Length  = fMMG4fiberLength;
  
  //PTS: should be made configurable
  fClad2RX = fClad1RX + 0.03*fMMG4fiberRX;
  fClad2RY = fClad1RY + 0.03*fMMG4fiberRY;
  fClad2Length  = fMMG4fiberLength;
  
  //  fWorldSizeX = fClad2RX   + fSiPMDist + fSiPMHalfL + 1.*cm;
  //  fWorldSizeY = fClad2RY   + fSiPMDist + fSiPMHalfL + 1.*cm;
  //  fWorldSizeZ = fMMG4fiberLength + fSiPMDist + fSiPMHalfL + 1.*cm;
  //
  fWorldSizeX = fPanelBase*2;
  fWorldSizeY = fPanelBase*2;
  fWorldSizeZ = fPanelBase*2;
  
  fCoupleRX = fPanelWidth;
  fCoupleRY = fPanelWidth;
  fCoupleLength = (fPanelBase-fMMG4fiberLength)/2.0;
  
  fClrfiberHalfL = fSiPMHalfL;
  
  fMirrorRmax = fClad2RY;
  
  //  fCoupleOrigin = fMMG4fiberOrigin + fMMG4fiberLength + fCoupleLength;
  fCoupleOrigin = fMMG4fiberLength+fCoupleLength;
  fMirrorOrigin = fMMG4fiberOrigin - fMMG4fiberLength - fMirrorZ;
  fSiPMOriginX  = std::sin(fSiPMTheta) * (fSiPMDist + fClrfiberLength);
  fSiPMOriginZ  = -fCoupleLength + std::cos(fSiPMTheta) * (fSiPMDist + fClrfiberLength);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4RotationMatrix
MMG4DetectorConstruction::StringToRotationMatrix(G4String rotation)
{
  // We apply successive rotations OF THE OBJECT around the FIXED
  // axes of the parent's local coordinates; rotations are applied
  // left-to-right (rotation="r1,r2,r3" => r1 then r2 then r3).
  
  G4RotationMatrix rot;
  
  unsigned int place = 0;
  
  while (place < rotation.size()) {
    
    G4double angle;
    char* p;
    
    const G4String tmpstring=rotation.substr(place+1);
    
    angle = strtod(tmpstring.c_str(),&p) * deg;
    
    if (!p || (*p != (char)',' && *p != (char)'\0')) {
      G4cerr << "Invalid rotation specification: " <<
      rotation.c_str() << G4endl;
      return rot;
    }
    
    G4RotationMatrix thisRotation;
    
    switch(rotation.substr(place,1).c_str()[0]) {
      case 'X': case 'x':
        thisRotation = G4RotationMatrix(CLHEP::HepRotationX(angle));
        break;
      case 'Y': case 'y':
        thisRotation = G4RotationMatrix(CLHEP::HepRotationY(angle));
        break;
      case 'Z': case 'z':
        thisRotation = G4RotationMatrix(CLHEP::HepRotationZ(angle));
        break;
      default:
        G4cerr << " Invalid rotation specification: "
        << rotation << G4endl;
        return rot;
    }
    
    rot = thisRotation * rot;
    place = rotation.find(',',place);
    if (place > rotation.size()) break;
    ++place;
  }
  
  return rot;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPhotonDetGeometry (G4String shape)
// Set the Geometry of the PhotonDet detector
// Pre:  shape must be either "Circle" and "Square"
{
  if (shape == "Circle" || shape == "Square" ) fSiPMShape = shape;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetNumberOfCladding(G4int num)
// Set the number of claddings
// Pre: 0 <= num <= 2
{
  fNumOfCladLayers = num;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetWLSLength (G4double length)
// Set the TOTAL length of the MMG4 fiber
{
  fMMG4fiberLength = length;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetWLSRadius (G4double radius)
// Set the Y radius of MMG4 fiber
{
  fMMG4fiberRY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetClad1Radius (G4double radius)
// Set the Y radius of Cladding 1
{
  fClad1RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetClad2Radius (G4double radius)
// Set the Y radius of Cladding 2
{
  fClad2RY = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPhotonDetHalfLength(G4double halfL)
// Set the half length of the PhotonDet detector
// The half length will be the radius if PhotonDet is circular
{
  fSiPMHalfL = halfL;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetGap (G4double gap)
// Set the distance between fiber end and PhotonDet
{ 
  fSiPMDist = gap;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPhotonDetAlignment(G4double theta)
// Set the Aligment of PhotonDet with respect to the z axis
// If theta is 0 deg, then the detector is perfectly aligned
// PhotonDet will be deviated by theta from z axis
// facing towards the center of the fiber
{
  fSiPMTheta = theta;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetSurfaceRoughness(G4double roughness)
// Set the Surface Roughness between Cladding 1 and MMG4 fiber
// Pre: 0 < roughness <= 1
{
  fSurfaceRoughness = roughness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetMirrorPolish(G4double polish)
// Set the Polish of the mirror, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fMirrorPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetMirrorReflectivity(G4double reflectivity)
// Set the Reflectivity of the mirror, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fMirrorReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPhotonDetPolish(G4double polish)
// Set the Polish of the PhotonDet, polish of 1 is a perfect mirror surface
// Pre: 0 < polish <= 1
{
  fSiPMPolish = polish;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPhotonDetReflectivity(G4double reflectivity)
// Set the Reflectivity of the PhotonDet, reflectivity of 1 is a perfect mirror
// Pre: 0 < reflectivity <= 1
{
  fSiPMReflectivity = reflectivity;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetMirror(G4bool flag)
// Toggle to place the mirror or not at one end (-z end) of the fiber
// True means place the mirror, false means otherwise
{
  fMirrorToggle = flag;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetXYRatio(G4double r)
// Set the ratio of the x and y radius of the ellipse (x/y)
// a ratio of 1 would produce a circle
{
  fXYRatio = r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPanelWidth (G4double width)
// Set the length of the scintillator panel
{
  fPanelWidth = width;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetPanelBase (G4double side)
// Set the side of the scintillator panel
{
  fPanelBase = side;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetHoleRadius (G4double radius)
// Set the radius of the fiber hole
{
  fHoleRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MMG4DetectorConstruction::SetCoatingThickness (G4double thick)
// Set thickness of the coating on the panel
{
  fCoatingThickness = thick;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void MMG4DetectorConstruction::SetCoatingRadius (G4double radius)
// Set inner radius of the corner panel coating
{
  fCoatingRadius = radius;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetWLSFiberLength() { return fMMG4fiberLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetPanelWidth() { return fPanelWidth; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetPanelBase() { return fPanelBase; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetHoleRadius() { return fHoleRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetHoleLength() { return fHoleLength; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetFiberRadius() { return GetWLSFiberRMax(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetCoatingThickness()
{ return fCoatingThickness; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4double MMG4DetectorConstruction::GetCoatingRadius() { return fCoatingRadius; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetWLSFiberEnd()
{
  return fMMG4fiberOrigin + fMMG4fiberLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetWLSFiberRMax()
{
  if (fNumOfCladLayers == 2) return fClad2RY;
  if (fNumOfCladLayers == 1) return fClad1RY;
  return fMMG4fiberRY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double MMG4DetectorConstruction::GetSurfaceRoughness()
{
  return fSurfaceRoughness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Return True if the fiber construction is ideal
G4bool MMG4DetectorConstruction::IsPerfectFiber()
{
  return     fSurfaceRoughness == 1. && fXYRatio == 1.
  && (!fMirrorToggle    ||
      (fMirrorPolish    == 1. && fMirrorReflectivity == 1.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* MMG4DetectorConstruction::FindMaterial(G4String name) {
  G4Material* material = G4Material::GetMaterial(name,true);
  return material;
}










//Do I have to define fGrooveDepth etc. here if it is already assigned in a macro?
//DO I USE GETTER FUNCTIONS WHEN I WANT TO ACCESS MACRO STUFF???
//WHATS UP WITH THE G4ROT????? Do I only use that for like G4tubs and such? anything that has to do with radians?
/*
fpRming = 100 * mm
fpRmaxg = 103 * mm  
fpDzg = fGrooveDepth/2 * mm 
fpDPhig = 0 * Degree
fpSPhig = 180 * Degree 
fGrooveDepth = 3*mm;
fFiberLength1 = 100*mm;
fpRminf = 0 * mm
fpRmaxf = 1.00 * mm 
fpRtorf = 100 * mm 
fpDPhif = 0 * Degree
fpSPhif = 180 * Degree
fSiPMbase = 3*mm
fSiPMy = 12.5*mm //random location in the 25 mm zone between edge and fiber
//idk Material of SiPM; I used the one from Pranava's previous code 
//Groove Extrusion
//put into corresponding scintillator logical volume defined previously
//Groove1 Extrusion
G4VSolid* Groove1Extrusion = new G4Box("Extrusion",fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove1logicExtrusion = new G4LogicalVolume(Groove1Extrusion, FindMaterial("G4_AIR"), "Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), -1*fFiberLength1/2, (fPanelLength/2 - 25*mm)),
                    Groove1logicExtrusion, "Extrusion", logicScintillator, false, 0);
//Groove2 SemiCircle Extrusion
G4VSolid* SemicircGrooveExtrusion = new G4Tubs("Semicircle Groove Extrusion", fpRming, fpRmaxg, fpDzg, fpSPhig, fpDPhig);
G4VLogicalVolume* SemicricGrooveLogExt = new G4LogicalVolume(SemicircGrooveExtrusion, FindMaterial("G4_AIR"), 
      "Semicircle Groove Extrusion Logic");
G4VPhysicalVolume* SemicircGroovePhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2),0,0), 
      SemicircGrooveLogExt, "Semicircle Groove Extrusion Physics", logicScintillator, false, 0);
//Groove3 Extrusion
G4VSolid* Groove3Extrusion = new G4Box("Extrusion",fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove3logicExtrusion =  new G4LogicalVolume(Groove3Extrusion,FindMaterial("G4_AIR"),"Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), (-1*fFiberLength1/2), -1*(fPanelLength/2 - 25*mm)),
                    Groove3logicExtrusion, "Extrusion",logicScintillator,false,0);
    
//Groove4/Connectorized Groove Extrusion
G4VSolid* Groove4Extrusion = new G4Box("Extrusion", 1.5*mm, 12.5*mm, 1.5*mm); //3mm groove depth for connector
G4LogicalVolume* Groove4logicExtrusion =  new G4LogicalVolume(Groove4Extrusion,FindMaterial("G4_AIR"),"Extrusion");
G4VPhysicalVolume* Groove1PhysExtrusion = new G4PVPlacement(0, 
                    G4ThreeVector((fPanelWidth/2 - 1.5*mm), -1*(fFiberLength1 + 12.5*mm), -1*(fPanelLength/2 - 25*mm)),
                    Groove4logicExtrusion, "Extrusion",logicScintillator,false,0);
//Groove Construction 
//put into corresponding groove extrusion logical volume for G4PVPlacement
//Groove1 Construction 
G4VSolid* Groove1 = new G4Box("Groove for Fiber 1", fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove1Log = new G4LogicalVolume(Groove1, FindMaterial("PMMA"), "Groove for Fiber 1", 0,0,0);
G4VPhysicalVolume* Groove1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), 
      (-1*fFiberLength1/2), (fPanelLength/2 - 25*mm)), Groove1Log, "Groove1Phys", Groove1logicExtrusion, false, 0);
//Groove2 SemiCircle Construction
G4Vsolid* SemicircGroove = new G4Tubs("Semicircle Groove", fpRming, fpRmaxg, fpDzg, fpSPhig, fpDPhig);
G4VLogicalVolume* SemicricGrooveLog = new G4LogicalVolume(SemicircGroove, FindMaterial("PMMA"), 
      "Semicircle Groove", 0,0,0);
G4VPhysicalVolume* SemicircGroovePhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2),0,0), 
      SemicircGrooveLog, "SemicircGroove Physics", SemicricGrooveLogExt, false, 0);
//Groove3 Construction
G4VSolid* Groove3 = new G4Box("Groove for Fiber 3", fGrooveDepth/2, fFiberLength1/2, fGrooveDepth/2);
G4LogicalVolume* Groove3Log = new G4LogicalVolume(Groove3, FindMaterial("PMMA"), "Groove for Fiber 3", 0,0,0);
G4VPhysicalVolume* Groove3Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth/2), 
       -1*fFiberLength1/2, -1*(fPanelLength/2 - 25*mm)), Groove3Log, "Groove3Phys", Groove3logicExtrusion, false, 0);
//Groove4/Connectorized Groove Construction
G4VSolid* Groove4 = new G4Box("Groove for Connector", 1.5*mm, 12.5*mm, 1.5*mm); //3mm groove depth for connector
G4LogicalVolume* Groove4Log =  new G4LogicalVolume(Groove4,FindMaterial("PMMA"),"Groove for Fiber 4");
G4VPhysicalVolume* Groove4Phys = new G4PVPlacement(0, 
                    G4ThreeVector((fPanelWidth/2 - 1.5*mm), -1*(fFiberLength1 + 12.5*mm), -1*(fPanelLength/2 - 25*mm)),
                    Groove4Log, "Connectorized Groove or Fiber 4 Phys",Groove4logicExtrusion,false,0);
//Fiber Extrusion
//put into corresponding grooves logical volume
//Fiber1 Holes inside Groove Construction logic volumes (?)...change fLogicHole
G4VSolid* Fiber1Hole = new G4Tubs("Fiber1 Hole",0.0*cm,fFiberRadius, fFiberLength1,0.*deg, 360.*deg);
G4LogicalVolume* Fiber1HoleLog = new G4LogicalVolume(Fiber1Hole, FindMaterial("G4_AIR"), "Fiber1 Logical Hole"); 
    fPhysiHole = new G4PVPlacement(g4rot, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
                                                 (-1*fFiberLength1/2),(fPanelLength/2 - 25*mm)),Fiber1HoleLog,"Fiber1 Hole",
                                   Groove1Log, //want to put the hole inside the groove1 constructed logical volume
                                   false, 0);
//Fiber2 Semicircle Hole
G4VSolid* SemicircSolidHole = new G4Torus("Semicircle Fiber Region Hole", fpRmin, fpRmaz, fpRtor, fDPhif, fpSPhif);
G4VLogicalVolume* SemicircLogHole = new G4LogicalVolume(SemicircSolidHole, FindMaterial("G4_AIR"), 
                                                    "Semicircle Fiber Region LogicHole");
G4VPhysicalVolume* SemicircPhysHole = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 0,0), 
                                                        SemicircLogHole, "SemicircPhys", SemicricGrooveLog, false, 0);
//Fiber3 Hole
G4VSolid* Fiber3Hole = new G4Tubs("Fiber3 Hole", 0.0*cm,fFiberRadius, fFiberLength1,0.*deg,360.*deg);
G4LogicalVolume* Fiber3HoleLog = new G4LogicalVolume(Fiber3Hole,FindMaterial("G4_AIR"),"Fiber3 Logical Hole");
    fPhysiHole = new G4PVPlacement(g4rot, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), (-1*fFiberLength1/2),
                                   -1*(fPanelLength/2 - 25*mm)),Fiber3HoleLog, "Fiber3 Hole",
                                   Groove3Log, false, 0);
//Fiber Construction
//put into corresponding Fiber hole logical volume
//Fiber1
G4VSolid* Fiber1Solid = new G4Tubs("Fiber region 1", fFiberRx, fFiberLength1, 0.0*rad, twopi*rad);
G4LogicalVolume* Fiber1Log = new G4LogicalVolume(Fiber1Solid, FindMaterial("PMMA"), "Fiber1Log");
G4vPhysicalVolume* Fiber1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
      -1*fFiberLength1/2, (fPanelLength/2 - 25*mm)), Fiber1Log, "Fiber1 Physical Volume", Fiber1HoleLog, false, 0);
//Fiber2, Fiber Semicircle
G4VSolid* SemicircSolid = new G4Torus("Semicircle Fiber Region", fpRmin, fpRmaz, fpRtor, fDPhif, fpSPhif);
G4VLogicalVolume* SemicircLog = new G4LogicalVolume(SemicircSolid, FindMaterial("PMMA"), "Semicircle Fiber Region", 
       0,0,0);
G4VPhysicalVolume* SemicircPhys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius),
       0,0), SemicircLog, "SemicircPhys", SemicircLogHole, false, 0);
//Fiber3
G4VSolid* Fiber1Solid = new G4Tubs("Fiber region 1", fFiberRx, fFiberLength1, 0.0*rad, twopi*rad);
G4LogicalVolume* Fiber1Log = new G4LogicalVolume(Fiber1Solid, FindMaterial("PMMA"), "Fiber1Log");
G4vPhysicalVolume* Fiber1Phys = new G4PVPlacement(0, G4ThreeVector((fPanelWidth/2 - fGrooveDepth + fFiberRadius), 
      -1*fFiberLength1/2, -1*(fPanelLength/2 - 25*mm)), Fiber1Log, "Fiber1 Physical Volume", Fiber3HoleLog, false, 0);
 //SiPM Extrusion
//Put SiPM as an upright square in the connectorized groove/Groove4
  G4VSolid* SiPMExtrusion = new G4Box("SiPM Solid Extrusion", fSiPMbase/2, fSiPMthickness/2, fSiPMbase/2);
  G4LogicalVolume* SiPMLogExtrusion = new G4LogicalVolume(SiPMExtrusion, FindMaterial("G4_AIR"), "Extrusion");
G4vPhysicalVolume* SiPMPhysExtrusion = new G4PVPlacement(0,
                    G4ThreeVector((fPanelWidth/2 - fSiPMbase/2), -1*(fFiberLength1 + fSiPMy),(fPanelLength/2 - 25*mm)),
                    SiPMLogExtrusion, "Extrusion", Groove4Log, false, 0);
//SiPM construction
//put into SiPM Extrusion Logical Volume
G4VSolid* SiPMsolid = new G4Box("SiPM Solid", fSiPMbase/2, fSiPMthickness/2, fSiPMbase/2);
G4VLogicalVolume* SiPMlog = new G4LogicalVolume(SiPMsolid, FindMaterial("G4_Al"), "SiPM Solid", 0,0,0);
G4VPhysicalVolume* SiPMPhys = new G4PVPlacement(0, 
                     G4ThreeVector((fPanelWidth/2 - fSiPMbase/2), -1*(fFiberLength1 + fSiPMy),(fPanelLength/2 - 25*mm)), 
                     SiPMlog, "SiPMPhys", SiPMLogExtrusion, false, 0); 
*/
