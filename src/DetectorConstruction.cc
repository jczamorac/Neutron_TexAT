// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Implements mandatory user class DetectorConstruction.
*/
// Local headers
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "MagneticField.hh"
#include "SensitiveDetector.hh"

// Default headers
#include <sstream>

// Geant4 headers
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Transform3D.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

using namespace std;
using namespace CLHEP;

SensitiveDetector *Sensitive = (SensitiveDetector *)0;

//--------------------------------------------------------------------------------------------------------------//

DetectorConstruction::DetectorConstruction()
{
  //Create a messanger (defines custom UI commands)
  messenger = new DetectorMessenger(this);

  //--------- Material definition ---------//
  DefineMaterials();

  //--------- Sizes of the principal geometrical components (solids)  ---------//
  ComputeParameters();
  ConstructSetup();
  magneticField = new MagneticField();
}

//--------------------------------------------------------------------------------------------------------------//

DetectorConstruction::~DetectorConstruction()
{
  delete messenger;
}

//--------------------------------------------------------------------------------------------------------------//

Detector::Detector(const G4int DetectorNumber, const G4double Height, const G4double Lenght, const G4double Width, G4int nStrips)
    : rDetector(DetectorNumber), rHeight(Height), rLenght(Lenght), rWidth(Width), rStrips(nStrips)
{
  // Getting Detector Name
  // e.g
  // Detector number is X so its name is going to be Detector_X and its strips SensorStripD0X
  ostringstream DetectorN, StripsN;

  DetectorN << "Detector " << rDetector;
  StripsN << "SensorStripD0" << rDetector;

  DetectorName = DetectorN.str().data();
  StripsName = StripsN.str().data();
}

//--------------------------------------------------------------------------------------------------------------//

Detector::~Detector()
{
}

//--------------------------------------------------------------------------------------------------------------//

void Detector::Rotate(G4double RotX, G4double RotY, G4double RotZ)
{
  // Angles in degree
  RotationMatrix->rotateX(RotX * deg);
  RotationMatrix->rotateY(RotY * deg);
  RotationMatrix->rotateZ(RotZ * deg);
}

//--------------------------------------------------------------------------------------------------------------//

G4VPhysicalVolume *Detector::Construct(G4LogicalVolume *LogicalMother)
{
  // Detector Material
  G4Material *Material = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

  // Creating a Physical Volume
  G4Box *DetectorBox = new G4Box(DetectorName, rLenght / 2, rHeight / 2, rWidth / 2);

  // Creating a Logical Volume
  G4LogicalVolume *LogVolume = new G4LogicalVolume(DetectorBox, Material, DetectorName);

  // Placing Detector
  G4VPhysicalVolume *PhysicalVolume = new G4PVPlacement(RotationMatrix, // Rotation Matrix
                                                        Position,       // Position
                                                        LogVolume,      // Respective Logical Volume
                                                        DetectorName,   // Name
                                                        LogicalMother,  // Logical Mother Volume
                                                        false,          // No Boolean Operation
                                                        rDetector,      // Copy Number
                                                        true);          // IDK

  // Strips parameters
  G4double StripHeight = rHeight / 2;
  G4double StripLenght = rLenght / (2 * rStrips);
  G4double StripWidth = rWidth / 2;

  // Creating Physical Volume
  G4Box *StripsBox = new G4Box(StripsName, StripLenght, StripHeight, StripWidth);

  // Creating LogicalVolume
  G4LogicalVolume *LogStripVolume = new G4LogicalVolume(StripsBox, Material, StripsName);

  //Placing Strips
  G4VPhysicalVolume *PhysicalStripsVolume = new G4PVReplica(StripsName,
                                                            LogStripVolume,
                                                            LogVolume,
                                                            kXAxis,
                                                            rStrips,
                                                            2.0 * StripLenght);

  // Set "User Limits"
  G4UserLimits *userLimits = new G4UserLimits(0.1 * mm);
  LogStripVolume->SetUserLimits(userLimits);

  // Setting Detector Sensitive
  if (!Sensitive)
  {
    Sensitive = new SensitiveDetector("/myDet/SiStripSD");
    //We register now the SD with the manager
    G4SDManager::GetSDMpointer()->AddNewDetector(Sensitive);
  }
  LogStripVolume->SetSensitiveDetector(Sensitive);

  // Color
  LogStripVolume->SetVisAttributes(new G4VisAttributes(Color));
  LogVolume->SetVisAttributes(new G4VisAttributes(Color));

  // Return Your Beaultiful Detector
  return PhysicalVolume;
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::DefineMaterials()
{

  G4double a, z, density; // z=mean number of protons;
  G4double temperature, pressure;
  G4int ncomponents, natoms;
  G4String name, symbol; // a=mass of a CLHEP::mole;

  // Define Elements
  a = 1.01 * CLHEP::g / CLHEP::mole;
  G4Element *elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1., a);
  density = 1.33e-11 * CLHEP::g / CLHEP::cm3;
  pressure = 1.0913e-10 * CLHEP::atmosphere;
  temperature = 200. * CLHEP::kelvin;
  H2 = new G4Material(name = "Hydrogen gas", density, ncomponents = 1,
                      kStateGas, temperature, pressure);
  H2->AddElement(elH, natoms = 2);

  // Define the AlN
  G4Element* elAl  = new G4Element(name="Aluminum",symbol="Al" , z= 13., a=26.98*g/mole);
  G4Element* elN  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a=14.00*g/mole);
    density     = 3.26*g/cm3;
  AlN = new G4Material(name="Aluminum Nitrode", density, ncomponents=2);
  AlN->AddElement(elAl, natoms=1);
  AlN->AddElement(elN, natoms=1);


  G4Element *elC = new G4Element(name = "Carbono", symbol = "C", z = 6., a = 12.01 * CLHEP::g / CLHEP::mole);
  G4Element *elHi = new G4Element(name = "Hidrogênio", symbol = "H", z = 1., a = 1.01 * CLHEP::g / CLHEP::mole);
  density = 0.98 * CLHEP::kg / CLHEP::m3;
  temperature = 200. * CLHEP::kelvin;
  pressure = 1. * CLHEP::atmosphere;
  CH4 = new G4Material(name = "Metano", density, ncomponents = 2, kStateGas, temperature, pressure);
  CH4->AddElement(elC, natoms = 1);
  CH4->AddElement(elHi, natoms = 4);

  // Define POLYETHYLENE deuterado
  G4Element *deuteron = new G4Element(name = "Deuteron", symbol = "D", z = 1., a = 2. * CLHEP::g / CLHEP::mole);
  CD2 = new G4Material(name = "CD2", 0.94 * g / cm3, ncomponents = 2, kStateSolid); /* , temperature, pressure); */
  CD2->AddElement(elC, natoms = 1);
  CD2->AddElement(elHi, natoms = 2);

  // Get Materials from NIST database
  G4NistManager *man = G4NistManager::Instance();
  man->SetVerbose(0);

  G4Element *Helio4 = new G4Element(name = "Helio4", symbol = "He4", z = 2., 4 * CLHEP::g / CLHEP::mole);
  G4double He4Density = 2.18e-5 * CLHEP::g / CLHEP::cm3;
  G4double He4Pressure = 0.13 * CLHEP::atmosphere;
  G4double He4Temperature = 293. * CLHEP::kelvin;
  He4 = new G4Material(name = "Helio4 gas", He4Density, ncomponents = 1,
                       kStateGas, He4Temperature, He4Pressure);
  He4->AddElement(Helio4, natoms = 1);

  G4Element *Argon40 = new G4Element(name = "Argon40", symbol = "Ar40", z = 18., 40 * CLHEP::g / CLHEP::mole);
  G4double Argon40Density = 6.6e-4 * CLHEP::g / CLHEP::cm3;
  G4double Argon40Pressure = 0.39 * CLHEP::atmosphere;
  G4double Argon40Temperature = 293. * CLHEP::kelvin;
  Ar40 = new G4Material(name = "Argon40 gas", Argon40Density, ncomponents = 1,
                       kStateGas, Argon40Temperature, Argon40Pressure);
  Ar40->AddElement(Argon40, natoms = 1);

  // Define NIST materials
  air = man->FindOrBuildMaterial("G4_AIR");
  silicon = man->FindOrBuildMaterial("G4_Si");
  vacuum = man->FindOrBuildMaterial("G4_Galactic");
  steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  tungsten = man->FindOrBuildMaterial("G4_W");
  tantalum = man->FindOrBuildMaterial("G4_Ta");
  sio2 = man->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  CH2 = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  alumi = man->FindOrBuildMaterial("G4_Al");
}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::ComputeParameters()
{
  // This function defines the defaults of the geometry construction
  // World size
  halfWorldLength = 20.0 * CLHEP::m;

  // Number of strips in the detector
  noOfSensorStrips = 16;

  // Detector Parameters
  Lengthy_neudet_t1 = 20*16* CLHEP::mm;
  Lengthx_neudet_t1 = 20*6* CLHEP::mm;
  Thickness_neudet_t1 = 20. * CLHEP::mm;
  Lengthy2_neudet_t1 = 20*4* CLHEP::mm;
  Thickness_neudet_t2 = 500. * CLHEP::um;

  // Solenoid 1 parameters
  Solenoid1_lenght = 100.0 * CLHEP::cm;
  Solenoid1_inner_diameter = 30.05 * CLHEP::cm;
  Solenoid1_outer_diameter = 100.0 * CLHEP::cm;
  Solenoid_pos[0] = G4ThreeVector(0., 0., 0.) * cm;

  // Solenoid 2 parameters
  Solenoid2_inner_diameter = 30.05 * CLHEP::cm;
  Solenoid2_outer_diameter = 100.0 * CLHEP::cm;
  Solenoid2_lenght = 100.0 * CLHEP::cm;
  Solenoid_pos[1] = G4ThreeVector(0., 0., 301.) * cm;

  // Solenoid 3 parameters
  Solenoid3_lenght = 150.0 * CLHEP::cm;
  Solenoid3_inner_diameter = 50.05 * CLHEP::cm;
  Solenoid3_outer_diameter = 70.0 * CLHEP::cm;
  Solenoid_pos[2] = G4ThreeVector(0., 0., 501.) * cm;

  // Magnetic field parameters
  Mag1_diameter = 30.0 * cm;
  Mag1_lenght = 140.0 * cm;
  Mag2_diameter = 30.0 * cm;
  Mag2_lenght = 140.0 * cm;
  Mag3_diameter = 50.0 * cm;
  Mag3_lenght = 90.0 * cm;

  // Target parameters
  Target_height = 1 * cm;
  Target_lenght = 1 * cm;

  // dummy parameters
  dummy_height = 100 * mm;
  dummy_lenght = 50*5 * mm;
  ang_det0 = 90.0;
  ang_det1 = 90;
  ang_det2 = -90;
  ang_det3 = -90;
  ang_det4 = 180;
  ang_det5 = 180;
  ang_det6 = 0;
  ang_det7 = 0;
  ang_det8 = 0;
  ang_det9 = 0;

  // Detector position
  G4double rPosition_z = -20.;
  G4double distaDet = 120*mm;
  G4double distaDet0 = 280*mm +10*mm;
  G4double coorx_D_0 = -(distaDet0)*cos((90-ang_det0)*3.141562/180.);
  G4double coory_D_0 = 0.*mm;
  G4double coorz_D_0 = -(distaDet0)*sin((90-ang_det0)*3.141562/180.);
  G4double distaDet1 = 280*mm +10*mm + 20*mm + 10*mm;
  G4double coorx_D_1 = -(distaDet1)*cos((90.- ang_det1)*3.141562/180.);
  G4double coory_D_1 = 0.*mm;
  G4double coorz_D_1 = -(distaDet1)*sin((90.- ang_det1)*3.141562/180.);
  G4double distaDet2 = 280*mm +10*mm;
  G4double coorx_D_2 = -(distaDet2)*cos((90. - ang_det2)*3.141562/180.);
  G4double coory_D_2 = 0.*mm;
  G4double coorz_D_2 = -(distaDet2)*sin((90.- ang_det2)*3.141562/180.);
  G4double distaDet3 = 280*mm +10*mm + 20*mm + 10*mm;
  G4double coorx_D_3 = -(distaDet3)*cos((90.- ang_det3)*3.141562/180.);
  G4double coory_D_3 = 0.*mm;
  G4double coorz_D_3 = -(distaDet3)*sin((90.- ang_det3)*3.141562/180.);
  G4double distaDet4 = 280*mm +10*mm;
  G4double coorx_D_4 = -(distaDet4)*cos((90.- ang_det4)*3.141562/180.);
  G4double coory_D_4 = 0.*mm;
  G4double coorz_D_4 = -(distaDet4)*sin((90.- ang_det4)*3.141562/180.);
  G4double distaDet5 = 280*mm +10*mm + 20*mm + 10*mm;
  G4double coorx_D_5 = -(distaDet5)*cos((90.- ang_det5)*3.141562/180.);
  G4double coory_D_5 = 0.*mm;
  G4double coorz_D_5 = -(distaDet5)*sin((90.- ang_det5)*3.141562/180.);
  G4double distaDet6 = 280*mm +10*mm;
  G4double coorx_D_6 = -(distaDet6)*cos((90.- ang_det6)*3.141562/180.) - 120*mm;
  G4double coory_D_6 = 0.*mm;
  G4double coorz_D_6 = -(distaDet6)*sin((90.- ang_det6)*3.141562/180.);
  G4double distaDet7 = 280*mm +10*mm + 20*mm + 10*mm;
  G4double coorx_D_7 = -(distaDet7)*cos((90.- ang_det7)*3.141562/180.) - 120*mm;
  G4double coory_D_7 = 0.*mm;
  G4double coorz_D_7 = -(distaDet7)*sin((90.- ang_det7)*3.141562/180.);
  G4double distaDet8 = 280*mm +10*mm;
  G4double coorx_D_8 = -(distaDet8)*cos((90.- ang_det8)*3.141562/180.) + 120*mm;
  G4double coory_D_8 = 0.*mm;
  G4double coorz_D_8 = -(distaDet8)*sin((90.- ang_det8)*3.141562/180.);
  G4double distaDet9 = 280*mm +10*mm + 20*mm + 10*mm;
  G4double coorx_D_9 = -(distaDet9)*cos((90.- ang_det9)*3.141562/180.) + 120*mm;
  G4double coory_D_9 = 0.*mm;
  G4double coorz_D_9 = -(distaDet9)*sin((90.- ang_det9)*3.141562/180.);

  // Rear detectors
  DetectorPosition[0] = G4ThreeVector(coorx_D_0, coory_D_0, coorz_D_0);
  DetectorPosition[1] = G4ThreeVector(coorx_D_1, coory_D_1, coorz_D_1);
  DetectorPosition[2] = G4ThreeVector(coorx_D_2, coory_D_2, coorz_D_2);
  DetectorPosition[3] = G4ThreeVector(coorx_D_3, coory_D_3, coorz_D_3);
  DetectorPosition[4] = G4ThreeVector(coorx_D_4, coory_D_4, coorz_D_4);
  DetectorPosition[5] = G4ThreeVector(coorx_D_5, coory_D_5, coorz_D_5);
  DetectorPosition[6] = G4ThreeVector(coorx_D_6, coory_D_6, coorz_D_6);
  DetectorPosition[7] = G4ThreeVector(coorx_D_7, coory_D_7, coorz_D_7);
  DetectorPosition[8] = G4ThreeVector(coorx_D_8, coory_D_8, coorz_D_8);
  DetectorPosition[9] = G4ThreeVector(coorx_D_9, coory_D_9, coorz_D_9);






}

//--------------------------------------------------------------------------------------------------------------//

void DetectorConstruction::ConstructSetup(void)
{
}

//--------------------------------------------------------------------------------------------------------------//

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  // Retrieving inputs
  Inputs *Inputs = &Inputs::GetInputs();

  G4NistManager *man = G4NistManager::Instance();
  man->SetVerbose(0);

  // Available colors
  G4Color green(0.0, 1.0, 0.0),
      red(1.0, 0.0, 0.0),
      yellow(1.0, 1.0, 0.0),
      cyan(0.2, 0.5, 0.5),
      blue(0.0, 0.0, 1.0);

  //---- This function is called by G4 when the detector has to be created -----//
  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------//

  //----------------------------------------------------//
  // World - Aqui o "mundo" é criado. O mundo é Lógico. //
  //----------------------------------------------------//

  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.0 * halfWorldLength);

  G4Box *solidWorld = new G4Box("world", halfWorldLength, halfWorldLength, halfWorldLength);
  logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World", 0, 0, 0);

  //  Must place the World Physical volume unrotated at (0,0,0).
  G4VPhysicalVolume *physiWorld = new G4PVPlacement(0,               // no rotation
                                                    G4ThreeVector(), // at (0,0,0)
                                                    logicWorld,      // its logical volume
                                                    "World",         // its name
                                                    0,               // its mother  volume
                                                    false,           // no boolean operations
                                                    0);              // copy number

  logicWorld->SetVisAttributes(G4VisAttributes::Invisible); // Turn the world invisible

  // ---------------------------------------------------------------------------------------------------- //


  //------- Building Target -------//

  //////////////
  //  Target  //
  //////////////

  // Target position
  G4ThreeVector Target_pos = Inputs->target_pos;

  // Target Material
  //G4Material *TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Au");
  //G4Material *TargetMaterial = G4NistManager::Instance()->FindOrBuildMaterial(Inputs->g4_material_name);
  //Inputs->TargetMaterial = TargetMaterial;

   //G4cout<<TargetMaterial<<"  "<<Inputs->width<<G4endl;



  G4ThreeVector chamber_out = G4ThreeVector(558/2*mm,368/2*mm,558/2*mm);
  G4ThreeVector chamber_thickness = G4ThreeVector(25.4/2*mm,25.4/2*mm,25.4/2*mm);
  G4ThreeVector chamber_in = chamber_out - chamber_thickness;

  // Creating Target
  Sol_Target = new G4Box("Sol_Target", chamber_in.x(),chamber_in.y(),chamber_in.z());

  // Creating dummy
  Sol_dummy = new G4Box("Sol_dummy",  dummy_lenght / 2, dummy_height/2,  0.5/ 2 * mm);
  G4ThreeVector dummy_pos = G4ThreeVector(0,0,194) *mm;

  G4SubtractionSolid* gas_wo_si =
      new G4SubtractionSolid("silicons",  Sol_Target, Sol_dummy,0,dummy_pos);




  Log_Target = new G4LogicalVolume(gas_wo_si, Ar40, "Log_Target");

  // Placing Target
  //Phys_Target = new G4PVPlacement(rotation, Target_pos, Log_Target, "Target", Log_Magnet2, false, 8, true);
  Phys_Target = new G4PVPlacement(0, Target_pos, Log_Target, "Target", logicWorld, false, 8, true);

  // Color
  Log_Target->SetVisAttributes(new G4VisAttributes(red));



//--- AT Chamber
G4Box * solid_in_bx = new G4Box("In_Bx",
                           chamber_in.x(),chamber_in.y(),chamber_in.z());


G4Box * solid_out_bx = new G4Box("Out_Bx",
                          chamber_out.x(),chamber_out.y(),chamber_out.z());




G4SubtractionSolid* solid_ch_bx =
        new G4SubtractionSolid("chamb_bx", solid_out_bx, solid_in_bx);

G4LogicalVolume * logic_ch_bx = new G4LogicalVolume(solid_ch_bx, // its solid
                                                        steel,     //its material
                                                      "Ch_bx");   //its name


 phys_ch_bx =
    new G4PVPlacement(0,
        G4ThreeVector(0,0,0),
              logic_ch_bx,
               "Cham_bx",
              logicWorld,
                false,
                0);


//---------------------dummy detectors
double angledet = 0;
// dummy position

//G4ThreeVector dummy_pos = Target_pos + G4ThreeVector(sin(angledet*3.1415/180)*89.8,0,cos(angledet*3.1415/180)*89.8) *mm;

// dummy Material
G4Material *dummyMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");


// Rotating dummy
G4RotationMatrix *rotationdummy = new G4RotationMatrix;
rotationdummy->rotateY(-1*angledet * CLHEP::deg);



Log_dummy = new G4LogicalVolume(Sol_dummy, dummyMaterial, "Log_dummy");

// Placing dummy
Phys_dummy = new G4PVPlacement(rotationdummy, dummy_pos, Log_dummy, "dummy", logicWorld, false, 8, true);

// Color
Log_dummy->SetVisAttributes(new G4VisAttributes(blue));

//------------------


  //-------------------------------------- Building All Detectors ----------------------------------------//

  ///////////////
  // Detectors //
  ///////////////

  Detector Detector_0(0, Lengthy_neudet_t1, Lengthx_neudet_t1,  Thickness_neudet_t1, noOfSensorStrips);
  Detector_0.Rotate(0., -1*ang_det0, 90);
  Detector_0.SetPosition(DetectorPosition[0]);
  Detector_0.SetColor(cyan);
  //Detector_0.Construct(Log_Magnet2); // This method requires a Logical Mother Volume
  Detector_0.Construct(logicWorld);

  Detector Detector_1(1, Lengthy_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_1.Rotate(0., ang_det1, 90.);
  Detector_1.SetPosition(DetectorPosition[1]);
  Detector_1.SetColor(cyan);
  //Detector_1.Construct(Log_Magnet2); // This method requires a Logical Mother Volume
  Detector_1.Construct(logicWorld);

  Detector Detector_2(2, Lengthy_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_2.Rotate(0., -1*ang_det2, 90.);
  Detector_2.SetPosition(DetectorPosition[2]);
  Detector_2.SetColor(cyan);
  //Detector_2.Construct(Log_Magnet2); // This method requires a Logical Mother Volume
  Detector_2.Construct(logicWorld);

  Detector Detector_3(3, Lengthy_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_3.Rotate(0., -1*ang_det3, 90.);
  Detector_3.SetPosition(DetectorPosition[3]);
  Detector_3.SetColor(cyan);
  //Detector_3.Construct(Log_Magnet2); // This method requires a Logical Mother Volume
  Detector_3.Construct(logicWorld);

  Detector Detector_4(4, Lengthy_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_4.Rotate(0., -1*ang_det4, 90.);
  Detector_4.SetPosition(DetectorPosition[4]);
  Detector_4.SetColor(cyan);
  Detector_4.Construct(logicWorld);

  Detector Detector_5(5, Lengthy_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_5.Rotate(0., -1*ang_det5, 90.);
  Detector_5.SetPosition(DetectorPosition[5]);
  Detector_5.SetColor(cyan);
  Detector_5.Construct(logicWorld);

  Detector Detector_6(6, Lengthy2_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_6.Rotate(0., -1*ang_det6, 90.);
  Detector_6.SetPosition(DetectorPosition[6]);
  Detector_6.SetColor(cyan);
  Detector_6.Construct(logicWorld);

  Detector Detector_7(7, Lengthy2_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_7.Rotate(0., -1*ang_det7, 90.);
  Detector_7.SetPosition(DetectorPosition[7]);
  Detector_7.SetColor(cyan);
  Detector_7.Construct(logicWorld);

  Detector Detector_8(8, Lengthy2_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_8.Rotate(0., -1*ang_det8, 90.);
  Detector_8.SetPosition(DetectorPosition[8]);
  Detector_8.SetColor(cyan);
  Detector_8.Construct(logicWorld);

  Detector Detector_9(9, Lengthy2_neudet_t1, Lengthx_neudet_t1, Thickness_neudet_t1, noOfSensorStrips);
  Detector_9.Rotate(0., -1*ang_det9, 90.);
  Detector_9.SetPosition(DetectorPosition[9]);
  Detector_9.SetColor(cyan);
  Detector_9.Construct(logicWorld);

  //--------------------------------------------------------------------------------------------------------------//

  // Return World!
  return physiWorld;
}

//--------------------------------------------------------------------------------------------------------------//

#include "G4RunManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

void DetectorConstruction::UpdateGeometry()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Reconstruct world with new parameters
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
//--------------------------------------------------------------------------------------------------------------//
