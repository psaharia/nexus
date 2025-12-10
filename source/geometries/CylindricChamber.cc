// ----------------------------------------------------------------------------
// nexus | CylindricChamber.cc
//
// General-purpose cylindric chamber.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "CylindricChamber.h"
#include "GenericPhotosensor.h"
#include "GenericWLSFiber.h"
#include "SensorSD.h"

#include "PmtR11410.h"
#include "NextNewKDB.h"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "UniformElectricDriftField.h"
#include "IonizationSD.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include "G4SystemOfUnits.hh"     
#include "G4UnitsTable.hh"  
#include <G4Tubs.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SDManager.hh>
#include <G4ApplicationState.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include "Randomize.hh"


namespace nexus {

  REGISTER_CLASS(CylindricChamber, GeometryBase)

  using namespace CLHEP;

  CylindricChamber::CylindricChamber():
  fVertexZPosition(5.5* cm),
  fFiberRadius(0.5 *mm),
  fMaxPhotonSteps_(2000),
  fUserLimits_(0), //Initialize pointer to null
    GeometryBase() //msg_(0)
  {
    msg_ = new G4GenericMessenger(this, "/source/geometries/CylindricChamber/",
                                  "Control commands for the cylindric chamber geometry");

    msg_->DeclareMethodWithUnit("set_vertex_z", "cm", &CylindricChamber::SetVertexZ,
                                "Set the Z position used by CylindricChamber::GenerateVertex.")
      .SetUnitCategory("Length") // Suggests appropriate units in help
      .SetDefaultValue("5.5")    // Default value in the specified unit (cm)
      .SetStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);
      // Available for states: G4State_PreInit (before run initialization)
      // and G4State_Idle (after a run, before a new one). This is important!

    msg_->DeclareMethodWithUnit("set_fiber_radius", "mm", &CylindricChamber::SetFiberRadius,
                                "Set the radius of the fiber's cross section for vertex generation.")
      .SetUnitCategory("Length") 
      .SetDefaultValue("0.5")    // Default value in the specified unit (mm)
      .SetStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);

    msg_->DeclareProperty("maxPhotonSteps", fMaxPhotonSteps_,
      "Maximum number of steps for optical photons inside this chamber.")
      .SetStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);
  }



  CylindricChamber::~CylindricChamber()
  {
    delete msg_;
  }



  void CylindricChamber::Construct()
  {
    // CHAMBER ///////////////////////////////////////////////////////

    /*const G4double chamber_diam   =  20. * cm;
    const G4double chamber_length = 100. * cm;
    const G4double chamber_thickn =   1. * cm;

    G4Tubs* chamber_solid =
      new G4Tubs("CHAMBER", 0., (chamber_diam/2. + chamber_thickn),
        (chamber_length/2. + chamber_thickn), 0., twopi);

    G4LogicalVolume* chamber_logic =
      new G4LogicalVolume(chamber_solid, materials::Steel(), "CHAMBER");

    this->SetLogicalVolume(chamber_logic);*/

    // GAS CYLINDER /////////////////////////////////////////////////

    const G4double chamber_diam   =  250. * cm;
    const G4double chamber_length = 500. * cm;

    G4Tubs* gas_solid =
      new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2., 0., twopi);

    G4Material* gxe = materials::GXe(10.*bar);
    gxe->SetMaterialPropertiesTable(opticalprops::Vacuum());

    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

    this->SetLogicalVolume(gas_logic);

    // G4UserLimits constructor parameters (in order):
    // max_step, max_track_length, max_time, max_nb_steps, min_ekin, min_range
    // We set max_nb_steps, and use DBL_MAX for others to signify "no limit"
    fUserLimits_ = new G4UserLimits();

    fUserLimits_->SetUserMaxTrackLength(20000000 * mm);

    //gas_logic->SetUserLimits(fUserLimits_);

    // GAS ///////////////////////////////////////////////////////////

    /*G4Tubs* gas_solid =
      new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2., 0., twopi);

    G4Material* gxe = materials::GXe(10.*bar);
    gxe->SetMaterialPropertiesTable(opticalprops::GXe(10.*bar, 303));

    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), gas_logic, "GAS",
		      chamber_logic, false, 0, true);*/

    /*
    // ACTIVE ////////////////////////////////////////////////////////

    const G4double active_diam   = chamber_diam;
    const G4double active_length = chamber_length/2.;

    G4Tubs* active_solid =
      new G4Tubs("ACTIVE", 0., active_diam/2., active_length/2., 0, twopi);

    G4LogicalVolume* active_logic =
      new G4LogicalVolume(active_solid, gxe, "ACTIVE");

    new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), active_logic, "ACTIVE",
		      gas_logic, false, 0, true);

    // Define this volume as an ionization sensitive detector
    IonizationSD* sensdet = new IonizationSD("/CYLINDRIC_CHAMBER/ACTIVE");
    active_logic->SetSensitiveDetector(sensdet);
    G4SDManager::GetSDMpointer()->AddNewDetector(sensdet);

    // Define an electric drift field for this volume
    UniformElectricDriftField* drift_field = new UniformElectricDriftField();
    drift_field->SetCathodePosition(-active_length/2.);
    drift_field->SetAnodePosition(active_length/2.);
    drift_field->SetDriftVelocity(1.*mm/microsecond);
    drift_field->SetTransverseDiffusion(1.*mm/sqrt(cm));
    drift_field->SetLongitudinalDiffusion(.5*mm/sqrt(cm));

    G4Region* drift_region = new G4Region("DRIFT_REGION");
    drift_region->SetUserInformation(drift_field);
    drift_region->AddRootLogicalVolume(active_logic);*/

    /*
    // EL GAP ////////////////////////////////////////////////////////

    const G4double elgap_diam   = active_diam;
    const G4double elgap_length = 1. * cm;

    G4Tubs* elgap_solid =
      new G4Tubs("EL_GAP", 0., elgap_diam/2., elgap_length/2., 0, twopi);

    G4LogicalVolume* elgap_logic =
      new G4LogicalVolume(elgap_solid, gxe, "EL_GAP");

    G4double pos_z = active_length/2. + elgap_length/2.;

    new G4PVPlacement(0, G4ThreeVector(0.,0.,pos_z), elgap_logic, "EL_GAP",
		      gas_logic, false, 0, true);

    // Define an EL field for this volume
    UniformElectricDriftField* el_field = new UniformElectricDriftField();
    el_field->SetCathodePosition(active_length/2.);
    el_field->SetAnodePosition(active_length/2. + elgap_length);
    el_field->SetDriftVelocity(5.*mm/microsecond);
    el_field->SetTransverseDiffusion(1.*mm/sqrt(cm));
    el_field->SetLongitudinalDiffusion(.5*mm/sqrt(cm));
    el_field->SetLightYield(1000./cm);

    G4Region* el_region = new G4Region("EL_REGION");
    el_region->SetUserInformation(el_field);
    el_region->AddRootLogicalVolume(elgap_logic);*/

    /*
    // PHOTOMULTIPLIER ///////////////////////////////////////////////

    PmtR11410 pmt_geom;
    pmt_geom.SetSensorDepth(0);
    pmt_geom.Construct();
    G4LogicalVolume* pmt_logic = pmt_geom.GetLogicalVolume();

    pos_z = -30. * cm;

    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_z), pmt_logic,
      "PMT", gas_logic, false, 0, true);*/


    // PHOTOSENSOR ///////////////////////////////////////////////////
    
    G4MaterialPropertiesTable* window_mat = opticalprops::Vacuum();
    //window_mat->SetMaterialPropertiesTable(opticalprops::Vacuum());

    G4MaterialPropertiesTable* photosensor_mpt = new G4MaterialPropertiesTable();
    G4double energy[]       = {0.2 * eV, 3.5 * eV, 3.6 * eV, 11.5 * eV};
    G4double reflectivity[] = {0.0     , 0.0     , 0.0     ,  0.0     };
    G4double efficiency[]   = {1.0     , 1.0     , 0.0     ,  0.0     };
    photosensor_mpt->AddProperty("REFLECTIVITY", energy, reflectivity, 4);
    photosensor_mpt->AddProperty("EFFICIENCY",   energy, efficiency,   4);

    G4MaterialPropertyVector* wind_rindex = 
      window_mat->GetProperty("RINDEX");
    photosensor = new GenericPhotosensor("photosensor", 2*mm); //originally 5cm 
    //photosensor.SetDimensions(5. * cm, 5. * cm, 0.5 * cm);
    //photosensor.SetSensitiveDepth(1);
    photosensor->SetWindowRefractiveIndex(wind_rindex);
    photosensor->SetOpticalProperties(photosensor_mpt);
    photosensor->SetSensorDepth(1);
    photosensor->SetVisibility(true);
    photosensor->Construct();

    G4LogicalVolume* photosensor_logic = photosensor->GetLogicalVolume();
    G4RotationMatrix* rotation = new G4RotationMatrix();
    //G4RotationMatrix* rotation = nullptr;
    rotation->rotateY(180 * deg);
    //rotation->rotateZ(180 * deg);
    //rotation->rotateX(0 * deg);

    // Place the photosensor inside the chamber
    //G4double pos_x = -30. * cm;

    new G4PVPlacement(rotation, G4ThreeVector(0 *mm, 0. *mm, 402.0 *mm), photosensor_logic,
                      photosensor_logic->GetName(), gas_logic, false, 0, true);

    // Sensitive detector
    SensorSD* pmtsd = new SensorSD("/PMT_R7378A/Pmt");
    pmtsd->SetDetectorVolumeDepth(2);
    pmtsd->SetTimeBinning(100.*nanosecond);
    G4SDManager::GetSDMpointer()->AddNewDetector(pmtsd);
    photosensor_logic->SetSensitiveDetector(pmtsd);
    

    //WLS FIBRE    ////////////////////////////////////////////////
    
    // Fiber dimensions
    const G4double fiber_diameter = 2. * fFiberRadius; //1 mm fFiberRadius is 0.5mm
    const G4double fiber_length   = 800 * mm; // 800 mm for BCF91 and 92 
    //const G4double fiber_length   = 3500 * mm; //for Y11 simulations 

    //G4String core_mat_ = "EJ280"; 
    //G4Material* core_mat = nullptr; 
    G4Material* core_mat = materials::PVT();
    core_mat->SetMaterialPropertiesTable(opticalprops::BCF91A());

    //G4Material* core_mat = "Y11"; 
    //G4Material* coating_material = materials::TPB();             

    GenericWLSFiber* wls_fiber = new GenericWLSFiber(
      "wls_fiber",        // Name
      true,               // Verbosity (for debugging)
      true,               // isRound
      fiber_diameter,     // Thickness (diameter for round fiber)
      fiber_length,       // Length
      false,               // Double cladding
      false,               // With coating
      core_mat,           // Core material
      0,   // Coating material
      true                // Visibility
      );

    
    wls_fiber->Construct();

    // Place the fiber inside the gas volume
    G4LogicalVolume* fiber_logic = wls_fiber->GetLogicalVolume();
    G4ThreeVector fiber_position(0.0, 0.0, 0.0 * mm);  //originally 0,0,5.5 for the 500nm simulations

    new G4PVPlacement(0, fiber_position, fiber_logic, "wls_fiber",
                  gas_logic, false, 0, true);

    fiber_logic->SetUserLimits(fUserLimits_);

    


    // SPECTROMETER FIBRE - PERPENDICULAR //////////////////////////////
    /*
    // Fiber dimensions
    const G4double spec_fiber_length = 1 * m; //1 meter og
    const G4double spec_fiber_diameter = 50 * micrometer; //50 micrometres og 

    G4Material* spec_fiber_mat = materials::FusedSilica();
    spec_fiber_mat->SetMaterialPropertiesTable(opticalprops::SpecFibre());

    GenericWLSFiber* spec_fiber = new GenericWLSFiber(
      "spec_fiber",     // Name
      true,               // Verbosity
      true,               // isRound
      spec_fiber_diameter,     // Thickness (diameter for round fiber)
      spec_fiber_length,// Length
      false,              // Double cladding
      false,              // With coating
      spec_fiber_mat,         // Core material
      0,                  // Coating material
      true                // Visibility
      );

    spec_fiber->Construct();

    // Place the perpendicular fiber inside the gas volume
    G4LogicalVolume* spec_fiber_logic = spec_fiber->GetLogicalVolume();
    G4ThreeVector spec_fiber_position(50.05 * cm, 0 * cm, 0 * cm);
    //G4RotationMatrix* rotation = new G4RotationMatrix();
    //rotation->rotateY(90 * deg);

    new G4PVPlacement(rotation, spec_fiber_position, spec_fiber_logic, "spec_fiber",
                      gas_logic, false, 0, true);
  }*/

    /*   
    // DICE BOARD ////////////////////////////////////////////////////

    NextNewKDB kdb_geom(5,5);
    kdb_geom.Construct();

    G4LogicalVolume* kdb_logic = kdb_geom.GetLogicalVolume();

    pos_z = active_length/2. + elgap_length + 5.0*mm;

    new G4PVPlacement(0, G4ThreeVector(0., 0., pos_z), kdb_logic,
      "KDB", gas_logic, false, 0, true);*/

  }


  
  G4ThreeVector nexus::CylindricChamber::GenerateVertex(const G4String& /*region*/) const
  {
    //return G4ThreeVector(0.,0.,5.5 *cm);
    const G4double safety_margin = 0.1 * mm; // Safety margin to avoid edge effects
    G4double effective_fiber_radius = fFiberRadius - safety_margin;

    if (effective_fiber_radius <= 0.) {
      G4ExceptionDescription ed;
      ed << "Effective fiber radius for vertex generation is non-positive ("
         << G4BestUnit(effective_fiber_radius, "Length")
         << "). This means the safety margin (" << G4BestUnit(safety_margin, "Length")
         << ") is too large for the fiber radius (" << G4BestUnit(fFiberRadius, "Length") << ").";
      G4Exception("CylindricChamber::GenerateVertex", "InvalidRadius",
                  FatalException, ed, "Adjust fFiberRadius or safety_margin.");
    }

    //G4double r_uniform = G4UniformRand(); // random number in [0,1)
    G4double r_uniform = 1 *mm;
    G4double theta_uniform = CLHEP::twopi * G4UniformRand(); // random angle between 0 and 2pi 

    //calculate actual radiu 'r' using sqrt for uniform area distribution
    G4double r = effective_fiber_radius * std::sqrt(r_uniform); 

    //convert polar coordinates to cartesian coordinates
    G4double x = r * std::cos(theta_uniform);
    G4double y = r * std::sin(theta_uniform);



    return G4ThreeVector(x, y, fVertexZPosition);
  }

  void CylindricChamber::SetVertexZ(G4double z)
  {
    // Set the Z position used by GenerateVertex
    fVertexZPosition = z;
    G4cout << "CylindricChamber: Vertex Z position set to " 
           << G4BestUnit(fVertexZPosition, "Length") << G4endl;
  }

  void CylindricChamber::SetFiberRadius(G4double radius)
  {
    // Set the radius of the fiber's cross section for vertex generation
    fFiberRadius = radius;
    G4cout << "CylindricChamber: Fiber radius set to " 
           << G4BestUnit(fFiberRadius, "Length") << G4endl;
  }

/*  void CylindricChamber::SetMaxPhotonSteps(G4int steps)
  {
    // Set the maximum number of steps for optical photons
    fMaxPhotonSteps_ = steps;
    G4cout << "CylindricChamber: Maximum photon steps set to " 
           << fMaxPhotonSteps_ << G4endl;
  }*/
  
}
// end namespace nexus
