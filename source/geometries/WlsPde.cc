// ----------------------------------------------------------------------------
// nexus | WLSpde.cc
//
// WLS fibre photon detection efficiency simulation.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "WlsPde.h"
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
#include "PmtR7378A.h"

namespace nexus {

  REGISTER_CLASS(WlsPde, GeometryBase)

  using namespace CLHEP;

  WlsPde::WlsPde():
  //fVertexXPosition(2.0 * mm),
    GeometryBase() //msg_(0)
  {
    msg_ = new G4GenericMessenger(this, "/source/geometries/WlsPde/",
                                  "Control commands for the WLS PDE geometry");

     //msg_->DeclareMethodWithUnit("set_vertex_x", "mm", &WlsPde::SetVertexX,
       //                         "Set the X position used by WlsPde::GenerateVertex.")
      //.SetUnitCategory("Length") // Suggests appropriate units in help
      //.SetDefaultValue("2.0")    // Default value in the specified unit (mm)
      //.SetStates(G4ApplicationState::G4State_PreInit, G4ApplicationState::G4State_Idle);
  }



  WlsPde::~WlsPde()
  {
    delete msg_;
  }



  void WlsPde::Construct()
  {
    // GAS VOLUME ///////////////////////////////////////////////////////

    const G4double chamber_diam   =  20. * cm;
    const G4double chamber_length = 50. * cm;

    G4Tubs* gas_solid =
      new G4Tubs("GAS", 0., chamber_diam/2., chamber_length/2., 0., twopi);

    G4Material* gxe = materials::GXe(10.*bar);
    gxe->SetMaterialPropertiesTable(opticalprops::Vacuum());

    G4LogicalVolume* gas_logic = new G4LogicalVolume(gas_solid, gxe, "GAS");

    this->SetLogicalVolume(gas_logic);


    // PHOTOSENSOR ///////////////////////////////////////////////////
  
    const G4double fiber_length   = 110 * mm; // 800 mm for BCF91 and 92 
    G4RotationMatrix sensor_rot;
    //sensor_rot.rotateY(pi);
    //G4ThreeVector sensor_pos = G4ThreeVector(0 *mm, 0. *mm, -402.0 *mm);
    G4ThreeVector sensor_pos = G4ThreeVector(0 *mm, 0. *mm, (-fiber_length/2. - 21.55) *mm);   
    PmtR7378A pmt;
    pmt.Construct();
    G4LogicalVolume* pmt_logic_ = pmt.GetLogicalVolume();
    new G4PVPlacement(G4Transform3D(sensor_rot, sensor_pos),
  		      pmt_logic_, "PMT", gas_logic, true, 0, true);


    //WLS FIBRE    ////////////////////////////////////////////////
    
    // Fiber dimensions
    const G4double FiberRadius = 0.5 * mm;
    const G4double fiber_diameter = 2. * FiberRadius; //1 mm FiberRadius is 0.5mm

    G4Material* core_mat = materials::PVT();
    core_mat->SetMaterialPropertiesTable(opticalprops::Y11());

    //G4Material* core_mat = "Y11"; 
    //G4Material* coating_material = materials::TPB();             

    GenericWLSFiber* wls_fiber = new GenericWLSFiber(
      "wls_fiber",        // Name
      true,               // Verbosity (for debugging)
      true,               // isRound
      fiber_diameter,     // Thickness (diameter for round fiber)
      fiber_length,       // Length
      true,               // Double cladding
      false,               // With coating
      core_mat,           // Core material
      0,   // Coating material
      true                // Visibility
      );

    
    wls_fiber->Construct();

    // Place the fiber inside the gas volume
    G4LogicalVolume* fiber_logic = wls_fiber->GetLogicalVolume();
  

    // Fiber 1
    G4ThreeVector fiber_position_1(0.0, 0.0, 0.0 * mm);
    new G4PVPlacement(0,                     // Rotation (null for no rotation)
                  fiber_position_1,      // Position
                  fiber_logic,           // The single logical volume
                  "wls_fiber",         // A unique name for this physical volume instance
                  gas_logic,             // Mother volume
                  false,                 // Not a replica
                  0,                     // Copy number
                  true);                 // Check for overlaps
    /*
    // Fiber 2
    G4ThreeVector fiber_position_2(0.0, 3 * mm, 0.0 * mm); 
    new G4PVPlacement(0,
                  fiber_position_2,
                  fiber_logic,
                  "wls_fiber",
                  gas_logic,
                  false,
                  1,                     // A unique copy number
                  true);

    // Fiber 3
    G4ThreeVector fiber_position_3(0.0, 6 * mm, 0 * mm); 
    new G4PVPlacement(0,
                  fiber_position_3,
                  fiber_logic,
                  "wls_fiber",
                  gas_logic,
                  false,
                  2,
                  true);

    // Fiber 4
    G4ThreeVector fiber_position_4(0.0, -3 * mm, 0 * mm); 
    new G4PVPlacement(0,
                  fiber_position_4,
                  fiber_logic,
                  "wls_fiber",
                  gas_logic,
                  false,
                  3,
                  true);
                  
    // Fiber 5
    G4ThreeVector fiber_position_5(0.0, -6 * mm, 0.0 * mm); 
    new G4PVPlacement(0,
                  fiber_position_5,
                  fiber_logic,
                  "wls_fiber",
                  gas_logic,
                  false,
                  4,
                  true);*/

  }


  
  G4ThreeVector nexus::WlsPde::GenerateVertex(const G4String& /*region*/) const
  {

    /*G4double r_max = 12.7 *mm; 

    // Generate the radius 'r' such that the points are uniform in area.
    // This requires r^2 to be uniformly distributed between 0 and r_max^2.
    // G4UniformRand() returns a random number between 0 and 1.
    G4double r_squared = r_max * r_max * G4UniformRand(); 
    G4double r = std::sqrt(r_squared);
    
    // Generate the angle 'theta' uniformly between 0 and 2pi
    G4double theta_uniform = CLHEP::twopi * G4UniformRand(); 

    // Calculate y and z coordinates using polar-to-Cartesian conversion
    // Note: The original code used cos for z and sin for y. We maintain that convention.
    G4double y = r * std::sin(theta_uniform);
    G4double z = r * std::cos(theta_uniform);

    G4cout << "Generated vertex at (x, y, z): (" 
           << G4BestUnit(fVertexXPosition, "Length") << ", "
           << G4BestUnit(y, "Length") << ", "
           << G4BestUnit(z, "Length") << ")" << G4endl;

    // The vertex position is (fVertexXPosition, y, z)*/
    //return G4ThreeVector(fVertexXPosition, y, z);
    //return G4ThreeVector(0.65 *mm, 0.0, 0.0);
    return G4ThreeVector(0.65 *mm, 0.0, 0.0);
  }

  /*void WlsPde::SetVertexX(G4double x)
  {
    // Set the X position used by GenerateVertex
    fVertexXPosition = x;
    G4cout << "WlsPde: Vertex X position set to " 
           << G4BestUnit(fVertexXPosition, "Length") << G4endl;
  }*/

}
// end namespace nexus
