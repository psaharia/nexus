// ----------------------------------------------------------------------------
// nexus | LEDGenerator.h
//
// This class is the primary generator of a number of optical photons with
// energy following the LED spectrum used for fibre measurements at IFIC. This 
// is copied and modified from the already existing framework
// in ScintillationGenerator.cc
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef LED_GENERATOR_H
#define LED_GENERATOR_H

#include <G4VPrimaryGenerator.hh>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4PhysicsOrderedFreeVector.hh>
#include <vector>

class G4GenericMessenger;
class G4Event;

namespace nexus {

  class GeometryBase;

  class LEDGenerator: public G4VPrimaryGenerator
  {
  public:
    /// Constructor
    LEDGenerator();
    /// Destructor
    ~LEDGenerator();

    /// This method is invoked at the beginning of the event. It sets
    /// a primary vertex (that is, a particle in a given position and time)
    /// in the event.
    void GeneratePrimaryVertex(G4Event*);

  private:

    void ComputeCumulativeDistribution(const G4PhysicsOrderedFreeVector&,
                                       G4PhysicsOrderedFreeVector&);
                                
    void LoadSpectrumFromFile(G4String filename);

    G4GenericMessenger* msg_;
    G4Navigator* geom_navigator_; ///< Geometry Navigator
    G4PhysicsOrderedFreeVector* led_spectrum_cdf_;
    G4double sc_max_;
    const GeometryBase* geom_; ///< Pointer to the detector geometry
    G4String spectrum_file_;

    G4String region_;
    G4int    nphotons_;

  };

} // end namespace nexus

#endif // __LEDGENERATOR__
