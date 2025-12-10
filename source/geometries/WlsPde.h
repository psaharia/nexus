// ----------------------------------------------------------------------------
// nexus | WLSpde.h
//
// WLS fibre photon detection efficiency simulation.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#ifndef WLS_PDE_H
#define WLS_PDE_H

#include "GeometryBase.h"
#include "G4ThreeVector.hh"
#include "G4UserLimits.hh"
class G4GenericMessenger;


namespace nexus {

    class GenericPhotosensor;
    class GenericWLSFiber;
    class WlsPde: public GeometryBase
    
  {
  public:
    /// Constructor
    WlsPde();
    /// Destructor
    ~WlsPde();

    /// Return vertex within region <region> of the chamber
    virtual G4ThreeVector GenerateVertex(const G4String& region) const;

    virtual void Construct();
    //void SetVertexX(G4double x);

  private:
    /// Messenger for the definition of control commands
    G4GenericMessenger* msg_;
    GenericPhotosensor* photosensor;
    //G4double fVertexXPosition;
    GenericWLSFiber* wls_fiber;
  };

} // end namespace nexus

#endif
