// ----------------------------------------------------------------------------
// nexus | DefaultTrackingAction.cc
//
// This class is the default tracking action of the NEXT simulation.
// It stores in memory the trajectories of all particles, except optical photons
// and ionization electrons, with the relevant tracking information that will be
// saved to the output file.
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "DefaultTrackingAction.h"

#include "Trajectory.h"
#include "TrajectoryMap.h"
#include "IonizationElectron.h"
#include "FactoryBase.h"

#include <G4Track.hh>
#include <G4TrackingManager.hh>
#include <G4Trajectory.hh>
#include <G4ParticleDefinition.hh>
#include <G4OpticalPhoton.hh>
#include <fstream>
#include <G4SystemOfUnits.hh>

using namespace nexus;

REGISTER_CLASS(DefaultTrackingAction, G4UserTrackingAction)

//write the output to a text file
static std::ofstream outputfile("/Users/psaharia/Pokhee-IFIC/Installs/nexus/outputs/PDE_comp/eff_conv_PDE_thickcheck.txt");

DefaultTrackingAction::DefaultTrackingAction() : G4UserTrackingAction()
{
  if (outputfile.is_open()){
    outputfile << "WLS conversion" << std::endl;
  }

}

DefaultTrackingAction::~DefaultTrackingAction()
{
  if (outputfile.is_open()){
    outputfile.close();
  }
  
}

void DefaultTrackingAction::PreUserTrackingAction(const G4Track *track)
{
  // Do nothing if the track is an optical photon or an ionization electron
  if (//track->GetDefinition() == G4OpticalPhoton::Definition() ||
      track->GetDefinition() == IonizationElectron::Definition())
  {
    fpTrackingManager->SetStoreTrajectory(false);
    return;
  }

  // Create a new trajectory associated to the track.
  // N.B. If the processesing of a track is interrupted to be resumed
  // later on (to process, for instance, its secondaries) more than
  // one trajectory associated to the track will be created, but
  // the event manager will merge them at some point.
  G4VTrajectory *trj = new Trajectory(track);

  // Set the trajectory in the tracking manager
  fpTrackingManager->SetStoreTrajectory(true);
  fpTrackingManager->SetTrajectory(trj);
}

void DefaultTrackingAction::PostUserTrackingAction(const G4Track *track)
{
  // Do nothing if the track is an optical photon or an ionization electron
  //if (//track->GetDefinition() == G4OpticalPhoton::Definition() ||
    //  track->GetDefinition() == IonizationElectron::Definition())
    //return;
  if (track->GetDefinition() != G4OpticalPhoton::Definition())
    return;

  //Get the process that ended the track
  G4String proc_name = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //check for WLS absorption
  if (proc_name == "OpWLS"){
    //write a specific marker and info for a WLS conversion event
    outputfile << track->GetTrackID() << std::endl;
  }
  return;

  Trajectory *trj = (Trajectory *)TrajectoryMap::Get(track->GetTrackID());

  // Do nothing if the track has no associated trajectory in the map
  if (!trj) return;

  // Record final time and position of the track
  trj->SetFinalPosition(track->GetPosition());
  trj->SetFinalTime(track->GetGlobalTime());
  trj->SetTrackLength(track->GetTrackLength());
  trj->SetFinalVolume(track->GetVolume()->GetName());
  trj->SetFinalMomentum(track->GetMomentum());

  // Record last process of the track
  //G4String proc_name = track->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //trj->SetFinalProcess(proc_name);

  //G4String next_volume = track->GetStep()->GetPostStepPoint()->GetPhysicalVolume()->GetName();

  //if (proc_name == "Transportation" && next_volume == "photosensor_SENSAREA") {
  //outputfile << track->GetTrackLength() / mm << std::endl;
  //}

  
}

