// ----------------------------------------------------------------------------
// nexus | AnalysisSteppingAction.cc
//
// This class allows the user to print the total number of photons detected by
// all kinds of photosensors at the end of the run.
// It also shows examples of information that can be accessed at the stepping
// level, so it is useful for debugging.
//
// The  NEXT Collaboration
// ----------------------------------------------------------------------------

#include "AnalysisSteppingAction.h"
#include "FactoryBase.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4ProcessManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4VPhysicalVolume.hh>
#include <G4SystemOfUnits.hh>

using namespace nexus;

REGISTER_CLASS(AnalysisSteppingAction, G4UserSteppingAction)

AnalysisSteppingAction::AnalysisSteppingAction(): G4UserSteppingAction(),
fWlsReemissionCount(0)
{
}



AnalysisSteppingAction::~AnalysisSteppingAction()
{
  G4double total_counts = 0;
  detectorCounts::iterator it = my_counts_.begin();
  while (it != my_counts_.end()) {
    G4cout << "Detector " << it->first << ": " << it->second << " counts" << G4endl;
    total_counts += it->second;
    it ++;
  }
  G4cout << "TOTAL COUNTS: " << total_counts << G4endl;
  /*
  //print the average track length 
  if (track_lengths_.size() > 0) {
    G4double sum_of_lengths = 0;
    for (G4double length : track_lengths_) {
      sum_of_lengths += length;
    }
    G4double average_length = sum_of_lengths / track_lengths_.size();
    G4cout << "Average Track Length: " << average_length / mm << " mm" << G4endl;
  }*/
  G4cout << "\n--------------------------------------------------" << G4endl;
  G4cout << "TOTAL WLS RE-EMITTED PHOTONS: " << fWlsReemissionCount << G4endl;
  G4cout << "--------------------------------------------------" << G4endl;
}



void AnalysisSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4ParticleDefinition* pdef = step->GetTrack()->GetDefinition();
  G4Track* track = step->GetTrack();

  const G4TrackVector* secondary = step->GetSecondary();
  if (secondary) {
      for (size_t i = 0; i < secondary->size(); ++i) {
          G4Track* sec_track = (*secondary)[i];
          const G4VProcess* creatorProcess = sec_track->GetCreatorProcess();
          //G4cout<<"creator process check" << creatorProcess << G4endl;

          // Check if the created particle is an optical photon AND its creator was WLS
          if (sec_track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() &&
              //creatorProcess && creatorProcess->GetProcessName() == "initStep" && 
              sec_track->GetParentID()>1)
          {
              fWlsReemissionCount++; // Increment the counter
          }
      }
  }
  /*
  if (track->GetTrackStatus() == fStopAndKill) {
    G4double final_track_length = track->GetTrackLength();
    //store this for further analysis
    track_lengths_.push_back(final_track_length);
    G4cout << "Track #" <<track->GetTrackID()
           << " terminated. Final length: " << final_track_length/ mm << " mm"
           << G4endl;
  }*/

  //Check whether the track is an optical photon
  if (pdef != G4OpticalPhoton::Definition()) return;

  /*
  // example of information one can access about optical photons

  G4Track* track = step->GetTrack();
  G4int pid = track->GetParentID();
  G4int tid = track->GetTrackID();
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  G4TouchableHandle touch1 = point1->GetTouchableHandle();
  G4TouchableHandle touch2 = point2->GetTouchableHandle();
  G4String vol1name = touch1->GetVolume()->GetName();
  G4String vol2name = touch2->GetVolume()->GetName();

  G4String proc_name = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4int copy_no = step->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(1);
  */

  // Retrieve the pointer to the optical boundary process.
  // We do this only once per run defining our local pointer as static.
  static G4OpBoundaryProcess* boundary = 0;

  if (!boundary) { // the pointer is not defined yet
    // Get the list of processes defined for the optical photon
    // and loop through it to find the optical boundary process.
    G4ProcessVector* pv = pdef->GetProcessManager()->GetProcessList();
    for (size_t i=0; i<pv->size(); i++) {
      if ((*pv)[i]->GetProcessName() == "OpBoundary") {
	boundary = (G4OpBoundaryProcess*) (*pv)[i];
	break;
      }
    }
  }

  if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
    if (boundary->GetStatus() == Detection ){
      G4String detector_name = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
      //G4cout << "##### Sensitive Volume: " << detector_name << G4endl;

      detectorCounts::iterator it = my_counts_.find(detector_name);
      if (it != my_counts_.end()) my_counts_[it->first] += 1;
      else my_counts_[detector_name] = 1;
    }
  }

  return;
}
