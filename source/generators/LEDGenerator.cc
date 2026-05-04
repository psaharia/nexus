// ----------------------------------------------------------------------------
// nexus | LEDGenerator.cc
//
// This class is the primary generator of a number of optical photons with
// energy following the LED spectrum used for fibre measurements at IFIC. This 
// is copied and modified from the already existing framework
// in ScintillationGenerator.cc
//
// The NEXT Collaboration
// ----------------------------------------------------------------------------

#include "LEDGenerator.h"

#include "DetectorConstruction.h"
#include "GeometryBase.h"
#include "OpticalMaterialProperties.h"
#include "FactoryBase.h"

#include <G4GenericMessenger.hh>
#include <G4ParticleDefinition.hh>
#include <G4RunManager.hh>
#include <G4ParticleTable.hh>
#include <G4PrimaryVertex.hh>
#include <G4Event.hh>
#include <G4RandomDirection.hh>
#include <G4OpticalPhoton.hh>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "CLHEP/Units/SystemOfUnits.h"

using namespace nexus;
using namespace CLHEP;

REGISTER_CLASS(LEDGenerator, G4VPrimaryGenerator)


LEDGenerator::LEDGenerator() :
  G4VPrimaryGenerator(), msg_(0), geom_(0), nphotons_(1), led_spectrum_cdf_(0)
{
  msg_ = new G4GenericMessenger(this, "/Generator/LEDGenerator/",
    "Control commands of LED spectrum generator.");

  msg_->DeclareProperty("region", region_,
                        "Region of the geometry where the vertex will be generated.");

  msg_->DeclareProperty("nphotons", nphotons_, "Number of photons");

  msg_->DeclareProperty("spectrum_file_", spectrum_file_, "Path to LED CSV spectrum file");

  
  DetectorConstruction* detconst =
    (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  geom_ = detconst->GetGeometry();
}

void LEDGenerator::LoadSpectrumFromFile(G4String filename)
{
  std::ifstream file(filename);
  if (!file.is_open()) {
    // This will print the actual string the code received
    G4String msg = "Could not find file: '" + filename + "'. Please check the path and filename.";
    G4Exception("[LEDGenerator]", "LoadSpectrumFromFile()", FatalException, msg.c_str());
  }

  std::vector<std::pair<G4double, G4double>> energy_intensity_list;
  std::string line;
  bool data_reached = false;
  
  // Conversion constant: Energy(eV) = 1239.84193 / Wavelength(nm)
  const G4double hc = 1239.84193 * eV * nm;

  while (std::getline(file, line)) {
    // 1. Look for the [Data] tag to start parsing
    if (line.find("[Data]") != std::string::npos) {
      data_reached = true;
      continue;
    }
    if (!data_reached || line.empty() || line.find("[EndOfFile]") != std::string::npos) continue;

    // 2. Parse Thorlabs semicolon format: Wavelength;Intensity
    std::replace(line.begin(), line.end(), ';', ' ');
    std::stringstream ss(line);
    G4double wavelength_nm, intensity;

    if (ss >> wavelength_nm >> intensity) {
      if (wavelength_nm <= 0) continue; 
      
      // We only care about positive intensities to build the CDF
      if (intensity < 0) intensity = 0;

      G4double energy = hc / (wavelength_nm * nm);
      energy_intensity_list.push_back({energy, intensity});
    }
  }
  file.close();

  if (energy_intensity_list.empty()) {
    G4Exception("[LEDGenerator]", "LoadSpectrumFromFile()",
                FatalException, "No valid data points found in CSV!");
  }

  // 3. IMPORTANT: Sort by energy ascending (CSV is descending energy)
  std::sort(energy_intensity_list.begin(), energy_intensity_list.end());

  G4PhysicsOrderedFreeVector led_pdf;
  for (auto const& point : energy_intensity_list) {
    led_pdf.InsertValues(point.first, point.second);
  }

  // 4. Compute Cumulative Distribution
  led_spectrum_cdf_ = new G4PhysicsOrderedFreeVector();
  ComputeCumulativeDistribution(led_pdf, *led_spectrum_cdf_);
  sc_max_ = led_spectrum_cdf_->GetMaxValue();
  
  G4cout << "LED Generator: Successfully loaded " << energy_intensity_list.size() 
         << " sampling points from Thorlabs CSV." << G4endl;
}

LEDGenerator::~LEDGenerator()
{
  delete msg_;
  if (led_spectrum_cdf_) delete led_spectrum_cdf_; 
}

void LEDGenerator::GeneratePrimaryVertex(G4Event* event)
{
  if (!led_spectrum_cdf_) {
    if (spectrum_file_ == "") {
        G4Exception("[LEDGenerator]", "GeneratePrimaryVertex()", FatalException, "No spectrum file provided in macro!");
    }
    LoadSpectrumFromFile(spectrum_file_);
  }

  G4ParticleDefinition* particle_definition = G4OpticalPhoton::Definition();
  // Generate an initial position for the particle using the geometry and set time to 0.
  G4ThreeVector position = geom_->GenerateVertex(region_);
  G4double time = 0.;

  // Energy is sampled from integral (like it is done in G4Scintillation)

  G4PrimaryVertex* vertex = new G4PrimaryVertex(position, time);

  for (G4int i = 0; i < nphotons_; i++)
  {
      // DIRECTION: Set to (1, 0, 0) to go perpendicularly into the fiber 
      // (assuming vertex is at -X as per your previous code)
      G4ThreeVector _momentum_direction(-1.0, 0.0, 0.0); 

      // ENERGY: Sample from our hard-coded LED spectrum
      G4double sc_value = G4UniformRand() * sc_max_;
      G4double pmod = led_spectrum_cdf_->GetEnergy(sc_value); // This returns Energy
      
      G4double px = pmod * _momentum_direction.x();
      G4double py = pmod * _momentum_direction.y();
      G4double pz = pmod * _momentum_direction.z();

      G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition, px, py, pz);
      
      // Randomize polarization (Standard for LED/unpolarized sources)
      G4ThreeVector polarization = G4RandomDirection();
      particle->SetPolarization(polarization);

      vertex->SetPrimary(particle);
  }
  event->AddPrimaryVertex(vertex);
  
}

void LEDGenerator::ComputeCumulativeDistribution(
  const G4PhysicsOrderedFreeVector& pdf, G4PhysicsOrderedFreeVector& cdf)
{
  G4double sum = 0.;
  cdf.InsertValues(pdf.Energy(0), sum);

  for (unsigned int i=1; i<pdf.GetVectorLength(); ++i) {
    G4double area =
      0.5 * (pdf.Energy(i) - pdf.Energy(i-1)) * (pdf[i] + pdf[i-1]);
    sum = sum + area;
    cdf.InsertValues(pdf.Energy(i), sum);
  }
}
