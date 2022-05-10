#include "ABPrimaryGeneratorAction.hh"
#include <stdio.h>
#include <math.h>
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "time.h"
#include "G4GenericMessenger.hh"
#include "G4PhysicalVolumeStore.hh"


double theta_source=0.5;
int nPrimeries=1;
int mode=1;

ABPrimaryGeneratorAction::ABPrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction(), fParticleGun(0){

  G4int n_particle=1;
  fParticleGun=new G4ParticleGun(n_particle);

  fMessenger=new G4GenericMessenger(this,"/ABALONE_Event/generator/","Primary generator controls");

  G4GenericMessenger::Command& modeCMD = fMessenger->DeclareProperty("mode",mode);
  modeCMD.SetParameterName("mode",true);
  modeCMD.SetRange("1:pe emission; 2:LYSO decay");
  modeCMD.SetDefaultValue("1");

  G4GenericMessenger::Command& angleCMD=fMessenger->DeclareProperty("angle",theta_source,"Angle of Event");
  angleCMD.SetParameterName("angle",true);
  angleCMD.SetRange("angle>0");
  angleCMD.SetDefaultValue("0.5");
  
  G4GenericMessenger::Command& particleCMD = fMessenger->DeclareProperty("number",nPrimeries);
  particleCMD.SetParameterName("number",true);
  particleCMD.SetRange("number<999");
  particleCMD.SetDefaultValue("1");

    // default particle kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle=particleTable->FindParticle(particleName="e-");
  fParticleGun->SetParticleDefinition(particle);
}


ABPrimaryGeneratorAction::~ABPrimaryGeneratorAction(){
  delete fParticleGun;
  delete fMessenger;
}



void ABPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  G4double r0=4.995*cm;

  if (mode==1){

  double phi_source=0.0; //2*pi*G4UniformRand()<->45/360
  double theta_source_rad=theta_source/360.0*2*3.1415926;
  G4double height=1.844*cm;
  G4double x_particle=(r0*sin(theta_source_rad)*cos(phi_source)), y_particle=(r0*sin(theta_source_rad)*sin(phi_source)), z_particle=r0*abs(cos(theta_source_rad))+height; //abs
  G4ThreeVector particle_pos=G4ThreeVector(x_particle, y_particle, z_particle);

  G4double x=x_particle, y=y_particle, z=z_particle;
  G4double x_dir=-1.*x/sqrt((std::pow(x, 2)+(std::pow(y, 2)+(std::pow(z, 2)))));
  G4double y_dir=-1.*y/sqrt((std::pow(x, 2)+(std::pow(y, 2)+(std::pow(z, 2)))));
  G4double z_dir=-1.*z/sqrt((std::pow(x, 2)+(std::pow(y, 2)+(std::pow(z, 2)))));
 
  fParticleGun->SetParticleEnergy(0.25*eV); // eqv to 10eV
  fParticleGun->SetParticlePosition(particle_pos); //particle_pos <-> temp_particle_pos

  G4double x_dir_prime, y_dir_prime, z_dir_prime;

  double x_dir_test , y_dir_test , z_dir_test; //for testing the dot product
  double x_dir_prime_test, y_dir_prime_test, z_dir_prime_test;


  double theta_rot=0., phi_rot=0., psi_rot=0.; //rotation angle around the axes
  double dot_product=-1.;

  if (nPrimeries!=0){
	  for (int j=0; j<nPrimeries;j++){
		x_dir_test=x_dir, y_dir_test=y_dir, z_dir_test=z_dir;
		while (dot_product < 0.5){ //around 84.26 degree angle

		  theta_rot=3.141515926*G4UniformRand(); //(3.141515926)*G4UniformRand();  //theta_source)/
		  phi_rot=(2. * 3.141515926)*G4UniformRand(); //(2 * 3.141515926) * G4UniformRand(); 
		  psi_rot=(2. * 3.141515926)*G4UniformRand(); //(2 * 3.141515926) * G4UniformRand();

		  x_dir_prime=(cos(psi_rot) * cos(phi_rot) - cos(theta_rot) * sin(phi_rot) * sin(psi_rot)) * x_dir + (cos(psi_rot) * sin(phi_rot) + cos(theta_rot) * cos(phi_rot) * sin(psi_rot)) * y_dir + (sin(psi_rot) * sin(theta_rot)) * z_dir;
		  y_dir_prime=(-sin(psi_rot) * cos(phi_rot) - cos(theta_rot) * sin(phi_rot) * cos(psi_rot)) * x_dir + (-sin(psi_rot) * sin(phi_rot) + cos(theta_rot) * cos(phi_rot) * cos(psi_rot)) * y_dir + (cos(psi_rot) * sin(theta_rot)) * z_dir;
		  z_dir_prime=sin(theta_rot) * sin(phi_rot) * x_dir - sin(theta_rot) * cos(phi_rot) * y_dir + cos(theta_rot) * z_dir;

		  x_dir_prime_test=x_dir_prime; y_dir_prime_test=y_dir_prime; z_dir_prime_test=z_dir_prime;

		  dot_product=x_dir_test*x_dir_prime_test+y_dir_test*y_dir_prime_test+z_dir_test*z_dir_prime_test; // prevent electrons going to the direction of photocathode
		}
		dot_product=-1.0;
    
		fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x_dir_prime, y_dir_prime, z_dir_prime)); // x_dir_prime, y_dir_prime, z_dir_prime<-> 0.,0.,-1 <-> x_dir, y_dir, z_dir
		fParticleGun->GeneratePrimaryVertex(anEvent);
    }
	}
  }

  if (mode==2){

    G4double xlo=-0.5*cm; G4double ylo=-0.5*cm; G4double zlo=-0.15004*cm;
    G4double xhi=0.5*cm; G4double yhi=0.5*cm; G4double zhi=-0.00004*cm;
    const G4String name="Scintillator";
    G4bool verbose=true;
    G4PhysicalVolumeStore * volumestore;
    G4VPhysicalVolume * pvol=volumestore->GetVolume(name,verbose);
    G4LogicalVolume * lvol=pvol->GetLogicalVolume();
    G4VSolid * solid=lvol->GetSolid();
    G4ThreeVector point;
    G4int maxtries=10000, itry=1;
    do{
      point.set(xlo+G4UniformRand()*(xhi-xlo), ylo+G4UniformRand()*(yhi-ylo), zlo+G4UniformRand()*(zhi-zlo));
    } while (!solid->Inside(point) && ++itry<maxtries);
    
    if (itry == maxtries)
      G4cerr << "Unable to find a point inside your volume!" << G4endl;

    fParticleGun->SetParticlePosition(point);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

