//
// ********************************************************************
// * License and Disclaimer    *
// **
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.   *
// **
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.  *
// **
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.   *
// ********************************************************************
//
//
/// \file ABDetectorConstruction.cc
/// \brief Implementation of the ABDetectorConstruction class

#include "ABDetectorConstruction.hh"
//#include "ABCalorimeterSD.hh"
#include <math.h>
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PropagatorInField.hh"
#include "G4UserLimits.hh"

#include "G4NistManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PhysicalConstants.hh"

#include "ABElectricField.hh" //Magnetic <-> Electric
#include "G4SDManager.hh"
#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
//#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
//#include "G4DormandPrince745.hh"

#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"


#include "G4Sipm.hh"
#include "MaterialFactory.hh"
#include "model/G4SipmModelFactory.hh"
#include "housing/G4SipmHousing.hh"
#include "housing/impl/HamamatsuCeramicHousing.hh"
#include "housing/impl/HamamatsuSmdHousing.hh"


#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <stdlib.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal ABElectricField* ABDetectorConstruction::fElectricField = 0;
G4ThreadLocal G4FieldManager* ABDetectorConstruction::fFieldMgr = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectricField* fEMfield;
G4EqMagElectricField* fEquation;
G4MagIntegratorStepper* fStepper;
G4FieldManager* fFieldMgr;
G4double fMinStep;
G4ChordFinder* fChordFinder;

ABDetectorConstruction::ABDetectorConstruction(std::string sipmModelName, std::string housingName) 
: G4VUserDetectorConstruction(),
  fMessenger(nullptr),
  fElectricLogical(nullptr),
  fScoringVolume(0),
  fNofLayers(-1)
{ 
    G4SipmModel* model = createSipmModel(sipmModelName);
    housing = createHousing(housingName, new G4Sipm(model));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABDetectorConstruction::~ABDetectorConstruction()
{ delete housing; }


G4SipmModel* ABDetectorConstruction::createSipmModel(std::string name) const {
	if (name == "generic") {
		return G4SipmModelFactory::getInstance()->createGenericSipmModel();
	}
	if (name == "HamamatsuS1036211100") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS1036211100();
	}
	if (name == "HamamatsuS1036233100") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS1036233100();
	}
	if (name == "HamamatsuS10985100") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS10985100();
	}
	if (name == "HamamatsuS12651050") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS12651050();
	}
	if (name == "HamamatsuS1036233050") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS1036233050();
	}
	if (name == "HamamatsuS12573100C") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS12573100C();
	}
	if (name == "HamamatsuS12573100X") {
		return G4SipmModelFactory::getInstance()->createHamamatsuS12573100X();
	}
	return G4SipmModelFactory::getInstance()->createConfigFileModel(name);
}

G4SipmHousing* ABDetectorConstruction::createHousing(std::string name, G4Sipm* sipm) const {
	if (name == "ceramic") {
		return new HamamatsuCeramicHousing(sipm);
	}
	if (name == "smd") {
		return new HamamatsuSmdHousing(sipm);
	}
	if (name == "default") {
		return new G4SipmHousing(sipm);
	}
	std::cout << "G4SipmDetectorConstruction::createHousingForName(name = \"" << name
			<< "\"): housing type does not exist." << std::endl;
	throw 1;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ABDetectorConstruction::Construct()
{  
	return DefineVolumes();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4VPhysicalVolume* ABDetectorConstruction::DefineVolumes() 
{
	fNofLayers = 1;
	G4NistManager* nist = G4NistManager::Instance();
	G4double env_sizeXY = 14 * cm, env_sizeZ = 14 * cm;
	
	G4bool checkOverlaps = true;

	////////////////////////////////
	////// Creating elements //////
	//////////////////////////////

	G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galatic");
	G4Material* alloy = nist->FindOrBuildMaterial("G4_Au");
	G4Material* In = nist->FindOrBuildMaterial("G4_In");
	G4Material* Cr = nist->FindOrBuildMaterial("G4_Cr");
	G4Material* Au = nist->FindOrBuildMaterial("G4_Au");
	G4Material* Ti = nist->FindOrBuildMaterial("G4_Ti");
	G4Material* quartz = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");

	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4double density;

	G4Material* fSilicon = new G4Material("Silicon", z = 14., a = 28.08 * g / mole,
		density = universe_mean_density); //2.3290 *g/cm3

	G4Material* fVacuum = new G4Material("Vacuum", z = 1., a = 1.01 * g / mole,
		density = universe_mean_density, kStateGas, 0.1 * kelvin,
		1.e-19 * pascal);
	
	G4Material* fCr = new G4Material("Chromium", z = 24., a = 51.996 * g / mole,
		density = 7.19*g/cm3);

	G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV,7.14 * eV };
	const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double);
	G4double vacuum_RIND[] = { 1.,1.,1. };
	G4MaterialPropertiesTable * vacuum_mt = new G4MaterialPropertiesTable();
	vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum);

	fVacuum->SetMaterialPropertiesTable(vacuum_mt);
	fCr->SetMaterialPropertiesTable(vacuum_mt);

	G4double silicon_Energy[] = { 2.0 * eV, 7.0 * eV,7.14 * eV };
	const G4int siliconnum = sizeof(silicon_Energy) / sizeof(G4double);
	G4double silicon_RIND[] = { 1.,1.,1. };
	G4MaterialPropertiesTable* silicon_mt = new G4MaterialPropertiesTable();
	silicon_mt->AddProperty("RINDEX", silicon_Energy, silicon_RIND, siliconnum);

	fSilicon->SetMaterialPropertiesTable(silicon_mt);

	//Scintillator material
	G4NistManager* man = G4NistManager::Instance();
	G4Element* Silicon = man->FindOrBuildElement("Si");
	G4Element* Lutetium = man->FindOrBuildElement("Lu");
	G4Element* Yttrium = man->FindOrBuildElement("Y");
	G4Element* Oxygen = man->FindOrBuildElement("O");
	G4Element* Cerium = man->FindOrBuildElement("Ce");
	G4Material* LYSO = new G4Material("LYSO", 7.1 * g / cm3, 5);
	LYSO->AddElement(Lutetium, 71.43 * perCent);
	LYSO->AddElement(Yttrium, 4.03 * perCent);
	LYSO->AddElement(Silicon, 6.37 * perCent);
	LYSO->AddElement(Oxygen, 18.14 * perCent);
	LYSO->AddElement(Cerium, 0.02 * perCent);
	LYSO->GetIonisation()->SetBirksConstant((0.0076 * g / (MeV * cm2)) / (7.4 * g / cm3));
	const G4int nEntries = 11;
	G4double PhotonEnergy[nEntries] = { 2.478 * eV, 2.53 * eV, 2.58 * eV, 2.636 * eV, 2.69 * eV, 2.75 * eV, 2.816 * eV, 2.88 * eV, 2.95 * eV, 3.022 * eV, 3.097 * eV };
	G4double RILyso[nEntries] = { 1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82,1.82 };// Is isotropic
	G4double ABLyso[nEntries] = { 100. * cm, 100. * cm, 100. * cm, 100. * cm, 100. * cm, 100. * cm, 100. * cm, 80. * cm, 50. * cm, 30. * cm, 2. * cm }; // To be updated
	G4double SFLyso[nEntries] = { 0., 0., 0., 0., 0., 0., 0.,0.1, 0.8, 0.1, 0. }; // To be updated
	G4double ReLyso[nEntries] = { 0.95, 0.95,0.95, 0.95, 0.95, 0.95, 0.95,0.95, 0.95, 0.95, 0.95 }; // To be updated
	G4MaterialPropertiesTable * MPTLyso = new G4MaterialPropertiesTable();
	MPTLyso->AddProperty("RINDEX", PhotonEnergy, RILyso, nEntries);
	MPTLyso->AddProperty("ABSLENGTH", PhotonEnergy, ABLyso, nEntries);
	MPTLyso->AddProperty("FASTCOMPONENT", PhotonEnergy, SFLyso, nEntries);
	MPTLyso->AddProperty("SLOWCOMPONENT", PhotonEnergy, SFLyso, nEntries);
	MPTLyso->AddConstProperty("SCINTILLATIONYIELD", 30000. / MeV);
	MPTLyso->AddConstProperty("RESOLUTIONSCALE", 1.04);//R = resolutionscale*sqrt(yield)
	MPTLyso->AddConstProperty("FASTTIMECONSTANT", 40. * ns);
	MPTLyso->AddConstProperty("SLOWTIMECONSTANT", 44. * ns);
	MPTLyso->AddConstProperty("YIELDRATIO", 1.0);
	LYSO->SetMaterialPropertiesTable(MPTLyso);
	
	/////////////////////////////////////
	////// 2. Creating geometries ///////
	////////////////////////////////////


	// 2.1 parameters
	G4double world_sizeXY =  env_sizeXY;
	G4double world_sizeZ =  env_sizeZ;

	// Hemisphere geometry description 
	G4double inner_radius = 5.0 * cm, outer_radius = 5.50164 * cm; 
	G4double phi_i = 0. * deg, phi_f = 180. * deg;
	G4double theta_i = 0. * deg, theta_f = 180. * deg;

	// sealing + photocathode
	G4double Glass_shell_diameter = 10.0 * cm;
	G4double sealing_Au_thickness = 50.0 * nm;
	G4double sealing_In_thickness = 0.005 * mm;
	G4double sealing_Cr_thickness = 30.0 * nm;
	G4double photo_cathode_thickness = 1.0 * nm; 
	G4double inner_radius_Au_sealing = (Glass_shell_diameter - (sealing_Au_thickness)) / 2.0;
	G4double outer_radius_Au_sealing = 5.0 * cm;
	G4double inner_radius_In_sealing = (Glass_shell_diameter - sealing_Au_thickness - sealing_In_thickness) / 2.0;
	G4double outer_radius_In_sealing = inner_radius_Au_sealing;
	G4double inner_radius_Cr_sealing = (Glass_shell_diameter - sealing_Au_thickness - sealing_In_thickness - sealing_Cr_thickness) / 2.0;
	G4double outer_radius_Cr_sealing = inner_radius_In_sealing;
	G4double inner_radius_photocathode = (Glass_shell_diameter - sealing_Au_thickness - sealing_In_thickness - sealing_Cr_thickness - photo_cathode_thickness) / 2.0;
	G4double outer_radius_photocathode = inner_radius_Cr_sealing;
	
	// Base disk construction
	G4double total_height = 1.844 * cm;
	G4double diameter_disk = 11.00328 * cm;
	
	//Base Disk (Upper part)
	G4double shell_width = 11.0328 - 2 * 4.0259 * cm;
	G4double inner_rad_up_1 = 2.01295 * cm, outter_rad_up_1 = 5.50164 * cm;  //lower part
	G4double inner_rad_up_2 = 4.0259 * cm, outter_rad_up_2 = 5.50164 * cm; // upper part
	G4double disk_height_up = 0.5 * 0.53082 * cm; //half height actually 0.286531343 = 0.5*(0.5* diameter_disk-shell_width-inner_rad_up_1)*tan(tilted_angle)								
	
												  // Base Disk (lower part 1)
	G4double inner_rad_low_1 = 0.2032 * cm, outter_rad_low_1 = 5.50164 * cm;  //lower part
	G4double inner_rad_low_2 = 2.01295 * cm, outter_rad_low_2 = 5.50164 * cm; // upper part
	G4double disk_height_low = 0.5 * 1.2827 * cm; //half height actually
	
	// Base Disk (lower part 2) most bottom thing
	G4double inner_rad_low2_1 = 0.2032 * cm, outter_rad_low2_1 = 5.50164 * cm;  //lower part
	G4double inner_rad_low2_2 = 0.2032 * cm, outter_rad_low2_2 = 5.50164 * cm; // upper part
	G4double disk_height_low2 = (0.5 * 0.03048) * cm; //half height actually height = 0.1cm 
	G4double scint_sealing_Au = 0.5 * 5. * nm;
	G4double scint_sealing_Cr = 0.5 * 35. * nm;
	G4double scintillator_thickness = 0.15 / 2. * cm; 
	G4double scintillator_dimension = 1.0 * cm;
	G4double scintillator_coating_thickness = 1.0 *cm;

	G4double sipm_thickness = 0.55 / 2. * mm; 

	G4double Ti_coating_thickness = 30. * micrometer;

	// Position of volumes
	G4ThreeVector pos1 = G4ThreeVector(0, 0, 1.844 * cm);
	G4ThreeVector pos_basedisk_lower2 = G4ThreeVector(0., 0., disk_height_low2);
	G4ThreeVector pos_basedisk_lower1 = G4ThreeVector(0., 0., disk_height_low2 * 2. + disk_height_low);
	G4ThreeVector pos_basedisk_upper = G4ThreeVector(0., 0., disk_height_low2 * 2. + disk_height_low * 2. + disk_height_up);
	G4ThreeVector pos_scintillator_Au = G4ThreeVector(0, 0, -1 * scint_sealing_Au);
	G4ThreeVector pos_scintillator_Cr = G4ThreeVector(0, 0, -1 * (2. * (scint_sealing_Au)+scint_sealing_Cr));
	G4ThreeVector pos_scintillator = G4ThreeVector(0, 0, -1 * (2. * (scint_sealing_Au + scint_sealing_Cr) + scintillator_thickness)); 
	G4ThreeVector pos_sipm = G4ThreeVector(0, 0, -1 * (2. * (scint_sealing_Au + scint_sealing_Cr + scintillator_thickness) + housing->getDz()/ 2.)); //
	
	// Rotation
	G4RotationMatrix* yRot = new G4RotationMatrix;
	yRot->rotateX(270. * deg);

	//Volume creation
	G4Box*solidWorld = new G4Box("World",0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);
	
	G4Sphere* glass_shell = new G4Sphere("Glass shell", inner_radius, outer_radius, phi_i, phi_f, theta_i, theta_f);
	G4Sphere* Shape_Au_sealing = new G4Sphere("Photocathode", inner_radius_Au_sealing, outer_radius_Au_sealing, phi_i, phi_f, theta_i, theta_f);
	G4Sphere* Shape_In_sealing = new G4Sphere("Photocathode", inner_radius_In_sealing, outer_radius_In_sealing, phi_i, phi_f, theta_i, theta_f);
	G4Sphere* Shape_Cr_sealing = new G4Sphere("Photocathode", inner_radius_Cr_sealing, outer_radius_Cr_sealing, phi_i, phi_f, theta_i, theta_f);
	G4Sphere* Shape_photocathode = new G4Sphere("Photocathode", inner_radius_photocathode, outer_radius_photocathode, phi_i, phi_f, theta_i, theta_f);
	
	G4Cons* Base_lower = new G4Cons("Base", inner_rad_low_1, outter_rad_low_1, inner_rad_low_2, outter_rad_low_2, disk_height_low, 0.0 * deg, 360. * deg);
	G4Cons* Base_lower_coating = new G4Cons("Base", inner_rad_low_1 - Ti_coating_thickness, inner_rad_low_1, inner_rad_low_2 - Ti_coating_thickness, inner_rad_low_2, disk_height_low, 0.0 * deg, 360. * deg);
	G4Cons* Base_lower2 = new G4Cons("Base", inner_rad_low2_1, outter_rad_low2_1, inner_rad_low2_2, outter_rad_low2_2, disk_height_low2, 0.0 * deg, 360. * deg);
	G4Cons* Base_upper = new G4Cons("Base", inner_rad_up_1, outter_rad_up_1, inner_rad_up_2, outter_rad_up_2, disk_height_up, 0.0 * deg, 360. * deg);
	
	G4Box* scint_sealing_Au_shape = new G4Box("scintillator", 1. * cm, 1. * cm, scint_sealing_Au); 
	G4Box* scint_sealing_Cr_shape = new G4Box("scintillator", 1. * cm, 1. * cm, scint_sealing_Cr); 
	G4Box* scintillator = new G4Box("scintillator", scintillator_dimension, scintillator_dimension, scintillator_thickness);

	G4Box* sipm = new G4Box("SiPM", 0.3 *cm, 0.3 * cm, sipm_thickness); 

	// Logical Volume creation
	fElectricLogical = new G4LogicalVolume(solidWorld,fVacuum,"World");
	G4LogicalVolume* glass_shell_logical = new G4LogicalVolume(glass_shell, quartz, "GlassShellLV");
	G4LogicalVolume* logicShape_Au_sealing = new G4LogicalVolume(Shape_Au_sealing, Au, "AU_sealingLV");
	G4LogicalVolume* logicShape_In_sealing = new G4LogicalVolume(Shape_In_sealing, In, "In_sealingLV");
	G4LogicalVolume* logicShape_Cr_sealing = new G4LogicalVolume(Shape_Cr_sealing, Cr, "Cr_sealingLV");
	G4LogicalVolume* logicShape_photocathode = new G4LogicalVolume(Shape_photocathode, alloy, "PhotocathodeLV");
	G4LogicalVolume* Base_upperLV = new G4LogicalVolume(Base_upper, quartz, "Base_logical");
	G4LogicalVolume* Base_lower1LV = new G4LogicalVolume(Base_lower, quartz, "Base_logical");
	G4LogicalVolume* Base_lower1_coating_LV = new G4LogicalVolume(Base_lower_coating, quartz, "Base_logical");
	G4LogicalVolume* Base_lower2LV = new G4LogicalVolume(Base_lower2, quartz, "Base_logical");
	G4LogicalVolume* logicShape_Sealing_Au = new G4LogicalVolume(scint_sealing_Au_shape, Au, "Scint_sealing_Au_LV");
	G4LogicalVolume* logicShape_Sealing_Cr = new G4LogicalVolume(scint_sealing_Cr_shape, fCr, "Scint_sealing_Cr_LV");
	G4LogicalVolume* logicalShape_scint = new G4LogicalVolume(scintillator, LYSO, "ScintillatorLV");

	G4LogicalVolume* logicShapeSiPM = new G4LogicalVolume(sipm, fSilicon, "SiPM_logical");

	//placement of volume
	G4VPhysicalVolume* physWorld =new G4PVPlacement(0,G4ThreeVector(),fElectricLogical,"World",0,false,0,checkOverlaps);
	new G4PVPlacement(yRot,pos1,glass_shell_logical,"Glass Shell",fElectricLogical,false,0,checkOverlaps);
	new G4PVPlacement(yRot,pos1,logicShape_Au_sealing,"Sealing_Au",fElectricLogical,false,0,checkOverlaps);
	new G4PVPlacement(yRot,pos1,logicShape_In_sealing,"Sealing_In",fElectricLogical,false,0,checkOverlaps);
	new G4PVPlacement(yRot,pos1,logicShape_Cr_sealing,"Photocathode",fElectricLogical,false,0,checkOverlaps);
	new G4PVPlacement(yRot,pos1,logicShape_photocathode,"Photocathode",fElectricLogical,false,0,checkOverlaps);	
	new G4PVPlacement(0, pos_basedisk_upper, Base_upperLV, "Base", fElectricLogical, false, 0, checkOverlaps);
	new G4PVPlacement(0, pos_basedisk_lower1, Base_lower1LV, "Base", fElectricLogical, false, 0, checkOverlaps);
	new G4PVPlacement(0,pos_basedisk_lower1,Base_lower1_coating_LV,"Base",fElectricLogical,false,0,checkOverlaps);
	G4VPhysicalVolume*bottom_base_plate_Phys = new G4PVPlacement(0,pos_basedisk_lower2,Base_lower2LV,"Base",fElectricLogical,false,0,checkOverlaps);
	new G4PVPlacement(0,pos_scintillator_Au,logicShape_Sealing_Au,"Scint_Sealing_Au",fElectricLogical,false,0,checkOverlaps);
	G4VPhysicalVolume* scintillator_Cr_Phys= new G4PVPlacement(0,pos_scintillator_Cr,logicShape_Sealing_Cr,"Scint_Sealing_Cr",fElectricLogical,false,0,checkOverlaps);
	G4VPhysicalVolume* scintillator_LYSO = new G4PVPlacement(0,pos_scintillator,logicalShape_scint,"Scintillator",fElectricLogical,false,0,checkOverlaps);
	housing->setPosition(pos_sipm);
    G4VPhysicalVolume* sipm_phys = housing->buildAndPlace(physWorld);

	///////////////////////////////////////
	// Optical property of Scintillator //
	/////////////////////////////////////

	G4OpticalSurface* opScintSurface_world = new G4OpticalSurface("Scintillator_Coating");
	opScintSurface_world->SetType(dielectric_dielectric); opScintSurface_world->SetFinish(polished); opScintSurface_world->SetModel(unified);
	G4double pp[] = { 2.0 * eV, 3.5 * eV }; const G4int num = sizeof(pp) / sizeof(G4double); G4double reflectivity[] = { 0.6, 0.6 }; G4double efficiency[] = { 0.0, 0.0 };

	G4MaterialPropertiesTable* scintWrapProperty = new G4MaterialPropertiesTable();
	scintWrapProperty->AddProperty("REFLECTIVITY", pp, reflectivity, num);
	scintWrapProperty->AddProperty("EFFICIENCY", pp, efficiency, num);
	opScintSurface_world->SetMaterialPropertiesTable(scintWrapProperty);

	G4LogicalBorderSurface* surface1= new G4LogicalBorderSurface("Scintillator_Coat_LBS", scintillator_LYSO, scintillator_Cr_Phys, opScintSurface_world);
    G4LogicalBorderSurface* surface1b= new G4LogicalBorderSurface("Scintillator_Coat_2_LBS", scintillator_Cr_Phys,scintillator_LYSO, opScintSurface_world);


	G4OpticalSurface* opScintSurface_other = new G4OpticalSurface("Scintillator_world");
	opScintSurface_other->SetType(dielectric_dielectric); opScintSurface_other->SetFinish(polished); opScintSurface_other->SetModel(unified);
	G4double pp_other[] = { 2.0 * eV, 3.5 * eV }; const G4int num_other = sizeof(pp) / sizeof(G4double); G4double reflectivity_other[] = { 0.2, 0.2 }; G4double efficiency_other[] = { 0.0, 0.0 };

	G4MaterialPropertiesTable* scintWrapProperty_other = new G4MaterialPropertiesTable();
	scintWrapProperty_other->AddProperty("REFLECTIVITY", pp_other, reflectivity_other, num_other);
	scintWrapProperty_other->AddProperty("EFFICIENCY", pp_other, efficiency_other, num_other);
	opScintSurface_other->SetMaterialPropertiesTable(scintWrapProperty_other);

	G4LogicalBorderSurface* surface2 = new G4LogicalBorderSurface("Scintillator_world_LBS", scintillator_LYSO, physWorld, opScintSurface_other);

    
	G4OpticalSurface* opScintSurface_sipm = new G4OpticalSurface("Scintillator_SiPM");
	opScintSurface_other->SetType(dielectric_dielectric); opScintSurface_other->SetFinish(polished); opScintSurface_other->SetModel(unified);
    G4double pp_sipm[] = { 2.0 * eV, 3.5 * eV }; const G4int num_sipm = sizeof(pp) / sizeof(G4double); G4double reflectivity_sipm[] = { 0.1, 0.1 }; G4double efficiency_sipm[] = { 0.0, 0.0 };
    
    G4MaterialPropertiesTable* scintWrapProperty_sipm = new G4MaterialPropertiesTable();
	scintWrapProperty_sipm->AddProperty("REFLECTIVITY", pp_sipm, reflectivity_sipm, num_sipm);
	scintWrapProperty_sipm->AddProperty("EFFICIENCY", pp_sipm, efficiency_sipm, num_sipm);
	opScintSurface_sipm->SetMaterialPropertiesTable(scintWrapProperty_sipm);

    G4LogicalBorderSurface* surface3 = new G4LogicalBorderSurface("Scintillator_SiPM_LBS", scintillator_LYSO, sipm_phys, opScintSurface_sipm);
	G4LogicalBorderSurface* surface3b = new G4LogicalBorderSurface("Scintillator_SiPM_2_LBS", sipm_phys, scintillator_LYSO, opScintSurface_sipm);
	
    
    /////////////////////////////////////////////////
	////////// set color for the geometries ////////
	///////////////////////////////////////////////

	auto visAttributes = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	visAttributes = new G4VisAttributes(G4Colour(0.905882, 1., 0.996078, 1.));
	glass_shell_logical->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.996078, 0.909804, 0.239216, 1));
	logicShape_photocathode->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.333333, 1, 0, 1));
	logicalShape_scint->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(1., 0.827451, 0.737255, 1.));
	Base_lower1LV->SetVisAttributes(visAttributes);
	Base_lower2LV->SetVisAttributes(visAttributes);
	Base_upperLV->SetVisAttributes(visAttributes);


	return physWorld;
}



void ABDetectorConstruction::ConstructSDandField()
{

	//auto SipmSD= new ABCalorimeterSD("SiPMSD", "AllHitsCollection", 1); // or set fNofLayers to 1?  
	//G4SDManager::GetSDMpointer()->AddNewDetector(SipmSD);
	//SetSensitiveDetector("SiPM_logical", SipmSD);

	// Electric field ----------------------------------------------------------

	fElectricField = new ABElectricField();
	G4EqMagElectricField*fEquation = new G4EqMagElectricField(fElectricField);  //fEMfield<->fElectricField
	std::cerr << "Electric Field Constructed..." << std::endl;

	G4int nvar = 8;
	fStepper = new G4ClassicalRK4(fEquation,nvar);

	//global Feild Manager
	fFieldMgr = new G4FieldManager();
	fFieldMgr->SetDetectorField(fElectricField); //fEMfield<->fElectricField

	fMinStep = 0.1* nm; // minimal step of 10 microns //1*mm  0.0005?!

	// just added 

	G4FieldManager* globalFieldManager;

	G4TransportationManager* transportMgr =
		G4TransportationManager::GetTransportationManager();

	double MaxTrackingStep = 0.00001;
	globalFieldManager = transportMgr->GetFieldManager();

	// Relative accuracy values:
	G4double minEps = 1.0e-5 * mm;  //   Minimum & value for smallest steps  (1e-5)
	G4double maxEps = 1.0e-4 *mm;  //   Maximum & value for largest steps  (1e-4)

	fFieldMgr->SetDeltaOneStep(0.5e-4 * mm);  // 0.5 micrometer

	// just added 

	G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep,
		fStepper,
		fStepper->GetNumberOfVariables());



	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldMgr->SetChordFinder(fChordFinder);
	fElectricLogical->SetFieldManager(fFieldMgr, true);
	
	transportMgr->GetPropagatorInField()->SetLargestAcceptableStep(10.*mm);

	// ** the problem of particle passing through volume is due to the step size 
	// problem each step the particle is losing all its energy
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
