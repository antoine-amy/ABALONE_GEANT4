//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file ABElectricField.cc
/// \brief Implementation of the ABElectricField class

#include "ABElectricField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <stdlib.h>
#include <typeinfo>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::string filename_datar = "../E_fields/x_Efield_25.csv";
std::string filename_dataz = "../E_fields/z_Efield_25.csv";
std::vector<std::vector<float> > data_r;
std::vector<std::vector<float> > data_z;
double r_min, r_max, z_min, z_max, interpolation_length_r, interpolation_length_z, seperation;
int count = 0;

void import_field1(std::string file_name, std::vector<std::vector<float> >& data) {
	
	std::string row;
	std::ifstream init_file(file_name.c_str());

	// open file
	if (init_file.is_open()) {
		std::cerr << "Opening Electric field file : "<< file_name << " ... " << std::endl;
		// loop in each line in the file
		int count_line = 0;
		while (getline(init_file, row)) {
			count_line++;
			std::istringstream iss(row);
			// initialize a vector to store the row elements
			std::vector<float> row;
			std::string token;
			// get each element by splitting the row string by commas
			while (std::getline(iss, token, ',')) {
				// convert the string to int (or use stof for floats)
				row.push_back((float)atof(token.c_str()));
			}
			data.push_back(row);
			//std::cerr << "Reading the " << count_line << " th line in file " << file_name <<" ... "<< std::endl;
		}
		init_file.close();
		std::cerr << "E filed loaded " << std::endl;
	}
	else
		std::cerr << "Unable to open file" << std::endl;
	return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABElectricField::ABElectricField(): G4ElectricField(), fMessenger(nullptr){ // E_x(0.1*volt/cm), E_y(0.1*volt/cm), E_z(0.1*volt/cm) //fEy useless
	std::cerr<<"Reading r-direction E field"<<std::endl;
	import_field1(filename_datar, data_r);
	std::cerr<<"Reading z-direction E field"<<std::endl;
	import_field1(filename_dataz, data_z);

	r_min=data_r[0][0]; r_max=data_r[0][1]; z_min=data_r[1][0]; z_max=data_r[1][1];
	interpolation_length_r=data_r[0][2]; interpolation_length_z=data_r[1][2];
	seperation=((r_max-r_min)/(interpolation_length_r-1.));

	DefineCommands();

	std::cerr << r_min << r_max << z_min << z_max << interpolation_length_r  << interpolation_length_z << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ABElectricField::~ABElectricField(){ 
  delete fMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABElectricField::GetFieldValue(const G4double point[4],double *eField) const{
	double x=(point[0])/cm, y=(point[1])/cm, z=(point[2])/cm; // in mm so we must convert it to cm / m
	double t=point[3];
	double r_particle=sqrt((std::pow(x, 2))+(std::pow(y, 2))+(std::pow(z, 2))); 
	count=count+1;

 	if (r_particle<=7.0){
		double sqrt_r=sqrt((std::pow(x, 2))+(std::pow(y, 2)));
		double index_z=(z-z_min)/seperation;
		double index_r=(sqrt_r-r_min)/seperation+1.0;

		double x11=floor(index_r); double y11=floor(index_z); double x22=x11+1.0; double y22=y11+1.0;

		int x1=x11; int y1=y11; int x2=x22; int y2=y22; 

		//double Er_11=data_r[y1][x1]; double Er_12=data_r[y2][x1]; double Er_21=data_r[y1][x2]; double Er_22=data_r[y2][x2];
		//double Ez_11=data_z[y1][x1]; double Ez_12=data_z[y2][x1]; double Ez_21=data_z[y1][x2]; double Ez_22=data_z[y2][x2];
		
		//double E_r=((y22-index_z)*((x22-index_r)*Er_11+(index_r-x11)*Er_21)+(index_z-y11)*((x22-index_r)*Er_12+(index_r-x11)*Er_22));
		//double E_zz=((y22-index_z)*((x22-index_r)*Ez_11+(index_r-x11)*Ez_21)+(index_z-y11)*((x22-index_r)*Ez_12+(index_r-x11)*Ez_22));

		//G4double E_x=E_r*x/sqrt_r *volt/m; G4double E_y=E_r*y/sqrt_r *volt/m; G4double E_z=E_zz*volt/m;
		
		eField[0] = 0.; //Bx
		eField[1] = 0.; //By
		eField[2] = 0.; //Bz
		eField[3] =0; //eField[3] = E_x; //E_x 
		eField[4] = 0; //eField[4] = E_y; //E_y
		eField[5] = 0; //eField[5] = E_z; //E_z
		eField[6] = 0.; //Useless
	}else{
		eField[0] = 0.; //Bx
		eField[1] = 0.; //By
		eField[2] = 0.; //Bz
		eField[3] = 0; //E_x
		eField[4] = 0; //E_y
		eField[5] = 0; //E_z
		eField[6] = 0.; //Useless
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ABElectricField::DefineCommands(){
  fMessenger=new G4GenericMessenger(this, "/AB/field/", "Field control"); // Define /AB/field command directory using generic messenger class
  
  auto& valueCmd=fMessenger->DeclareMethodWithUnit("value","tesla", &ABElectricField::SetField, "Set field strength."); // fieldValue command 
  
  valueCmd.SetParameterName("field", true);
  valueCmd.SetDefaultValue("1.");
}
