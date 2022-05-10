#include "ParticleSourceMessenger.hh"

#include <boost/program_options.hpp>

#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>
#include <G4OpticalPhoton.hh>

#include "G4UiMessengerUtil.hh"
//#include "ProgramOptionsUtil.hh"

const std::string ParticleSourceMessenger::MACRO_FILE_NAME = "ps.mac";

ParticleSourceMessenger::ParticleSourceMessenger() :
		G4UImessenger() {
	setDefaultValues();
	// Create directories.
	new G4UIdirectory("/ps/");
	new G4UIdirectory("/ps/energy/");
	new G4UIdirectory("/ps/angle/");
	new G4UIdirectory("/ps/plane/");
	// Create basic commands.
	verboseCmd = G4UiMessengerUtil::createCmd(this, "/ps/", "verbose", verbose);
	nParticlesCmd = G4UiMessengerUtil::createCmd(this, "/ps/", "nParticles", nParticles);
	polarCmd = G4UiMessengerUtil::createCmd(this, "/ps/", "polar", polar, "deg");
	// Create time commands.
	tMinCmd = G4UiMessengerUtil::createCmd(this, "/ps/", "tMin", tMin, "s");
	tMaxCmd = G4UiMessengerUtil::createCmd(this, "/ps/", "tMax", tMax, "s");
	tInputCmd = G4UiMessengerUtil::createCmd(this, "/ps/time/", "input", tInput);
	// Create energy commands.
	eMinCmd = G4UiMessengerUtil::createCmd(this, "/ps/energy/", "eMin", eMin, "eV");
	eMaxCmd = G4UiMessengerUtil::createCmd(this, "/ps/energy/", "eMax", eMax, "eV");
	eInputCmd = G4UiMessengerUtil::createCmd(this, "/ps/energy/", "input", eInput);
	// Create angle commands.
	phiMinCmd = G4UiMessengerUtil::createCmd(this, "/ps/angle/", "phiMin", phiMin, "deg");
	phiMaxCmd = G4UiMessengerUtil::createCmd(this, "/ps/angle/", "phiMax", phiMax, "deg");
	thetaMinCmd = G4UiMessengerUtil::createCmd(this, "/ps/angle/", "thetaMin", thetaMin, "deg");
	thetaMaxCmd = G4UiMessengerUtil::createCmd(this, "/ps/angle/", "thetaMax", thetaMax, "deg");
	// Create plane commands.
	surfaceNormalCmd = G4UiMessengerUtil::createCmd(this, "/ps/plane/", "surfaceNormal", surfaceNormal, "mm");
	aCmd = G4UiMessengerUtil::createCmd(this, "/ps/plane/", "a", a, "mm");
	bCmd = G4UiMessengerUtil::createCmd(this, "/ps/plane/", "b", b, "mm");
	posCmd = G4UiMessengerUtil::createCmd(this, "/ps/plane/", "pos", pos, "mm");
	// Set guidance.
	polarCmd->SetGuidance("Set linear polarization angle w.r.t. (k,n) plane. "
			"Negative values will trigger a random polarization.");
	// Set parameter boundaries.
	verboseCmd->SetParameterName("verbose", false);
	verboseCmd->SetRange("verbose >= 0 && verbose <= 2");
	nParticlesCmd->SetParameterName("nParticles", false);
	nParticlesCmd->SetRange("nParticles > 0");
	// Initialize from macro.
	G4UiMessengerUtil::executeMacro(MACRO_FILE_NAME, true);
}

ParticleSourceMessenger::~ParticleSourceMessenger() {
	delete verboseCmd;
	delete nParticlesCmd;
	delete polarCmd;
	delete tMinCmd;
	delete tMaxCmd;
	delete tInputCmd;
	delete eMinCmd;
	delete eMaxCmd;
	delete eInputCmd;
	delete phiMinCmd;
	delete phiMaxCmd;
	delete thetaMinCmd;
	delete thetaMaxCmd;
	delete surfaceNormalCmd;
	delete aCmd;
	delete bCmd;
	delete posCmd;
}

void ParticleSourceMessenger::setDefaultValues() {
	verbose = 0;
	nParticles = 1;
	polar = -360 * deg;
	tMin = 0;
	tMax = 0;
	eMin = 1 * eV;
	eMax = 10 * eV;
	phiMin = 0;
	phiMax = 360 * deg;
	thetaMin = 0;
	thetaMax = 90 * deg;
	surfaceNormal = G4ThreeVector(0., 0., -1.);
	a = 1 * mm;
	b = 1 * mm;
	pos = G4ThreeVector(0, 0, 5 * cm);
}

void ParticleSourceMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {

}

template<typename T> void setNewValueIfArgumentGiven(const boost::program_options::variables_map& vm, std::string key,
		T* value) {
	if (vm.count(key)) {
		*value = vm[key].as<T>();
	}
}

void setNewValueIfArgumentGiven(const boost::program_options::variables_map& vm, std::string key, double* value,
		double unit = 1.) {
	if (vm.count(key)) {
		*value = vm[key].as<double>() * unit;
	}
}

void setNewValueIfArgumentGiven(const boost::program_options::variables_map& vm, std::string key, G4String* value) {
	if (vm.count(key)) {
		*value = vm[key].as<std::string>();
	}
}

void ParticleSourceMessenger::parseProgramOptions(int argc, char** argv) {
	
}

ParticleSourceMessenger* ParticleSourceMessenger::getInstance() {
	static ParticleSourceMessenger* instance = new ParticleSourceMessenger;
	return instance;
}

int ParticleSourceMessenger::getVerbose() const {
	return verbose;
}

void ParticleSourceMessenger::setVerbose(int _verbose) {
	verbose = _verbose;
}

int ParticleSourceMessenger::getNParticles() const {
	return nParticles;
}

void ParticleSourceMessenger::setNParticles(int _nParticles) {
	nParticles = _nParticles;
}

double ParticleSourceMessenger::getPolar() const {
	return polar;
}

void ParticleSourceMessenger::setPolar(double _polar) {
	polar = _polar;
}

double ParticleSourceMessenger::getEMax() const {
	return eMax;
}

void ParticleSourceMessenger::setEMax(double _eMax) {
	eMax = _eMax;
}

double ParticleSourceMessenger::getEMin() const {
	return eMin;
}

void ParticleSourceMessenger::setEMin(double _eMin) {
	eMin = _eMin;
}

double ParticleSourceMessenger::getPhiMax() const {
	return phiMax;
}

void ParticleSourceMessenger::setPhiMax(double _phiMax) {
	phiMax = _phiMax;
}

double ParticleSourceMessenger::getPhiMin() const {
	return phiMin;
}

void ParticleSourceMessenger::setPhiMin(double _phiMin) {
	phiMin = _phiMin;
}

double ParticleSourceMessenger::getThetaMax() const {
	return thetaMax;
}

void ParticleSourceMessenger::setThetaMax(double _thetaMax) {
	thetaMax = _thetaMax;
}

double ParticleSourceMessenger::getThetaMin() const {
	return thetaMin;
}

void ParticleSourceMessenger::setThetaMin(double _thetaMin) {
	thetaMin = _thetaMin;
}

G4ThreeVector ParticleSourceMessenger::getPos() const {
	return pos;
}

void ParticleSourceMessenger::setPos(G4ThreeVector _pos) {
	pos = _pos;
}

G4ThreeVector ParticleSourceMessenger::getSurfaceNormal() const {
	return surfaceNormal;
}

void ParticleSourceMessenger::setSurfaceNormal(G4ThreeVector _surfaceNormal) {
	surfaceNormal = _surfaceNormal;
}

double ParticleSourceMessenger::getA() const {
	return a;
}

void ParticleSourceMessenger::setA(double _a) {
	a = _a;
}

double ParticleSourceMessenger::getB() const {
	return b;
}

void ParticleSourceMessenger::setB(double _b) {
	b = _b;
}

double ParticleSourceMessenger::getTMin() const {
	return tMin;
}

void ParticleSourceMessenger::setTMin(double _tMin) {
	tMin = _tMin;
}

double ParticleSourceMessenger::getTMax() const {
	return tMax;
}

void ParticleSourceMessenger::setTMax(double _tMax) {
	tMax = _tMax;
}

G4String ParticleSourceMessenger::getTInput() const {
	return tInput;
}

G4String ParticleSourceMessenger::getEInput() const {
	return eInput;
}
