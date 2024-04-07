// Local headers
#include "BeamTestPrimaryGeneratorAction.hh"
#include "BeamTestPrimaryGeneratorMessenger.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <stdlib.h>

// Geant4 headers
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "DetectorConstruction.hh"

#define pi 3.141592
#define N 465
#define N_inel 129
using namespace CLHEP;
using namespace std;

G4ParticleGun *BeamTestPrimaryGeneratorAction::particleGun(0);

// ----------------------------------------------------------------------------- //
BeamTestPrimaryGeneratorAction::BeamTestPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  gunMessenger = new BeamTestPrimaryGeneratorMessenger(this);
}

// ----------------------------------------------------------------------------- //

BeamTestPrimaryGeneratorAction::~BeamTestPrimaryGeneratorAction()
{
  delete particleGun;
  delete gunMessenger;
}

// ----------------------------------------------------------------------------- //

  double BeamTestPrimaryGeneratorAction::omega(double x, double y, double z){
    return sqrt(x*x + y*y + z*z -2*x*y -2*y*z -2*x*z);

  }

  void BeamTestPrimaryGeneratorAction::KinElast(double mbeam, double mtar, double Kb){

    double inte = 0;
    double ran_F = 0;
    double thetacm = 0;
    double phiela =0;
    double th_ini = 3.0*pi/180;
    double th_fin = 6.0*pi/180;

    ran_F = G4UniformRand();
    thetacm = 180/pi*pow(1/(ran_F*(pow(1/th_fin,3)) + pow(1/th_ini,3)*(1-ran_F) ) ,0.333);

    //G4cout<<thetacm<<G4endl;

    double mbg1 = mbeam;
    double mbg2 = mtar;
    double mbg3 = mbg1;
    double mbg4 = mbg2;
    double Et1 = Kb + mbg1;
    double Et2 = mbg2;
    double s = pow(mbg1,2) + pow(mbg2,2) +2*mbg2*Et1;
    double a =4.*mbg2*s;
    double b =(pow(mbg3,2)-pow(mbg4,2)+s)*(pow(mbg1,2)-pow(mbg2,2)-s);
    double c = omega(s,pow(mbg1,2),pow(mbg2,2))*omega(s,pow(mbg3,2),pow(mbg4,2));
    double Et3 = (c*cos((thetacm)*3.1415/180)-b)/a;
	  double Et4 = (Et1 + mbg2 - Et3);
	  double K3 = Et3 - mbg3;
    //double K3 = 2.0;
	  double K4 = Et4 - mbg4;
    double t =  pow(mbg2,2) + pow(mbg4,2) - 2*mbg2*Et4;
	  double u =  pow(mbg2,2) + pow(mbg3,2) - 2*mbg2*Et3;
    double thetalab = acos(((s-pow(mbg1,2)-pow(mbg2,2))*(pow(mbg2,2)+pow(mbg3,2)-u)+2.*pow(mbg2,2)*(pow(mbg2,2)+pow(mbg4,2)-s-u))/(omega(s,pow(mbg1,2),pow(mbg2,2))*omega(u,pow(mbg2,2),pow(mbg3,2))));
    //SetThetaEls(thetalab*180/3.1415);
    //std::cout << " el angulo  : "<< thetalab*180/3.1415<< '\n';
    //std::cout << " la energia  : "<< K3<< '\n';
    //std::cout << " angulo cm  : "<< thetacm<< '\n';
    double Pbg3 = sqrt(pow((K3 + mbg2),2) - pow(mbg2,2))/1000.0;
    double pzela = Pbg3*cos(thetalab);
    double pxela = Pbg3*sin(thetalab)*cos(phiela);
    double pyela = Pbg3*sin(thetalab)*sin(phiela);
    //TVector3 Pelas(pxela,pyela,pzela);
    //SetPels(Pelas);

    SetThprod(thetacm);
    SetKprod(K3);

  }

void BeamTestPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  Inputs *Inputs = &Inputs::GetInputs();

  fposition = Inputs->primary_pos;
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4IonTable *ionTable = G4IonTable::GetIonTable();
  G4String particleName;

  G4ParticleDefinition *part = ionTable->GetIon(Inputs->primary_Z, Inputs->primary_A, 0.);;
  //G4double Energy = Inputs->primary_energy * MeV;

  particleGun->SetParticleDefinition(part);
  //particleGun->SetParticleEnergy(Energy);

  // ---------------------- Primary beam as a dot ------------------------- //

  // G4double phi = 2 * pi * G4UniformRand();
  // G4double theta;
  // do
  // {
  //   theta =
  //       0.066667 * pi * (G4UniformRand() - 0.5)
  //       // pi/2                                   //90°
  //       // pi/3                                   //60°
  //       // pi/6                                   //30°
  //       // (pi/6) + ( (pi/3) * G4UniformRand() )  //30° - 90°
  //       // (pi/6) + ( 0.261799 * G4UniformRand() )//30° - 45°
  //       // (pi/6) * G4UniformRand()               //0°  - 30°
  //       // (pi/48)*G4UniformRand()                //0°  - 15°
  //       // (pi/3) + ( (pi/6) * G4UniformRand() )  //60° - 90°
  //       // 0.0872665 * G4UniformRand()            //0°  - 5°
  //       // (1 * G4UniformRand() + 2) * pi / 180 //2°  - 6°

  //       ;
  // }

  // while (abs(theta) < 0.0349);

  // particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)));

  // ----------------------------------------------------------------------- //

  //------------- Primary beam as a Gaussian -------------//

  G4double x, y, z;
  G4double sigma_beam;
  G4double x0, y0, z0;

  x = 0. * mm;
  y = 0. * mm;
  z = 0. * mm;

  //sigma_beam = 2 * mm / 2.35; // alta intensidad 3mm-sigma
  sigma_beam = 0.0001 * mm / 2.35; // alta intensidad 3mm-sigma

  y0 = 0. * mm;
  x0 = 0. * mm;
  z0 = -250. * mm;

  //z0 = -115. * cm;
  //z0 = 148. * cm;

  G4double a = 0.5 * mm, b = 1.0 * mm;

  for (;;)
  {
    x = G4RandGauss::shoot(x0, sigma_beam);
    y = G4RandGauss::shoot(y0, sigma_beam);
    z = z0;

    if ( sqrt(x*x +y*y) <  3 *  sigma_beam )
    {
      break;
    }
    /*
    if ((((double)(x / a) * (x / a)) + (double)(y / b) * (y / b) <= 1))
    {
      break;
    }
    */
  }

  z = z0;
  G4ThreeVector fposition = G4ThreeVector(x, y, z);


  //KinElast(11.009305404*931.5, 26.981538578*931.5, Inputs->primary_energy);
   //KinElast(12.0*931.5, 6.01512288741*931.5, Inputs->primary_energy);

  //G4double theta = th_prod*pi/180;
  G4double theta = 0.0*pi/180;
  G4double phi = 2 * pi * G4UniformRand();

  //G4double Energy = K_prod * MeV;
  G4double Energy = Inputs->primary_energy * MeV;
  particleGun->SetParticleEnergy(G4RandGauss::shoot(Energy, 3.0e-2*Energy));

/*
  do
  {
    //theta = 0.066667*pi*(G4UniformRand()-0.5);
  } while (abs(theta) < 0.0349);
*/
  //-----------------------------------------------------------------------------//
  //theta = (180+45)*pi/180;
  //theta = (180+130)*pi/180;
  //theta = (180-90)*pi/180;
  //phi = 0*pi/180;

  //particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)));
  particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0, 1));
  particleGun->SetParticlePosition(fposition);
  particleGun->GeneratePrimaryVertex(anEvent);
}
