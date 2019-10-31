///// Fireball burst GRB model
// This calculates broken spectra for neutrinos from prompt and afterglow GRBs
// given redshift, gamma fluence and T90.
// Full theoretical details: 
// https://arxiv.org/pdf/astro-ph/9802280.pdf

// Relevant past searches
// https://arxiv.org/pdf/astro-ph/0605480.pdf - RICE
// https://arxiv.org/pdf/1102.3206.pdf - ANITA-II
/////
#include <iostream>
#include <libgen.h>
#include <string>
#include <map>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TMath.h" 
#include "TObjString.h" 
#include "TTimeStamp.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMarker.h"
#include "SkyMap.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TSpline.h"

void computeGRBSpectra()
{
  TFile *file = new TFile("./data/GRBIceCube.root"); 
  TTree *tree = (TTree*) file->Get("iceCubeTree");
  UInt_t entries = tree->GetEntries();

  // Root vars
  float RA = 0.;
  float dec = 0.;
  std::string *name = new std::string;
  int unixTriggerTime = 0;
  float T90 = 0.;
  float redshift = 0.;
  float fluence = 0.;
  
  // New vars
  int a_tmin;
  int a_tmax;
  int passed = 0;

  // Neutrino flux spectra
  const int grbNumberLimit = 25; // all time count = 6326, A4 = 25, A3 = 18
  
  // Graphing
  const int newPts=1000;
  double minAllowedNuE = 1e+4;
  double maxAllowedNuE = 1e+12; // in GeV
  int resolution = 1000;
  TGraph *stitchedSpectrum[grbNumberLimit];
  TGraph *stitchedSpectrumShock[grbNumberLimit];
  TGraph *stitchedSpectrumShockLog[grbNumberLimit];
  TGraph *spectrum[grbNumberLimit];
  double stitchedEnergy[newPts];
  double stitchedEnergyLog[newPts];
  double stitchedFlux[newPts];
  TGraph *spectrumShock[grbNumberLimit];
  double stitchedEnergyShock[newPts];
  double stitchedEnergyShockLog[newPts];
  double stitchedFluxShock[newPts];
  TH1F *histEnergyCutoff = new TH1F("histEnergyCutoff","histEnergyCutoff",100,6.,13.);
  
  // Cosmology
  double omegaRad = 9.236e-5;
  double omegaM = 0.315;
  double omegaLambda = 0.679;
  double H0 = 67.4 * pow(10,3); // We want this in ms^-1/Mpc, not the standard kms^-1/Mpc
  double c = 299792458; // ms^-1, as we we calculate c/H0
  double EGammaIso = 0.;
  double LGammaIso = 0.;
  double dL = 0.; // we will calculate in Mpc
  double MpcToCm = 3.08e+24;

  // GRB model vars
  double energyBreakGamma = 0.;
  double energyBreakNeutrino = 0.;
  double energyBreakNeutrinoSynchrotron = 0.;
  double r = 0.;
  double magField = 0.;
  double fluxPrefactor = 0.;
  
  // Model assumptions
  double lorentzFactor = 300.;
  double tV = 0.01;
  double eB = 0.1;
  double eE = 0.1;
  double fPi = 0.2;
  double alpha = 1;
  double beta = 2;

  /// Afterglow specific calcs
  double nEx = 1.; // cm^-3
  double EKinIso = 0.;
  double shockLorentzFactor = 0.;
  double energyCGamma = 0.;
  double fPiAG = 0.;
  double energyBreakNeutrinoAG = 0.;
  double protonGammaSigma = 5e-28; // cm^2
  double shockMagneticField = 0.;
  double shockRadius = 0.;
  double massProton = 0.0015; // erg/c^2
  double magFieldShock = 0.;
  double maxShockNuE = 0.;
  double shockPrefactor = 0.;

  // A4 for now
  a_tmin = 1480707643;
  a_tmax = 1483004563;

  // A3
  //a_tmin = 1418938406;
  //a_tmax = 1420777814;
  
  // A2 for comparison
  //a_tmin = 1230422400;
  //a_tmax = 1232236800;

  // others
  //a_tmin = 1117756800;
  //a_tmax = 1126310400;
  
  // Branches
  tree->SetBranchAddress("RA",&RA);
  tree->SetBranchAddress("dec",&dec);
  tree->SetBranchAddress("name",&name);
  tree->SetBranchAddress("unixTriggerTime",&unixTriggerTime);
  tree->SetBranchAddress("T90",&T90);
  tree->SetBranchAddress("redshift",&redshift);
  tree->SetBranchAddress("fluence",&fluence);
  
  TF1 *flux1[grbNumberLimit]; TF1* afterGlowFlux1[grbNumberLimit];
  TF1 *flux2[grbNumberLimit]; TF1* afterGlowFlux2[grbNumberLimit];
  TF1 *flux3[grbNumberLimit]; TF1* afterGlowFlux3[grbNumberLimit];

  double quasiDiffuseFlux[newPts] = {};
  double quasiDiffuseFluxShock[newPts] = {};
  TGraph* grQuasiDiffuse;
  TGraph* grQuasiDiffuseLog;
  TGraph* grQuasiDiffuseShock;
  TGraph* grQuasiDiffuseShockLog;
  /*
  TCanvas *can = new TCanvas("c","c",2000,1000);
  can->SetLogx();
  can->SetLogy();
  */
  for(unsigned int i = 0; i < 6326; i++) // do 2 for now but should be entries
    {
      tree->GetEntry(i);

      if(passed == grbNumberLimit){break;} // if you only want to look at a specified amount of GRBs
      if( (unixTriggerTime < a_tmin) || (unixTriggerTime > a_tmax) ){continue;}
      //if(abs(dec) > 30){continue;}

      //if(*name != "GRB090112B" && *name != "GRB090107A"){continue;} // just for comparison
      //if(*name != "GRB090112B"){continue;} // just for comparison
      
      // GRB spectra calculation starts
      // Some T90s are not recorded, set to 20s.
      if(T90==-999){T90 = 20;}
      // Some redshifts are not calculated, set to 2.
      if(redshift==-999){redshift = 2;}
      // Just set the fluence to this for now
      if(fluence==-999){fluence=1e-6;}
      
      // test variables
      /*
      if(passed==0)
	{
          redshift = 3.35;
	  T90 = 20;
	  fluence = 5.1e-7;
	}
      else
	{
	  redshift = 3.97;
	  T90 = 155;
	  fluence = 4.4e-6; 
	}
      */
      //
            
      // Cosmology stuff
      // Assume flat Universe
      // Curvature is k = 0, omegaK = 0
      
      TF1 Ez("E(z)", "1/(pow([0] * pow((1+x),3) + [1] + [2] * pow((1+x),4),0.5))", 0, redshift);
      //TF1 Ez("E(z)", "1/(pow([0] * pow((1+x),3) + [1],0.5))", 0, redshift);
      Ez.SetParameter(0,omegaM); Ez.SetParName(0,"omegaM");
      Ez.SetParameter(1,omegaLambda); Ez.SetParName(1,"omegaLambda");
      Ez.SetParameter(2,omegaRad); Ez.SetParName(2,"omegaRad");
      ROOT::Math::WrappedTF1 wEz(Ez);
      ROOT::Math::GaussIntegrator ig;
      ig.SetFunction(wEz);
      ig.SetRelTolerance(0.001);
      double integral =  ig.Integral(0, redshift);
      
      // Now calc lum dist
      dL = (1+redshift) * c/H0 * integral; // Mpc
      
      // Calc gamma-ray bolometric energy
      EGammaIso = 4 * TMath::Pi() * dL*dL * fluence * 1/(1+redshift); // Mpc^2 * erg cm^-2
      EGammaIso*=pow(MpcToCm,2); // erg
      
      // Calc GRB luminosity
      LGammaIso = 4 * TMath::Pi() * dL*dL * fluence *1/T90; // Mpc^2 * erg cm^-2 s^-1
      LGammaIso*=pow(MpcToCm,2); // erg s^-1
      
      // Use Band fit as a broken power-law
      // Calculate break energy for each spec (Ghirlanda relationship)
      energyBreakGamma = 300/(1 + redshift) * pow((EGammaIso/(1e53)),0.56); // in keV
      
      // We have photon break energy, now calc for neutrinos.
      energyBreakNeutrino = 0.015 * pow((lorentzFactor/(1+redshift)),2) * pow((energyBreakGamma * 1e-6),-1); // in GeV

      //r = 2 * lorentzFactor * lorentzFactor * c * tV; // m
      //magField = pow(((2*eB * LGammaIso)/(eE * r * r * lorentzFactor * lorentzFactor * c)),0.5)
      magField = 5e+4 * pow(((LGammaIso)/(1e+52)),0.5) * pow(((lorentzFactor)/(300.)),-3) * pow(((tV)/(0.01)),-1); // in Gauss

      double randomFactor = 1;
      // Synchrotron break energy
      energyBreakNeutrinoSynchrotron = 1e11 * lorentzFactor/(4*(1+redshift)) * pow(magField,-1) * randomFactor; // in GeV
      // 
      // Also different cosomological parameters
      
      // Calculate broken spectrum
      fluxPrefactor = 1.56e-6 * fPi/0.2 * fluence/(1e-6) * pow(((T90)/(10)),-1); //in GeV cm^-2 s^-1
      // prefactor is correct

      //int color = (passed*3+29);
      int color = passed+1;
      if(color==5){color=106;} // dont use yellow
      if(passed ==8){color = 1;}

      flux1[passed] = new TF1("flux1", "[0] * pow((x/[1]),([2]-1))", minAllowedNuE, energyBreakNeutrino);
      // params and conditionals
      flux1[passed]->SetParameter(0,fluxPrefactor); flux1[passed]->SetParName(0,"fluxPrefactor");
      flux1[passed]->SetParameter(1,energyBreakNeutrino); flux1[passed]->SetParName(1,"energyBreakNeutrino");
      flux1[passed]->SetParameter(2,beta); flux1[passed]->SetParName(2,"beta");
      flux1[passed]->SetParameter(4,alpha); flux1[passed]->SetParName(4,"alpha");
      flux1[passed]->SetParameter(5,energyBreakNeutrinoSynchrotron); flux1[passed]->SetParName(5,"energyBreakNeutrinoSynchrotron");
      //
      flux1[passed]->SetLineColor(color);
      flux1[passed]->SetLineWidth(3);
      flux1[passed]->SetLineStyle(1);
      
      flux2[passed] = new TF1("flux2", "[0] * pow((x/[1]),[2]-1)", energyBreakNeutrino, energyBreakNeutrinoSynchrotron);
      flux2[passed]->SetParameter(0,fluxPrefactor); flux2[passed]->SetParName(0,"fluxPrefactor");
      flux2[passed]->SetParameter(1,energyBreakNeutrino); flux2[passed]->SetParName(1,"energyBreakNeutrino");
      flux2[passed]->SetParameter(2,alpha); flux2[passed]->SetParName(2,"alpha"); 
      flux2[passed]->SetLineColor(color);
      flux2[passed]->SetLineWidth(3);
      flux2[passed]->SetLineStyle(1);

      // We use the tail of this spectrum as its between e18 and e21 eV
      flux3[passed] = new TF1("flux3", "[0] * pow(([1]/[2]),([3]-1)) * pow((x/[1]),-2)", energyBreakNeutrinoSynchrotron, maxAllowedNuE);
      flux3[passed]->SetParameter(0,fluxPrefactor); flux3[passed]->SetParName(0,"fluxPrefactor");
      flux3[passed]->SetParameter(1,energyBreakNeutrinoSynchrotron); flux3[passed]->SetParName(1,"energyBreakNeutrinoSynchrotron");
      flux3[passed]->SetParameter(2,energyBreakNeutrino); flux3[passed]->SetParName(2,"energyBreakNeutrino");
      flux3[passed]->SetParameter(3,alpha); flux3[passed]->SetParName(3,"alpha");
      flux3[passed]->SetLineColor(color);
      flux3[passed]->SetLineWidth(3);
      flux3[passed]->SetLineStyle(1);
      
      flux1[passed]->SetNpx(resolution); flux2[passed]->SetNpx(resolution); flux3[passed]->SetNpx(resolution);

      if(passed==0)
	{
	  //can->DrawFrame(minAllowedNuE,1e-12,maxAllowedNuE,1e2,"Predicted broken flux spectra & quasi-diffuse spectrum for GRBs during the ANITA-4 flight;Neutrino energy (GeV);E^{2} dF/dE (GeV cm^{-2})");
          //can->DrawFrame(1e9,1e-12,1e12,1e1,"Quasi-diffuse spectra for GRBs;Neutrino energy (GeV);E^{2} dF/dE (GeV cm^{-2})");
	}
      
      //flux1[passed]->Draw("same");
      //flux2[passed]->Draw("same");
      //flux3[passed]->Draw("same");

      // Dealing with broken spectra is annoying.
      // For flux 1 from minAllowedNuE to energyBreakNeutrino
      const int pts = 1000;
      const int stitchedPts = 3*1000;
      // Equals steps from minAllowedNuE to energyBreakNeutrino in log scale
      
      // Broken spec 1
      double logMin = TMath::Log10(minAllowedNuE);
      double logEbn = TMath::Log10(energyBreakNeutrino);
      double divisionStep1 = (logEbn - logMin)/double(pts);
      double currentEnergy[stitchedPts] = {};
      double currentFlux[stitchedPts] = {};

      // Broken spec 2
      double logEbs = TMath::Log10(energyBreakNeutrinoSynchrotron);
      double divisionStep2 = (logEbs - logEbn)/double(pts);
      
      // Broken spec 3
      double logMax = TMath::Log10(maxAllowedNuE);
      double divisionStep3 = (logMax - logEbs)/double(pts);
      
      // For each spectrum _piece_:
      for(int j = 0; j < pts; j++)
	{
          //Stitch broken spectrum pieces together
	  currentEnergy[j] = pow(10,logMin + divisionStep1*j); currentFlux[j] = flux1[passed]->Eval(currentEnergy[j]);
	  currentEnergy[pts+j] = pow(10,logEbn + divisionStep2*j); currentFlux[pts+j] = flux2[passed]->Eval(currentEnergy[pts+j]);
	  currentEnergy[2*pts+j] = pow(10,logEbs + divisionStep3*j); currentFlux[2*pts+j] = flux3[passed]->Eval(currentEnergy[2*pts+j]);
          //cout << currentEnergy[j] << endl;
	}

      // Finally, a non-broken spectrum :)
      
      stitchedSpectrum[passed] = new TGraph(stitchedPts,currentEnergy,currentFlux);
      stitchedSpectrum[passed]->SetLineWidth(3);
      stitchedSpectrum[passed]->SetLineStyle(5);
      //stitchedSpectrum[passed]->Draw();
      
      // Equally spaced stitched spectra
      double divisionStepGlobal= (logMax - logMin)/double(newPts);
      
      // Now re-evaluate the stitched spectrum so we get equal spacing in our global energy range
      // Then we can sum stitched spectra
      for(int k = 0; k < newPts; k++)
	{
	  stitchedEnergy[k] = pow(10,logMin + divisionStepGlobal*k);
          stitchedEnergyLog[k] = logMin + divisionStepGlobal*k;
          stitchedFlux[k] = stitchedSpectrum[passed]->Eval(stitchedEnergy[k],0,"");
          //cout << "energy = " << stitchedEnergy[k] << ", flux = " << stitchedFlux[k] << endl;

          quasiDiffuseFlux[k]+=stitchedFlux[k];
        }

      spectrum[passed] = new TGraph(newPts, stitchedEnergy, stitchedFlux);
      spectrum[passed]->SetLineColor(5);
      spectrum[passed]->SetLineStyle(9);
      //spectrum[passed]->Draw();

      grQuasiDiffuse = new TGraph(newPts, stitchedEnergy,quasiDiffuseFlux);
      grQuasiDiffuseLog = new TGraph(newPts, stitchedEnergyLog,quasiDiffuseFlux);
      //

      //////////////
      // Afterglow
      //////////////

      // Calculate total jet kinetic energy and bulk Lorentz factor of GRB afterglow
      EKinIso = EGammaIso/eE; //erg
      
      // Radii, depths
      shockRadius = pow(((3*EKinIso)/(4*TMath::Pi() * nEx * massProton * c*c * lorentzFactor*lorentzFactor)),(1./3.)); // in cm
      shockRadius*=100; // in meters, to compare with internal rad
      if(shockRadius < r)
        {
          cout << "Afterglow shock radius < burst radius. Something might be wrong." << endl;
        }

      shockRadius/=100;
      magFieldShock = pow(8 * TMath::Pi() * eB * nEx * massProton,0.5);
      
      shockLorentzFactor = 195 * pow(((EKinIso)/(pow(10,54))),0.125) * pow(T90/10,-0.375) * pow(nEx/1.,-0.125);

      // Peak sync gamma energy radiated by electrons in B field:
      energyCGamma = 0.4 * pow(((EKinIso)/(pow(10,54))),(-25./24.)) * pow(((lorentzFactor)/(300)),(4./3.)) * pow(((T90)/(10)),(9./8.)) * pow(((nEx)/(1)),(-11./24.));
      
      // Calculate proton-to-pion conversion factor (not as obvious in afterglow stage
      fPiAG = 0.2;
      //fPiAG = 0.2 * pow(((EKinIso)/(pow(10,54))),(33./24.)) * pow(T90/10,(-9./8.)) * pow(nEx/1,(9./8.));
      //cout << "fPi = " << fPiAG << endl;

      // convert to GeV
      energyCGamma/=1e+9;
      // break energies
      energyBreakNeutrinoAG = 0.015 * shockLorentzFactor * 1/(1+redshift) * pow(energyCGamma,-1); // GeV
      //cout << "Ebreak = " << energyBreakNeutrinoAG << endl;
      // max shock neutrino energy
      maxShockNuE = fPiAG/(4*(1+redshift)) * eE * shockLorentzFactor * magFieldShock * shockRadius;
      
      //cout << "Emax = " << maxShockNuE << endl;      

      shockPrefactor = 1.56e-6 * fPiAG/0.2 * fluence/(1e-6) * pow(((T90)/(10)),-1); //in GeV cm^-2 s^-1

      // Afterglow fluence
      afterGlowFlux1[passed] = new TF1("afterGlowFlux1", "[0] * (x/[1])", minAllowedNuE, energyBreakNeutrinoAG);
      // params and conditionals
      afterGlowFlux1[passed]->SetParameter(0,shockPrefactor); afterGlowFlux1[passed]->SetParName(0,"shockPrefactor");
      afterGlowFlux1[passed]->SetParameter(1,energyBreakNeutrinoAG); afterGlowFlux1[passed]->SetParName(1,"energyBreakNeutrinoAG");
      afterGlowFlux1[passed]->SetLineColor(color);
      afterGlowFlux1[passed]->SetLineWidth(3);
      afterGlowFlux1[passed]->SetLineStyle(10);
      //afterGlowFlux1[passed]->Draw("same");

      cout << "______Afterglow______ " << endl;
      cout << "fPiAG = " << fPiAG << endl;
      cout << "prefactor = " << shockPrefactor << endl;
      cout << "energy break = " << energyBreakNeutrinoAG << endl;
      cout << "maxShockNuE = " << maxShockNuE << endl;
      
      // normal case
      if(maxShockNuE > energyBreakNeutrinoAG)
        {
          afterGlowFlux2[passed] = new TF1("afterGlowFlux2", "[0] * pow((x/[1]),0.5)", energyBreakNeutrinoAG, maxShockNuE);
          // params and conditionals
          afterGlowFlux2[passed]->SetParameter(0,shockPrefactor); afterGlowFlux2[passed]->SetParName(0,"shockPrefactor");
          afterGlowFlux2[passed]->SetParameter(1,energyBreakNeutrinoAG); afterGlowFlux2[passed]->SetParName(1,"energyBreakNeutrinoAG");
        }
      // sometimes, maxShock < energyBreakNeutrino, RICE paper sees same thing.
      // thus anything > eneryBreakNeutrinoAG will be zero. just set maxShockNuE = energyBreakNeutrinoAG+energyBreakNeutrinoAG*0.01
      else
        {
          maxShockNuE = energyBreakNeutrinoAG+energyBreakNeutrinoAG*0.01;
          afterGlowFlux2[passed] = new TF1("afterGlowFlux2", "x*0", energyBreakNeutrinoAG, maxShockNuE);
        }
      afterGlowFlux2[passed]->SetLineColor(color);
      afterGlowFlux2[passed]->SetLineWidth(3);
      afterGlowFlux2[passed]->SetLineStyle(10);
      //afterGlowFlux2[passed]->Draw("same");

      // blank spectrum from theoretical max to max on chart for plotting
      afterGlowFlux3[passed] = new TF1("afterGlowFlux2", "x*0", maxShockNuE, maxAllowedNuE);

      const int stitchedPtsShock = 3*1000;
      // Equals steps from minAllowedNuE to energyBreakNeutrino in log scale
      
      // Broken spec 1: min to break
      double logEbShock = TMath::Log10(energyBreakNeutrinoAG);
      double divisionStepShock1 = (logEbShock - logMin)/double(pts);
      double currentEnergyShock[stitchedPtsShock] = {};
      double currentFluxShock[stitchedPtsShock] = {};

      // Broken spec 2: break to max shock
      double logMaxShock = TMath::Log10(maxShockNuE);
      double divisionStepShock2 = (logMaxShock - logEbShock)/double(pts);
      
      histEnergyCutoff->Fill(logMaxShock);
      
      // Broken spec 3: max shock to actual max
      double divisionStepShock3 = (logMax - logMaxShock)/double(pts);

      // For each spectrum _piece_:
      for(int l = 0; l < pts; l++)
	{
          //Stitch broken spectrum pieces together
	  currentEnergyShock[l] = pow(10,logMin + divisionStepShock1*l); currentFluxShock[l] = afterGlowFlux1[passed]->Eval(currentEnergyShock[l]);
	  currentEnergyShock[pts+l] = pow(10,logEbShock + divisionStepShock2*l); currentFluxShock[pts+l] = afterGlowFlux2[passed]->Eval(currentEnergyShock[pts+l]);
          currentEnergyShock[2*pts+l] = pow(10,logMaxShock + divisionStepShock3*l); currentFluxShock[2*pts+l] = afterGlowFlux3[passed]->Eval(currentEnergyShock[2*pts+l]);
          //cout << currentEnergyShock[l] << endl;
          //cout << currentFluxShock[l] << endl;
	}

      stitchedSpectrumShock[passed] = new TGraph(stitchedPtsShock,currentEnergyShock,currentFluxShock);
      stitchedSpectrumShock[passed]->SetLineWidth(3);
      stitchedSpectrumShock[passed]->SetLineStyle(1);
      
      // Now re-evaluate the stitched spectrum so we get equal spacing in our global energy range
      // Then we can sum stitched spectra
      for(int k = 0; k < newPts; k++)
	{
	  stitchedEnergyShock[k] = pow(10,logMin + divisionStepGlobal*k);
          stitchedEnergyShockLog[k] = logMin + divisionStepGlobal*k;
          stitchedFluxShock[k] = stitchedSpectrumShock[passed]->Eval(stitchedEnergyShock[k],0,"");
          cout << "index = " << k << ", energy = " << stitchedEnergyShock[k] << ", flux = " << stitchedFluxShock[k] << endl;

          quasiDiffuseFluxShock[k]+=stitchedFluxShock[k];
        }

      spectrumShock[passed] = new TGraph(newPts, stitchedEnergyShock, stitchedFluxShock);

      grQuasiDiffuseShock = new TGraph(newPts, stitchedEnergyShock,quasiDiffuseFluxShock);
      grQuasiDiffuseShockLog = new TGraph(newPts, stitchedEnergyShockLog,quasiDiffuseFluxShock);
      
      //// Calculations done and spectra produced, let's summarize:
      cout << "_____Calculation summary_____" << endl;
      cout << setw(10) << left <<  "name =" << *name << endl;
      cout << setw(10) << left <<  "dec =" << dec << endl;
      cout << setw(10) << left <<  "passed # = " << passed << endl;
      cout << setw(10) << left <<  "line color = " << color << endl;
      cout << setw(10) << left <<  "z =" << redshift << endl;
      cout << setw(10) << left <<  "T90 =" << T90 << " s" << endl;
      cout << setw(10) << left <<  "F =" << fluence << " erg cm^-2" << endl;
      cout << setw(10) << left <<  "dL =" << dL << " Mpc" << endl;
      cout << setw(10) << left <<  "Egiso =" << EGammaIso << " erg" << endl;
      cout << setw(10) << left <<  "Lgiso =" << LGammaIso << " erg s^-1" << endl;
      cout << setw(10) << left <<  "Egb =" << energyBreakGamma << " keV" << endl;
      cout << setw(10) << left <<  "Enub =" << energyBreakNeutrino << " GeV" << endl;
      cout << setw(10) << left <<  "B =" << magField << " G" << endl;
      cout << setw(10) << left <<  "EnubS =" << energyBreakNeutrinoSynchrotron << " GeV" << endl;
      cout << endl;

      passed++;
    }

  cout << "Integral = " << grQuasiDiffuseLog->Integral(625,999); // indices of 1e9 -> 1e12
  cout << "Integral ag = " << grQuasiDiffuseShockLog->Integral(625,999); // indices of 1e9 -> 1e12
  cout << "ratio = " << grQuasiDiffuseShockLog->Integral(625,999)/grQuasiDiffuseLog->Integral(625,999);
  // Prompt
  // Quasi diffuse graph from summed stitched graphs
  
  //grQuasiDiffuse->SetLineColor(2);
  //grQuasiDiffuse->SetLineWidth(9);
  //grQuasiDiffuse->Draw("same");

  // Afterglow
  //grQuasiDiffuseShock->SetLineColor(1);
  //grQuasiDiffuseShock->SetLineWidth(9);
  //grQuasiDiffuseShock->Draw("same");

  // Logs
  TCanvas *c2 = new TCanvas();
  //c2->SetLogy();
  
  grQuasiDiffuseShockLog->SetLineColor(1);
  grQuasiDiffuseShockLog->SetLineWidth(9);
  grQuasiDiffuseShockLog->GetYaxis()->SetRangeUser(0.00000000006,0.00025);
  grQuasiDiffuseShockLog->GetXaxis()->SetRangeUser(9.,12.);
  grQuasiDiffuseShockLog->GetXaxis()->SetTitle("Log(Energy [GeV]))");
  grQuasiDiffuseShockLog->GetYaxis()->SetTitle("E^{2} dF/dE (GeV cm^{-2})");
  grQuasiDiffuseShockLog->SetTitle("");
  grQuasiDiffuseShockLog->Draw();
  
  grQuasiDiffuseLog->SetLineColor(2);
  grQuasiDiffuseLog->SetLineWidth(9);
  grQuasiDiffuseLog->Draw("same");

  //TCanvas *c2 = new TCanvas();
  //histEnergyCutoff->GetXaxis()->SetTitle("Log(E [GeV])");
  //histEnergyCutoff->GetYaxis()->SetTitle("Counts");
  //histEnergyCutoff->Draw();
  
  return;
  
}
