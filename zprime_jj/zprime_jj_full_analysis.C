#include "pmatmatch.h"
#include "butil.h"

void zprime_jj_full_analysis_H(bool useArgs, char const *args) {
  
  ///////////////////////////
  // GENERAL FILE SETTINGS //
  ///////////////////////////

  TString directory = "data/";

  // ntuple file name
  TString fileName = "zprime_jj_ntup.root";
  auto recoTreeName = "MyLCTuple;1";
  auto truthTreeName = "JetHistogramGenJetTuple;1";
  
  TString saveDir = "zprime_jj_small";

  //////////////////////////////
  // ANALYSIS CONFIG SETTINGS //
  //////////////////////////////

  // Max pT a reconstructed particle can be from truth pT
  float maxERecoErr = 0.3;

  // Max R a reconstructed particle can be
  float maxRRecoErr = 0.1;

  // delta R that any semi-major reconstructed particles can be from the most major one
  float coneAcceptR = 0.13;

  // minimum weight that a match can be
  float matchWeightCutoff = 0.8;
  
  // If true, limits the reconstructed events used to those which only had one reconstructed particle
  bool limit1Reco = false;

  // Number of particles
  int maxN = 40000;

  // Low pT range
  int lowE = 0;
  // High pT range
  int highE = 5000;

  // Amount that the graphs show below/above the min/max energy
  int lowEErr = 200;
  int highEErr = 2000;

  // number of bins for res fits
  int fitDivs = 20;

  // number of bins for resolution plots
  int resBins = 100;

  // number of bins for showing particle distributions
  int partBins = 40;

  // Delta R testing values
  // double deltaRVals[10] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1};

  /////////////////////////////
  // GRAPH CREATION SETTINGS //
  /////////////////////////////

  bool eH = false;              // Efficiency graphs are made (Eff_PA, Eff_pt)
  bool pH = false;              // Phase graphs contain correlation of energy, polar angle, and energy correction factor (Ephase, K2Corr, K2Corr sub graphs, K2Corr2D, KEtrue)
  bool rH = false;              // Resolution graphs of delta R and energy (deltaR, Eres)
  bool rFH = false;             // Resolution histograms that have gauss fits applied which are used to create some phase plots, mainly used for debugging (ERE_*, ERT_*)
  bool ESH = false;             // Energy, polar angle, azimuth, and p_T for truth, linked, and reconstructed particles (allTruth_*, linkedTruth_*, allReco_*)
  bool ECH = false;             // Energy, polar angle, azimuth, and p_T for truth, linked, and reconstructed particles, but compiled into one graph for each attribute (All_E, All_PA, All_A)
  bool ECF = false;             // Energy correction fit phase plot (KEres_err_0, KTres_err_0)
  bool ECA = false;             // Energy correction applied to reconstructed energies (seperate variable because this causes another loop to run) (CorrTermE, CorrTermT, CorrE, CorrERes)
  
  ////////////////////////////

  if(useArgs) {
     eH = charToBool(args[0]);
     pH = charToBool(args[1]);
     rH = charToBool(args[2]);
     rFH = charToBool(args[3]);
     ESH = charToBool(args[4]);
     ECH = charToBool(args[5]);
     ECF = charToBool(args[6]);
     ECA = charToBool(args[7]);
  }

  TString fileDir = directory + fileName;
  //TFile *myFile = new TFile("ntuple.root");
  TFile *myFile = new TFile(fileDir);
  TTree *myRecoTree = (TTree*)myFile->Get(recoTreeName);
  TTree *myTruthTree = (TTree*)myFile->Get(truthTreeName);

  //opens the file to be read
  auto openFile = TFile::Open(fileDir);
    
  TCanvas *c = new TCanvas();
  
  // Histograms

  TH1F *photonTruth_E_0 = new TH1F("pt_E_0", "All Truth p_{T}", 20, lowE - lowEErr, highE + highEErr); //All true energies
  TH1F *recoPhotonTruth_E_0 = new TH1F("prt_E_0", "Reco True p_{T}", 20, lowE - lowEErr, highE + highEErr); //True energies of reconstructed photons
  TH1F *photonReco_E_0 = new TH1F("pr_E_0", "Reco p_{T}", 20, lowE - lowEErr, highE + highEErr); //Reconstructed energies of linked photons

  TH1F *photonTruth_pt_0 = new TH1F("pt_pt_0", "All Truth p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //All true transverse momenta
  TH1F *recoPhotonTruth_pt_0 = new TH1F("prt_pt_0", "Reco True p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //True transverse momentum of reconstructed photons
  TH1F *photonReco_pt_0 = new TH1F("pr_pt_0", "Reco p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //Reconstructed transverse momentum of linked photons

  TH1F *photonTruth_PA_0 = new TH1F("pt_PA_0", "All Truth Polar Angle", 20, 0, 3.2); //All true polar anglesx
  TH1F *recoPhotonTruth_PA_0 = new TH1F("prt_PA_0", "Reco True Polar Angle", 20, 0, 3.2); //True polar angle of reconstructed photons
  TH1F *photonReco_PA_0 = new TH1F("pr_PA_0", "Reco Polar Angle", 20, 0, 3.2); //Reconstructed polar angle of linked photons

  TH1F *photonTruth_A_0 = new TH1F("pt_A_0", "All Truth Azimuth", 20, -1.6, 1.6); //All true azimuths
  TH1F *recoPhotonTruth_A_0 = new TH1F("prt_A_0", "Reco True Azimuth", 20, -1.6, 1.6); //True azimuth of reconstructed photons
  TH1F *photonReco_A_0 = new TH1F("pr_A_0", "Reco Azimuth", 20, -1.6, 1.6); //Reconstructed azimuth of linked photons

  TH1F *deltaR_0 = new TH1F("dR_0", "Delta R", 20, 0 , 0.002); // delta R between mc truth and best match reco
  TH1F *Eres_0 = new TH1F("Er_0", "Energy resolution", resBins, -0.1, 0.1); //Resolution of reconstructed energy
  TH2F *Ephase_0 = new TH2F("Epha_0", "True vs Reco energy phase plot", 100, -lowEErr, highE + highEErr, 100, lowE - lowEErr, highE + highEErr); //Phase plot \
of true energy w.r.t reconstructed energy
  TH3F *K2Corr_0 = new TH3F("K2c_0", "Correction w.r.t. polar angle and energy", 20, lowE, highE, 20, 0, 3.2, 20, 0.94, 0.98); //Correction term plot w.r.t. polar angle and energy
  TGraph2D *K2Corr2D_0 = new TGraph2D(); //2D surface of the above scatter plot
  TH2F *KEtrue_0 = new TH2F("KEt_0", "Correction factor as function of energy", 100, lowE, highE, 100, 0.94, 1.00);
  TH1F *CorrE_Bins_0 = new TH1F("CEB_0", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr); //Corrected photon energies
  TH1F *CorrERes_Bins_0 = new TH1F("CERB_0", "Corrected energy resolution -- Bins", resBins, -0.1, 0.1); //Resolution of corrected photon energies
  TH1F *CorrE_Fit_0 = new TH1F("CEF_0", "Corrected energy -- Fit", 20, lowE - lowEErr, highE + highEErr); //Corrected photon energies
  TH1F *CorrERes_Fit_0 = new TH1F("CERF_0", "Corrected energy resolution -- Fit", resBins, -0.1, 0.1); //Resolution of corrected photon energies
  
  THStack *hs;
  TF2 *Corr2DFitFunc_0;
  TF2 *Corr2DFitFunc_1;
  TF2 *Corr2DFitFunc_2;

  // Root file reader stuff
  TTreeReader myRecoReader(recoTreeName, openFile);
  TTreeReader myTruthReader(truthTreeName, openFile);
  
  TTreeReaderArray<int> njet_RRA(myRecoReader, "njet");
  TTreeReaderArray<Float_t> jmox_RRA(myRecoReader, "jmox");
  TTreeReaderArray<Float_t> jmoy_RRA(myRecoReader, "jmoy");
  TTreeReaderArray<Float_t> jmoz_RRA(myRecoReader, "jmoz");
  TTreeReaderArray<Float_t> jene_RRA(myRecoReader, "jene");
  
  TTreeReaderArray<int> njet_TRA(myTruthReader, "njet");
  TTreeReaderArray<Float_t> jmox_TRA(myTruthReader, "jmox");
  TTreeReaderArray<Float_t> jmoy_TRA(myTruthReader, "jmoy");
  TTreeReaderArray<Float_t> jmoz_TRA(myTruthReader, "jmoz");
  TTreeReaderArray<Float_t> jene_TRA(myTruthReader, "jene");
  
  int r2fTemp, r2tTemp;
  TVector3 mcmoTemp, rcmoTemp, coneTemp, hcalTemp;
  Float_t r2wTemp, dPATemp, dATemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp, weightTemp;
  Float_t r2wMax = 0;
  int r2wMax_Index = -1;
  int arrCount = 0;
  int evtCount = 0;
  
  double corrEAvgs[fitDivs];
  double corrTAvgs[fitDivs];
  double corrESTD[fitDivs];
  double corrTSTD[fitDivs];
  
  TH1F* sigFitEHists[fitDivs];
  TH1F* sigFitTHists[fitDivs];
  TH1F* corrSigFitEHists[fitDivs];
  TH1F* corrSigFitTHists[fitDivs];
  
  TH1F* recoTypesE[5];
  TH1F* recoTypesT[5];
  TH1F* mcRecoTypesE[5];
  TH1F* mcRecoTypesT[5];

  double corrEBins[fitDivs];
  int corrEBinN[fitDivs];
  double corrTBins[fitDivs];
  int corrTBinN[fitDivs];
  double recoEnergyArr[maxN];
  double trueEnergyArr[maxN];
  double recoPAArr[maxN];

  cout << "Test0" << endl;
  
  for(int i = 0; i < fitDivs; i++) {
    corrEBins[i] = 0;
    corrEBinN[i] = 0;
    corrTBins[i] = 0;
    corrTBinN[i] = 0;
  }

  if(ECF || rFH) {
    for(int i = 0; i < fitDivs; i++) {
      char s[4];
      char outE[100];
      char outT[100];
      char outCorrE[100];
      char outCorrT[100];
      
      sprintf(s, "%d", i);
      
      strcpy(outE, "sfEH");
      strcat(outE, s);
      strcat(outE, "_0");
      strcpy(outT, "sfTH");
      strcat(outT, s);
      strcat(outT, "_0");

      strcpy(outCorrE, "csfEH");
      strcat(outCorrE, s);
      strcat(outCorrE, "_0");
      strcpy(outCorrT, "csfTH");
      strcat(outCorrT, s);
      strcat(outCorrT, "_0");
      
      sigFitEHists[i] = new TH1F(outE, "Energy resolution", 40, -0.1, 0.1);
      sigFitTHists[i] = new TH1F(outT, "Energy resolution", 40, -0.1, 0.1);
      corrSigFitEHists[i] = new TH1F(outCorrE, "Energy resolution", 40, -0.1, 0.1);
      corrSigFitTHists[i] = new TH1F(outCorrT, "Energy resolution", 40, -0.1, 0.1);
    }
  }

  cout << "Test1" << endl;

  // run 0
  while(myRecoReader.Next()) {
    myTruthReader.Next();

    // cout << "Test1.01" << endl;

    int njet_R = njet_RRA.At(0);
    int njet_T = njet_TRA.At(0);

    // cout << "Test1.02" << endl;

    TVector3 recoVec3[njet_R];
    TVector3 truthVec3[njet_T];

    // cout << "Test1.04" << endl;

    for(int i = 0; i < njet_R; i++) {
      recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
    }
    for(int i = 0; i < njet_T; i++) {
      truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
    }

    // cout << "Test1.06" << endl;

    int** r = build2DArrI(2, njet_T);

    // cout << "Test1.07" << endl;
    
    r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);

    // cout << "Test1.08" << endl;
    
    // int r[][] = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
    
    int matchArrR[njet_T];
    int matchArrT[njet_T];
    for(int i = 0; i < njet_T; i++) {
      matchArrR[i] = r[0][i];
      matchArrT[i] = r[1][i];
    }

    if(evtCount == 99) {
      for(int i = 0; i < njet_T; i++) {
	cout << "R: " << matchArrR[i] << " T: " << matchArrT[i] << endl;
      }
    }
    
    // cout << "Test 1.5" << endl;

    cout << "njet_T loop: " << njet_T << " matchCount: " << r[2][0] << endl;
    for(int i = 0; i < njet_T; i++) {
      
      bool matched = contains(matchArrT, njet_T, i);
      int recoMatch;

      mcmoTemp = truthVec3[i];
      
      if(matched) {
	recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	
        rcmoTemp = recoVec3[recoMatch];

        photonReco_E_0->Fill(rcmoTemp.Mag());
        photonReco_pt_0->Fill(rcmoTemp.Perp());
        photonReco_PA_0->Fill(rcmoTemp.Theta());
        photonReco_A_0->Fill(rcmoTemp.Phi());

        deltaR_0->Fill(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ));
        Ephase_0->Fill(rcmoTemp.Mag(), mcmoTemp.Mag());
        K2Corr_0->Fill(rcmoTemp.Mag(), rcmoTemp.Theta(), mcmoTemp.Mag() / rcmoTemp.Mag());
        K2Corr2D_0->SetPoint(arrCount, rcmoTemp.Mag(), rcmoTemp.Theta(), mcmoTemp.Mag() / rcmoTemp.Mag());

        KEtrue_0->Fill(mcmoTemp.Mag(), mcmoTemp.Mag() / rcmoTemp.Mag());

        if(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ) < maxRRecoErr) {

          int RecoEnergyBin = (int) floor((rcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
          int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
          int TrueEnergyBin = (int) floor((mcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
          int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));

          while(RecoEnergyBin < 0) { RecoEnergyBin++; }
          while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
          while(RecoPABin < 0) { RecoPABin++; }
          while(RecoPABin >= fitDivs) { RecoPABin--; }

          if(ECF || rFH) {

            sigFitEHists[ RecoEnergyBin ]->Fill((rcmoTemp.Mag() - mcmoTemp.Mag()) / mcmoTemp.Mag());
            sigFitTHists[ RecoPABin ]->Fill((rcmoTemp.Mag() - mcmoTemp.Mag()) / mcmoTemp.Mag());
          }

          if( (rcmoTemp - mcmoTemp).Mag() / mcmoTemp.Mag() < maxERecoErr) {

            if(ECA) {

              corrEBins[ RecoEnergyBin ] += (mcmoTemp.Mag() / rcmoTemp.Mag());
              corrEBinN[ RecoEnergyBin ]++;
              corrTBins[ RecoPABin ] += (mcmoTemp.Mag() / rcmoTemp.Mag());
              corrTBinN[ RecoPABin ]++;
            }

          }

	  if(ESH || ECH || eH) {
	    recoPhotonTruth_E_0->Fill(mcmoTemp.Mag());
	    recoPhotonTruth_pt_0->Fill(mcmoTemp.Perp());
	    recoPhotonTruth_PA_0->Fill(mcmoTemp.Theta());
	    recoPhotonTruth_A_0->Fill(mcmoTemp.Phi());
          }

        }
	
	Eres_0->Fill((rcmoTemp.Mag() - mcmoTemp.Mag()) / mcmoTemp.Mag());
	
	recoEnergyArr[arrCount] = rcmoTemp.Mag();
	trueEnergyArr[arrCount] = mcmoTemp.Mag();
	recoPAArr[arrCount] = rcmoTemp.Theta();
	arrCount++;
      }

      photonTruth_E_0->Fill(mcmoTemp.Mag());
      photonTruth_pt_0->Fill(mcmoTemp.Perp());
      photonTruth_PA_0->Fill(mcmoTemp.Theta());
      photonTruth_A_0->Fill(mcmoTemp.Phi());
    }
    evtCount++;
  }

  cout << evtCount << endl;

  cout << "Test2" << endl;
  
  
  if(ECA) {

    cout << "Test2.01" << endl;

    // restart run
    myRecoReader.Restart();
    cout << "Test2.02" << endl;
    myTruthReader.Restart();
    cout << "Test2.03" << endl;

    myRecoReader.Next();
    myTruthReader.Next();
    
    int njet_R = njet_RRA.At(0);
    cout << "Test2.04" << endl;
    int njet_T = njet_TRA.At(0);
    cout << "Test2.05" << endl;
    
    TVector3 recoVec3[njet_R];
    cout << "Test2.06" << endl;
    TVector3 truthVec3[njet_T];
    cout << "Test2.07" << endl;
    
    for(int i = 0; i < njet_R; i++) {
      cout << "Test2.08" << endl;
      recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
      cout << "Test2.09" << endl;
    }
    for(int i = 0; i < njet_T; i++) {
      cout << "Test2.010" << endl;
      truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
      cout << "Test2.011" << endl;
    }
    
    int** r = build2DArrI(2, njet_T);
    cout << "Test2.012" << endl;
    r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
    cout << "Test2.013" << endl;

    int matchArrR[njet_T];
    int matchArrT[njet_T];
    for(int i = 0; i < njet_T; i++) {
      matchArrR[i] = r[0][i];
      matchArrT[i] = r[1][i];
    }

    for(int i = 0; i < njet_T; i++) {

      bool matched = contains(matchArrT, njet_T, i);
      int recoMatch;

      mcmoTemp = truthVec3[i];
      
      if(matched) {
        recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];

        rcmoTemp = recoVec3[recoMatch];
	
	int RecoEnergyBin = (int) floor((rcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
	int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	int TrueEnergyBin = (int) floor((mcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
	int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	
	while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	while(RecoPABin < 0) { RecoPABin++; }
	while(RecoPABin >= fitDivs) { RecoPABin--; }
	
	if( (rcmoTemp - mcmoTemp).Mag() / mcmoTemp.Mag() < maxERecoErr) {
	  if(ECA) {
	    if( pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2) < 1 ) {
	      corrESTD[ RecoEnergyBin ] += pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2);
	    } else {
	      cout << "BAD: " << pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2) << endl;
	    }
	    if( pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2) < 1 ) {
	      corrTSTD[ RecoPABin ] += pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2);
	    } else {
	      cout << "BAD: " << pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2) << endl;
	    }	      
	  }
	}
      }
    }

    cout << "Test2.5" << endl;

    // run 1
    while(myRecoReader.Next()) {
      myRecoReader.Next();
      
      int njet_R = njet_RRA.At(0);
      int njet_T = njet_TRA.At(0);

      TVector3 recoVec3[njet_R];
      TVector3 truthVec3[njet_T];

      for(int i = 0; i < njet_R; i++) {
	recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
      }
      for(int i = 0; i < njet_T; i++) {
	truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
      }
      
      int** r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
      
      int matchArrR[njet_T];
      int matchArrT[njet_T];
      for(int i = 0; i < njet_T; i++) {
	matchArrR[i] = r[0][i];
	matchArrT[i] = r[1][i];
      }

      for(int i = 0; i < njet_T; i++) {

	bool matched = contains(matchArrT, njet_T, i);
	int recoMatch;

	mcmoTemp = truthVec3[i];

	if(matched) {
	  recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];

	  rcmoTemp = recoVec3[recoMatch];

	  int RecoEnergyBin = (int) floor((rcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
	  int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	  int TrueEnergyBin = (int) floor((mcmoTemp.Mag() - lowE) * fitDivs / (highE - lowE));
	  int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));

	  while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	  while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	  while(RecoPABin < 0) { RecoPABin++; }
	  while(RecoPABin >= fitDivs) { RecoPABin--; }

	  if( (rcmoTemp - mcmoTemp).Mag() / mcmoTemp.Mag() < maxERecoErr) {
	    if(ECA) {
	      if( pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2) < 1 ) {
		corrESTD[ RecoEnergyBin ] += pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2);
	      } else {
		cout << "BAD: " << pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrEBins[ RecoEnergyBin ] / corrEBinN[ RecoEnergyBin ]), 2) << endl;
	      }
	      if( pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2) < 1 ) {
		corrTSTD[ RecoPABin ] += pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2);
	      } else {
		cout << "BAD: " << pow((mcmoTemp.Mag() / rcmoTemp.Mag()) - (corrTBins[ RecoPABin ] / corrTBinN[ RecoPABin ]), 2) << endl;
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "Test3" << endl;

  if(ECA) {

    double recoEnergyArr_Bins[arrCount];
    double recoEnergyArr_Fit[arrCount];

    // Arrays for plotting TGraphErrors
    Double_t xE[fitDivs];
    Double_t exE[fitDivs];
    Double_t eyE[fitDivs];

    Double_t xT[fitDivs];
    Double_t exT[fitDivs];
    Double_t eyT[fitDivs];

    double stepE = (highE - lowE) / fitDivs;
    double currE = lowE + (stepE / 2);

    double stepT = (TMath::Pi() * 8 ) / (9 * fitDivs);
    double currT = (TMath::Pi() / 18) + (stepT / 2);

    // Prints out the bins used to find the weights
    for(int i = 0; i < fitDivs; i++) {
      cout << corrEBins[i] << " " << corrEBinN[i] << " " << corrESTD[i] << endl;
      cout << corrTBins[i] << " " << corrTBinN[i] << " " << corrTSTD[i] << endl;
    }

    // Set up arrays for plotting TGraphErrors for the correction terms (for energy and polar angle)
    for(int i = 0; i < fitDivs; i++) {
      if(corrEBinN[i] != 0) {
	corrEBins[i] = corrEBins[i] / corrEBinN[i];
	corrESTD[i] = corrESTD[i] / corrEBinN[i];
      }
      if(corrTBinN[i] != 0) {
	corrTBins[i] = corrTBins[i] / corrTBinN[i];
	corrTSTD[i] = corrTSTD[i] / corrTBinN[i];
      }
      
      xE[i] = currE;
      exE[i] = stepE / 2;
      eyE[i] = 0;

      xT[i] = currT;
      exT[i] = stepT / 2;
      eyT[i] = 0;
      
      currE += stepE;
      currT += stepT;
    }

    // 1D fit functions for the correction factors
    TF1 *CorrTermTFit_0 = new TF1("CorrTermTFit_0", "([0]*cos(( x - [1]) * [2]) + [3])", (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    TF1 *CorrTermTFit_1 = new TF1("CorrTermTFit_1", "([0]*cos(( x - (3.1415926535897932384626433832795028 / 2) ) * [1]) + [2])", atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    TF1 *CorrTermTFit_2 = new TF1("CorrTermTFit_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - x ) * [2]) + [3])", TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    //2D fit functions for the correction factors 
    // TF2 *Corr2DFitFunc_0 = new TF2("Corr2DFitFunc_0", "([0]*cos(( y - [1]) * [2]) + [3]) + sqrt((([4]/sqrt(x)) * ([4]/sqrt(x))) + (([5] / x) * ([5] / x)) + ([6] * [6]))", lowE, highE, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    // TF2 *Corr2DFitFunc_1 = new TF2("Corr2DFitFunc_1", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) + [2]) + sqrt((([3]/sqrt(x)) * ([3]/sqrt(x))) + (([4] / x) * ([4] / x)) + ([5] * [5]))", lowE, highE, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    // TF2 *Corr2DFitFunc_2 = new TF2("Corr2DFitFunc_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) + [3]) + sqrt((([4]/sqrt(x)) * ([4]/sqrt(x))) + (([5] / x) * ([5] / x)) + ([6] * [6]))", lowE, highE, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_0 = new TF2("Corr2DFitFunc_0", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_1 = new TF2("Corr2DFitFunc_1", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(TMath::E(), -x*[3]) + [4]", lowE + 100, highE - 100, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_2 = new TF2("Corr2DFitFunc_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    // TF1 *spQuadFit = new TF1("spQuadFit", "sqrt(([0] / sqrt(x))**2 + ([1] / x)**2 + ([2])**2)", lowE, highE);
    TF1 *spQuadFit = new TF1("spQuadFit", "[0]*pow(TMath::E(), -x * [1]) + [2]", lowE + 100, highE - 100);
    // TF1 *spQuadFit = new TF1("spQuadFit", "[0] + ([1]*x) + ([2]*x**2) + ([3]*x**3) + ([4]*x**4) + ([5]*x**5)");
    
    // Set initial 1D fit parameters
    CorrTermTFit_0->SetParameters(0.07, 0.2, 1.6, 0.89);
    CorrTermTFit_1->SetParameters(-0.0065, 3.4, 0.92);
    CorrTermTFit_2->SetParameters(0.07, 0.2, 1.6, 0.89);

    // Set initial 2D fit parameters
    Corr2DFitFunc_0->SetParameters(0.07, 0.2, 1.6, 0.0157202, 0.000649653, 0.95403);
    Corr2DFitFunc_1->SetParameters(-0.0065, 3.4, 0.0157202, 0.000649653, 0.95403);
    Corr2DFitFunc_2->SetParameters(0.07, 0.2, 1.6, 0.0157202, 0.000649653, 0.95403);

    spQuadFit->SetParameters(0.016, 0.001, 0.9552);
    // spQuadFit->SetParameters(0.976, -0.0000095, 0.00000000104, 0.0, 0.0, 0.0);

    // Fit 2D graph
    K2Corr2D_0->Fit("Corr2DFitFunc_0", "RW0");
    K2Corr2D_0->Fit("Corr2DFitFunc_1", "RW0");
    K2Corr2D_0->Fit("Corr2DFitFunc_2", "RW0");

    // Apply energy corrections
    for(int i = 0; i < arrCount; i++) {
      
      // Energy correction for bins
      int energyBin = (int) floor((recoEnergyArr[i] - lowE) * fitDivs / (highE - lowE));
      int PABin = (int) floor((recoPAArr[i] - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
      while(energyBin < 0) { energyBin++; }
      while(energyBin >= fitDivs) { energyBin--; }
      while(PABin < 0) { PABin++; }
      while(PABin >= fitDivs) { PABin--; }

      // Apply bin corrections
      recoEnergyArr_Bins[i] = recoEnergyArr[i] * ( ( corrEBins[ energyBin ] + corrTBins[ PABin ] ) / 2 );
      CorrE_Bins_0->Fill(recoEnergyArr_Bins[i]);
      CorrERes_Bins_0->Fill((recoEnergyArr_Bins[i] - trueEnergyArr[i]) / trueEnergyArr[i]);

      // Apply fit correction
      if( recoPAArr[i] < atan(1500.0/2210.0) ) {
	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_0->Eval(recoEnergyArr[i], recoPAArr[i]);
      } else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_1->Eval(recoEnergyArr[i], recoPAArr[i]);
      } else {
	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_2->Eval(recoEnergyArr[i], recoPAArr[i]);
      }
      CorrE_Fit_0->Fill(recoEnergyArr_Fit[i]);
      CorrERes_Fit_0->Fill((recoEnergyArr_Fit[i] - trueEnergyArr[i]) / trueEnergyArr[i]);
    }
    
    // Draw and save graphs    
    TGraphErrors *CorrTermE_0 = new TGraphErrors(fitDivs, xE, corrEBins, exE, corrESTD);
    
    CorrTermE_0->Draw();
    CorrTermE_0->SetTitle("Correction Terms w.r.t. Energy");
    CorrTermE_0->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    CorrTermE_0->GetYaxis()->SetTitle("Correction Factor");

    CorrTermE_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    CorrTermE_0->Fit("spQuadFit", "RW");
    spQuadFit->Draw("SAME");

    c->SaveAs(saveDir + "/phaseHists/CorrTermE.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *CorrTermT_0 = new TGraphErrors(fitDivs, xT, corrTBins, exT, corrTSTD);
    
    CorrTermT_0->Draw();
    CorrTermT_0->SetTitle("Correction Terms w.r.t. Polar Angle");
    CorrTermT_0->GetXaxis()->SetTitle("Polar Angle (Rads)");
    CorrTermT_0->GetYaxis()->SetTitle("Correction Factor");

    CorrTermT_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    CorrTermT_0->Fit("CorrTermTFit_0", "RW");
    CorrTermT_0->Fit("CorrTermTFit_1", "RW");
    CorrTermT_0->Fit("CorrTermTFit_2", "RW");
    
    CorrTermTFit_0->Draw("SAME");
    CorrTermTFit_1->Draw("SAME");
    CorrTermTFit_2->Draw("SAME");
    
    c->SaveAs(saveDir + "/phaseHists/CorrTermT.png");
    c->Close();
    c = new TCanvas();

    CorrE_Bins_0->Draw();
    CorrE_Bins_0->SetTitle("Corrected reconstructed energy -- Bins");
    CorrE_Bins_0->GetXaxis()->SetTitle("Corrected Energy (GeV)");
    CorrE_Bins_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/CorrE_Bins.png");
    c->Close();
    c = new TCanvas();

    CorrERes_Bins_0->Draw();
    CorrERes_Bins_0->SetTitle("Corrected reconstructed energy resolution -- Bins");
    CorrERes_Bins_0->GetXaxis()->SetTitle("Energy resolution (GeV)");
    CorrERes_Bins_0->GetYaxis()->SetTitle("Count");

    CorrERes_Bins_0->Fit("gaus");

    c->SaveAs(saveDir + "/resolutionHists/CorrERes_Bins.png");
    c->Close();
    c = new TCanvas();

    CorrE_Fit_0->Draw();
    CorrE_Fit_0->SetTitle("Corrected reconstructed energy -- Fit");
    CorrE_Fit_0->GetXaxis()->SetTitle("Corrected Energy (GeV)");
    CorrE_Fit_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/CorrE_Fit.png");
    c->Close();
    c = new TCanvas();

    CorrERes_Fit_0->Draw();
    CorrERes_Fit_0->SetTitle("Corrected reconstructed energy resolution -- Fit");
    CorrERes_Fit_0->GetXaxis()->SetTitle("Energy resolution (GeV)");
    CorrERes_Fit_0->GetYaxis()->SetTitle("Count");

    CorrERes_Fit_0->Fit("gaus");

    c->SaveAs(saveDir + "/resolutionHists/CorrERes_Fit.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_0 = new TGraph2D((TH2D*)Corr2DFitFunc_0->GetHistogram());

    C2DFFplot_0->GetYaxis()->SetRange(0, 4);

    C2DFFplot_0->Draw("surf1");
    C2DFFplot_0->SetTitle("Reco Energy Correction Factors 2D Fit -- 0-34.166 deg");
    C2DFFplot_0->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_0->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_0->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_0.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_1 = new TGraph2D((TH2D*)Corr2DFitFunc_1->GetHistogram());

    C2DFFplot_1->GetYaxis()->SetRange(0, 4);

    C2DFFplot_1->Draw("surf1");
    C2DFFplot_1->SetTitle("Reco Energy Correction Factors 2D Fit -- 34.167-145.833 deg");
    C2DFFplot_1->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_1->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_1->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_1.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_2 = new TGraph2D((TH2D*)Corr2DFitFunc_2->GetHistogram());

    C2DFFplot_2->GetYaxis()->SetRange(0, 4);

    C2DFFplot_2->Draw("surf1");
    C2DFFplot_2->SetTitle("Reco Energy Correction Factors 2D Fit -- 145.834-180 deg");
    C2DFFplot_2->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_2->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_2->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_2.png");
    c->Close();
    c = new TCanvas();    
  }

  cout << "Test4" << endl;

  if(ECF || rFH) {

    myRecoReader.Restart();
    myTruthReader.Restart();

    int njet_R = njet_RRA.At(0);
    int njet_T = njet_TRA.At(0);

    TVector3 recoVec3[njet_R];
    TVector3 truthVec3[njet_T];

    for(int i = 0; i < njet_R; i++) {
      recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
    }
    for(int i = 0; i < njet_T; i++) {
      truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
    }

    int** r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);

    int matchArrR[njet_T];
    int matchArrT[njet_T];
    for(int i = 0; i < njet_T; i++) {
      matchArrR[i] = r[0][i];
      matchArrT[i] = r[1][i];
    }
    
    for(int i = 0; i < njet_T; i++) {

      bool matched = contains(matchArrT, njet_T, i);
      int recoMatch;

      mcmoTemp = truthVec3[i];

      if(matched) {
        recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];

        rcmoTemp = recoVec3[recoMatch];

	Double_t tempE;
	// correct energy
	if( rcmoTemp.Theta() < atan(1500.0/2210.0) ) {
	  tempE = rcmoTemp.Mag() * Corr2DFitFunc_0->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	} else if( rcmoTemp.Theta() < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  tempE = rcmoTemp.Mag() * Corr2DFitFunc_1->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	} else {
	  tempE = rcmoTemp.Mag() * Corr2DFitFunc_2->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	}

	int RecoEnergyBin = (int) floor((tempE - lowE) * fitDivs / (highE - lowE));
	int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));

	while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	while(RecoPABin < 0) { RecoPABin++; }
	while(RecoPABin >= fitDivs) { RecoPABin--; }

	// fill in hist array
	corrSigFitEHists[ RecoEnergyBin ]->Fill((tempE - mcmoTemp.Mag()) / mcmoTemp.Mag());
	corrSigFitTHists[ RecoPABin ]->Fill((tempE - mcmoTemp.Mag()) / mcmoTemp.Mag());

      }
    }

    while(myRecoReader.Next()) {
      myTruthReader.Next();
      
      int njet_R = njet_RRA.At(0);
      int njet_T = njet_TRA.At(0);

      TVector3 recoVec3[njet_R];
      TVector3 truthVec3[njet_T];

      for(int i = 0; i < njet_R; i++) {
	recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
      }
      for(int i = 0; i < njet_T; i++) {
	truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
      }
      
      int** r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);

      int matchArrR[njet_T];
      int matchArrT[njet_T];
      for(int i = 0; i < njet_T; i++) {
	matchArrR[i] = r[0][i];
	matchArrT[i] = r[1][i];
      }
      
      for(int i = 0; i < njet_T; i++) {
	
	bool matched = contains(matchArrT, njet_T, i);
	int recoMatch;

	mcmoTemp = truthVec3[i];

	if(matched) {
	  recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];

	  rcmoTemp = recoVec3[recoMatch];

	  Double_t tempE;
	  // correct energy
	  if( rcmoTemp.Theta() < atan(1500.0/2210.0) ) {
	    tempE = rcmoTemp.Mag() * Corr2DFitFunc_0->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	  } else if( rcmoTemp.Theta() < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	    tempE = rcmoTemp.Mag() * Corr2DFitFunc_1->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	  } else {
	    tempE = rcmoTemp.Mag() * Corr2DFitFunc_2->Eval(rcmoTemp.Mag(), rcmoTemp.Theta());
	  }

	  int RecoEnergyBin = (int) floor((tempE - lowE) * fitDivs / (highE - lowE));
	  int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));

	  while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	  while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	  while(RecoPABin < 0) { RecoPABin++; }
	  while(RecoPABin >= fitDivs) { RecoPABin--; }

	  // fill in hist array
	  corrSigFitEHists[ RecoEnergyBin ]->Fill((tempE - mcmoTemp.Mag()) / mcmoTemp.Mag());
	  corrSigFitTHists[ RecoPABin ]->Fill((tempE - mcmoTemp.Mag()) / mcmoTemp.Mag());

	}
      }
    }
    
    Double_t pxE[fitDivs];
    Double_t pyE[fitDivs];
    Double_t exE[fitDivs];
    Double_t eyE[fitDivs];
    
    Double_t pxT[fitDivs];
    Double_t pyT[fitDivs];
    Double_t exT[fitDivs];
    Double_t eyT[fitDivs];

    Double_t pyCE[fitDivs];
    Double_t eyCE[fitDivs];

    Double_t pyCT[fitDivs];
    Double_t eyCT[fitDivs];
      
    double stepE = (highE - lowE) / fitDivs;
    double currE = lowE + (stepE / 2);
      
    double stepT = (TMath::Pi() * 8 ) / (9 * fitDivs);
    double currT = (TMath::Pi() / 18) + (stepT / 2);
      
    for(int i = 0; i < fitDivs; i++) {
      TF1 *iGaussFit = new TF1("iGaussFit", "gaus", -1, 1);
	
      char s[4];
      char outE[100];
      char outT[100];
      char outCE[100];
      char outCT[100];
	
      if(rFH) {
	sprintf(s, "%d", i);
	strcpy(outE, saveDir);
	strcat(outE, "/seperateRes/ERE_");
	strcat(outE, s);
	strcat(outE, ".png");
	strcpy(outT, saveDir);
	strcat(outT, "/seperateRes/ERT_");
	strcat(outT, s);
	strcat(outT, ".png");
	strcpy(outCE, saveDir);
        strcat(outCE, "/seperateRes/CERE_");
        strcat(outCE, s);
        strcat(outCE, ".png");
        strcpy(outCT, saveDir);
        strcat(outCT, "/seperateRes/CERT_");
        strcat(outCT, s);
        strcat(outCT, ".png");
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);
      
      sigFitEHists[i]->Fit("iGaussFit");
	
      if(rFH) {
	sigFitEHists[i]->Draw();
	sigFitEHists[i]->SetTitle("Energy resolution - E");
	sigFitEHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
	sigFitEHists[i]->GetYaxis()->SetTitle("Count");

	c->SaveAs(outE);
	c->Close();
	c = new TCanvas();
      }

      if(ECF) {
	pxE[i] = currE;
	pyE[i] = iGaussFit->GetParameter(2);
	exE[i] = stepE / 2;
	eyE[i] = iGaussFit->GetParError(2);
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);

      sigFitTHists[i]->Fit("iGaussFit");

      if(rFH) {
	sigFitTHists[i]->Draw();
	sigFitTHists[i]->SetTitle("Energy resolution - PA");
	sigFitTHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
	sigFitTHists[i]->GetYaxis()->SetTitle("Count");

	c->SaveAs(outT);
	c->Close();
	c = new TCanvas();
      }

      if(ECF) {
	pxT[i] = currT;
	pyT[i] = iGaussFit->GetParameter(2);
	exT[i] = stepT / 2;
	eyT[i] = iGaussFit->GetParError(2);
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);

      corrSigFitEHists[i]->Fit("iGaussFit");

      if(rFH) {
        corrSigFitEHists[i]->Draw();
        corrSigFitEHists[i]->SetTitle("Corrected Energy resolution - E");
        corrSigFitEHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
        corrSigFitEHists[i]->GetYaxis()->SetTitle("Count");

        c->SaveAs(outCE);
        c->Close();
        c = new TCanvas();
      }

      if(ECF) {
        pxE[i] = currE;
        pyCE[i] = iGaussFit->GetParameter(2);
        exE[i] = stepE / 2;
        eyCE[i] = iGaussFit->GetParError(2);
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);

      corrSigFitTHists[i]->Fit("iGaussFit");

      if(rFH) {
        corrSigFitTHists[i]->Draw();
        corrSigFitTHists[i]->SetTitle("Corrected Energy resolution - PA");
        corrSigFitTHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
        corrSigFitTHists[i]->GetYaxis()->SetTitle("Count");

        c->SaveAs(outCT);
        c->Close();
        c = new TCanvas();
      }

      if(ECF) {
        pxT[i] = currT;
        pyCT[i] = iGaussFit->GetParameter(2);
        exT[i] = stepT / 2;
        eyCT[i] = iGaussFit->GetParError(2);
      }
  
      currE += stepE;
      currT += stepT;
    }

    TF1 *CorrTermTFit_0 = new TF1("CorrTermTFit_0", "([0]*cos(( x - [1]) * [2]) + [3])", (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    TF1 *CorrTermTFit_1 = new TF1("CorrTermTFit_1", "([0]*cos(( x - (3.1415926535897932384626433832795028 / 2) ) * [1]) + [2])", atan(1500.0/2210.0) + 0.05, TMath::Pi() - atan(1500.0/2210.0) - 0.05);
    TF1 *CorrTermTFit_2 = new TF1("CorrTermTFit_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - x ) * [2]) + [3])", TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    TF1 *spQuadFit = new TF1("spQuadFit", "sqrt(([0] / sqrt(x))**2 + ([1] / x)**2 + ([2])**2)");

    // Set initial 1D fit parameters
    CorrTermTFit_0->SetParameters(0.005, 5.3, 0.7, 0.0097);
    CorrTermTFit_1->SetParameters(0.0035, 3.2, 0.0105);
    CorrTermTFit_2->SetParameters(0.005, 5.3, 0.7, 0.0097);

    // Set initial quad fit parameters
    spQuadFit->SetParameters(0.015, 0.001, 0.9552);
    
    TGraphErrors *KEres_err_0 = new TGraphErrors(fitDivs, pxE, pyE, exE, eyE);

    KEres_err_0->Draw();
    KEres_err_0->SetTitle("Error plot of Gauss fit sigma vs true energy");
    KEres_err_0->GetXaxis()->SetTitle("E_{True} (GeV)");
    KEres_err_0->GetYaxis()->SetTitle("Sigma (GeV)");

    KEres_err_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    KEres_err_0->Fit("spQuadFit");

    c->SaveAs(saveDir + "/phaseHists/KEres_err.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *KTres_err_0 = new TGraphErrors(fitDivs, pxT, pyT, exT, eyT);

    KTres_err_0->Draw();
    KTres_err_0->SetTitle("Error plot of Gauss fit sigma vs true polar angle");
    KTres_err_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    KTres_err_0->GetYaxis()->SetTitle("Sigma (GeV)");

    KTres_err_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    KTres_err_0->Fit("CorrTermTFit_0", "RW");
    KTres_err_0->Fit("CorrTermTFit_1", "RW");
    KTres_err_0->Fit("CorrTermTFit_2", "RW");

    CorrTermTFit_0->Draw("SAME");
    CorrTermTFit_1->Draw("SAME");
    CorrTermTFit_2->Draw("SAME");

    c->SaveAs(saveDir + "/phaseHists/KTres_err.png");
    c->Close();
    c = new TCanvas();

    // Set initial 1D fit parameters
    CorrTermTFit_0->SetParameters(0.005, 5.3, 0.7, 0.0097);
    CorrTermTFit_1->SetParameters(0.0035, 3.2, 0.0105);
    CorrTermTFit_2->SetParameters(0.005, 5.3, 0.7, 0.0097);

    // Set initial quad fit parameters
    spQuadFit->SetParameters(0.015, 0.001, 0.9552);

    TGraphErrors *CorrKEres_err_0 = new TGraphErrors(fitDivs, pxE, pyCE, exE, eyCE);

    CorrKEres_err_0->Draw();
    CorrKEres_err_0->SetTitle("Error plot of Corrected Gauss fit sigma vs true energy");
    CorrKEres_err_0->GetXaxis()->SetTitle("E_{True} (GeV)");
    CorrKEres_err_0->GetYaxis()->SetTitle("Sigma (GeV)");
    
    CorrKEres_err_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    CorrKEres_err_0->Fit("spQuadFit");

    c->SaveAs(saveDir + "/phaseHists/CorrKEres_err.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *CorrKTres_err_0 = new TGraphErrors(fitDivs, pxT, pyCT, exT, eyCT);

    CorrKTres_err_0->Draw();
    CorrKTres_err_0->SetTitle("Error plot of Corrected Gauss fit sigma vs true polar angle");
    CorrKTres_err_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    CorrKTres_err_0->GetYaxis()->SetTitle("Sigma (GeV)");
    
    CorrKTres_err_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    CorrKTres_err_0->Fit("CorrTermTFit_0", "RW");
    CorrKTres_err_0->Fit("CorrTermTFit_1", "RW");
    CorrKTres_err_0->Fit("CorrTermTFit_2", "RW");

    CorrTermTFit_0->Draw("SAME");
    CorrTermTFit_1->Draw("SAME");
    CorrTermTFit_2->Draw("SAME");

    c->SaveAs(saveDir + "/phaseHists/CorrKTres_err.png");
    c->Close();
    c = new TCanvas();
  }

  cout << "Test5" << endl;

  ///////////////////
  // STACKED HISTS //
  ///////////////////

  gStyle->SetPalette(kRust);

  if(ECH) {

    // Energy
    hs = new THStack("hs", "Energy of photons;p_{T} (GeV);Count");

    c->SetLogy();

    photonTruth_E_0->SetLineColor(kBlue);
    recoPhotonTruth_E_0->SetLineColor(kRed);
    photonReco_E_0->SetLineColor(kBlack);

    hs->Add(photonTruth_E_0);
    hs->Add(recoPhotonTruth_E_0);
    hs->Add(photonReco_E_0);

    hs->Draw("nostack");

    gPad->SetGrid(1, 0);
    gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");

    c->SaveAs(saveDir + "/All_E.png");
    c->Close();
    c = new TCanvas();

    // pT
    hs = new THStack("hs", "Transverse momentum of photons;p_{T} (GeV / c);Count");
  
    c->SetLogy();

    photonTruth_pt_0->SetLineColor(kBlue);
    recoPhotonTruth_pt_0->SetLineColor(kRed);
    photonReco_pt_0->SetLineColor(kBlack);

    hs->Add(photonTruth_pt_0);
    hs->Add(recoPhotonTruth_pt_0);
    hs->Add(photonReco_pt_0);

    hs->Draw("nostack");

    gPad->SetGrid(1, 0);
    gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");

    c->SaveAs(saveDir + "/All_pT.png");
    c->Close();
    c = new TCanvas();

    // Polar angle
    hs = new THStack("hs", "Polar angle of photons;Angle (Rads);Count");

    c->SetLogy();

    photonTruth_PA_0->SetLineColor(kBlue);
    recoPhotonTruth_PA_0->SetLineColor(kRed);
    photonReco_PA_0->SetLineColor(kBlack);

    hs->Add(photonTruth_PA_0);
    hs->Add(recoPhotonTruth_PA_0);
    hs->Add(photonReco_PA_0);

    hs->Draw("nostack");

    gPad->SetGrid(1, 0);
    gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");

    c->SaveAs(saveDir + "/All_PA.png");
    c->Close();
    c = new TCanvas();

    // Azimuth
    hs = new THStack("hs", "Azimuth of photons;Azimuth (Rads);Count");

    c->SetLogy();

    photonTruth_A_0->SetLineColor(kBlue);
    recoPhotonTruth_A_0->SetLineColor(kRed);
    photonReco_A_0->SetLineColor(kBlack);

    hs->Add(photonTruth_A_0);
    hs->Add(recoPhotonTruth_A_0);
    hs->Add(photonReco_A_0);

    hs->Draw("nostack");

    gPad->SetGrid(1, 0);
    gPad->BuildLegend(0.75, 0.75, 0.9, 0.9, "");

    c->SaveAs(saveDir + "/All_A.png");
    c->Close();
    c = new TCanvas();
  }

  ///////////////////////
  // SIMPLE HISTOGRAMS //
  ///////////////////////

  if(ESH) {

    // ENERGY

    photonTruth_E_0->Draw();
    photonTruth_E_0->SetTitle("All truth energies");
    photonTruth_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    photonTruth_E_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allTruth_E.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_E_0->Draw();
    recoPhotonTruth_E_0->SetTitle("Linked truth energies");
    recoPhotonTruth_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    recoPhotonTruth_E_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_E.png");
    c->Close();
    c = new TCanvas();

    photonReco_E_0->Draw();
    photonReco_E_0->SetTitle("Reconstructed energies");
    photonReco_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    photonReco_E_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allReco_E.png");
    c->Close();
    c = new TCanvas();
    
    // TRANSVERSE MOMENTUM
    
    photonTruth_pt_0->Draw();
    photonTruth_pt_0->SetTitle("All truth p_{T}");
    photonTruth_pt_0->GetXaxis()->SetTitle("p_{T} (GeV)");
    photonTruth_pt_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allTruth_pt.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_pt_0->Draw();
    recoPhotonTruth_pt_0->SetTitle("Linked truth p_{T}");
    recoPhotonTruth_pt_0->GetXaxis()->SetTitle("p_{T} (GeV)");
    recoPhotonTruth_pt_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_pt.png");
    c->Close();
    c = new TCanvas();

    photonReco_pt_0->Draw();
    photonReco_pt_0->SetTitle("Reconstructed p_{T}");
    photonReco_pt_0->GetXaxis()->SetTitle("p_{T} (GeV)");
    photonReco_pt_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allReco_pt.png");
    c->Close();
    c = new TCanvas();

    // POLAR ANGLE

    photonTruth_PA_0->Draw();
    photonTruth_PA_0->SetTitle("All truth polar angles");
    photonTruth_PA_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    photonTruth_PA_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allTruth_PA.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_PA_0->Draw();
    recoPhotonTruth_PA_0->SetTitle("Linked truth polar angles");
    recoPhotonTruth_PA_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    recoPhotonTruth_PA_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_PA.png");
    c->Close();
    c = new TCanvas();

    photonReco_PA_0->Draw();
    photonReco_PA_0->SetTitle("Reconstructed polar angles");
    photonReco_PA_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    photonReco_PA_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allReco_PA.png");
    c->Close();
    c = new TCanvas();

    // AZIMUTH

    photonTruth_A_0->Draw();
    photonTruth_A_0->SetTitle("All truth azimuths");
    photonTruth_A_0->GetXaxis()->SetTitle("Azimuth (Rads)");
    photonTruth_A_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allTruth_A.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_A_0->Draw();
    recoPhotonTruth_A_0->SetTitle("Linked truth azimuths");
    recoPhotonTruth_A_0->GetXaxis()->SetTitle("Azimuth (Rads)");
    recoPhotonTruth_A_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_A.png");
    c->Close();
    c = new TCanvas();

    photonReco_A_0->Draw();
    photonReco_A_0->SetTitle("Reconstructed azimuths");
    photonReco_A_0->GetXaxis()->SetTitle("Azimuths (Rads)");
    photonReco_A_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/allReco_A.png");
    c->Close();
    c = new TCanvas();
  }

  // OTHERS

  if(rH) {
    deltaR_0->Draw();
    deltaR_0->SetTitle("Delta R");
    deltaR_0->GetXaxis()->SetTitle("Delta R (Rads)");
    deltaR_0->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/resolutionHists/deltaR.png");
    c->Close();
    c = new TCanvas();

    // matchDRphoton_0->Draw();
    // matchDRphoton_0->SetTitle("Delta R between best match photon reco particles");
    // matchDRphoton_0->GetXaxis()->SetTitle("Delta R (Rads)");
    // matchDRphoton_0->GetYaxis()->SetTitle("Count");
    // 
    // c->SaveAs(saveDir + "/resolutionHists/matchDRphoton.png");
    // c->Close();
    // c = new TCanvas();
    // 
    // matchDRneutron_0->Draw();
    // matchDRneutron_0->SetTitle("Delta R between best match neutron reco particles");
    // matchDRneutron_0->GetXaxis()->SetTitle("Delta R (Rads)");
    // matchDRneutron_0->GetYaxis()->SetTitle("Count");
    // 
    // c->SaveAs(saveDir + "/resolutionHists/matchDRneutron.png");
    // c->Close();
    // c = new TCanvas();
    // 
    // matchDRelectron_0->Draw();
    // matchDRelectron_0->SetTitle("Delta R between best match electron reco particles");
    // matchDRelectron_0->GetXaxis()->SetTitle("Delta R (Rads)");
    // matchDRelectron_0->GetYaxis()->SetTitle("Count");
    // 
    // c->SaveAs(saveDir + "/resolutionHists/matchDRelectron.png");
    // c->Close();
    // c = new TCanvas();
    // 
    // matchDRpion_0->Draw();
    // matchDRpion_0->SetTitle("Delta R between best match pion reco particles");
    // matchDRpion_0->GetXaxis()->SetTitle("Delta R (Rads)");
    // matchDRpion_0->GetYaxis()->SetTitle("Count");
    // 
    // c->SaveAs(saveDir + "/resolutionHists/matchDRpion.png");
    // c->Close();
    // c = new TCanvas();
    // 
    // hs = new THStack("hs", "Delta R between best match and other reco particles;Delta R (Rads);Count");
    // 
    // matchDRphoton_0->SetLineColor(1);
    // matchDRneutron_0->SetLineColor(2);
    // matchDRelectron_0->SetLineColor(4);
    // matchDRpion_0->SetLineColor(8);
    // 
    // hs->Add(matchDRphoton_0);
    // hs->Add(matchDRneutron_0);
    // hs->Add(matchDRelectron_0);
    // hs->Add(matchDRpion_0);
    // 
    // hs->Draw();
    // 
    // gPad->SetGrid(1, 0);
    // gPad->BuildLegend(0.7, 0.7, 0.9, 0.9, "");
    // 
    // c->SaveAs(saveDir + "/resolutionHists/matchDRStack.png");
    // c->Close();
    // c = new TCanvas();
    
    Eres_0->Draw();
    Eres_0->SetTitle("Energy Resolution");
    Eres_0->GetXaxis()->SetTitle("Energy Resolution");
    Eres_0->GetYaxis()->SetTitle("Count");

    Eres_0->Fit("gaus");

    c->SaveAs(saveDir + "/resolutionHists/Eres.png");
    c->Close();
    c = new TCanvas();
  }

  if(pH) {
    gStyle->SetPalette(kRainBow);

    Ephase_0->Draw("colz");
    Ephase_0->SetTitle("Reco vs Truth Energy Phase Plot");
    Ephase_0->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    Ephase_0->GetYaxis()->SetTitle("True Energy (GeV)");

    Ephase_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/Ephase.png");
    c->Close();
    c = new TCanvas();

    K2Corr_0->SetMarkerSize(10);

    K2Corr_0->Draw("surf1");
    K2Corr_0->SetTitle("Energy, Polar Angle, and Energy Correction Phase Plot");
    K2Corr_0->GetXaxis()->SetTitle("E_{Reco} (GeV)");
    K2Corr_0->GetYaxis()->SetTitle("Polar Angle_{Reco} (Rads)");
    K2Corr_0->GetZaxis()->SetTitle("Correction term");

    K2Corr_0->GetXaxis()->SetTitleOffset(2.0);
    K2Corr_0->GetYaxis()->SetTitleOffset(2.0);
    K2Corr_0->GetZaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    K2Corr_0->SetMarkerSize(10);

    c->SaveAs(saveDir + "/phaseHists/K2Corr.png");
    c->Close();
    c = new TCanvas();

    // K2Corr axes
    K2Corr_0->ProjectionX()->Draw();
    K2Corr_0->SetTitle("Energy distribution");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_aX.png");
    c->Close();
    c = new TCanvas();

    K2Corr_0->ProjectionY()->Draw();
    K2Corr_0->SetTitle("Polar Angle distribution");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_aY.png");
    c->Close();
    c = new TCanvas();

    K2Corr_0->ProjectionZ()->Draw();
    K2Corr_0->SetTitle("Energy Correction distribution");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_aZ.png");
    c->Close();
    c = new TCanvas();

    // K2Corr profiles

    K2Corr_0->Project3D("yx")->Draw("colz");
    K2Corr_0->SetTitle("Energy, Polar Angle Phase Plot");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_Pxy.png");
    c->Close();
    c = new TCanvas();
    K2Corr_0->Project3D("zx")->Draw("colz");
    K2Corr_0->SetTitle("Energy, Energy Correction Phase Plot");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_Pzx.png");
    c->Close();
    c = new TCanvas();

    K2Corr_0->Project3D("zy")->Draw("colz");
    K2Corr_0->SetTitle("Polar Angle, and Energy Correction Phase Plot");
  
    K2Corr_0->SetStats(0);

    c->SaveAs(saveDir + "/phaseHists/K2Corr_Pzy.png");
    c->Close();
    c = new TCanvas();

    K2Corr2D_0->SetNpx(20);
    K2Corr2D_0->SetNpy(20);

    K2Corr2D_0->SetMinimum(0.94);
    K2Corr2D_0->SetMaximum(1.00);
  
    K2Corr2D_0->Draw("surf1");
    K2Corr2D_0->SetTitle("E_{True}, Polar Angle, and Energy Correction Phase Plot");
    K2Corr2D_0->GetXaxis()->SetTitle("E_{Reco} (GeV)");
    K2Corr2D_0->GetYaxis()->SetTitle("Polar Angle_{Reco} (Rads)");

    K2Corr2D_0->GetXaxis()->SetRangeUser(0, 5000);
    K2Corr2D_0->GetYaxis()->SetRangeUser(0, 3.2);
    K2Corr2D_0->GetZaxis()->SetRangeUser(0.94, 1.00);

    K2Corr2D_0->GetXaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetYaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetZaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/K2Corr2D.png");
    c->Close();
    c = new TCanvas();

    KEtrue_0->Draw();
    KEtrue_0->SetTitle("Correction factor as function of true energy");
    KEtrue_0->GetXaxis()->SetTitle("E_{True} (GeV)");
    KEtrue_0->GetYaxis()->SetTitle("Correction term");

    c->SaveAs(saveDir + "/phaseHists/KEtrue.png");
    c->Close();
    c = new TCanvas();    
  }

  ///////////////////////////
  // EFFICIENCY HISTOGRAMS //
  ///////////////////////////

  if(eH) {
    
    TEfficiency *E_Eff = new TEfficiency("E_Eff", "Efficiency of photon energy reconstruction;E (GeV);Efficiency", 20, lowE, highE);
    TFile* eff_File = new TFile("effFile.root", "recreate");

    if(TEfficiency::CheckConsistency(*recoPhotonTruth_E_0, *photonTruth_E_0)) {
      E_Eff->SetPassedHistogram(*recoPhotonTruth_E_0, "f");
      E_Eff->SetTotalHistogram(*photonTruth_E_0, "f");
      eff_File->Write();
    }
    
    TEfficiency* pt_Eff = new TEfficiency("pt_Eff", "Efficiency of photon p_{T} reconstruction;p_{T} (GeV / c);Efficiency", 20, lowE, highE);
    
    if(TEfficiency::CheckConsistency(*recoPhotonTruth_pt_0, *photonTruth_pt_0)) {
      pt_Eff->SetPassedHistogram(*recoPhotonTruth_pt_0, "f");
      pt_Eff->SetTotalHistogram(*photonTruth_pt_0, "f");
      eff_File->Write();
    }

    TEfficiency* PA_Eff = new TEfficiency("PA_Eff", "Efficiency of photon polar angle reconstruction;Polar angle (Rads);Efficiency", 20, -1.6, 1.6);
    eff_File = new TFile("effFile.root", "recreate");

    if(TEfficiency::CheckConsistency(*recoPhotonTruth_PA_0, *photonTruth_PA_0)) {
      PA_Eff->SetPassedHistogram(*recoPhotonTruth_PA_0, "f");
      PA_Eff->SetTotalHistogram(*photonTruth_PA_0, "f");
      eff_File->Write();
    }
    
    E_Eff->SetLineColor(4);
    pt_Eff->SetLineColor(4);
    PA_Eff->SetLineColor(4);

    //draws an efficiency graph for the energy
    TMultiGraph *mg = new TMultiGraph();
    c = new TCanvas();

    mg->Add(E_Eff->CreateGraph());

    mg->SetTitle("Energy reconstruction efficiency");
    mg->GetXaxis()->SetTitle("Energy (GeV)");
    mg->GetYaxis()->SetTitle("Efficiency");

    mg->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    mg->Draw("aZ");

    c->SaveAs(saveDir + "/efficiency/Eff_E.png");
    c->Close();

    //draws an efficiency graph for the transverse momentum
    mg = new TMultiGraph();
    c = new TCanvas();

    mg->Add(pt_Eff->CreateGraph());
    
    mg->SetTitle("Transverse momentum reconstruction efficiency");
    mg->GetXaxis()->SetTitle("p_{T} (GeV / c)");
    mg->GetYaxis()->SetTitle("Efficiency");

    mg->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    mg->Draw("aZ");

    c->SaveAs(saveDir + "/efficiency/Eff_pt.png");
    c->Close();

    //draws an efficiency graph for the polar angle
    mg = new TMultiGraph();
    c = new TCanvas();

    mg->Add(PA_Eff->CreateGraph());

    mg->SetTitle("Polar angle reconstruction efficiency");
    mg->GetXaxis()->SetTitle("Polar angle (rads)");
    mg->GetYaxis()->SetTitle("Efficiency");

    mg->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    mg->Draw("aZ");

    c->SaveAs(saveDir + "/efficiency/Eff_PA.png");
    c->Close();
  }
}

void zprime_jj_full_analysis() {
  zprime_jj_full_analysis_H(false, "00000000");
}

void zprime_jj_full_analysis(char const *args) {
  zprime_jj_full_analysis_H(true, args);
}

// int** matchJets(int njet_R, int njet_T, TVector3 recoVec3[], TVector3 truthVec3[]) {
  
//   int matchArrR[njet_T];
//   int matchArrT[njet_T];
//   int matchCount = 0;

//   int minAxisLen = min(njet_R, njet_T);

//   Double_t matchArrW[njet_R][njet_T];

//   for(int i = 0; i < njet_R; i++) {
//     for(int j = 0; j < njet_T; j++) {
//       matchArrW[i][j] = 1 - ( abs((recoVec3[j] - truthVec3[i]).Mag() / truthVec3[i].Mag()) / 2 ) + ( sqrt( pow(recoVec3[j].Eta() - truthVec3[i].Eta(), 2) + pow(recoVec3[j].Phi() - truthVec3[i].Phi(), 2) ) / (2 * TMath::Pi()) );
//     }
//   }

//   int banR[njet_R];
//   int banT[njet_T];
//   int banCount = 0;
//   for(int l = 0; l < minAxisLen; l++) {
//     for(int i = 0; i < njet_T; i++) {
//       Double_t maxWR = 0;
//       Double_t maxWT = 0;
//       int maxWRIndex = -1;
      
//       for(int j = 0; i < njet_R; j++) {
// 	if(!contains(banR, banCount, j)) {
// 	  if(matchArrW[j, i] > maxWR) {
// 	    maxWR = matchArrW[j, i];
// 	    maxWRIndex = j;
// 	  }
// 	}
//       }
//       if(maxWRIndex != -1) {
// 	for(int j = 0; i < njet_T; j++) {
// 	  if(!contains(banT, banCount, j)) {
// 	    if(matchArrW[maxWRIndex, j] > maxWT) {
// 	      maxWT = matchArrW[maxWRIndex, j];
// 	    }
// 	  }
// 	}
//       }
//       if((maxWT == maxWR) && (maxWT >= matchWeightCutoff)) {
// 	matchArrR[matchCount] = maxWRIndex;
// 	matchArrT[matchCount] = i;
// 	matchCount++;
// 	banR[banLen] = maxWRIndex;
// 	banT[banLen] = i;
// 	banCount++;
//       }
//     }
//     if(banCount == minAxisLen) {
//       break;
//     }
//   }

//   int ret[][] = {matchArrR, matchArrT};
//   return ret;
// }

// {  
  // int redoTArr[njet_T];
  // int redoRArr[njet_T];
  // int redoCount = 0;
  // 
  // for(int i = 0; i < njet_T; i++) {
  //   Double_t weight = 0;
  //   Double_t weightMax = 0;
  //   int maxWeightIndex = -1;
  //   for(int j = 0; j < njet_R; j++) {
  //     weight = 1 - ( abs((recoVec3[j] - truthVec3[i]).Mag() / truthVec3[i].Mag()) / 2 ) + ( sqrt( pow(recoVec3[j].Eta() - truthVec3[i].Eta(), 2) + pow(recoVec3[j].Phi() - truthVec3[i].Phi(), 2) ) / (2 * TMath::Pi()) );
  //     if( ( weight > weightMax ) > matchWeightCutoff ) {
  // 	weightMax = weight;
  // 	maxWeightIndex = j;
  //     }
  //   }
  //   if(maxWeightIndex >= 0) {
  //     for(int j = 0; j < i; j++) {
  // 	if(matchArrR[j] == maxWeightIndex) {
  // 	  int fixIndex;
  // 	  int banIndex;
  // 	  if(matchArrW[j] < weightMax) {
  // 	    redoTArr[redoCount] = matchArrT[j];
  // 	    redoRArr[redoCount] = matchArrR[j];
  // 	    redoCount ++;
  // 	    matchArrR[j] = maxWeightIndex;
  // 	    matchArrT[j] = i;
  // 	    matchArrW[j] = weightMax;
  // 	  } else {
  // 	    redoTArr[redoCount] = i;
  // 	    redoCount ++;
  // 	  }
  // 	} else if(j+1 == i) {
  // 	  matchArrR[matchCount] = maxWeightIndex;
  // 	  matchArrT[matchCount] = i;
  // 	  matchArrW[matchCount] = weightMax;
  // 	  matchCount++;
  // 	}
  //     }
  //   }
  // }
// }

// int** matchJetsRedo(int njet_R, int njet_T, intTVector3  recoVec3[], TTreeReaderArray truthVec3[], int fixIndex, int banIndex[], int rArr[][]) {
//   int matchArrR[njet_T] = rArr[0];
//   int matchArrT[njet_T] = rArr[1];
//   int matchArrW[njet_T] = rArr[2];
//   
//   Double_t weight = 0;
//   Double_t weightMax = 0;
//   int maxWeightIndex = -1;
//   for(int j = 0; i < njet_R; i++) {
//     weight = 1 - ( abs((recoVec3[j] - truthVec3[fixIndex]).Mag() / truthVec3[fixIndex].Mag()) / 2 ) + ( sqrt( pow(recoVec3[j].Eta() - truthVec3[fixIndex].Eta(), 2) + pow(recoVec3[j].Phi() - truthVec3[fixIndex].Phi(), 2) ) / (2 * TMath::Pi()) );
//     if( ( weight > weightMax ) > matchWeightCutoff ) {
//       weightMax = weight;
//       maxWeightIndex = j;
//     }
//   }
//   if(maxWeightIndex >= 0) {
//     for(int j = 0; j < i; j++) {
//       if(matchArrR[j] == maxWeightIndex) {
// 	if(matchArrW[j] < weightMax) {
// 	  redoTArr[redoCount] = matchArrT[j];
// 	  redoRArr[redoCount] = matchArrR[j];
// 	  redoCount ++;
// 	  matchArrR[j] = maxWeightIndex;
// 	  matchArrT[j] = i;
// 	  matchArrW[j] = weightMax;
// 	} else {
// 	  redoTArr[redoCount] = i;
// 	  redoCount ++;
// 	}
//       } else if(j+1 == i) {
// 	matchArrR[matchCount] = maxWeightIndex;
// 	matchArrT[matchCount] = i;
// 	matchArrW[matchCount] = weightMax;
// 	matchCount++;
//       }
//     }
//   }
// }
