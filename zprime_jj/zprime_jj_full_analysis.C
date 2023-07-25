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
  float maxERecoErr = 6;

  // Max R a reconstructed particle can be
  float maxRRecoErr = 0.1;

  // delta R that any semi-major reconstructed particles can be from the most major one
  float coneAcceptR = 0.13;

  // minimum weight that a match can be
  float matchWeightCutoff = 0.96;
  
  // If true, limits the reconstructed events used to those which only had one reconstructed particle
  bool limit1Reco = false;

  // Number of particles
  int maxN = 40000;

  // Low pT range
  int lowE = 0;
  // High pT range
  int highE = 5100;

  // Amount that the graphs show below/above the min/max energy
  int lowEErr = 200;
  int highEErr = 1000;

  // bimodal distribution fit
  double bimodSplit = 900;

  // minimum true energy cutoff for eliminating bad jets
  double minTruthE = 0;

  // number of bins for res fits
  int fitDivs = 50;

  // number of bins for resolution plots
  int resBins = 100;

  // number of bins for showing particle distributions
  int partBins = 40;

  // number of bins used for sig fits
  int sigFitBins = 40;

  // IVM bins
  int IVMBins = 30;

  // scatter correction factor, energy, polar angle min Z
  double K2CorrMin = 0.00;
  // scatter correction factor, energy, polar angle max Z
  double K2CorrMax = 3.00;

  // energy resolution plot min
  double resMin = -1.0;
  // energy resolution plot max
  double resMax = 1.0;

  // IVM min
  double IVMMin = 0;
  // IVM max
  double IVMMax = 12000;

  // Correction term distribution hist minimum
  double CorrTermMin = 0.5;
  // Correction term distribution hist maximum
  double CorrTermMax = 3.5;

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
  bool IVM = false;             // Invarient mass
  bool BMD = false;             // Bimod fits for plots
  bool EP = false;              // Extra (usually debug/spam) plots for certain categories above
  
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
    IVM = charToBool(args[8]);
    BMD = charToBool(args[9]);
    EP = charToBool(args[10]);
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
  TH1F *recoPhotonReco_E_0 = new TH1F("rPR_E_0", "Reco Matched Energy", 20, lowE - lowEErr, highE + highEErr);

  TH1F *photonTruth_pt_0 = new TH1F("pt_pt_0", "All Truth p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //All true transverse momenta
  TH1F *recoPhotonTruth_pt_0 = new TH1F("prt_pt_0", "Reco True p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //True transverse momentum of reconstructed photons
  TH1F *photonReco_pt_0 = new TH1F("pr_pt_0", "Reco p_{T}", 20, (lowE * sin(TMath::Pi() / 18)) - lowEErr, highE + highEErr); //Reconstructed transverse momentum of linked photons

  TH1F *photonTruth_PA_0 = new TH1F("pt_PA_0", "All Truth Polar Angle", 20, 0, 3.2); //All true polar anglesx
  TH1F *recoPhotonTruth_PA_0 = new TH1F("prt_PA_0", "Reco True Polar Angle", 20, 0, 3.2); //True polar angle of reconstructed photons
  TH1F *photonReco_PA_0 = new TH1F("pr_PA_0", "Reco Polar Angle", 20, 0, 3.2); //Reconstructed polar angle of linked photons

  TH1F *photonTruth_A_0 = new TH1F("pt_A_0", "All Truth Azimuth", 20, -1.6, 1.6); //All true azimuths
  TH1F *recoPhotonTruth_A_0 = new TH1F("prt_A_0", "Reco True Azimuth", 20, -1.6, 1.6); //True azimuth of reconstructed photons
  TH1F *photonReco_A_0 = new TH1F("pr_A_0", "Reco Azimuth", 20, -1.6, 1.6); //Reconstructed azimuth of linked photons

  TH1F *deltaR_0 = new TH1F("dR_0", "Delta R", resBins, 0 , 0.05); // delta R between mc truth and best match reco
  TH1F *deltaRBimod0 = new TH1F("DRB0", "Delta R", resBins, 0, 0.002);
  TH1F *deltaRBimod1 = new TH1F("DRB1", "Delta R", resBins, 0, 0.002);
  TH1F *ERes_0 = new TH1F("ER_0", "Energy accuracy", resBins, resMin, resMax); //Resolution of reconstructed energy
  TH1F *EResBimod0 = new TH1F("ERB0", "Energy resolution", resBins, resMin, resMax);
  TH1F *EResBimod1 = new TH1F("ERB1", "Energy resolution", resBins, resMin, resMax);
  TH2F *Ephase_0 = new TH2F("Epha_0", "True vs Reco energy phase plot", 100, -lowEErr, highE + highEErr, 100, lowE - lowEErr, highE + highEErr); //Phase plot of true energy w.r.t reconstructed energy
  TH3F *K2Corr_0 = new TH3F("K2c_0", "Correction w.r.t. polar angle and energy", 20, lowE, highE, 20, 0, 3.2, 20, K2CorrMin, K2CorrMax); //Correction term plot w.r.t. polar angle and energy
  TGraph2D *K2Corr2D_0 = new TGraph2D(); //2D surface of the above scatter plot
  TGraph2DErrors *CorrTerm2D_0;
  TH2F *KEtrue_0 = new TH2F("KEt_0", "Correction factor as function of energy", 100, lowE, highE, 100, K2CorrMin, K2CorrMax);
  TH1F *CorrE_Bins_0 = new TH1F("CEB_0", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr); //Corrected photon energies
  TH1F *CorrERes_Bins_0 = new TH1F("CERB_0", "Corrected energy resolution -- Bins", resBins, resMin, resMax); //Resolution of corrected photon energies
  TH1F *CorrE_Fit_0 = new TH1F("CEF_0", "Corrected energy -- Fit", 20, lowE - lowEErr, highE + highEErr); //Corrected photon energies
  TH1F *CorrERes_Fit_0 = new TH1F("CERF_0", "Corrected energy resolution -- Fit", resBins, resMin, resMax); //Resolution of corrected photon energies

  TH1F *CorrE_BinsBimod0 = new TH1F("CEBB0", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr);
  TH1F *CorrE_BinsBimod1 = new TH1F("CEBB1", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr);
  TH1F *CorrE_FitBimod0 = new TH1F("CEFB0", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr);
  TH1F *CorrE_FitBimod1 = new TH1F("CEFB1", "Corrected energy -- Bins", 20, lowE - lowEErr, highE + highEErr);
  
  TH1F *CorrERes_BinsBimod0 = new TH1F("CREBB0", "Corrected energy resolution -- Bins", resBins, resMin, resMax);
  TH1F *CorrERes_BinsBimod1 = new TH1F("CREBB1", "Corrected energy resolution -- Bins", resBins, resMin, resMax);
  TH1F *CorrERes_FitBimod0 = new TH1F("CREFB0", "Corrected energy resolution -- Fit", resBins, resMin, resMax);
  TH1F *CorrERes_FitBimod1 = new TH1F("CREFB1", "Corrected energy resolution -- Fit", resBins, resMin, resMax);

  TH1F *TruthIVM = new TH1F("TIVM", "Truth Invarient Mass", IVMBins, IVMMin, IVMMax);
  TH1F *RecoIVM = new TH1F("RIVM", "Reconstructed Invarient Mass", IVMBins, IVMMin, IVMMax);
  TH1F *CorrIVM = new TH1F("CIVM", "Corrected Reconstructed Invarient Mass", IVMBins, IVMMin, IVMMax);

  // TH1F *CorrETermsMeans = new TH1F("CETM" "Correction Term Fit Means", fitDivs, )
  
  THStack *hs;
  TF2 *Corr2DFitFunc_00, *Corr2DFitFunc_01, *Corr2DFitFunc_02, *Corr2DFitFunc_10, *Corr2DFitFunc_11, *Corr2DFitFunc_12, *Corr2DFitFunc_20, *Corr2DFitFunc_21, *Corr2DFitFunc_22, *Corr2DFitFunc_30, *Corr2DFitFunc_31, *Corr2DFitFunc_32;
  TF2 *Corr2DFitFunc_0, *Corr2DFitFunc_1, *Corr2DFitFunc_2;

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
  
  double** corr2DBins;
  double** corr2DBinN;
  double** corr2DSTD;
  corr2DBins = build2DArrD(fitDivs, fitDivs);
  corr2DBinN = build2DArrD(fitDivs, fitDivs);
  corr2DSTD = build2DArrD(fitDivs, fitDivs);
  
  TH1F* sigFitEHists[fitDivs];
  TH1F* sigFitTHists[fitDivs];
  TH1F* corrSigFitEHists[fitDivs];
  TH1F* corrSigFitTHists[fitDivs];

  TH1F* sigFitEHistsBimod0[fitDivs];
  TH1F* sigFitEHistsBimod1[fitDivs];
  TH1F* sigFitTHistsBimod0[fitDivs];
  TH1F* sigFitTHistsBimod1[fitDivs];
  
  TH1F* recoTypesE[5];
  TH1F* recoTypesT[5];
  TH1F* mcRecoTypesE[5];
  TH1F* mcRecoTypesT[5];

  TH1F* CorrETerms[fitDivs];
  TH1F* CorrTTerms[fitDivs];

  double corrEBins[fitDivs];
  int corrEBinN[fitDivs];
  double corrTBins[fitDivs];
  int corrTBinN[fitDivs];
  double corrESTD[fitDivs];
  double corrTSTD[fitDivs];

  double resEBins[fitDivs];
  int resEBinN[fitDivs];
  double resTBins[fitDivs];
  int resTBinN[fitDivs];
  double resESTD[fitDivs];
  double resTSTD[fitDivs];
  
  double recoEnergyArr[maxN];
  double trueEnergyArr[maxN];
  double recoPAArr[maxN];

  double recoIVMEArr0[maxN];
  double recoIVMEArr1[maxN];
  double recoIVMTArr0[maxN];
  double recoIVMTArr1[maxN];

  double recoIVMxArr0[maxN];
  double recoIVMxArr1[maxN];
  double recoIVMyArr0[maxN];
  double recoIVMyArr1[maxN];
  double recoIVMzArr0[maxN];
  double recoIVMzArr1[maxN];
  int recoIVMCount = 0;

  cout << "Test0" << endl;
  
  for(int i = 0; i < fitDivs; i++) {
    corrEBins[i] = 0;
    corrEBinN[i] = 0;
    corrESTD[i] = 0;
    corrTBins[i] = 0;
    corrTBinN[i] = 0;
    corrTSTD[i] = 0;

    resEBins[i] = 0;
    resEBinN[i] = 0;
    resTBins[i] = 0;
    resTBinN[i] = 0;
    resESTD[i] = 0;
    resTSTD[i] = 0;
    
    // memset(corr2DBins[i], 0.0, sizeof(*corr2DBins[i]));
    // memset(corr2DBinN[i], 0, sizeof(*corr2DBinN[i]));
    // memset(corr2DSTD[i], 0.0, sizeof(*corr2DSTD[i]));
    
    for(int j = 0; j < fitDivs; j++) {
      corr2DSTD[i][j] = 0;
      corr2DBins[i][j] = 0;
      corr2DBinN[i][j] = 0;
    }
  }
  
  cout << "Test0.01" << endl;
  
  if(ECF || rFH) {
    for(int i = 0; i < fitDivs; i++) {
      char s[4];
      
      char outE[100] = "";
      char outT[100] = "";
      char outCorrERes[100] = "";
      char outCorrTRes[100] = "";
      
      char outEBimod0[100] = "";
      char outEBimod1[100] = "";
      char outTBimod0[100] = "";
      char outTBimod1[100] = "";

      char outCorrETerms[100] = "";
      char outCorrTTerms[100] = "";
      
      sprintf(s, "%d", i);
      
      strcpy(outE, "sfEH");
      strcat(outE, s);
      strcpy(outT, "sfTH");
      strcat(outT, s);
      
      strcpy(outCorrERes, "csfEH");
      strcat(outCorrERes, s);
      strcpy(outCorrTRes, "csfTH");
      strcat(outCorrTRes, s);

      strcpy(outCorrETerms, "CET");
      strcat(outCorrETerms, s);
      strcpy(outCorrTTerms, "CTT");
      strcat(outCorrTTerms, s);
      
      sigFitEHists[i] = new TH1F(outE, "Energy resolution", sigFitBins, resMin, resMax);
      sigFitTHists[i] = new TH1F(outT, "Energy resolution", sigFitBins, resMin, resMax);
      corrSigFitEHists[i] = new TH1F(outCorrERes, "Energy resolution", sigFitBins, resMin, resMax);
      corrSigFitTHists[i] = new TH1F(outCorrTRes, "Energy resolution", sigFitBins, resMin, resMax);

      CorrETerms[i] = new TH1F(outCorrETerms, "Energy Correction Factors - Energy", fitDivs, CorrTermMin, CorrTermMax);
      CorrTTerms[i] = new TH1F(outCorrTTerms, "Energy Correction Factors - Polar Angle", fitDivs, CorrTermMin, CorrTermMax);
      
      if(BMD) {
	strcat(outEBimod0, outE);
        strcat(outEBimod0, "_0");
        strcat(outEBimod1, outE);
        strcat(outEBimod1, "_1");
        strcat(outTBimod0, outT);
        strcat(outTBimod0, "_0");
        strcat(outTBimod1, outT);
        strcat(outTBimod1, "_1");
	
	sigFitEHistsBimod0[i] = new TH1F(outEBimod0, "Energy resolution", sigFitBins, resMin, resMax);
	sigFitEHistsBimod1[i] = new TH1F(outEBimod1, "Energy resolution", sigFitBins, resMin, resMax);
	sigFitTHistsBimod0[i] = new TH1F(outTBimod0, "Energy resolution", sigFitBins, resMin, resMax);
	sigFitTHistsBimod1[i] = new TH1F(outTBimod1, "Energy resolution", sigFitBins, resMin, resMax);
      }
    }
  }
  
  // run 0
  while(myRecoReader.Next()) {
    myTruthReader.Next();
    
    int njet_R = njet_RRA.At(0);
    int njet_T = njet_TRA.At(0);
    
    TVector3 recoVec3[njet_R];
    TVector3 truthVec3[njet_T];
    
    for(int i = 0; i < njet_R; i++) {
      recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
      photonReco_E_0->Fill(jene_RRA.At(i));
      photonReco_pt_0->Fill(recoVec3[i].Perp());
      photonReco_PA_0->Fill(recoVec3[i].Theta());
      photonReco_A_0->Fill(recoVec3[i].Phi());
    }
    for(int i = 0; i < njet_T; i++) {
      truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
    }
    
    int** r = build2DArrI(2, njet_T);
    
    r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
    
    int matchArrR[njet_T];
    int matchArrT[njet_T];
    for(int i = 0; i < njet_T; i++) {
      matchArrR[i] = r[0][i];
      matchArrT[i] = r[1][i];
    }

    double addE = pow(jene_TRA.At(0) + jene_TRA.At(1), 2);
    double jvx = pow(jmox_TRA.At(0) + jmox_TRA.At(1), 2);
    double jvy = pow(jmoy_TRA.At(0) + jmoy_TRA.At(1), 2);
    double jvz = pow(jmoz_TRA.At(0) + jmoz_TRA.At(1), 2);

    TruthIVM->Fill( sqrt(addE - jvx - jvy - jvz) );

    if(contains(matchArrT, njet_T, 0) && contains(matchArrT, njet_T, 1)) {
      int R0 = matchArrR[indexOf(matchArrT, sizeof(matchArrT), 0)];
      int R1 = matchArrR[indexOf(matchArrT, sizeof(matchArrT), 1)];

      addE = pow(jene_RRA.At(0) + jene_RRA.At(1), 2);
      jvx = pow(jmox_RRA.At(0) + jmox_RRA.At(1), 2);
      jvy = pow(jmoy_RRA.At(0) + jmoy_RRA.At(1), 2);
      jvz = pow(jmoz_RRA.At(0) + jmoz_RRA.At(1), 2);
      
      RecoIVM->Fill(sqrt( addE - jvx - jvy - jvz) );
    }
    
    for(int i = 0; i < njet_T; i++) {
      
      bool matched = contains(matchArrT, njet_T, i);
      int recoMatch;
      
      mcmoTemp = truthVec3[i];
      
      if(minTruthE <= jene_TRA.At(i)) {
	
	if(matched) {
	  recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	  
	  rcmoTemp = recoVec3[recoMatch];

	  recoPhotonTruth_E_0->Fill(jene_TRA.At(i));
	  recoPhotonTruth_pt_0->Fill(mcmoTemp.Perp());
	  recoPhotonTruth_PA_0->Fill(mcmoTemp.Theta());
	  recoPhotonTruth_A_0->Fill(mcmoTemp.Phi());

	  recoPhotonReco_E_0->Fill(jene_RRA.At(recoMatch));
	  
	  if(contains(matchArrT, njet_T, 0) && contains(matchArrT, njet_T, 1)) {
	    if(i == 0) {
	      recoIVMEArr0[recoIVMCount] = jene_RRA.At(recoMatch);
	      recoIVMTArr0[recoIVMCount] = rcmoTemp.Theta();

	      recoIVMxArr0[recoIVMCount] = jmox_RRA.At(recoMatch);
	      recoIVMyArr0[recoIVMCount] = jmoy_RRA.At(recoMatch);
	      recoIVMzArr0[recoIVMCount] = jmoz_RRA.At(recoMatch);
	    } else if(i == 1) {
	      recoIVMEArr1[recoIVMCount] = jene_RRA.At(recoMatch);
              recoIVMTArr1[recoIVMCount] = rcmoTemp.Theta();

	      recoIVMxArr1[recoIVMCount] = jmox_RRA.At(recoMatch);
              recoIVMyArr1[recoIVMCount] = jmoy_RRA.At(recoMatch);
              recoIVMzArr1[recoIVMCount] = jmoz_RRA.At(recoMatch);

	      recoIVMCount++;
	    }
	  }
	  
	  deltaR_0->Fill(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ));
	  Ephase_0->Fill(jene_RRA.At(recoMatch), jene_TRA.At(i));
	  K2Corr_0->Fill(jene_RRA.At(recoMatch), rcmoTemp.Theta(), jene_TRA.At(i) / jene_RRA.At(recoMatch));
	  K2Corr2D_0->SetPoint(arrCount, jene_RRA.At(recoMatch), rcmoTemp.Theta(), jene_TRA.At(i) / jene_RRA.At(recoMatch));
	  
	  if(BMD) {
	    if(jene_RRA.At(recoMatch) < bimodSplit) {
	      deltaRBimod0->Fill(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ));
	    } else {
	      deltaRBimod1->Fill(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ));
	    }
	  }
	  
	  KEtrue_0->Fill(jene_TRA.At(i), jene_TRA.At(i) / jene_RRA.At(recoMatch));
	  
	  if(sqrt( pow(rcmoTemp.Eta() - mcmoTemp.Eta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ) < maxRRecoErr) {
	    
	    int RecoEnergyBin = (int) floor((jene_RRA.At(recoMatch) - lowE) * fitDivs / (highE - lowE));
	    int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	    int TrueEnergyBin = (int) floor((jene_TRA.At(i) - lowE) * fitDivs / (highE - lowE));
	    int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	    
	    while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	    while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	    while(RecoPABin < 0) { RecoPABin++; }
	    while(RecoPABin >= fitDivs) { RecoPABin--; }
	    
	    if(ECF || rFH) {
	      
	      sigFitEHists[ RecoEnergyBin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
	      sigFitTHists[ RecoPABin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
	      
	      if(BMD) {
		if(jene_RRA.At(recoMatch) < bimodSplit) {
		  sigFitEHistsBimod0[ TrueEnergyBin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
		  sigFitTHistsBimod0[ TruePABin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
		} else {
		  sigFitEHistsBimod1[ TrueEnergyBin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
		  sigFitTHistsBimod1[ TruePABin ]->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
		}
	      }
	    }

	    if( (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) < maxERecoErr) && (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) > (1.0/maxERecoErr)) ) {

	      evtCount++;
	      
	      if(ECA) {
		corrEBins[ RecoEnergyBin ] += (jene_TRA.At(i) / jene_RRA.At(recoMatch));
		corrEBinN[ RecoEnergyBin ]++;
		corrTBins[ RecoPABin ] += (jene_TRA.At(i) / jene_RRA.At(recoMatch));
		corrTBinN[ RecoPABin ]++;

		CorrETerms[ RecoEnergyBin ]->Fill(jene_TRA.At(i) / jene_RRA.At(recoMatch));
		CorrTTerms[ RecoPABin ]->Fill(jene_TRA.At(i) / jene_RRA.At(recoMatch));

		resEBins[ TrueEnergyBin ] += (jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i);
		resEBinN[ TrueEnergyBin ]++;
		resTBins[ TruePABin ] += (jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i);
		resTBinN[ TruePABin ]++;
		
		corr2DBins[ RecoEnergyBin ][ RecoPABin ] += (jene_TRA.At(i) / jene_RRA.At(recoMatch));
		corr2DBinN[ RecoEnergyBin ][ RecoPABin ]++;
	      }
	 
	    }
	  }
	  
	  ERes_0->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
	  if(BMD) {
	    if(jene_RRA.At(recoMatch) < bimodSplit) {
	      EResBimod0->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
	    } else {
	      EResBimod1->Fill((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i));
	    }
	  }
	  
	  recoEnergyArr[arrCount] = jene_RRA.At(recoMatch);
	  trueEnergyArr[arrCount] = jene_TRA.At(i);
	  recoPAArr[arrCount] = rcmoTemp.Theta();
	  arrCount++;
	}
	
	photonTruth_E_0->Fill(jene_TRA.At(i));
	photonTruth_pt_0->Fill(mcmoTemp.Perp());
	photonTruth_PA_0->Fill(mcmoTemp.Theta());
	photonTruth_A_0->Fill(mcmoTemp.Phi());
      }
    }
    // evtCount++;
  }
  
  cout << "EVTCOUNT1: " << evtCount << endl;

  //cout << "Test2" << endl;
  
  if(ECA) {
    
    // restart run
    myRecoReader.Restart();
    myTruthReader.Restart();
    
    myRecoReader.Next();
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
    
    int** r = build2DArrI(2, njet_T);
    r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
    
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

      if(minTruthE <= jene_TRA.At(i)) {
	
	if(matched) {
	  recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	  
	  rcmoTemp = recoVec3[recoMatch];
	  
	  int RecoEnergyBin = (int) floor((jene_RRA.At(recoMatch) - lowE) * fitDivs / (highE - lowE));
	  int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	  int TrueEnergyBin = (int) floor((jene_TRA.At(i) - lowE) * fitDivs / (highE - lowE));
	  int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	  
	  while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	  while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	  while(RecoPABin < 0) { RecoPABin++; }
	  while(RecoPABin >= fitDivs) { RecoPABin--; }
	  
	  if( (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) < maxERecoErr) && (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) > (1.0/maxERecoErr)) ) {
	    if(ECA) {
	      if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2) < 3 ) {
		corrESTD[ RecoEnergyBin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2);
	      } else {
		cout << "BAD: nEres " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " E STD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
	      }
	      
	      if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2) < 3 ) {
		corrTSTD[ RecoPABin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2);
	      } else {
		cout << "BAD: nEres: " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " PA STD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
	      }

	      if( pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2) < 3 ) {
		resESTD[ TrueEnergyBin ] += pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2);
	      } else {
		cout << "BAD: resE " << ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) << " E STD: " << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2) << endl;
	      }
		  
	      if( pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2) < 3 ) {
                resTSTD[ TruePABin ] += pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2);
              } else {
                cout << "BAD: resT " << ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) << " T STD: " << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2) << endl;
              }
	      
	      if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2) < 3 ) {
		corr2DSTD[ RecoEnergyBin ][ RecoPABin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2);
	      } else {
		cout << "BAD: nEres " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " 2D STD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
	      }
	    }
	  }
	}
      }
    }
    
    // cout << "Test2.5" << endl;
    evtCount = 0; 
    // run 1
      
    while(myRecoReader.Next()) {
      myTruthReader.Next();
      
      int njet_R = njet_RRA.At(0);
      int njet_T = njet_TRA.At(0);

      TVector3 recoVec3[njet_R];
      TVector3 truthVec3[njet_T];

      for(int i = 0; i < njet_R; i++) {
	recoVec3[i].SetXYZ(jmox_RRA.At(i), jmoy_RRA.At(i), jmoz_RRA.At(i));
	// cout << "recoVec3: " << recoVec3[i].Mag() << endl;
      }
      for(int i = 0; i < njet_T; i++) {
	truthVec3[i].SetXYZ(jmox_TRA.At(i), jmoy_TRA.At(i), jmoz_TRA.At(i));
	// cout << "truthVec3: " << truthVec3[i].Mag() << endl;
      }

      int** r = build2DArrI(2, njet_T);

      r = matchJets(njet_R, njet_T, recoVec3, truthVec3, matchWeightCutoff);
      
      int matchArrR[njet_T];
      int matchArrT[njet_T];
      for(int i = 0; i < njet_T; i++) {
	matchArrR[i] = r[0][i];
	matchArrT[i] = r[1][i];
	// cout << "MatchR: " << matchArrR[i] << " MatchT: " << matchArrT[i] << endl;
      }

      for(int i = 0; i < njet_T; i++) {

	// evtCount++;

	bool matched = contains(matchArrT, njet_T, i);
	int recoMatch;

	mcmoTemp = truthVec3[i];

	if(minTruthE <= jene_TRA.At(i)) {

	  // evtCount++;
	  
	  if(matched) {

	    evtCount++;
	    
	    recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	    
	    rcmoTemp = recoVec3[recoMatch];
	    
	    int RecoEnergyBin = (int) floor((jene_RRA.At(recoMatch) - lowE) * fitDivs / (highE - lowE));
	    int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	    int TrueEnergyBin = (int) floor((jene_TRA.At(i) - lowE) * fitDivs / (highE - lowE));
	    int TruePABin = (int) floor((mcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	    
	    while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	    while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	    while(RecoPABin < 0) { RecoPABin++; }
	    while(RecoPABin >= fitDivs) { RecoPABin--; }
	    
	    if( (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) < maxERecoErr) && (abs(jene_TRA.At(i) / jene_RRA.At(recoMatch)) > (1.0/maxERecoErr)) ) {
	      if(ECA) {
		if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2) < 99 ) {
		  corrESTD[ RecoEnergyBin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2);
		  // cout << "1:" << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2) << " ";
		} else {
		  cout << "BAD: nEres " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " CorrESTD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrEBins[ RecoEnergyBin ] / (double)corrEBinN[ RecoEnergyBin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
		}
		
		if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2) < 99 ) {
		  corrTSTD[ RecoPABin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2);
		  // cout << "2:" << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2);
		} else {
		  cout << "BAD: nEres " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " CorrTSTD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corrTBins[ RecoPABin ] / (double)corrTBinN[ RecoPABin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
		}

		if( pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2) < 99 ) {
		  resESTD[ TrueEnergyBin ] += pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2);
		  // cout << "3:" << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2);
		} else {
		  cout << "BAD: resE " << ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) << " E STD: " << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resEBins[ TrueEnergyBin ] / (double)resEBinN[ TrueEnergyBin ]), 2) << endl;
		}
		
		if( pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2) < 99 ) {
		  resTSTD[ TruePABin ] += pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2);
		  // cout << "4:" << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2);
		} else {
		  cout << "BAD: resT " << ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) << " T STD: " << pow( ((jene_RRA.At(recoMatch) - jene_TRA.At(i)) / jene_TRA.At(i)) - (resTBins[ TruePABin ] / (double)resTBinN[ TruePABin ]), 2) << endl;
		}
		
		if( pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2) < 99 ) {
		  corr2DSTD[ RecoEnergyBin ][ RecoPABin ] += pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2);
		  // cout << "5:" << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2);
		} else {
		  cout << "BAD: nEres: " << (rcmoTemp - mcmoTemp).Mag() / jene_TRA.At(i) << " Corr2DSTD: " << pow((jene_TRA.At(i) / jene_RRA.At(recoMatch)) - (corr2DBins[ RecoEnergyBin ][ RecoPABin ] / (double)corr2DBinN[ RecoEnergyBin ][ RecoPABin ]), 2) << " corr: " << jene_TRA.At(i) / jene_RRA.At(recoMatch) << endl;
		}
	      }
	    }
	  }
	}
      }
    }

    cout << "EVTCOUNT: " << evtCount << endl;
  }
    
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

    int FDS = fitDivs * fitDivs;

    double x2D[FDS];
    double y2D[FDS];
    double z2D[FDS];
    double ex2D[FDS];
    double ey2D[FDS];
    double ez2D[FDS];

    double stepE = (highE - lowE) / (double)fitDivs;
    double halfStepE = stepE / 2.0;
    double currE = lowE + (stepE / 2.0);

    double stepT = (TMath::Pi() * 8.0 ) / (double)(9 * fitDivs);
    double halfStepT = stepT / 2;
    double currT = (TMath::Pi() / 18.0) + ((double)stepT / 2.0);

    double curr2DTBase = currT;
    double curr2DT = curr2DTBase;

    // Prints out the bins used to find the weights
    for(int i = 0; i < fitDivs; i++) {
      cout << corrEBins[i] << " " << corrEBinN[i] << " " << corrESTD[i] << endl;
      cout << corrTBins[i] << " " << corrTBinN[i] << " " << corrTSTD[i] << endl;
      cout << resEBins[i] << " " << resEBinN[i] << " " << resESTD[i] << endl;
      cout << resTBins[i] << " " << resTBinN[i] << " " << resTSTD[i] << endl;
    }

    for(int i = 0; i < fitDivs; i++) {
      for(int j = 0; j < fitDivs; j++) {
	// cout << i << " " << j << " " << corr2DBins[i][j] << " " << corr2DBinN[i][j] << " " << corr2DSTD[i][j] << endl;
      }
    }

    // Set up arrays for plotting TGraphErrors for the correction terms (for energy and polar angle)
    for(int i = 0; i < fitDivs; i++) {
      if(corrEBinN[i] != 0) {
	corrEBins[i] = corrEBins[i] / (double)corrEBinN[i];
	corrESTD[i] = sqrt(corrESTD[i]) / (double)corrEBinN[i];
      }
      if(corrTBinN[i] != 0) {
	corrTBins[i] = corrTBins[i] / (double)corrTBinN[i];
	corrTSTD[i] = sqrt(corrTSTD[i]) / (double)corrTBinN[i];
      }
      if(resEBinN[i] != 0) {
        resEBins[i] = resEBins[i] / (double)resEBinN[i];
        resESTD[i] = sqrt(resESTD[i]) / (double)resEBinN[i];
      }
      if(resTBinN[i] != 0) {
        resTBins[i] = resTBins[i] / (double)resTBinN[i];
        resTSTD[i] = sqrt(resTSTD[i]) / (double)resTBinN[i];
      }
      
      xE[i] = currE;
      exE[i] = halfStepE;
      eyE[i] = 0;

      xT[i] = currT;
      exT[i] = halfStepT;
      eyT[i] = 0;

      for(int j = 0; j < fitDivs; j++) {
        if(corr2DBinN[i][j] != 0) {
          corr2DBins[i][j] = corr2DBins[i][j] / (double)corr2DBinN[i][j];
          corr2DSTD[i][j] = corr2DSTD[i][j] / sqrt((double)corr2DBinN[i][j]);
        }
	int index = (fitDivs * i) + j;
	x2D[index] = currE;
	y2D[index] = curr2DTBase + (stepT * j);
	z2D[index] = corr2DBins[i][j];
	ex2D[index] = halfStepE;
        ey2D[index] = halfStepT;
	ez2D[index] = corr2DSTD[i][j];
      }
      
      currE += stepE;
      currT += stepT;
    }

    for(int i = 0; i < fitDivs; i++) {
      cout << corrESTD[i] << " ";
    }
    cout << endl;
    for(int i = 0; i < fitDivs; i++) {
      cout << corrTSTD[i] << " ";
    }
    cout << endl;
    
    for(int i = 0; i < fitDivs; i++) {
      cout << resESTD[i] << " ";
    }
    cout << endl;
    for(int i = 0; i < fitDivs; i++) {
      cout << resTSTD[i] << " ";
    }
    cout << endl;

    // 1D fit functions for the energy correction factors
    TF1 *CorrTermEFit_0 = new TF1("CorrTermEFit_0", "([0]*pow(e, -x*[1])) + [2]", 0.0, 550.0);
    TF1 *CorrTermEFit_1 = new TF1("CorrTermEFit_1", "([0]*pow(e, -x*[1])) + [2]", 550, 900);
    TF1 *CorrTermEFit_2 = new TF1("CorrTermEFit_2", "([0]*pow(e, -x*[1])) + [2]", 900.0, 2700.0);
    TF1 *CorrTermEFit_3 = new TF1("CorrTermEFit_3", "[0]*x + [1]", 2700.0, 5000.0);
    
    // 1D fit functions for the polar angle correction factors
    TF1 *CorrTermTFit_0 = new TF1("CorrTermTFit_0", "([0]*cos(( x - [1]) * [2]) + [3])", (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    TF1 *CorrTermTFit_1 = new TF1("CorrTermTFit_1", "([0]*cos(( x - (3.1415926535897932384626433832795028 / 2) ) * [1]) + [2])", atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    TF1 *CorrTermTFit_2 = new TF1("CorrTermTFit_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - x ) * [2]) + [3])", TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    //2D fit functions for the correction factors 
    // TF2 *Corr2DFitFunc_0 = new TF2("Corr2DFitFunc_0", "([0]*cos(( y - [1]) * [2]) + [3]) + sqrt((([4]/sqrt(x)) * ([4]/sqrt(x))) + (([5] / x) * ([5] / x)) + ([6] * [6]))", lowE, highE, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    // TF2 *Corr2DFitFunc_1 = new TF2("Corr2DFitFunc_1", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) + [2]) + sqrt((([3]/sqrt(x)) * ([3]/sqrt(x))) + (([4] / x) * ([4] / x)) + ([5] * [5]))", lowE, highE, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    // TF2 *Corr2DFitFunc_2 = new TF2("Corr2DFitFunc_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) + [3]) + sqrt((([4]/sqrt(x)) * ([4]/sqrt(x))) + (([5] / x) * ([5] / x)) + ([6] * [6]))", lowE, highE, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);
    
    // Corr2DFitFunc_0 = new TF2("Corr2DFitFunc_0", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    // Corr2DFitFunc_1 = new TF2("Corr2DFitFunc_1", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(TMath::E(), -x*[3]) + [4]", lowE + 100, highE - 100, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    // Corr2DFitFunc_2 = new TF2("Corr2DFitFunc_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_0 = new TF2("Corr2DFitFunc_0", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_1 = new TF2("Corr2DFitFunc_1", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(TMath::E(), -x*[3]) + [4]", lowE + 100, highE - 100, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_2 = new TF2("Corr2DFitFunc_2", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(TMath::E(), -x*[4]) + [5]", lowE + 100, highE - 100, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_00 = new TF2("Corr2DFitFunc_00", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", lowE, 550.0, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_01 = new TF2("Corr2DFitFunc_01", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(e, -x*[3]) + [4]", lowE, 550.0, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_02 = new TF2("Corr2DFitFunc_02", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", lowE, 550.0, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_10 = new TF2("Corr2DFitFunc_10", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", 550.0, 900.0, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_11 = new TF2("Corr2DFitFunc_11", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(e, -x*[3]) + [4]", 550.0, 900.0, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_12 = new TF2("Corr2DFitFunc_12", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", 550.0, 900.0, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_20 = new TF2("Corr2DFitFunc_20", "([0]*cos(( y - [1]) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", 990.0, 2700.0, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_21 = new TF2("Corr2DFitFunc_21", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*pow(e, -x*[3]) + [4]", 900.0, 2700.0, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_22 = new TF2("Corr2DFitFunc_22", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*pow(e, -x*[4]) + [5]", 900.0, 2700.0, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    Corr2DFitFunc_30 = new TF2("Corr2DFitFunc_30", "([0]*cos(( y - [1]) * [2]) ) + [3]*x + [4]", 2700.0, highE, (10.0 * TMath::Pi() / 180.0) + 0, atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_31 = new TF2("Corr2DFitFunc_31", "([0]*cos(( y - (3.1415926535897932384626433832795028 / 2) ) * [1]) ) + [2]*x + [3]", 2700.0, highE, atan(1500.0/2210.0) + 0, TMath::Pi() - atan(1500.0/2210.0) - 0);
    Corr2DFitFunc_32 = new TF2("Corr2DFitFunc_32", "([0]*cos(( 3.1415926535897932384626433832795028 - [1] - y ) * [2]) ) + [3]*x + [4]", 2700.0, highE, TMath::Pi() - atan(1500.0/2210.0) + 0, TMath::Pi() - (10.0 * TMath::Pi() / 180.0) - 0);

    // TF1 *spQuadFit = new TF1("spQuadFit", "sqrt(([0] / sqrt(x))**2 + ([1] / x)**2 + ([2])**2)", lowE, highE);
    TF1 *spQuadFit = new TF1("spQuadFit", "[0]*pow(TMath::E(), -x * [1]) + [2]", lowE + 100, highE - 100);

    // Set initial 1D energy fit parameters
    CorrTermEFit_0->SetParameters(0.3, -0.0015, 1.7);
    CorrTermEFit_1->SetParameters(0.3, 0.0025, 1.5);
    CorrTermEFit_2->SetParameters(2.8, 0.0003, 0.25);
    CorrTermEFit_3->SetParameters(-0.0002, 2.04);
    
    // Set initial 1D polar angle fit parameters
    CorrTermTFit_0->SetParameters(0.07, 0.2, 1.6, 0.89);
    CorrTermTFit_1->SetParameters(-0.0065, 3.4, 0.92);
    CorrTermTFit_2->SetParameters(0.07, 0.2, 1.6, 0.89);

    // Set initial 2D fit parameters
    // Corr2DFitFunc_0->SetParameters(0.434114, 0.566458, 1.10895, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_1->SetParameters(-0.236795, 0.823947, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_2->SetParameters(0.45647, 0.610311, 1.20088, -0.000120153, -0.00162367, 1.30008);

    // Set initial 2D fit parameters
    // Corr2DFitFunc_0->SetParameters(1.21372, 0.389307, 3.156, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_1->SetParameters(0.0686716, 4.31098, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_2->SetParameters(3.85189, 0.45284, 1.35512, -0.000120153, -0.00162367, 1.30008);
    
    spQuadFit->SetParameters(0.016, 0.001, 0.9552);

    CorrTerm2D_0 = new TGraph2DErrors(FDS, x2D, y2D, z2D, ex2D, ey2D, ez2D);

    CorrTerm2D_0->SetMinimum(K2CorrMin);
    CorrTerm2D_0->SetMaximum(K2CorrMax);
    
    CorrTerm2D_0->SetMarkerStyle(20);
    CorrTerm2D_0->Draw("pcol");
    CorrTerm2D_0->SetTitle("Correction Terms w.r.t. Energy And Polar Angle");
    CorrTerm2D_0->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    CorrTerm2D_0->GetYaxis()->SetTitle("Polar Angle (Rads)");
    CorrTerm2D_0->GetZaxis()->SetTitle("Correction Factor");

    K2Corr2D_0->GetXaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetYaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetZaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    // Set initial 2D fit parameters
    // Corr2DFitFunc_0->SetParameters(0.434114, 0.566458, 1.10895, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_1->SetParameters(-0.236795, 0.823947, -0.000120153, -0.00162367, 1.30008);
    // Corr2DFitFunc_2->SetParameters(0.45647, 0.610311, 1.20088, -0.000120153, -0.00162367, 1.30008);

    Corr2DFitFunc_0->SetParameters(1.21372, 0.389307, 3.156, -0.000120153, -0.00162367, 1.30008);
    Corr2DFitFunc_1->SetParameters(0.0686716, 4.31098, -0.000120153, -0.00162367, 1.30008);
    Corr2DFitFunc_2->SetParameters(3.85189, 0.45284, 1.35512, -0.000120153, -0.00162367, 1.30008);

    Corr2DFitFunc_00->SetParameters(2.09325, 0.437224, 1.69391, 0.786146, 0.00752194, 1.81893);
    Corr2DFitFunc_01->SetParameters(0.027288, 5.09138, 0.786146, 0.00752194, 1.81893);
    Corr2DFitFunc_02->SetParameters(0.204818, 0.467179, 5.94263, 0.786146, 0.00752194, 1.81893);
    Corr2DFitFunc_10->SetParameters(2.09325, 0.437224, 1.69391, 0.000040101, -0.011017, 1.76968);
    Corr2DFitFunc_11->SetParameters(0.027288, 5.09138, 0.000040101, -0.011017, 1.76968);
    Corr2DFitFunc_12->SetParameters(0.204818, 0.467179, 5.94263, 0.000040101, -0.011017, 1.76968);
    Corr2DFitFunc_20->SetParameters(2.09325, 0.437224, 1.69391, 2.52256, 0.000453489, 0.734863);
    Corr2DFitFunc_21->SetParameters(0.027288, 5.09138, 2.52256, 0.000453489, 0.734863);
    Corr2DFitFunc_22->SetParameters(0.204818, 0.467179, 5.94263, 2.52256, 0.000453489, 0.734863);
    Corr2DFitFunc_30->SetParameters(2.09325, 0.437224, 1.69391, -0.000228342, 0.00778522);
    Corr2DFitFunc_31->SetParameters(0.027288, 5.09138, -0.000228342, 0.00778522);
    Corr2DFitFunc_32->SetParameters(0.204818, 0.467179, 5.94263, -0.000228342, 0.00778522);
    
    // Fit 2D graph
    CorrTerm2D_0->Fit("Corr2DFitFunc_0", "R WW 0");
    CorrTerm2D_0->Fit("Corr2DFitFunc_1", "R WW 0");
    CorrTerm2D_0->Fit("Corr2DFitFunc_2", "R WW 0");

    cout << "00" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_00", "R WW 0");
    cout << "01" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_01", "R WW 0");
    cout << "02" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_02", "R WW 0");
    cout << "10" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_10", "R WW 0");
    cout << "11" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_11", "R WW 0");
    cout << "12" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_12", "R WW 0");
    cout << "20" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_20", "R WW 0");
    cout << "21" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_21", "R WW 0");
    cout << "22" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_22", "R WW 0");
    cout << "30" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_30", "R WW 0");
    cout << "31" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_31", "R WW 0");
    cout << "32" << endl;
    CorrTerm2D_0->Fit("Corr2DFitFunc_32", "R WW 0");

    CorrTerm2D_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D.png");
    c->Close();
    c = new TCanvas();

    CorrTerm2D_0->Project("X")->Draw();

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D_X.png");
    c->Close();
    c = new TCanvas();

    CorrTerm2D_0->Project("Y")->Draw();

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D_Y.png");
    c->Close();
    c = new TCanvas();

    CorrTerm2D_0->Project("XY")->Draw();

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D_XY.png");
    c->Close();
    c = new TCanvas();

    CorrTerm2D_0->Project("ZX")->Draw();

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D_ZX.png");
    c->Close();
    c = new TCanvas();

    CorrTerm2D_0->Project("ZY")->Draw();

    c->SaveAs(saveDir + "/phaseHists/CorrTerm2D_ZY.png");
    c->Close();
    c = new TCanvas();

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
      // if( recoPAArr[i] < atan(1500.0/2210.0) ) {
      // 	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_0->Eval(recoEnergyArr[i], recoPAArr[i]);
      // } else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
      // 	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_1->Eval(recoEnergyArr[i], recoPAArr[i]);
      // } else {
      // 	recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_2->Eval(recoEnergyArr[i], recoPAArr[i]);
      // }
      if( recoEnergyArr[i] < 550 ) {
	if( recoPAArr[i] < atan(1500.0/2210.0) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_00->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_01->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_02->Eval(recoEnergyArr[i], recoPAArr[i]);
	}
      } else if( recoEnergyArr[i] < 900 ) {
	if( recoPAArr[i] < atan(1500.0/2210.0) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_10->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_11->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_12->Eval(recoEnergyArr[i], recoPAArr[i]);
	}
      } else if( recoEnergyArr[i] < 2700 ) {
	if( recoPAArr[i] < atan(1500.0/2210.0) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_20->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_21->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_22->Eval(recoEnergyArr[i], recoPAArr[i]);
	}
      } else {
	if( recoPAArr[i] < atan(1500.0/2210.0) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_30->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else if( recoPAArr[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_31->Eval(recoEnergyArr[i], recoPAArr[i]);
	} else {
	  recoEnergyArr_Fit[i] = recoEnergyArr[i] * Corr2DFitFunc_32->Eval(recoEnergyArr[i], recoPAArr[i]);
	}
      }
      CorrE_Fit_0->Fill(recoEnergyArr_Fit[i]);
      CorrERes_Fit_0->Fill((recoEnergyArr_Fit[i] - trueEnergyArr[i]) / trueEnergyArr[i]);

      if(BMD) {
	if(recoEnergyArr_Bins[i] < bimodSplit) {
	  CorrE_BinsBimod0->Fill(recoEnergyArr_Bins[i]);
	  CorrE_FitBimod0->Fill(recoEnergyArr_Fit[i]);
	  CorrERes_BinsBimod0->Fill((recoEnergyArr_Bins[i] - trueEnergyArr[i]) / trueEnergyArr[i]);
	  CorrERes_FitBimod0->Fill((recoEnergyArr_Fit[i] - trueEnergyArr[i]) / trueEnergyArr[i]);
	} else {
	  CorrE_BinsBimod1->Fill(recoEnergyArr_Bins[i]);
	  CorrE_FitBimod1->Fill(recoEnergyArr_Fit[i]);
	  CorrERes_BinsBimod1->Fill((recoEnergyArr_Bins[i] - trueEnergyArr[i]) / trueEnergyArr[i]);
	  CorrERes_FitBimod1->Fill((recoEnergyArr_Fit[i] - trueEnergyArr[i]) / trueEnergyArr[i]);
	}
      }
    }

    for(int i = 0; i < recoIVMCount; i++) {
      
      if( recoIVMEArr0[i] < 550 ) {
        if( recoIVMTArr0[i] < atan(1500.0/2210.0) ) {
          recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_00->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else if( recoIVMTArr0[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_00->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_00->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        }
      } else if( recoIVMEArr0[i] < 900 ) {
        if( recoIVMTArr0[i] < atan(1500.0/2210.0) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_10->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else if( recoIVMTArr0[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_11->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_12->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        }
      } else if( recoIVMEArr0[i] < 2700 ) {
        if( recoIVMTArr0[i] < atan(1500.0/2210.0) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_20->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else if( recoIVMTArr0[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_21->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_22->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        }
      } else {
        if( recoIVMTArr0[i] < atan(1500.0/2210.0) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_30->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else if( recoIVMTArr0[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_31->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        } else {
	  recoIVMEArr0[i] = recoIVMEArr0[i] * Corr2DFitFunc_32->Eval(recoIVMEArr0[i], recoIVMTArr0[i]);
        }
      }

      if( recoIVMEArr1[i] < 550 ) {
        if( recoIVMTArr1[i] < atan(1500.0/2210.0) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_00->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else if( recoIVMTArr1[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_00->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_00->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        }
      } else if( recoIVMEArr1[i] < 900 ) {
        if( recoIVMTArr1[i] < atan(1500.0/2210.0) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_10->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else if( recoIVMTArr1[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_11->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_12->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        }
      } else if( recoIVMEArr1[i] < 2700 ) {
        if( recoIVMTArr1[i] < atan(1500.0/2210.0) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_20->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else if( recoIVMTArr1[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_21->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_22->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        }
      } else {
        if( recoIVMTArr1[i] < atan(1500.0/2210.0) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_30->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else if( recoIVMTArr1[i] < (TMath::Pi() - atan(1500.0/2210.0)) ) {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_31->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        } else {
          recoIVMEArr1[i] = recoIVMEArr1[i] * Corr2DFitFunc_32->Eval(recoIVMEArr1[i], recoIVMTArr1[i]);
        }
      }
      
      double addE = pow(recoIVMEArr0[i] + recoIVMEArr1[i], 2);
      double jvx = pow(recoIVMxArr0[i] + recoIVMxArr1[i], 2);
      double jvy = pow(recoIVMyArr0[i] + recoIVMyArr1[i], 2);
      double jvz = pow(recoIVMzArr0[i] + recoIVMzArr1[i], 2);
      
      CorrIVM->Fill( sqrt(addE - jvx - jvy - jvz) );
    }

    double pyCEM[fitDivs];
    double eyCEM[fitDivs];
    
    double pyCTM[fitDivs];
    double eyCTM[fitDivs];

    double pyCES[fitDivs];
    double eyCES[fitDivs];

    double pyCTS[fitDivs];
    double eyCTS[fitDivs];

    for(int i = 0; i < fitDivs; i++) {
      TF1 *iGaussFit = new TF1("iGaussFit", "gaus", -1, 1);
      
      char s[4];
      char outE[100];
      char outT[100];

      if(EP) {
	sprintf(s, "%d", i);
	strcpy(outE, saveDir);
	strcat(outE, "/seperateCorr/CET_");
	strcat(outE, s);
	strcat(outE, ".png");
	strcpy(outT, saveDir);
	strcat(outT, "/seperateCorr/CTT_");
	strcat(outT, s);
	strcat(outT, ".png");
      }
      
      iGaussFit->SetParameters(200, 2.0, 0.05);
      
      CorrETerms[i]->Fit("iGaussFit");
      
      pyCES[i] = iGaussFit->GetParameter(2);
      eyCES[i] = iGaussFit->GetParError(2);
      
      pyCEM[i] = iGaussFit->GetParameter(1);
      eyCEM[i] = iGaussFit->GetParError(1);
      
      if(EP) {
        CorrETerms[i]->Draw();
        CorrETerms[i]->SetTitle("Energy Correction Terms - E");
        CorrETerms[i]->GetXaxis()->SetTitle("E_{True} / E_{Reco}");
        CorrETerms[i]->GetYaxis()->SetTitle("Count");
	
	CorrETerms[i]->Fit("gaus");
	
        c->SaveAs(outE);
        c->Close();
        c = new TCanvas();
      }

      iGaussFit->SetParameters(200, 2.0, 0.05);
      
      CorrTTerms[i]->Fit("iGaussFit");
      
      pyCTS[i] = iGaussFit->GetParameter(2);
      eyCTS[i] = iGaussFit->GetParError(2);
      
      pyCTM[i] = iGaussFit->GetParameter(1);
      eyCTM[i] = iGaussFit->GetParError(1);

      if(EP) {
        CorrTTerms[i]->Draw();
        CorrTTerms[i]->SetTitle("Energy Correction Terms - PA");
        CorrTTerms[i]->GetXaxis()->SetTitle("E_{True} / E_{Reco}");
        CorrTTerms[i]->GetYaxis()->SetTitle("Count");
	
	CorrTTerms[i]->Fit("gaus");
	
        c->SaveAs(outT);
        c->Close();
        c = new TCanvas();
      }
    }
      
    // Draw and save graphs

    TGraphErrors *CorrEFitMeans = new TGraphErrors(fitDivs, xT, pyCEM, exT, eyCEM);

    CorrEFitMeans->Draw();
    CorrEFitMeans->SetTitle("Energy Correction Fit Means -- Energy");
    CorrEFitMeans->GetXaxis()->SetTitle("Energy_{True} (GeV)");
    CorrEFitMeans->GetYaxis()->SetTitle("Correction Factor");

    CorrEFitMeans->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/CEFitMeans.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *CorrTFitMeans = new TGraphErrors(fitDivs, xT, pyCTM, exT, eyCTM);

    CorrTFitMeans->Draw();
    CorrTFitMeans->SetTitle("Energy Correction Fit Means -- Polar Angle");
    CorrTFitMeans->GetXaxis()->SetTitle("Polar Angle_{True} (GeV)");
    CorrTFitMeans->GetYaxis()->SetTitle("Correction Factor");

    CorrTFitMeans->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/CTFitMeans.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *CorrEFitSigmas = new TGraphErrors(fitDivs, xT, pyCES, exT, eyCES);

    CorrEFitSigmas->Draw();
    CorrEFitSigmas->SetTitle("Energy Correction Fit Sigmas -- Energy");
    CorrEFitSigmas->GetXaxis()->SetTitle("Energy_{True} (GeV)");
    CorrEFitSigmas->GetYaxis()->SetTitle("Correction Factor");

    CorrEFitSigmas->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/CEFitSigmas.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *CorrTFitSigmas = new TGraphErrors(fitDivs, xT, pyCTS, exT, eyCTS);

    CorrTFitSigmas->Draw();
    CorrTFitSigmas->SetTitle("Energy Correction Fit Sigmas -- Polar Angle");
    CorrTFitSigmas->GetXaxis()->SetTitle("Polar Angle_{True} (GeV)");
    CorrTFitSigmas->GetYaxis()->SetTitle("Correction Factor");

    CorrTFitSigmas->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/CTFitSigmas.png");
    c->Close();
    c = new TCanvas();

    CorrIVM->Draw();
    CorrIVM->SetTitle("Corrected Reconstructed Invariant Mass");
    CorrIVM->GetXaxis()->SetTitle("Invariant Mass (GeV/ c^2)");
    CorrIVM->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/CorrIVM.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *TruthERes_0 = new TGraphErrors(fitDivs, xE, resEBins, exE, resESTD);

    TruthERes_0->Draw();
    TruthERes_0->SetTitle("Energy Resolution w.r.t. Energy");
    TruthERes_0->GetXaxis()->SetTitle("Energy_{True} (GeV)");
    TruthERes_0->GetYaxis()->SetTitle("Correction Factor");

    TruthERes_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/TruthERes.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *TruthTRes_0 = new TGraphErrors(fitDivs, xT, resTBins, exT, resTSTD);

    TruthTRes_0->Draw();
    TruthTRes_0->SetTitle("Energy Resolution w.r.t. Polar Angle");
    TruthTRes_0->GetXaxis()->SetTitle("Polar Angle_{True} (GeV)");
    TruthTRes_0->GetYaxis()->SetTitle("Correction Factor");

    TruthTRes_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/TruthTRes.png");
    c->Close();
    c = new TCanvas();
    
    TGraphErrors *CorrTermE_0 = new TGraphErrors(fitDivs, xE, corrEBins, exE, corrESTD);
    
    CorrTermE_0->Draw();
    CorrTermE_0->SetTitle("Correction Terms w.r.t. Energy");
    CorrTermE_0->GetXaxis()->SetTitle("Reconstructed Energy (GeV)");
    CorrTermE_0->GetYaxis()->SetTitle("Correction Factor");

    CorrTermE_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    CorrTermE_0->Fit("CorrTermEFit_0", "WW R");
    CorrTermE_0->Fit("CorrTermEFit_1", "WW R");
    CorrTermE_0->Fit("CorrTermEFit_2", "WW R");
    CorrTermE_0->Fit("CorrTermEFit_3", "WW R");

    CorrTermEFit_0->Draw("SAME");
    CorrTermEFit_1->Draw("SAME");
    CorrTermEFit_2->Draw("SAME");
    
    // CorrTermE_0->Fit("spQuadFit", "RW");
    // spQuadFit->Draw("SAME");

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
    
    CorrTermT_0->Fit("CorrTermTFit_0", "WW R");
    CorrTermT_0->Fit("CorrTermTFit_1", "WW R");
    CorrTermT_0->Fit("CorrTermTFit_2", "WW R");
    
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

    if(BMD) {
      CorrE_BinsBimod0->SetLineColor(kMagenta);
      CorrE_BinsBimod1->SetLineColor(kBlack);
      
      CorrE_BinsBimod0->Draw("SAME");
      CorrE_BinsBimod1->Draw("SAME");
      
      CorrE_BinsBimod0->Fit("gaus");
      CorrE_BinsBimod1->Fit("gaus");
    }

    c->SaveAs(saveDir + "/simpleHists/CorrE_Bins.png");
    c->Close();
    c = new TCanvas();

    CorrERes_Bins_0->Draw();
    CorrERes_Bins_0->SetTitle("Corrected reconstructed energy resolution -- Bins");
    CorrERes_Bins_0->GetXaxis()->SetTitle("Energy resolution (GeV)");
    CorrERes_Bins_0->GetYaxis()->SetTitle("Count");

    CorrERes_Bins_0->Fit("gaus");

    if(BMD) {
      CorrERes_BinsBimod0->SetLineColor(kMagenta);
      CorrERes_BinsBimod1->SetLineColor(kBlack);
      
      CorrERes_BinsBimod0->Draw("SAME");
      CorrERes_BinsBimod1->Draw("SAME");
      
      CorrERes_BinsBimod0->Fit("gaus");
      CorrERes_BinsBimod1->Fit("gaus");
    }

    c->SaveAs(saveDir + "/resolutionHists/CorrERes_Bins.png");
    c->Close();
    c = new TCanvas();

    CorrE_Fit_0->Draw();
    CorrE_Fit_0->SetTitle("Corrected reconstructed energy -- Fit");
    CorrE_Fit_0->GetXaxis()->SetTitle("Corrected Energy (GeV)");
    CorrE_Fit_0->GetYaxis()->SetTitle("Count");

    if(BMD) {
      CorrE_FitBimod0->SetLineColor(kMagenta);
      CorrE_FitBimod1->SetLineColor(kBlack);
      
      CorrE_FitBimod0->Draw("SAME");
      CorrE_FitBimod1->Draw("SAME");
      
      // CorrE_FitBimod0->Fit("gaus");
      // CorrE_FitBimod1->Fit("gaus");
    }

    c->SaveAs(saveDir + "/simpleHists/CorrE_Fit.png");
    c->Close();
    c = new TCanvas();

    CorrERes_Fit_0->Draw();
    CorrERes_Fit_0->SetTitle("Corrected reconstructed energy accuracy -- Fit");
    CorrERes_Fit_0->GetXaxis()->SetTitle("(E_{Reco} - E_{Truth}) / E_{Truth}");
    CorrERes_Fit_0->GetYaxis()->SetTitle("Count");

    // CorrERes_Fit_0->Fit("gaus");

    if(BMD) {
      CorrERes_FitBimod0->SetLineColor(kMagenta);
      CorrERes_FitBimod1->SetLineColor(kBlack);
      
      CorrERes_FitBimod0->Draw("SAME");
      CorrERes_FitBimod1->Draw("SAME");
      
      // CorrERes_FitBimod0->Fit("gaus");
      // CorrERes_FitBimod1->Fit("gaus");
    }

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

    
    TGraph2D *C2DFFplot_00 = new TGraph2D((TH2D*)Corr2DFitFunc_00->GetHistogram());

    C2DFFplot_00->GetYaxis()->SetRange(0, 4);
    C2DFFplot_00->Draw("surf1");
    C2DFFplot_00->SetTitle("Reco Energy Correction Factors 2D Fit -- 0-550 GeV, 0-34.166 deg");
    C2DFFplot_00->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_00->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_00->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_00.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_01 = new TGraph2D((TH2D*)Corr2DFitFunc_01->GetHistogram());

    C2DFFplot_01->GetYaxis()->SetRange(0, 4);
    C2DFFplot_01->Draw("surf1");
    C2DFFplot_01->SetTitle("Reco Energy Correction Factors 2D Fit -- 0-550 GeV, 34.167-145.833 deg");
    C2DFFplot_01->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_01->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_01->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_01.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_02 = new TGraph2D((TH2D*)Corr2DFitFunc_02->GetHistogram());

    C2DFFplot_02->GetYaxis()->SetRange(0, 4);
    C2DFFplot_02->Draw("surf1");
    C2DFFplot_02->SetTitle("Reco Energy Correction Factors 2D Fit -- 0-550 GeV, 145.834-180 deg");
    C2DFFplot_02->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_02->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_02->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_02.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_10 = new TGraph2D((TH2D*)Corr2DFitFunc_10->GetHistogram());

    C2DFFplot_10->GetYaxis()->SetRange(0, 4);
    C2DFFplot_10->Draw("surf1");
    C2DFFplot_10->SetTitle("Reco Energy Correction Factors 2D Fit -- 550-900 GeV, 0-34.166 deg");
    C2DFFplot_10->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_10->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_10->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_10.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_11 = new TGraph2D((TH2D*)Corr2DFitFunc_11->GetHistogram());

    C2DFFplot_11->GetYaxis()->SetRange(0, 4);
    C2DFFplot_11->Draw("surf1");
    C2DFFplot_11->SetTitle("Reco Energy Correction Factors 2D Fit -- 550-900 GeV, 34.167-145.833 deg");
    C2DFFplot_11->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_11->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_11->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_11.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_12 = new TGraph2D((TH2D*)Corr2DFitFunc_12->GetHistogram());

    C2DFFplot_12->GetYaxis()->SetRange(0, 4);
    C2DFFplot_12->Draw("surf1");
    C2DFFplot_12->SetTitle("Reco Energy Correction Factors 2D Fit -- 550-900 GeV, 145.834-180 deg");
    C2DFFplot_12->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_12->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_12->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_12.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_20 = new TGraph2D((TH2D*)Corr2DFitFunc_20->GetHistogram());

    C2DFFplot_20->GetYaxis()->SetRange(0, 4);
    C2DFFplot_20->Draw("surf1");
    C2DFFplot_20->SetTitle("Reco Energy Correction Factors 2D Fit -- 900-2700 GeV, 0-34.166 deg");
    C2DFFplot_20->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_20->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_20->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_20.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_21 = new TGraph2D((TH2D*)Corr2DFitFunc_21->GetHistogram());

    C2DFFplot_21->GetYaxis()->SetRange(0, 4);
    C2DFFplot_21->Draw("surf1");
    C2DFFplot_21->SetTitle("Reco Energy Correction Factors 2D Fit -- 900-2700 GeV, 34.167-145.833 deg");
    C2DFFplot_21->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_21->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_21->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_21.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_22 = new TGraph2D((TH2D*)Corr2DFitFunc_22->GetHistogram());

    C2DFFplot_22->GetYaxis()->SetRange(0, 4);
    C2DFFplot_22->Draw("surf1");
    C2DFFplot_22->SetTitle("Reco Energy Correction Factors 2D Fit -- 900-2700 GeV, 145.834-180 deg");
    C2DFFplot_22->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_22->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_22->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_22.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_30 = new TGraph2D((TH2D*)Corr2DFitFunc_30->GetHistogram());

    C2DFFplot_30->GetYaxis()->SetRange(0, 4);
    C2DFFplot_30->Draw("surf1");
    C2DFFplot_30->SetTitle("Reco Energy Correction Factors 2D Fit -- 2700-5000 GeV, 0-34.166 deg");
    C2DFFplot_30->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_30->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_30->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_30.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_31 = new TGraph2D((TH2D*)Corr2DFitFunc_31->GetHistogram());

    C2DFFplot_31->GetYaxis()->SetRange(0, 4);
    C2DFFplot_31->Draw("surf1");
    C2DFFplot_31->SetTitle("Reco Energy Correction Factors 2D Fit -- 2700-5000 GeV, 34.167-145.833 deg");
    C2DFFplot_31->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_31->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_31->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_31.png");
    c->Close();
    c = new TCanvas();

    TGraph2D *C2DFFplot_32 = new TGraph2D((TH2D*)Corr2DFitFunc_32->GetHistogram());

    C2DFFplot_32->GetYaxis()->SetRange(0, 4);
    C2DFFplot_32->Draw("surf1");
    C2DFFplot_32->SetTitle("Reco Energy Correction Factors 2D Fit -- 2700-5000 GeV, 145.834-180 deg");
    C2DFFplot_32->GetXaxis()->SetTitle("Reco Energy (GeV)");
    C2DFFplot_32->GetYaxis()->SetTitle("Reco Polar Angle (Rads)");
    C2DFFplot_32->GetZaxis()->SetTitle("Correction Factor");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2DFit_32.png");
    c->Close();
    c = new TCanvas();
  }

  // cout << "Test4" << endl;

  if(ECF || rFH || ECA) {

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

      if(minTruthE <= jene_TRA.At(i)) {
	
	if(matched) {
	  recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	  
	  rcmoTemp = recoVec3[recoMatch];
	  
	  Double_t tempE;
	  // correct energy
	  if( rcmoTemp.Theta() < atan(1500.0/2210.0) ) {
	    tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_0->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	  } else if( rcmoTemp.Theta() < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	    tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_1->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	  } else {
	    tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_2->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	  }
	  
	  int RecoEnergyBin = (int) floor((tempE - lowE) * fitDivs / (highE - lowE));
	  int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	  
	  while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	  while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	  while(RecoPABin < 0) { RecoPABin++; }
	  while(RecoPABin >= fitDivs) { RecoPABin--; }
	  
	  // fill in hist array
	  corrSigFitEHists[ RecoEnergyBin ]->Fill((tempE - jene_TRA.At(i)) / jene_TRA.At(i));
	  corrSigFitTHists[ RecoPABin ]->Fill((tempE - jene_TRA.At(i)) / jene_TRA.At(i));
	  
	}
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

	if(minTruthE <= jene_TRA.At(i)) {
	  
	  if(matched) {
	    recoMatch = matchArrR[indexOf(matchArrT, sizeof(matchArrT), i)];
	    
	    rcmoTemp = recoVec3[recoMatch];
	    
	    Double_t tempE;
	    // correct energy
	    if( rcmoTemp.Theta() < atan(1500.0/2210.0) ) {
	      tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_0->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	    } else if( rcmoTemp.Theta() < (TMath::Pi() - atan(1500.0/2210.0)) ) {
	      tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_1->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	    } else {
	      tempE = jene_RRA.At(recoMatch) * Corr2DFitFunc_2->Eval(jene_RRA.At(recoMatch), rcmoTemp.Theta());
	    }
	    
	    int RecoEnergyBin = (int) floor((tempE - lowE) * fitDivs / (highE - lowE));
	    int RecoPABin = (int) floor((rcmoTemp.Theta() - (TMath::Pi() / 18)) * fitDivs / (TMath::Pi() * 8 / 9));
	    
	    while(RecoEnergyBin < 0) { RecoEnergyBin++; }
	    while(RecoEnergyBin >= fitDivs) { RecoEnergyBin--; }
	    while(RecoPABin < 0) { RecoPABin++; }
	    while(RecoPABin >= fitDivs) { RecoPABin--; }
	    
	    // fill in hist array
	    corrSigFitEHists[ RecoEnergyBin ]->Fill((tempE - jene_TRA.At(i)) / jene_TRA.At(i));
	    corrSigFitTHists[ RecoPABin ]->Fill((tempE - jene_TRA.At(i)) / jene_TRA.At(i));

	  }
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

    double pyBE[fitDivs];
    double eyBE[fitDivs];

    double pyBT[fitDivs];
    double eyBT[fitDivs];
      
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

      iGaussFit->SetParameters(200, -0.3, 0.02);
      
      sigFitEHists[i]->Fit("iGaussFit");

      if(ECF) {
        pxE[i] = currE;
        pyE[i] = iGaussFit->GetParameter(2);
        exE[i] = stepE / 2;
        eyE[i] = iGaussFit->GetParError(2);

	pyBE[i] = 1 / (1 + iGaussFit->GetParameter(1));
	eyBE[i] = iGaussFit->GetParError(1);
      }
	
      if(rFH) {
	sigFitEHists[i]->Draw();
	sigFitEHists[i]->SetTitle("Energy resolution - E");
	sigFitEHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
	sigFitEHists[i]->GetYaxis()->SetTitle("Count");

	if(BMD) {
	  sigFitEHistsBimod0[i]->SetLineColor(kMagenta);
	  sigFitEHistsBimod1[i]->SetLineColor(kBlack);
	  
	  sigFitEHistsBimod0[i]->Draw("SAME");
	  sigFitEHistsBimod1[i]->Draw("SAME");
	  
	  sigFitEHistsBimod0[i]->Fit("gaus");
	  sigFitEHistsBimod1[i]->Fit("gaus");
	}

	c->SaveAs(outE);
	c->Close();
	c = new TCanvas();
      }

      iGaussFit->SetParameters(200, -0.03, 0.02);

      sigFitTHists[i]->Fit("iGaussFit");

      if(ECF) {
        pxT[i] = currT;
        pyT[i] = iGaussFit->GetParameter(2);
        exT[i] = stepT / 2;
        eyT[i] = iGaussFit->GetParError(2);

	pyBT[i] = 1 / (1 + iGaussFit->GetParameter(1));
	eyBT[i] = iGaussFit->GetParError(1);
      }

      if(rFH) {
	sigFitTHists[i]->Draw();
	sigFitTHists[i]->SetTitle("Energy resolution - PA");
	sigFitTHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
	sigFitTHists[i]->GetYaxis()->SetTitle("Count");

	if(BMD) {
	  sigFitTHistsBimod0[i]->SetLineColor(kMagenta);
	  sigFitTHistsBimod1[i]->SetLineColor(kBlack);
	  
	  sigFitTHistsBimod0[i]->Draw("SAME");
	  sigFitTHistsBimod1[i]->Draw("SAME");
	  
	  sigFitTHistsBimod0[i]->Fit("gaus");
	  sigFitTHistsBimod1[i]->Fit("gaus");
	}
	
	c->SaveAs(outT);
	c->Close();
	c = new TCanvas();
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);

      corrSigFitEHists[i]->Fit("iGaussFit");

      if(ECF) {
        pxE[i] = currE;
        pyCE[i] = iGaussFit->GetParameter(2);
        exE[i] = stepE / 2;
        eyCE[i] = iGaussFit->GetParError(2);
      }
      
      if(rFH) {
        corrSigFitEHists[i]->Draw();
        corrSigFitEHists[i]->SetTitle("Corrected Energy resolution - E");
        corrSigFitEHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
        corrSigFitEHists[i]->GetYaxis()->SetTitle("Count");

        c->SaveAs(outCE);
        c->Close();
        c = new TCanvas();
      }

      iGaussFit->SetParameters(50, 0.001, 0.02);

      corrSigFitTHists[i]->Fit("iGaussFit");

      if(ECF) {
        pxT[i] = currT;
        pyCT[i] = iGaussFit->GetParameter(2);
        exT[i] = stepT / 2;
        eyCT[i] = iGaussFit->GetParError(2);
      }
      
      if(rFH) {
        corrSigFitTHists[i]->Draw();
        corrSigFitTHists[i]->SetTitle("Corrected Energy resolution - PA");
        corrSigFitTHists[i]->GetXaxis()->SetTitle("E_{Reco} - E_{True}");
        corrSigFitTHists[i]->GetYaxis()->SetTitle("Count");

        c->SaveAs(outCT);
        c->Close();
        c = new TCanvas();
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
    KEres_err_0->GetXaxis()->SetTitle("E_{Reco} (GeV)");
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
    KTres_err_0->GetXaxis()->SetTitle("Polar angle_{Reco} (Rads)");
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

    TF1 *expFit = new TF1("expFit", "([0]*(e**(x*[1]))) + [2]", lowE, highE);

    // Set initial 1D fit parameters
    CorrTermTFit_0->SetParameters(0.005, 5.3, 0.7, 0.0097);
    CorrTermTFit_1->SetParameters(0.0035, 3.2, 0.0105);
    CorrTermTFit_2->SetParameters(0.005, 5.3, 0.7, 0.0097);

    // Set initial quad fit parameters
    expFit->SetParameters(1.0, -0.0006, 1.2);
    
    TGraphErrors *KERes_Bins = new TGraphErrors(fitDivs, pxE, pyBE, exE, eyBE);

    KERes_Bins->Draw();
    KERes_Bins->SetTitle("Error Plot Of Gauss Fit Correction Factors");
    KERes_Bins->GetXaxis()->SetTitle("E_{Reco} (GeV)");
    KERes_Bins->GetYaxis()->SetTitle("Correction Factor");
    
    KERes_Bins->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    KERes_Bins->Fit("expFit", "WMI", "SAME", lowE, highE);

    // expFit->Draw("SAME");

    c->SaveAs(saveDir + "/phaseHists/KERes_Bins.png");
    c->Close();
    c = new TCanvas();

    TGraphErrors *KTRes_Bins = new TGraphErrors(fitDivs, pxT, pyBT, exT, eyBT);

    KTRes_Bins->Draw();
    KTRes_Bins->SetTitle("Error Plot Of Gauss Fit Correction Factors");
    KTRes_Bins->GetXaxis()->SetTitle("Polar angle_{Reco} (Rads)");
    KTRes_Bins->GetYaxis()->SetTitle("Correction Factor");
    
    KTRes_Bins->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    KTRes_Bins->Fit("CorrTermTFit_0", "RW");
    KTRes_Bins->Fit("CorrTermTFit_1", "RW");
    KTRes_Bins->Fit("CorrTermTFit_2", "RW");

    CorrTermTFit_0->Draw("SAME");
    CorrTermTFit_1->Draw("SAME");
    CorrTermTFit_2->Draw("SAME");

    c->SaveAs(saveDir + "/phaseHists/KTRes_Bins.png");
    c->Close();
    c = new TCanvas();
  }

  // cout << "Test5" << endl;

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

    photonTruth_E_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/simpleHists/allTruth_E.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_E_0->Draw();
    recoPhotonTruth_E_0->SetTitle("Linked truth energies");
    recoPhotonTruth_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    recoPhotonTruth_E_0->GetYaxis()->SetTitle("Count");

    recoPhotonTruth_E_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_E.png");
    c->Close();
    c = new TCanvas();

    photonReco_E_0->Draw();
    photonReco_E_0->SetTitle("Reconstructed energies");
    photonReco_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    photonReco_E_0->GetYaxis()->SetTitle("Count");

    photonReco_E_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/simpleHists/allReco_E.png");
    c->Close();
    c = new TCanvas();

    recoPhotonReco_E_0->Draw();
    recoPhotonReco_E_0->SetTitle("Linked reconstructed energies");
    recoPhotonReco_E_0->GetXaxis()->SetTitle("Energy (GeV)");
    recoPhotonReco_E_0->GetYaxis()->SetTitle("Count");

    recoPhotonReco_E_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/simpleHists/LinkedReco_E.png");
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

    photonTruth_PA_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/simpleHists/allTruth_PA.png");
    c->Close();
    c = new TCanvas();

    recoPhotonTruth_PA_0->Draw();
    recoPhotonTruth_PA_0->SetTitle("Linked truth polar angles");
    recoPhotonTruth_PA_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    recoPhotonTruth_PA_0->GetYaxis()->SetTitle("Count");

    recoPhotonTruth_PA_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);
    
    c->SaveAs(saveDir + "/simpleHists/LinkedTruth_PA.png");
    c->Close();
    c = new TCanvas();

    photonReco_PA_0->Draw();
    photonReco_PA_0->SetTitle("Reconstructed polar angles");
    photonReco_PA_0->GetXaxis()->SetTitle("Polar angle (Rads)");
    photonReco_PA_0->GetYaxis()->SetTitle("Count");

    photonReco_PA_0->GetYaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

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

    if(BMD) {
      deltaRBimod0->Draw("SAME");
      deltaRBimod1->Draw("SAME");
      
      deltaRBimod0->Fit("gaus");
      deltaRBimod1->Fit("gaus");
      
      deltaRBimod0->SetLineColor(kMagenta);
      deltaRBimod1->SetLineColor(kBlack);
    }
      
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
    
    ERes_0->Draw();
    ERes_0->SetTitle("Energy Accuracy");
    ERes_0->GetXaxis()->SetTitle("(E_{Reco} - E_{Truth}) / E_{Truth}");
    ERes_0->GetYaxis()->SetTitle("Count");

    // ERes_0->Fit("gaus");

    if(BMD) {
      EResBimod0->Draw("SAME");
      EResBimod1->Draw("SAME");
      
      // EResBimod0->Fit("gaus");
      // EResBimod1->Fit("gaus");
      
      EResBimod0->SetLineColor(kMagenta);
      EResBimod1->SetLineColor(kBlack);
    }

    c->SaveAs(saveDir + "/resolutionHists/ERes.png");
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

    K2Corr2D_0->SetMinimum(K2CorrMin);
    K2Corr2D_0->SetMaximum(K2CorrMax);
  
    K2Corr2D_0->Draw("surf1");
    K2Corr2D_0->SetTitle("E_{True}, Polar Angle, and Energy Correction Phase Plot");
    K2Corr2D_0->GetXaxis()->SetTitle("E_{Reco} (GeV)");
    K2Corr2D_0->GetYaxis()->SetTitle("Polar Angle_{Reco} (Rads)");

    K2Corr2D_0->GetXaxis()->SetRangeUser(0, 5000);
    K2Corr2D_0->GetYaxis()->SetRangeUser(0, 3.2);
    K2Corr2D_0->GetZaxis()->SetRangeUser(K2CorrMin, K2CorrMax);

    K2Corr2D_0->GetXaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetYaxis()->SetTitleOffset(2.0);
    K2Corr2D_0->GetZaxis()->SetTitleOffset(2.0);
    c->SetLeftMargin(0.2);

    c->SaveAs(saveDir + "/phaseHists/K2Corr2D.png");
    c->Close();
    c = new TCanvas();

    K2Corr2D_0->GetHistogram()->Draw("surf1");

    c->SaveAs(saveDir + "/phaseHists/K2Corr2D-Delaunay.png");
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

  if(IVM) {

    TruthIVM->Draw();
    TruthIVM->SetTitle("Truth Invariant Mass");
    TruthIVM->GetXaxis()->SetTitle("Invariant Mass (GeV/ c^2)");
    TruthIVM->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/TruthIVM.png");
    c->Close();
    c = new TCanvas();

    RecoIVM->Draw();
    RecoIVM->SetTitle("Reco Invariant Mass");
    RecoIVM->GetXaxis()->SetTitle("Invariant Mass (GeV/ c^2)");
    RecoIVM->GetYaxis()->SetTitle("Count");

    c->SaveAs(saveDir + "/simpleHists/RecoIVM.png");
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
  zprime_jj_full_analysis_H(false, "00000000000");
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
