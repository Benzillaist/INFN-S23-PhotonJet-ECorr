void photon_energy()
{
  /////////////////////////////
  //GENERAL ANALYSIS SETTINGS//
  /////////////////////////////

  TString directory = "data/";

  TString fileName = "photonGun_ntup_3000-5000.root";
  auto treeName = "MyLCTuple";

  TString saveDir = "photonGun_3000-5000";

  ////////////////////////////
  //ANALYSIS CONFIG SETTINGS//
  ////////////////////////////

  // Max pT a reconstructed particle can be from truth pT
  float maxpTRecoErr = 0.1;

  // Number of particles
  int maxN = 1000;

  ////////////////////////////

  TString fileDir = directory + fileName;
  //TFile *myFile = new TFile("ntuple.root");
  TFile *myFile = new TFile(fileDir);
  TTree *myTree = (TTree*)myFile->Get(treeName);

  TCanvas *c = new TCanvas();

  //opens the file to be read
  auto openFile = TFile::Open(fileDir);

  // Histograms

  TH1F *photonTruth_pt_0 = new TH1F("pt_pt_0", "All Truth p_{T}", 20, 3000, 5000); //All true transverse momenta
  TH1F *recoPhotonTruth_pt_0 = new TH1F("prt_pt_0", "Reco True p_{T}", 20, 3000, 5000); //True transverse momentum of reconstructed photons
  TH1F *photonReco_pt_0 = new TH1F("pr_pt_0", "Reco p_{T}", 20, 3000, 5000); //Reconstructed transverse momentum of linked photons

  TH1F *photonTruth_PA_0 = new TH1F("pt_PA_0", "All Truth Polar Angle", 20, 0, 1.6); //All true polar angles
  TH1F *recoPhotonTruth_PA_0 = new TH1F("prt_PA_0", "Reco True Polar Angle", 20, 0, 1.6); //True polar angle of reconstructed photons
  TH1F *photonReco_PA_0 = new TH1F("pr_PA_0", "Reco Polar Angle", 20, 0, 1.6); //Reconstructed polar angle of linked photons

  TH1F *photonTruth_A_0 = new TH1F("pt_A_0", "All Truth Azimuth", 20, -1.6, 1.6); //All true azimuths
  TH1F *recoPhotonTruth_A_0 = new TH1F("prt_A_0", "Reco True Azimuth", 20, -1.6, 1.6); //True azimuth of reconstructed photons
  TH1F *photonReco_A_0 = new TH1F("pr_A_0", "Reco Azimuth", 20, -1.6, 1.6); //Reconstructed azimuth of linked photons

  TH1F *deltaR_0 = new TH1F("dR_0", "Delta R", 20, 0 , 0.01);
  TH1F *Eres_0 = new TH1F("Er_0", "Energy resolution", 40, -0.1, 0.1); //Resolution of reconstructed energy
  TH2F *Ephase_0 = new TH2F("Epha_0", "True vs Reco energy phase plot", 100, -500, 30000, 100, 2500, 30000); //Phase plot of true energy w.r.t reconstructed energy
  TH3F *K2Corr_0 = new TH3F("K2c_0", "Correction w.r.t polar angle and p_{T}", 20, 3000, 5000, 20, 0, 3.2, 20, 0.5, 1.5); //Correction term plot w.r.t. polar angle and transverse momentum
  TGraph2D *K2Corr2D_0 = new TGraph2D();

  

  // Root file reader stuff
  TTreeReader myReader("MyLCTuple", openFile);
  
  TTreeReaderArray<int> r2f_RA(myReader, "r2f"); //link reco particle index number
  TTreeReaderArray<int> r2t_RA(myReader, "r2t"); //link truth particle index number
  TTreeReaderArray<Float_t> r2w_RA(myReader, "r2w"); //link weight
  TTreeReaderArray<Float_t> mcmox_RA(myReader, "mcmox"); //truth x momentum
  TTreeReaderArray<Float_t> mcmoy_RA(myReader, "mcmoy"); //truth y momentum
  TTreeReaderArray<Float_t> mcmoz_RA(myReader, "mcmoz"); //truth z momentum
  TTreeReaderArray<Float_t> rcmox_RA(myReader, "rcmox"); //reco x momentum
  TTreeReaderArray<Float_t> rcmoy_RA(myReader, "rcmoy"); //reco y momentum
  TTreeReaderArray<Float_t> rcmoz_RA(myReader, "rcmoz"); //reco z momentum
  TTreeReaderArray<int> rctyp_RA(myReader, "rctyp"); //reconstructed type (11 = electron, 22 = photon, 2112 = neutron)
  TTreeReaderArray<int> mcpdg_RA(myReader, "mcpdg"); //truth type (see list above)
  TTreeReaderArray<int> mcgst_RA(myReader, "mcgst"); //is this particle a generating particle or not

  bool photonMatch = false;
  int r2fTemp, r2tTemp;
  TVector3 mcmoTemp, rcmoTemp;
  Float_t r2wTemp, dPATemp, dATemp;
  Float_t mcmoxTemp, mcmoyTemp, mcmozTemp, rcmoxTemp, rcmoyTemp, rcmozTemp, weightTemp;
  Float_t r2wMax = 0;
    //std::numeric_limits<Float_t>::infinity();
  int r2wMax_Index = -1;
  int arrCount = 0;

  // TArrayF K2Corr_X = new TArrayF(maxN);
  // TArrayF K2Corr_Y = new TArrayF(maxN);
  // TArrayF K2Corr_Z = new TArrayF(maxN);
  
  while(myReader.Next()) {

    mcmoxTemp = mcmox_RA.At(0);
    mcmoyTemp = mcmoy_RA.At(0);
    mcmozTemp = mcmoz_RA.At(0);
    
    mcmoTemp.SetXYZ(mcmoxTemp, mcmoyTemp, mcmozTemp);
    
    for(int i = 0; i < rcmox_RA.GetSize(); i++) {
      rcmoxTemp = rcmox_RA.At(i);
      rcmoyTemp = rcmoy_RA.At(i);
      rcmozTemp = rcmoz_RA.At(i);

      weightTemp = mcmoTemp.Mag() / sqrt( pow(mcmoxTemp - rcmoxTemp, 2) + pow(mcmoyTemp - rcmoyTemp, 2) +pow(mcmozTemp - rcmozTemp, 2) );

      if(rctyp_RA.At(i) == 22) {
	photonMatch = true;
      }

      if(weightTemp > r2wMax) {
        r2wMax = weightTemp;
        r2wMax_Index = i;
      }
    }

    if(r2wMax > 0) {
      rcmoTemp.SetXYZ(rcmox_RA.At(r2wMax_Index), rcmoy_RA.At(r2wMax_Index), rcmoz_RA.At(r2wMax_Index));
      
      photonReco_pt_0->Fill(rcmoTemp.Perp());
      photonReco_PA_0->Fill(rcmoTemp.Theta());
      photonReco_A_0->Fill(rcmoTemp.Phi());

      deltaR_0->Fill(sqrt( pow(rcmoTemp.Theta() - mcmoTemp.Theta(), 2) + pow(rcmoTemp.Phi() - mcmoTemp.Phi(), 2) ));
      Eres_0->Fill((rcmoTemp.Mag() - mcmoTemp.Mag()) / mcmoTemp.Mag());
      Ephase_0->Fill(rcmoTemp.Mag(), mcmoTemp.Mag());
      K2Corr_0->Fill(mcmoTemp.Perp(), mcmoTemp.Theta(), mcmoTemp.Mag() / rcmoTemp.Mag());
      K2Corr2D_0->SetPoint(arrCount, mcmoTemp.Perp(), mcmoTemp.Theta(), mcmoTemp.Mag() / rcmoTemp.Mag());

      arrCount++;

      if((r2wMax > (1 / maxpTRecoErr)) && photonMatch) {
	rcmoTemp.SetXYZ(rcmox_RA.At(r2fTemp), rcmoy_RA.At(r2fTemp), rcmoz_RA.At(r2fTemp));

	recoPhotonTruth_pt_0->Fill(mcmoTemp.Perp());
	recoPhotonTruth_PA_0->Fill(mcmoTemp.Theta());
	recoPhotonTruth_A_0->Fill(mcmoTemp.Phi());
      }
    } else {
      // Eres_0->Fill(-1);
    }

    photonTruth_pt_0->Fill(mcmoTemp.Perp());
    photonTruth_PA_0->Fill(mcmoTemp.Theta());
    photonTruth_A_0->Fill(mcmoTemp.Phi());
    
    r2wMax = 0;
    r2wMax_Index = -1;
    photonMatch = false;
  }

  gStyle->SetPalette(kRust);

  THStack *hs = new THStack("hs", "Transverse momentum of photons;p_{T} (GeV);Count");
  
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

  ///////////////////////
  // SIMPLE HISTOGRAMS //
  ///////////////////////

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

  // OTHERS
  
  deltaR_0->Draw();
  deltaR_0->SetTitle("Delta R");
  deltaR_0->GetXaxis()->SetTitle("Delta R (Rads)");
  deltaR_0->GetYaxis()->SetTitle("Count");

  c->SaveAs(saveDir + "/resolutionHists/deltaR.png");
  c->Close();
  c = new TCanvas();

  Eres_0->Draw();
  Eres_0->SetTitle("Energy Resolution");
  Eres_0->GetXaxis()->SetTitle("Energy Resolution");
  Eres_0->GetYaxis()->SetTitle("Count");

  Eres_0->Fit("gaus");

  Eres_0->FitPanel();

  c->SaveAs(saveDir + "/resolutionHists/Eres.png");
  c->Close();
  c = new TCanvas();

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
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");
  K2Corr_0->GetZaxis()->SetTitle("Correction term");

  K2Corr_0->SetMarkerSize(10);

  // K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr.png");
  c->Close();
  c = new TCanvas();

  // K2Corr axes
  K2Corr_0->ProjectionX()->Draw();
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - X hist");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_aX.png");
  c->Close();
  c = new TCanvas();

  K2Corr_0->ProjectionY()->Draw();
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - Y hist");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_aY.png");
  c->Close();
  c = new TCanvas();

  K2Corr_0->ProjectionZ()->Draw();
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - Z hist");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_aZ.png");
  c->Close();
  c = new TCanvas();

  // K2Corr profiles

  K2Corr_0->Project3D("xy")->Draw("colz");
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - xy profile");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_Pxy.png");
  c->Close();
  c = new TCanvas();

  K2Corr_0->Project3D("zx")->Draw("colz");
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - zx profile");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_Pzx.png");
  c->Close();
  c = new TCanvas();

  K2Corr_0->Project3D("zy")->Draw("colz");
  K2Corr_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot - zy profile");
  K2Corr_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr_0->GetYaxis()->SetTitle("Polar Angle (Rads)");

  K2Corr_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr_Pzy.png");
  c->Close();
  c = new TCanvas();

  // TGraph2D :

  K2Corr2D_0->SetMinimum(0.9);
  K2Corr2D_0->SetMaximum(1.5);
  
  K2Corr2D_0->Draw("surf1");
  K2Corr2D_0->SetTitle("p_{T}, Polar Angle, and Energy Correction Phase Plot");
  K2Corr2D_0->GetXaxis()->SetTitle("p_{T} (GeV)");
  K2Corr2D_0->GetYaxis()->SetTitle("Polar Angle (Rads)");
  
  //K2Corr2D_0->SetStats(0);

  c->SaveAs(saveDir + "/phaseHists/K2Corr2D.png");
  c->Close();
  c = new TCanvas();
  

  ///////////////////////////
  // EFFICIENCY HISTOGRAMS //
  ///////////////////////////
  
  TEfficiency* pt_Eff = new TEfficiency("pt_Eff", "Efficiency of photon p_{T} reconstruction;p_{T} (GeV);Efficiency", 20, 100, 1000);
  TFile* eff_File = new TFile("effFile.root", "recreate");

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

  pt_Eff->SetLineColor(4);
  PA_Eff->SetLineColor(4);

  //draws an efficiency graph for the transverse momentum
  TMultiGraph *mg = new TMultiGraph();
  c = new TCanvas();

  mg->Add(pt_Eff->CreateGraph());

  mg->SetTitle("Transverse momentum reconstruction efficiency");
  mg->GetXaxis()->SetTitle("p_{T} (GeV)");
  mg->GetYaxis()->SetTitle("Efficiency");

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

  mg->Draw("aZ");

  c->SaveAs(saveDir + "/efficiency/Eff_PA.png");
  c->Close();
}
