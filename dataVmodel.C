#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    can->SetRightMargin(0.35);
}

void dataVmodel() {
    //import histograms
    TFile * dataHists = TFile::Open("EESIplots.root");
    TFile * modelHists = TFile::Open("toymodelplots.root");
    TFile * mixedeventHists = TFile::Open("MixedEventplots.root");
    TFile * goldHists = TFile::Open("goldEESIplots.root");


    dataHists->ls();
    auto * mDataPiPT = (TH1F*)dataHists->Get("mPiPT");
    auto * mDataPiEta = (TH1F*)dataHists->Get("mPiEta");
    auto * mDataPiAngle = (TH1F*)dataHists->Get("mPiAngle");

    auto * mDataMass = (TH1F*)dataHists->Get("mMass");
    auto * mDataPairPT = (TH1F*)dataHists->Get("mPperp");

    auto * mDataPairPhi = (TH1F*)dataHists->Get("mPairPhi");

    auto * mDataPhivsPT = (TH2F*)dataHists->Get("mPhivsPT");
    auto * mDataCos2PhivsPT = (TH2F*)dataHists->Get("mCos2phivsPT");
    auto * mDataCos4PhivsPT = (TH2F*)dataHists->Get("mCos4phivsPT");
    
    auto * mDataPhivsMass = (TH2F*)dataHists->Get("mPhivsMass");
    auto * mDataCos2PhivsMass = (TH2F*)dataHists->Get("mCos2phivsMass");
    auto * mDataCos4PhivsMass = (TH2F*)dataHists->Get("mCos4phivsMass");

    auto * mUDataCos2phivsPTmoments = (TProfile*)dataHists->Get("mPTcos2phimoments");

    goldHists->ls();
    auto * mAuDataCos2phivsPTmoments = (TProfile*)goldHists->Get("mPTcos2phimoments");

    modelHists->ls();
    auto * mToyMass = (TH1F*)modelHists->Get("mReconstructedM");
    auto * mToyPairPT = (TH1F*)modelHists->Get("mPairPT");

    auto * mToyPairPhi = (TH1F*)modelHists->Get("mPairPhi");
    auto * mToyPairPhiwMu = (TH1F*)modelHists->Get("mPairPhiwMu");
    auto * mToyNoisyPairPhi = (TH1F*)modelHists->Get("mNoisyPairPhi");
    auto * mToyCutPairPhi = (TH1F*)modelHists->Get("mCutPairPhi");

    auto * mToyPhivsPT = (TH2F*)modelHists->Get("mPhivsPT");
    auto * mToyCos2PhivsPT = (TH2F*)modelHists->Get("mCos2cutflatvsPT");
    //auto * mToyCos2PhivsPT = (TH2F*)modelHists->Get("mCos2decayphivsPT");
    auto * mToyCos4PhivsPT = (TH2F*)modelHists->Get("mCos4phivsPT");

    auto * mToyPhivsMass = (TH2F*)modelHists->Get("mPhivsMass");
    auto * mToyCos2PhivsMass = (TH2F*)modelHists->Get("mCos2phivsMass");
    auto * mToyCos4PhivsMass = (TH2F*)modelHists->Get("mCos4phivsMass");

    mixedeventHists->ls();
    auto * mDataLikePiPT = (TH1F*)mixedeventHists->Get("mPiDataLikePT");
    auto * mDataLikePiEta = (TH1F*)mixedeventHists->Get("mPiDataLikeEta");
    auto * mDataLikePiAngle = (TH1F*)mixedeventHists->Get("mPiDataLikeAngle");

    auto * mDataLikePT = (TH1F*)mixedeventHists->Get("mDataLikePT");

    auto * mMixedMass = (TH1F*)mixedeventHists->Get("mMass");
    auto * mMixedPairPT = (TH1F*)mixedeventHists->Get("mPperp");

    auto * mMixedPairPhi = (TH1F*)mixedeventHists->Get("mPairPhi");
    auto * mLikePairPhi = (TH1F*)mixedeventHists->Get("mLikePhi");

    auto * mMixedCos2PhivsPT = (TH2F*)mixedeventHists->Get("mCos2PhivsPT");
    auto * mLikeCos2PhivsPT = (TH2F*)mixedeventHists->Get("mLikeCos2PhivsPT");
    auto * mDataLikeCos2PhivsPT = (TH2F*)mixedeventHists->Get("mDataLikeCos2PhivsPT");

    auto * mMixedCos2PhivsMass = (TH2F*)mixedeventHists->Get("mCos2PhivsMass");
    
    //new TFile
    TFile * fo = new TFile( "dataVmodelplots.root", "RECREATE" );

    //analysis
/*    mDataPiPT->Scale(1.5/mDataPiPT->Integral("width"));
    mDataPiEta->Scale(5/mDataPiEta->Integral("width"));
    mDataPiAngle->Scale(6.18/mDataPiAngle->Integral("width"));

    mDataLikePiPT->Scale(1.5/mDataLikePiPT->Integral("width"));
    mDataLikePiEta->Scale(5/mDataLikePiEta->Integral("width"));
    mDataLikePiAngle->Scale(6.18/mDataLikePiAngle->Integral("width"));

    mToyMass->Scale(1/mToyMass->Integral());
    mDataMass->Scale(1/mDataMass->Integral());
    auto * mDataMinusToyMass = (TH1F*)mDataMass->Clone("mDataMinusToyMass");
    mDataMinusToyMass->Add(mToyMass, -1);

    mToyPairPT->Scale(1/mToyPairPT->Integral("width"));
    mDataPairPT->Scale(1/mDataPairPT->Integral("width"));
    mDataLikePT->Scale(1/mDataLikePT->Integral("width"));
    auto * mDataMinusToyPairPT = (TH1F*)mDataPairPT->Clone("mDataMinusToyPairPT");
    mDataMinusToyPairPT->Add(mToyPairPT, -1);

    mToyPairPhi->Scale(2*3.1415/(mToyPairPhi->Integral("width")));
    mDataPairPhi->Scale(2*3.1415/(mDataPairPhi->GetBinWidth(1) * mDataPairPhi->Integral()));
    auto * mDataMinusToyPhi = (TH1F*)mDataPairPhi->Clone("mDataMinusToyPhi");
    mDataMinusToyPhi->Add(mToyPairPhi, -1);

    mToyPairPhiwMu->Scale(2*3.1415/(mToyPairPhiwMu->GetBinWidth(1) * mToyPairPhiwMu->Integral()));
    auto * mDataMinusToyPhiwMu = (TH1F*)mDataPairPhi->Clone("mDataMinusToyPhiwMu");
    mDataMinusToyPhiwMu->Add(mToyPairPhiwMu, -1);

    mToyNoisyPairPhi->Scale(2*3.1415/(mToyNoisyPairPhi->GetBinWidth(1) * mToyNoisyPairPhi->Integral()));
    auto * mDataMinusToyNoisyPhi = (TH1F*)mDataPairPhi->Clone("mDataMinusToyNoisyPhi");
    mDataMinusToyNoisyPhi->Add(mToyNoisyPairPhi, -1);

    mToyCutPairPhi->Scale(2*3.1415/(mToyCutPairPhi->GetBinWidth(1) * mToyCutPairPhi->Integral()));
    auto * mDataMinusToyCutPhi = (TH1F*)mDataPairPhi->Clone("mDataMinusToyCutPhi");
    mDataMinusToyCutPhi->Add(mToyCutPairPhi, -1);

    mToyPhivsPT->Scale((2*3.1415*0.25)/mToyPhivsPT->Integral(1,-1,1,-1, "width"));
    mDataPhivsPT->Scale((2*3.1415*0.25)/mDataPhivsPT->Integral(1,-1,1,-1, "width"));
    auto * mDataMinusToyPhivsPT = (TH2F*)mDataPhivsPT->Clone("mDataMinusToyPhivsPT");
    mDataMinusToyPhivsPT->Add(mToyPhivsPT, -1);

    mToyCos2PhivsPT->Scale((2*3.1415*0.25)/mToyCos2PhivsPT->Integral(1,-1,1,-1, "width"));
    mDataCos2PhivsPT->Scale((2*3.1415*0.25)/mDataCos2PhivsPT->Integral(1,-1,1,-1, "width"));
    auto * mDataMinusToyCos2PhivsPT = (TH2F*)mDataCos2PhivsPT->Clone("mDataMinusToyCos2PhivsPT");
    mDataMinusToyCos2PhivsPT->Add(mToyCos2PhivsPT, -1);

    mToyPhivsMass->Scale((2*3.1415*0.25)/mToyPhivsMass->Integral(1,-1,1,-1, "width"));
    mDataPhivsMass->Scale((2*3.1415*0.25)/mDataPhivsMass->Integral(1,-1,1,-1, "width"));
    auto * mDataMinusToyPhivsMass = (TH2F*)mDataPhivsMass->Clone("mDataMinusToyPhivsMass");
    mDataMinusToyPhivsMass->Add(mToyPhivsMass, -1);

    mToyCos2PhivsMass->Scale((2*3.1415*0.25)/mToyCos2PhivsMass->Integral(1,-1,1,-1, "width"));
    mDataCos2PhivsMass->Scale((2*3.1415*0.25)/mDataCos2PhivsMass->Integral(1,-1,1,-1, "width"));
    auto * mDataMinusToyCos2PhivsMass = (TH2F*)mDataCos2PhivsMass->Clone("mDataMinusToyCos2PhivsMass");
    mDataMinusToyCos2PhivsMass->Add(mToyCos2PhivsMass, -1);

auto * mToyPTCos2PhiMoments = mToyCos2PhivsPT->ProfileY("mToyPTCos2PhiMoments", 1, -1);
auto * mDataPTCos2PhiMoments = mDataCos2PhivsPT->ProfileY("mDataPTCos2PhiMoments", 1, -1);
auto * mMixedPTCos2PhiMoments = mMixedCos2PhivsPT->ProfileY("mMixedPTCos2PhiMoments", 1, -1);
auto * mLikePTCos2PhiMoments = mLikeCos2PhivsPT->ProfileY("mLikePTCos2PhiMoments", 1, -1);
auto * mDataLikePTCos2PhiMoments = mDataLikeCos2PhivsPT->ProfileY("mDataLikePTCos2PhiMoments", 1, -1);
auto * mDatabyToyPTn2Moments = mDataCos2PhivsPT->ProfileY("mDatabyToyn2Moments", 1, -1);
mDatabyToyPTn2Moments->Divide(mToyPTCos2PhiMoments); 

auto * mToyMassCos2PhiMoments = mToyCos2PhivsMass->ProfileY("mToyMassCos2PhiMoments", 1, -1);
auto * mDataMassCos2PhiMoments = mDataCos2PhivsMass->ProfileY("mDataMassCos2PhiMoments", 1, -1);
auto * mMixedMassCos2PhiMoments = mMixedCos2PhivsMass->ProfileY("mMixedMassCos2PhiMoments", 1, -1);
*/

//make plots
fo->cd();
/*
makeCanvas();
mToyMass->SetLineColor(kBlack);
mToyMass->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyMass->Draw();
gPad->Print( "plots/dataVmodel/Mass/plot_mNormToyMass.pdf" );
gPad->Print( "plots/dataVmodel/Mass/plot_mNormToyMass.png" );

makeCanvas();
mDataMass->SetLineColor(kBlack);
mDataMass->SetTitle("Data mass distribution; Mass (GeV); counts");
mDataMass->Draw();
gPad->Print( "plots/dataVmodel/Mass/plot_mNormDataMass.pdf" );
gPad->Print( "plots/dataVmodel/Mass/plot_mNormDataMass.png" );

makeCanvas();
mDataMinusToyMass->SetLineColor(kBlack);
mDataMinusToyMass->SetTitle("Data-toy model mass distribution; Mass (GeV); counts");
mDataMinusToyMass->Draw();
gPad->Print( "plots/dataVmodel/Mass/plot_mDataMinusToyMass.pdf" );
gPad->Print( "plots/dataVmodel/Mass/plot_mDataminusToyMass.png" );

makeCanvas();
mToyPairPT->SetLineColor(kBlack);
mToyPairPT->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyPairPT->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyPairPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyPairPT.png" );

makeCanvas();
mDataPairPT->SetLineColor(kBlack);
mDataPairPT->SetTitle("Data pair P_{T} distribution; P_{T} (GeV/c); counts");
mDataPairPT->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataPairPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataPairPT.png" );

makeCanvas();
mDataMinusToyPairPT->SetLineColor(kBlack);
mDataMinusToyPairPT->SetTitle("Data-toy model pair P_{T} distribution; P_{T} (GeV/c); counts");
mDataMinusToyPairPT->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mDataMinusToyPairPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDataminusToyPairPT.png" );

makeCanvas();
mToyPairPhi->SetLineColor(kBlack);
mToyPairPhi->SetTitle("Toy model #phi distribution; #phi (rad); counts");
mToyPairPhiwMu->SetLineColor(kRed);
mToyNoisyPairPhi->SetLineColor(kBlue);
mToyCutPairPhi->SetLineColor(kOrange);
mToyPairPhi->Draw();
mToyPairPhiwMu->Draw("SAME");
mToyNoisyPairPhi->Draw("SAME");
mToyCutPairPhi->Draw("SAME");
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->AddEntry(mToyPairPhi,"No #pi->#mu decay, no noise, no cut");
legend->AddEntry(mToyPairPhiwMu,"With #pi->#mu decay, no noise, no cut");
legend->AddEntry(mToyNoisyPairPhi,"With #pi->#mu decay, with noise, no cut");
legend->AddEntry(mToyCutPairPhi,"With #pi->#mu decay, with noise, with cut");
legend->Draw();
gPad->Print( "plots/dataVmodel/plot_mNormToyPhi.pdf" );
gPad->Print( "plots/dataVmodel/plot_mNormToyPhi.png" );

makeCanvas();
mDataPairPhi->SetLineColor(kBlack);
mDataPairPhi->SetTitle("Data #phi distribution; #phi (rad); counts");
mDataPairPhi->Draw();
gPad->Print( "plots/dataVmodel/plot_mNormDataPhi.pdf" );
gPad->Print( "plots/dataVmodel/plot_mNormDataPhi.png" );

makeCanvas();
mDataMinusToyPhi->SetLineColor(kBlack);
mDataMinusToyPhi->SetTitle("Data-toy model #phi distribution; #phi (rad); counts");
mDataMinusToyPhiwMu->SetLineColor(kRed);
mDataMinusToyNoisyPhi->SetLineColor(kBlue);
mDataMinusToyCutPhi->SetLineColor(kOrange);
mDataMinusToyCutPhi->Draw();
mDataMinusToyPhiwMu->Draw("SAME");
mDataMinusToyNoisyPhi->Draw("SAME");
mDataMinusToyPhi->Draw("SAME");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->AddEntry(mDataMinusToyPhi,"No #pi->#mu decay, no noise, no cut");
legend2->AddEntry(mDataMinusToyPhiwMu,"With #pi->#mu decay, no noise, no cut");
legend2->AddEntry(mDataMinusToyNoisyPhi,"With #pi->#mu decay, with noise, no cut");
legend2->AddEntry(mDataMinusToyCutPhi,"With #pi->#mu decay, with noise, with cut");
legend2->Draw();
gPad->Print( "plots/dataVmodel/plot_mDataMinusToyPhi.pdf" );
gPad->Print( "plots/dataVmodel/plot_mDataminusToyPhi.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataPhivsPT->SetTitle("Data #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV/c); counts");
mDataPhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataPhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mToyPhivsPT->SetTitle("Toy model #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV/c); counts");
mToyPhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyPhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataMinusToyPhivsPT->SetTitle("Data-toy model #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV/c); counts");
mDataMinusToyPhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mDataMinusToyPhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDataMinusToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataCos2PhivsPT->SetTitle("Data cos(2#phi) vs. P_{T} (scaled to avg=1); cos(2#phi); P_{T} (GeV/c); counts");
mDataCos2PhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataCos2PhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormDataCos2PhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mToyCos2PhivsPT->SetTitle("Toy model cos(2#phi) vs. P_{T} (scaled to avg=1); cos(2#phi); P_{T} (GeV/c); counts");
mToyCos2PhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyCos2PhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mNormToyCos2PhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataMinusToyCos2PhivsPT->SetTitle("Data-toy model cos(2#phi) vs. P_{T}; cos(2#phi); P_{T} (GeV/c); counts");
mDataMinusToyCos2PhivsPT->Draw("colz");
gPad->Print( "plots/dataVmodel/PT/plot_mDataMinusToyCos2PhivsPT.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDataMinusToyCos2PhivsPT.png" );

makeCanvas();
mToyPTCos2PhiMoments->SetLineColor(kGreen);
mDataPTCos2PhiMoments->SetLineColor(kBlack);
mToyPTCos2PhiMoments->SetTitle("Strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV/c); 2<cos(2#phi)>");
mToyPTCos2PhiMoments->SetMaximum(0.5);
mToyPTCos2PhiMoments->SetMinimum(-0.5);
mToyPTCos2PhiMoments->Draw();
mDataPTCos2PhiMoments->Draw("SAME");
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->AddEntry(mDataPTCos2PhiMoments,"Run 12 U+U data");
legend3->AddEntry(mToyPTCos2PhiMoments,"Toy model");
legend3->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mDataVToyPTCos2PhiMoments.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDataVToyPTCos2PhiMoments.png" );

makeCanvas();
mDatabyToyPTn2Moments->SetLineColor(kBlack);
mDatabyToyPTn2Moments->SetTitle("Data/model cos2#phi signal strength; P_{T} (GeV/c); data 2<cos2#phi>/toy model 2<cos2#phi>");
mDatabyToyPTn2Moments->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mDatabyToyPTn2Moments.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDatabyToyPTn2Moments.png" );

makeCanvas();
mToyMassCos2PhiMoments->SetLineColor(kGreen);
mDataMassCos2PhiMoments->SetLineColor(kBlack);
mToyMassCos2PhiMoments->SetTitle("Strength of cos(2#phi) signal vs. pair combined mass; mass (GeV/c); 2<cos(2#phi)>");
mToyMassCos2PhiMoments->SetMaximum(0.7);
mToyMassCos2PhiMoments->SetMinimum(-0.5);
mToyMassCos2PhiMoments->Draw();
mDataMassCos2PhiMoments->Draw("SAME");
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->AddEntry(mDataMassCos2PhiMoments,"Run 12 U+U data");
legend4->AddEntry(mToyMassCos2PhiMoments,"Toy model");
legend4->Draw();
gPad->Print( "plots/dataVmodel/Massplot_mDataVToyMassCos2PhiMoments.pdf" );
gPad->Print( "plots/dataVmodel/Mass/plot_mDataVToyMassCos2PhiMoments.png" );

makeCanvas();
mMixedPTCos2PhiMoments->SetLineColor(kGreen);
mDataPTCos2PhiMoments->SetLineColor(kBlack);
//mLikePTCos2PhiMoments->SetLineColor(kMagenta);
mMixedPTCos2PhiMoments->SetTitle("Mixed parent P_{T} < original parent P_{T} (0.65 < M_{#rho^{0}} < 0.75 GeV); P_{T} (GeV/c); 2<cos(2#phi)>");
mMixedPTCos2PhiMoments->SetMaximum(0.5);
mMixedPTCos2PhiMoments->Draw();
mDataPTCos2PhiMoments->Draw("SAME");
//mLikePTCos2PhiMoments->Draw("SAME");
auto legend5 = new TLegend(0.65,0.1,0.95,0.4);
legend5->AddEntry(mDataPTCos2PhiMoments,"Run 12 U+U data");
legend5->AddEntry(mMixedPTCos2PhiMoments,"Mixed event like sign source");
//legend5->AddEntry(mLikePTCos2PhiMoments,"Mixed event like sign pairs");
legend5->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mDataVMixedPTCos2PhiMoments.pdf" );
gPad->Print( "plots/dataVmodel/PT/plot_mDataVMixedPTCos2PhiMoments.png" );

makeCanvas();
mDataPiPT->SetLineColor(kBlack);
mDataLikePiPT->SetLineColor(kMagenta);
mDataPiPT->Draw();
mDataLikePiPT->Draw("SAME");
auto legend6 = new TLegend(0.65, 0.1, 0.95, 0.4);
legend6->AddEntry(mDataPiPT, "#pi from unlike sign pair");
legend6->AddEntry(mDataLikePiPT, "#pi from like sign pair");
legend6->Draw();
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePiPT.png" );
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePiPT.pdf" );

makeCanvas();
mDataPiEta->SetLineColor(kBlack);
mDataLikePiEta->SetLineColor(kMagenta);
mDataLikePiEta->Draw();
mDataPiEta->Draw("SAME");
auto legend7 = new TLegend(0.65, 0.1, 0.95, 0.4);
legend7->AddEntry(mDataPiEta, "#pi from unlike sign pair");
legend7->AddEntry(mDataLikePiEta, "#pi from like sign pair");
legend7->Draw();
gPad->Print( "plots/data/Rapidity/plot_mLikevsUnlikePiEta.png" );
gPad->Print( "plots/data/Rapidity/plot_mLikevsUnlikePiEta.pdf" );

makeCanvas();
mDataPiAngle->SetLineColor(kBlack);
mDataLikePiAngle->SetLineColor(kMagenta);
mDataLikePiAngle->Draw();
mDataPiAngle->Draw("SAME");
auto legend8 = new TLegend(0.65, 0.1, 0.95, 0.4);
legend8->AddEntry(mDataPiAngle, "#pi from unlike sign pair");
legend8->AddEntry(mDataLikePiAngle, "#pi from like sign pair");
legend8->Draw();
gPad->Print( "plots/data/Angle/plot_mLikevsUnlikePiAngle.png" );
gPad->Print( "plots/data/Angle/plot_mLikevsUnlikePiAngle.pdf" );

makeCanvas();
mDataLikePTCos2PhiMoments->SetLineColor(kGreen);
mDataPTCos2PhiMoments->SetLineColor(kBlack);
mDataPTCos2PhiMoments->SetTitle("Strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV/c); 2<cos(2#phi)>");
mDataPTCos2PhiMoments->Draw();
mDataLikePTCos2PhiMoments->Draw("SAME");
auto legend9 = new TLegend(0.65,0.1,0.95,0.4);
legend9->AddEntry(mDataPTCos2PhiMoments,"Data unlike sign pairs");
legend9->AddEntry(mDataLikePTCos2PhiMoments,"Data like sign pairs");
legend9->Draw();
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePTCos2PhiMoments.pdf" );
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePTCos2PhiMoments.png" );

makeCanvas();
mDataLikePT->SetLineColor(kGreen);
mDataPairPT->SetLineColor(kBlack);
mDataPairPT->SetTitle("Pair P_{T}; P_{T} (GeV/c); counts");
mDataPairPT->Draw("HIST C");
mDataLikePT->Draw("HIST SAME C");
auto legend10 = new TLegend(0.65,0.1,0.95,0.4);
legend10->AddEntry(mDataPairPT,"Data unlike sign pairs");
legend10->AddEntry(mDataLikePT,"Data like sign pairs");
legend10->Draw();
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePT.pdf" );
gPad->Print( "plots/data/PT/plot_mLikevsUnlikePT.png" );

makeCanvas();
mMixedPTCos2PhiMoments->SetLineColor(kGreen);
mDataPTCos2PhiMoments->SetLineColor(kBlack);
//mLikePTCos2PhiMoments->SetLineColor(kMagenta);
mMixedMassCos2PhiMoments->SetTitle("PM_{T}(#rho^{0}) < 0.6 GeV/c); P_{T} (GeV/c); 2<cos(2#phi)>");
mMixedMassCos2PhiMoments->SetMaximum(0.5);
mMixedMassCos2PhiMoments->Draw();
mDataMassCos2PhiMoments->Draw("SAME");
auto legend11 = new TLegend(0.65,0.1,0.95,0.4);
legend11->AddEntry(mDataMassCos2PhiMoments,"Run 12 U+U data");
legend11->AddEntry(mMixedMassCos2PhiMoments,"Mixed event like sign source");
legend11->Draw();
gPad->Print( "plots/dataVmodel/Mass/plot_mDataVMixedMassCos2PhiMoments.pdf" );
gPad->Print( "plots/dataVmodel/Mass/plot_mDataVMixedMassCos2PhiMoments.png" );
*/

makeCanvas();
mUDataCos2phivsPTmoments->SetLineColor(kGreen-2);
mAuDataCos2phivsPTmoments->SetLineColor(kBlack);
mAuDataCos2phivsPTmoments->SetTitle("AuAu vs. UU 2<cos2#phi> vs. Pair P_{T} (0.65 < M_{#pi#pi} < 0.9 GeV); Pair P_{T}; 2<cos2#phi>");
mUDataCos2phivsPTmoments->Draw();
mAuDataCos2phivsPTmoments->Draw();
auto legend12 = new TLegend(0.65,0.1,0.95,0.4);
legend12->AddEntry(mUDataCos2phivsPTmoments,"U+U");
legend12->AddEntry(mAuDataCos2phivsPTmoments,"Au+Au");
legend12->Draw();
gPad->Print( "plots/dataVmodel/PT/plot_mDataAuvsU.pdf");

}
