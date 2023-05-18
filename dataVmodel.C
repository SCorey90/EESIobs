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

    dataHists->ls();
    auto * mDataMass = (TH1F*)dataHists->Get("mMass");
    auto * mDataPairPT = (TH1F*)dataHists->Get("mPperp");

    auto * mDataPairPhi = (TH1F*)dataHists->Get("mPairPhi");

    auto * mDataPhivsPT = (TH2F*)dataHists->Get("mPhivsPT");
    auto * mDataCos2phivsPT = (TH2F*)dataHists->Get("mCos2phivsPT");
    auto * mDataCos4phivsPT = (TH2F*)dataHists->Get("mCos4phivsPT");
    
    modelHists->ls();
    auto * mToyMass = (TH1F*)modelHists->Get("mReconstructedM");
    auto * mToyPairPT = (TH1F*)modelHists->Get("mPairPT");

    auto * mToyPairPhi = (TH1F*)modelHists->Get("mPairPhi");
    auto * mToyPairPhiwMu = (TH1F*)modelHists->Get("mPairPhiwMu");
    auto * mToyNoisyPairPhi = (TH1F*)modelHists->Get("mNoisyPairPhi");
    auto * mToyCutPairPhi = (TH1F*)modelHists->Get("mCutPairPhi");

    auto * mToyPhivsPT = (TH2F*)modelHists->Get("mPhivsPT");
    auto * mToyCos2phivsPT = (TH2F*)modelHists->Get("mCos2phivsPT");
    auto * mToyCos4phivsPT = (TH2F*)modelHists->Get("mCos4phivsPT");

    //new TFile
    TFile * fo = new TFile( "dataVmodelplots.root", "RECREATE" );

    //analysis
    mToyMass->Scale(1/mToyMass->Integral());
    mDataMass->Scale(1/mDataMass->Integral());
    auto * mDataMinusToyMass = (TH1F*)mDataMass->Clone("mDataMinusToyMass");
    mDataMinusToyMass->Add(mToyMass, -1);

    mToyPairPT->Scale(1/mToyPairPT->Integral());
    mDataPairPT->Scale(1/mDataPairPT->Integral());
    auto * mDataMinusToyPairPT = (TH1F*)mDataPairPT->Clone("mDataMinusToyPairPT");
    mDataMinusToyPairPT->Add(mToyPairPT, -1);

    mToyPairPhi->Scale(2*3.1415/mToyPairPhi->GetEntries());
    mDataPairPhi->Scale(2*3.1415/mDataPairPhi->GetEntries());
    auto * mDataMinusToyPhi = (TH1F*)mDataPairPhi->Clone("mDataMinusToyPhi");
    mDataMinusToyPhi->Add(mToyPairPhi, -1);
    
//make plots
fo->cd();

makeCanvas();
mToyMass->SetLineColor(kBlack);
mToyMass->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyMass->Draw();
gPad->Print( "plot_mNormToyMass.pdf" );
gPad->Print( "plot_mNormToyMass.png" );

makeCanvas();
mDataMass->SetLineColor(kBlack);
mDataMass->SetTitle("Data mass distribution; Mass (GeV); counts");
mDataMass->Draw();
gPad->Print( "plot_mNormDataMass.pdf" );
gPad->Print( "plot_mNormDataMass.png" );

makeCanvas();
mDataMinusToyMass->SetLineColor(kBlack);
mDataMinusToyMass->SetTitle("Data-toy model mass distribution; Mass (GeV); counts");
mDataMinusToyMass->Draw();
gPad->Print( "plot_mDataMinusToyMass.pdf" );
gPad->Print( "plot_mDataminusToyMass.png" );

makeCanvas();
mToyPairPT->SetLineColor(kBlack);
mToyPairPT->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyPairPT->Draw();
gPad->Print( "plot_mNormToyPairPT.pdf" );
gPad->Print( "plot_mNormToyPairPT.png" );

makeCanvas();
mDataPairPT->SetLineColor(kBlack);
mDataPairPT->SetTitle("Data pair P_{T} distribution; P_{T} (GeV); counts");
mDataPairPT->Draw();
gPad->Print( "plot_mNormDataPairPT.pdf" );
gPad->Print( "plot_mNormDataPairPT.png" );

makeCanvas();
mDataMinusToyPairPT->SetLineColor(kBlack);
mDataMinusToyPairPT->SetTitle("Data-toy model pair P_{T} distribution; P_{T} (GeV); counts");
mDataMinusToyPairPT->Draw();
gPad->Print( "plot_mDataMinusToyPairPT.pdf" );
gPad->Print( "plot_mDataminusToyPairPT.png" );

makeCanvas();
mToyPairPhi->SetLineColor(kBlack);
mToyPairPhi->SetTitle("Toy model #phi distribution; #phi (rad); counts");
mToyPairPhi->Draw();
gPad->Print( "plot_mNormToyPhi.pdf" );
gPad->Print( "plot_mNormToyPhi.png" );

makeCanvas();
mDataPairPhi->SetLineColor(kBlack);
mDataPairPhi->SetTitle("Data #phi distribution; #phi (rad); counts");
mDataPairPhi->Draw();
gPad->Print( "plot_mNormDataPhi.pdf" );
gPad->Print( "plot_mNormDataPhi.png" );

makeCanvas();
mDataMinusToyPhi->SetLineColor(kBlack);
mDataMinusToyPhi->SetTitle("Data-toy model #phi distribution; #phi (rad); counts");
mDataMinusToyPhi->Draw();
gPad->Print( "plot_mDataMinusToyPhi.pdf" );
gPad->Print( "plot_mDataminusToyPhi.png" );

}

