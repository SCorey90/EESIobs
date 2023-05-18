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
    auto * mDataPairPhi = (TH1F*)dataHists->Get("mPairPhi");

    auto * mDataPhivsPT = (TH2F*)dataHists->Get("mPhivsPT");
    auto * mDataCos2phivsPT = (TH2F*)dataHists->Get("mCos2phivsPT");
    auto * mDataCos4phivsPT = (TH2F*)dataHists->Get("mCos4phivsPT");
    
    modelHists->ls();
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
    mToyPairPhi->Scale(1/mToyPairPhi->GetMean());
    mDataPairPhi->Scale(1/mDataPairPhi->GetMean());
    
    auto * mDataMinusToyPhi = (TH1F*)mDataPairPhi->Clone("mDataMinusToyPhi");
    mDataMinusToyPhi->Add(mToyPairPhi, -1);
    
//make plots
fo->cd();

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

