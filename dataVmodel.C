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
    auto * mDataCos2PhivsPT = (TH2F*)dataHists->Get("mCos2phivsPT");
    auto * mDataCos4PhivsPT = (TH2F*)dataHists->Get("mCos4phivsPT");
    
    auto * mDataPhivsMass = (TH2F*)dataHists->Get("mPhivsMass");
    auto * mDataCos2PhivsMass = (TH2F*)dataHists->Get("mCos2phivsMass");
    auto * mDataCos4PhivsMass = (TH2F*)dataHists->Get("mCos4phivsMass");

    modelHists->ls();
    auto * mToyMass = (TH1F*)modelHists->Get("mReconstructedM");
    auto * mToyPairPT = (TH1F*)modelHists->Get("mPairPT");

    auto * mToyPairPhi = (TH1F*)modelHists->Get("mPairPhi");
    auto * mToyPairPhiwMu = (TH1F*)modelHists->Get("mPairPhiwMu");
    auto * mToyNoisyPairPhi = (TH1F*)modelHists->Get("mNoisyPairPhi");
    auto * mToyCutPairPhi = (TH1F*)modelHists->Get("mCutPairPhi");

    auto * mToyPhivsPT = (TH2F*)modelHists->Get("mPhivsPT");
    auto * mToyCos2PhivsPT = (TH2F*)modelHists->Get("mCos2phivsPT");
    auto * mToyCos4PhivsPT = (TH2F*)modelHists->Get("mCos4phivsPT");

    auto * mToyPhivsMass = (TH2F*)modelHists->Get("mPhivsMass");
    auto * mToyCos2PhivsMass = (TH2F*)modelHists->Get("mCos2phivsMass");
    auto * mToyCos4PhivsMass = (TH2F*)modelHists->Get("mCos4phivsMass");

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
auto * mDataMinusToyPTCos2PhiMoments = mDataMinusToyCos2PhivsPT->ProfileY("mDataMinusToyPTCos2PhiMoments", 1, -1);

auto * mToyMassCos2PhiMoments = mToyCos2PhivsMass->ProfileY("mToyMassCos2PhiMoments", 1, -1);
auto * mDataMassCos2PhiMoments = mDataCos2PhivsMass->ProfileY("mDataMassCos2PhiMoments", 1, -1);
auto * mDataMinusToyMassCos2PhiMoments = mDataMinusToyCos2PhivsMass->ProfileY("mDataMinusToyMassCos2PhiMoments", 1, -1);

//make plots
fo->cd();

makeCanvas();
mToyMass->SetLineColor(kBlack);
mToyMass->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyMass->Draw();
gPad->Print( "plots/plot_mNormToyMass.pdf" );
gPad->Print( "plots/plot_mNormToyMass.png" );

makeCanvas();
mDataMass->SetLineColor(kBlack);
mDataMass->SetTitle("Data mass distribution; Mass (GeV); counts");
mDataMass->Draw();
gPad->Print( "plots/plot_mNormDataMass.pdf" );
gPad->Print( "plots/plot_mNormDataMass.png" );

makeCanvas();
mDataMinusToyMass->SetLineColor(kBlack);
mDataMinusToyMass->SetTitle("Data-toy model mass distribution; Mass (GeV); counts");
mDataMinusToyMass->Draw();
gPad->Print( "plots/plot_mDataMinusToyMass.pdf" );
gPad->Print( "plots/plot_mDataminusToyMass.png" );

makeCanvas();
mToyPairPT->SetLineColor(kBlack);
mToyPairPT->SetTitle("Toy model mass distribution; Mass (GeV); counts");
mToyPairPT->Draw();
gPad->Print( "plots/plot_mNormToyPairPT.pdf" );
gPad->Print( "plots/plot_mNormToyPairPT.png" );

makeCanvas();
mDataPairPT->SetLineColor(kBlack);
mDataPairPT->SetTitle("Data pair P_{T} distribution; P_{T} (GeV); counts");
mDataPairPT->Draw();
gPad->Print( "plots/plot_mNormDataPairPT.pdf" );
gPad->Print( "plots/plot_mNormDataPairPT.png" );

makeCanvas();
mDataMinusToyPairPT->SetLineColor(kBlack);
mDataMinusToyPairPT->SetTitle("Data-toy model pair P_{T} distribution; P_{T} (GeV); counts");
mDataMinusToyPairPT->Draw();
gPad->Print( "plots/plot_mDataMinusToyPairPT.pdf" );
gPad->Print( "plots/plot_mDataminusToyPairPT.png" );

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
gPad->Print( "plots/plot_mNormToyPhi.pdf" );
gPad->Print( "plots/plot_mNormToyPhi.png" );

makeCanvas();
mDataPairPhi->SetLineColor(kBlack);
mDataPairPhi->SetTitle("Data #phi distribution; #phi (rad); counts");
mDataPairPhi->Draw();
gPad->Print( "plots/plot_mNormDataPhi.pdf" );
gPad->Print( "plots/plot_mNormDataPhi.png" );

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
gPad->Print( "plots/plot_mDataMinusToyPhi.pdf" );
gPad->Print( "plots/plot_mDataminusToyPhi.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataPhivsPT->SetTitle("Data #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV); counts");
mDataPhivsPT->Draw("colz");
gPad->Print( "plots/plot_mNormDataPhivsPT.pdf" );
gPad->Print( "plots/plot_mNormDataPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mToyPhivsPT->SetTitle("Toy model #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV); counts");
mToyPhivsPT->Draw("colz");
gPad->Print( "plots/plot_mNormToyPhivsPT.pdf" );
gPad->Print( "plots/plot_mNormToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataMinusToyPhivsPT->SetTitle("Data-toy model #phi vs. P_{T} (scaled to avg=1); #phi (rad); P_{T} (GeV); counts");
mDataMinusToyPhivsPT->Draw("colz");
gPad->Print( "plots/plot_mDataMinusToyPhivsPT.pdf" );
gPad->Print( "plots/plot_mDataMinusToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataCos2PhivsPT->SetTitle("Data cos(2#phi) vs. P_{T} (scaled to avg=1); cos(2#phi); P_{T} (GeV); counts");
mDataCos2PhivsPT->Draw("colz");
gPad->Print( "plots/plot_mNormDataCos2PhivsPT.pdf" );
gPad->Print( "plots/plot_mNormDataCos2PhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mToyCos2PhivsPT->SetTitle("Toy model cos(2#phi) vs. P_{T} (scaled to avg=1); cos(2#phi); P_{T} (GeV); counts");
mToyCos2PhivsPT->Draw("colz");
gPad->Print( "plots/plot_mNormToyCos2PhivsPT.pdf" );
gPad->Print( "plots/plot_mNormToyCos2PhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mDataMinusToyCos2PhivsPT->SetTitle("Data-toy model cos(2#phi) vs. P_{T}; cos(2#phi); P_{T} (GeV); counts");
mDataMinusToyCos2PhivsPT->Draw("colz");
gPad->Print( "plots/plot_mDataMinusToyCos2PhivsPT.pdf" );
gPad->Print( "plots/plot_mDataMinusToyCos2PhivsPT.png" );

makeCanvas();
mToyPTCos2PhiMoments->SetLineColor(kBlack);
mToyPTCos2PhiMoments->SetTitle("Toy model strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mToyPTCos2PhiMoments->Draw();

makeCanvas();
mDataPTCos2PhiMoments->SetLineColor(kBlack);
mDataPTCos2PhiMoments->SetTitle("Data strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mDataPTCos2PhiMoments->Draw();

makeCanvas();
mDataMinusToyPTCos2PhiMoments->SetLineColor(kBlack);
mDataMinusToyPTCos2PhiMoments->SetMinimum(-0.5);
mDataMinusToyPTCos2PhiMoments->SetMaximum(0.5);
mDataMinusToyPTCos2PhiMoments->SetTitle("Data minus toy model strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mDataMinusToyPTCos2PhiMoments->Draw();
gPad->Print( "plots/plot_mDataMinusToyPTCos2phiMoments.pdf" );
gPad->Print( "plots/plot_mDataMinusToyPTCos2phiMoments.png" );

makeCanvas();
mToyMassCos2PhiMoments->SetLineColor(kBlack);
mToyMassCos2PhiMoments->SetTitle("Toy model strength of cos(2#phi) signal vs. mass; mass (GeV); 2<cos(2#phi)>");
mToyMassCos2PhiMoments->Draw();

makeCanvas();
mDataMassCos2PhiMoments->SetLineColor(kBlack);
mDataMassCos2PhiMoments->SetTitle("Data strength of cos(2#phi) signal vs. mass; mass (GeV); 2<cos(2#phi)>");
mDataMassCos2PhiMoments->Draw();

makeCanvas();
mDataMinusToyMassCos2PhiMoments->SetLineColor(kBlack);
mDataMinusToyMassCos2PhiMoments->SetMinimum(-0.5);
mDataMinusToyMassCos2PhiMoments->SetMaximum(0.5);
mDataMinusToyMassCos2PhiMoments->SetTitle("Data minus toy model strength of cos(2#phi) signal vs. mass; mass (GeV); 2<cos(2#phi)>");
mDataMinusToyMassCos2PhiMoments->Draw();
gPad->Print( "plots/plot_mDataMinusToyMassCos2phiMoments.pdf" );
gPad->Print( "plots/plot_mDataMinusToyMassCos2phiMoments.png" );

}

