#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include <vector>

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    can->SetRightMargin(0.35);
}

void MixedAnalysis() {

    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "MixedAnalysisplots.root", "RECREATE" );

    TFile * hists = TFile::Open("MixedEventplotsv2.root");

    hists->ls();
    auto mDataCos1phivsPTvsMass = (TH3F*)hists->Get("mDataCos1phivsPTvsMass");
    auto mDataCos2phivsPTvsMass = (TH3F*)hists->Get("mDataCos2phivsPTvsMass");
    auto mDataCos3phivsPTvsMass = (TH3F*)hists->Get("mDataCos3phivsPTvsMass");
    auto mDataCos4phivsPTvsMass = (TH3F*)hists->Get("mDataCos4phivsPTvsMass");

    auto mMixedCos1phivsPTvsMass = (TH3F*)hists->Get("mMixedCos1phivsPTvsMass");
    auto mMixedCos2phivsPTvsMass = (TH3F*)hists->Get("mMixedCos2phivsPTvsMass");
    auto mMixedCos3phivsPTvsMass = (TH3F*)hists->Get("mMixedCos3phivsPTvsMass");
    auto mMixedCos4phivsPTvsMass = (TH3F*)hists->Get("mMixedCos4phivsPTvsMass");
    
    //cuts
    double MASS_MIN = 0.76;
    double MASS_MAX = 0.78;

    double PT_MIN = 0;
    double PT_MAX = 0.06;

    //make plots
    auto * mPperp = mDataCos2phivsPTvsMass->ProjectionY("mPperp", 1, -1, 1, -1  );
    auto * mMass = mDataCos2phivsPTvsMass->ProjectionZ("mMass", 1, -1, 1, -1  );

    auto * mMixedMass = mMixedCos2phivsPTvsMass->ProjectionZ("mMixedMass", 1, -1, 1, -1  );

    (mDataCos1phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mDataCos1phivsPT = (TH2D*)mDataCos1phivsPTvsMass->Project3D("yx");
    auto * mDataCos1phivsPTprofile = mDataCos1phivsPT->ProfileY("mDataCos1phivsPTprofile", 1, -1);
    
    (mDataCos2phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mDataCos2phivsPT = (TH2D*)mDataCos2phivsPTvsMass->Project3D("yx");
    auto * mDataCos2phivsPTprofile = mDataCos2phivsPT->ProfileY("mDataCos2phivsPTprofile", 1, -1);
    
    (mDataCos3phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mDataCos3phivsPT = (TH2D*)mDataCos3phivsPTvsMass->Project3D("yx");
    auto * mDataCos3phivsPTprofile = mDataCos3phivsPT->ProfileY("mDataCos3phivsPTprofile", 1, -1);
    
    (mDataCos4phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mDataCos4phivsPT = (TH2D*)mDataCos4phivsPTvsMass->Project3D("yx");
    auto * mDataCos4phivsPTprofile = mDataCos4phivsPT->ProfileY("mDataCos4phivsPTprofile", 1, -1);


    (mMixedCos1phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mMixedCos1phivsPT = (TH2D*)mMixedCos1phivsPTvsMass->Project3D("yx");
    auto * mMixedCos1phivsPTprofile = mMixedCos1phivsPT->ProfileY("mMixedCos1phivsPTprofile", 1, -1);

    (mMixedCos2phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mMixedCos2phivsPT = (TH2D*)mMixedCos2phivsPTvsMass->Project3D("yx");
    auto * mMixedCos2phivsPTprofile = mMixedCos2phivsPT->ProfileY("mMixedCos2phivsPTprofile", 1, -1);

    (mMixedCos3phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mMixedCos3phivsPT = (TH2D*)mMixedCos3phivsPTvsMass->Project3D("yx");
    auto * mMixedCos3phivsPTprofile = mMixedCos3phivsPT->ProfileY("mMixedCos3phivsPTprofile", 1, -1);

    (mMixedCos4phivsPTvsMass->GetZaxis())->SetRange( mMass->FindBin( MASS_MIN ), mMass->FindBin( MASS_MAX ) );
    auto * mMixedCos4phivsPT = (TH2D*)mMixedCos4phivsPTvsMass->Project3D("yx");
    auto * mMixedCos4phivsPTprofile = mMixedCos4phivsPT->ProfileY("mMixedCos4phivsPTprofile", 1, -1);

auto * mCorrectedCos2phivsPTprofile = new TH1F("mCorrectedCos2phivsPTprofile", "Corrected  (0.76 < M_{#pi#pi} < 0.78 GeV); P_{T} (GeV/c); 2<cos2#phi>", 50, 0, 0.5);
for (int i = 0; i < (mCorrectedCos2phivsPTprofile->GetNbinsX()) - 1; i++) {
    double gamma_2 = mDataCos2phivsPTprofile->GetBinContent( i+1 )/2;
    double omega_2 = mMixedCos2phivsPTprofile->GetBinContent( mMixedCos2phivsPTprofile->FindBin( mDataCos2phivsPTprofile->GetBinCenter( i+1 ) ) )/2;
    double delta_gamma_2 =  mDataCos2phivsPTprofile->GetBinError( i+1 );
    double delta_omega_2 =  mMixedCos2phivsPTprofile->GetBinError( mMixedCos2phivsPTprofile->FindBin( mDataCos2phivsPTprofile->GetBinCenter( i+1 ) ) );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * gamma_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedCos2phivsPTprofile->SetBinContent( i+1, 2*alpha_2 );
    mCorrectedCos2phivsPTprofile->SetBinError( i+1, delta_alpha_2);
}

fo->cd();
/*
makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plots/mixedevent/Mass/plot_mMass.pdf" );
gPad->Print( "plots/mixedevent/Mass/plot_mMass.png" );

makeCanvas();
mMixedMass->SetLineColor(kBlack);
mMixedMass->Draw();
gPad->Print( "plots/mixedevent/Mass/plot_mMixedMass.pdf" );
gPad->Print( "plots/mixedevent/Mass/plot_mMixedMass.png" );

makeCanvas();
mDataCos1phivsPTprofile->SetLineColor(kGreen-2);
mMixedCos1phivsPTprofile->SetLineColor(kBlack);
mDataCos1phivsPTprofile->SetLineWidth(2);
mMixedCos1phivsPTprofile->SetLineWidth(2);
mDataCos1phivsPTprofile->SetTitle("2<cos1#phi> (0.75 < M_{#pi#pi} < 0.76 GeV); Pair P_{T} (GeV/c); 2<cos1#phi>");
mDataCos1phivsPTprofile->Draw();
mMixedCos1phivsPTprofile->Draw("SAME");
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->SetHeader("Legend","C"); // option "C" allows to center the header
legend->AddEntry(mDataCos1phivsPTprofile,"Data");
legend->AddEntry(mMixedCos1phivsPTprofile,"Mixed Event Like Sign Source");
legend->Draw();
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos1PhivsPTmoments.png" );
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos1PhivsPTmoments.pdf" );
*/
makeCanvas();
mDataCos2phivsPTprofile->SetLineColor(kGreen-2);
mMixedCos2phivsPTprofile->SetLineColor(kBlack);
mCorrectedCos2phivsPTprofile->SetLineColor(kBlue);
mDataCos2phivsPTprofile->SetLineWidth(2);
mMixedCos2phivsPTprofile->SetLineWidth(2);
mCorrectedCos2phivsPTprofile->SetLineWidth(2);
mDataCos2phivsPTprofile->SetTitle("2<cos2#phi> (0.76 < M_{#pi#pi} < 0.78 GeV); Pair P_{T} (GeV/c); 2<cos2#phi>");
mDataCos2phivsPTprofile->Draw();
mMixedCos2phivsPTprofile->Draw("Same");
mCorrectedCos2phivsPTprofile->Draw("Same");
TLine *l2 = new TLine(0,0,0.5,0);
l2->SetLineColor(kBlack);
l2->SetLineStyle(kDashed);
l2->Draw("Same");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->SetHeader("Legend","C"); // option "C" allows to center the header
legend2->AddEntry(mDataCos2phivsPTprofile,"Data");
legend2->AddEntry(mMixedCos2phivsPTprofile,"Mixed Event Like Sign");
legend2->AddEntry(mCorrectedCos2phivsPTprofile,"Corrected");
legend2->Draw();
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos2PhivsPTmoments.png" );
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos2PhivsPTmoments.pdf" );

/*
makeCanvas();
mDataCos3phivsPTprofile->SetLineColor(kGreen-2);
mMixedCos3phivsPTprofile->SetLineColor(kBlack);
mDataCos3phivsPTprofile->SetLineWidth(2);
mMixedCos3phivsPTprofile->SetLineWidth(2);
mDataCos3phivsPTprofile->SetTitle("2<cos3#phi> (0.75 < M_{#pi#pi} < 0.76 GeV); Pair P_{T} (GeV/c); 2<cos3#phi>");
mDataCos3phivsPTprofile->Draw();
mMixedCos3phivsPTprofile->Draw("SAME");
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->SetHeader("Legend","C"); // option "C" allows to center the header
legend3->AddEntry(mDataCos3phivsPTprofile,"Data");
legend3->AddEntry(mMixedCos3phivsPTprofile,"Mixed Event Like Sign Source");
legend3->Draw();
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos3PhivsPTmoments.png" );
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos3PhivsPTmoments.pdf" );

makeCanvas();
mDataCos4phivsPTprofile->SetLineColor(kGreen-2);
mMixedCos4phivsPTprofile->SetLineColor(kBlack);
mDataCos4phivsPTprofile->SetLineWidth(2);
mMixedCos4phivsPTprofile->SetLineWidth(2);
mDataCos4phivsPTprofile->SetTitle("2<cos4#phi> (0.75 < M_{#pi#pi} < 0.76 GeV); Pair P_{T} (GeV/c); 2<cos4#phi>");
mMixedCos4phivsPTprofile->Draw();
mDataCos4phivsPTprofile->Draw("SAME");
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->SetHeader("Legend","C"); // option "C" allows to center the header
legend4->AddEntry(mDataCos4phivsPTprofile,"Data");
legend4->AddEntry(mMixedCos4phivsPTprofile,"Mixed Event Like Sign Source");
legend4->Draw();
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos4PhivsPTmoments.png" );
gPad->Print("plots/mixedevent/PT/plot_mDatavMixedCos4PhivsPTmoments.pdf" );
*/
fo->Write();
}
