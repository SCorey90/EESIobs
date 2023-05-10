#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "FemtoPairFormat.h"
#include "TRandom3.h"

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    can->SetRightMargin(0.35);
}

std::random_device global_rng;
TRandom3 rng(global_rng());

TLorentzVector *random_LV( double ptmin, double ptmax, double etamin, double etamax, double phimin, double phimax, double Mmin, double Mmax) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = rng.Uniform(ptmin, ptmax);
    double eta = rng.Uniform(etamin, etamax);
    double phi = rng.Uniform(phimin, phimax);
    double mass = rng.Uniform(Mmin, Mmax);
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
}

//Two body decay, but only returns one lv; the other is lv0-lv1
TLorentzVector *decay_LV( TLorentzVector *lv0, double daughter_mass ) {
    TLorentzVector *lv1 = new TLorentzVector;
    double phi = rng.Uniform(-3.141592, 3.141592);
    double theta = rng.Uniform(0, 3.141592);
    double absP1 = sqrt((lv0->M()*lv0->M()/4) - daughter_mass*daughter_mass);
    lv1->SetPxPyPzE( (absP1*sin(theta)*cos(phi)), (absP1*sin(theta)*sin(phi)), (absP1*cos(theta)), (sqrt(absP1*absP1 + daughter_mass*daughter_mass)) );
    lv1->Boost(lv0->BoostVector());
    return lv1;
}

void toymodel() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile *fo = new TFile( "toymodelplots.root", "RECREATE" );

    auto *mRhoM = new TH1F("mRhoM", "#rho^{0} mass; mass (GeV); counts", 1000, 0, 1.5);
    auto *mReconstructedM = new TH1F("mReconstructedM", "Sum of #pi^{+}#pi^{-} mass; mass (GeV); counts", 1000, 0, 1.5);
    auto *mRhoPT = new TH1F("mRhoPT", "#rho^{0} P_{T}; P_{T} (GeV); counts", 100, 0, 1.5);
    auto *mPairPT = new TH1F("mPairPT", "#pi^{+}#pi^{-} pair P_{T}; P_{T} (GeV); counts", 100, 0, 1.5);

    double m_pi = 0.139;
    for (int i = 0; i < 10000000; i++) {
        TLorentzVector *lv0 = random_LV(0, 1.5, -6, 6, -3.141592, 3.141592, 2*m_pi, 1.5);
        TLorentzVector *lv1 = decay_LV( lv0, m_pi);

        TLorentzVector rho, pi1, pi2, gauspi1, gauspi2, lvRecon;
        rho.SetPtEtaPhiM( lv0->Pt(), lv0->Eta(), lv0->Phi(), lv0->M() );
        pi1.SetPxPyPzE( lv1->Px(), lv1->Py(), lv1->Pz(), lv1->E() );
        pi2 = rho - pi1;
        gauspi1.SetPtEtaPhiM( pi1.Pt() + rng.Gaus( 0, 0.1 ), pi1.Eta(), pi1.Phi(), pi1.M() );
        gauspi2.SetPtEtaPhiM( pi2.Pt() + rng.Gaus( 0, 0.1 ), pi2.Eta(), pi2.Phi(), pi2.M() );
        lvRecon = gauspi1 + gauspi2;
        mRhoM->Fill(rho.M());
        mRhoPT->Fill(rho.Pt());
        if ( lvRecon.Pt() > 0.7 ) { 
            mReconstructedM->Fill(lvRecon.M()); 
            mPairPT->Fill(lvRecon.Pt());
        }
    }

fo -> cd();

makeCanvas();
mRhoM->SetLineColor(kBlack);
mRhoM->Draw();
gPad->Print( "plot_mToyRhoM.pdf" );
gPad->Print( "plot_mToyRhoM.png" );

makeCanvas();
mReconstructedM->SetLineColor(kBlack);
mReconstructedM->Draw();
gPad->Print( "plot_mToyReconstructedM.pdf" );
gPad->Print( "plot_mToyReconstructedM.png" );

makeCanvas();
mRhoPT->SetLineColor(kBlack);
mRhoPT->Draw();
gPad->Print( "plot_mToyRhoPT.pdf" );
gPad->Print( "plot_mToyRhoPT.png" );

makeCanvas();
mPairPT->SetLineColor(kBlack);
mPairPT->Draw();
gPad->Print( "plot_mToyPairPT.pdf" );
gPad->Print( "plot_mToyPairPT.png" );

fo->Write();
}
