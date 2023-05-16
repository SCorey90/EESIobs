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
#include "my_functions.h"

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
    double theta = acos(rng.Uniform(-1, 1));
    double absP1 = sqrt((lv0->M()*lv0->M()/4) - daughter_mass*daughter_mass);
    lv1->SetPxPyPzE( (absP1*sin(theta)*cos(phi)), (absP1*sin(theta)*sin(phi)), (absP1*cos(theta)), (sqrt(absP1*absP1 + daughter_mass*daughter_mass)) );
    lv1->Boost(lv0->BoostVector());
    return lv1;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    TLorentzVector lvMinus = lv1 - lv2;
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double absPperp = pow((Px*Px)+(Py*Py), 0.5);
    double absQperp = pow((Qx*Qx)+(Qy*Qy), 0.5);
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double PdotQ = (Px*Qx) + (Py*Qy);
    double cosphi = (Px*Qx + Py*Qy) / (absPperp*absQperp);
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){ 
        return PairPhi - 3.141592;
    } else { 
        return 3.141592 - PairPhi;
    }
}
    

void toymodel() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile *fo = new TFile( "toymodelplots.root", "RECREATE" );

    auto *mRhoM = new TH1F("mRhoM", "#rho^{0} mass; mass (GeV); counts", 1000, 0, 1.5);
    auto *mReconstructedM = new TH1F("mReconstructedM", "Sum of #pi^{+}#pi^{-} mass; mass (GeV); counts", 1000, 0, 1.5);
    auto *mRhoPT = new TH1F("mRhoPT", "#rho^{0} P_{T}; #rho^{0} P_{T} (GeV); counts", 100, 0, 0.2);
    auto *mPairPT = new TH1F("mPairPT", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV); counts", 100, 0, 0.2);
    auto *mNoisyPairPT = new TH1F("mNoisyPairPT", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV); counts", 100, 0, 0.2);
    auto *mCutPairPT = new TH1F("mCutPairPT", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV); counts", 100, 0, 0.2);

    auto * mPairPhi = new TH1F("mPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 500, -3.13, 3.13);
    auto * mNoisyPairPhi = new TH1F("mNoisyPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 500, -3.13, 3.13);
    auto * mCutPairPhi = new TH1F("mCutPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 500, -3.13, 3.13);

    auto * mPhivsPT = new TH2F("mPhivsPT", "model #pi^{#pm} #phi distribution vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV); Counts", 100, -3.14, 3.14, 100, 0, 0.25);
    auto * mPolarPhivsPT = new TH2F("mPolarPhivsPT", "model #pi^{#pm} #phi distribution vs. parent P_{T}; P_{T}*cos(#phi) (GeV); P_{T}*sin(#phi); Counts", 100, -0.025, 0.025, 100, -0.025, 0.025);
    auto * mExagPolarPhivsPT = new TH2F("mExagPolarPhivsPT", "model #pi^{#pm} #phi distribution with 20% P_{T} smear vs. parent P_{T}; P_{T}*cos(#phi) (GeV); P_{T}*sin(#phi); Counts", 100, -0.25, 0.25, 100, -0.25, 0.25);
    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "model cos2#phi distribution vs P_{T}; 2cos2#phi; Parent P_{T} (GeV); Counts", 100, -2, 2, 100, 0, 0.25);

    auto * mGaus = new TH2F("mGaus", "Gaussian Testing; P_{T}; Gaus(0, 0.1*P_{T})", 100, 0, 0.2, 100, -0.2, 0.2);
    auto * mLowPTPhi = new TH1F("mLowPTPhi", "#phi hist at p_{T}<10 MeV; #phi (rad); counts", 100, -3.5, 3.5); 
    auto * mMidPTPhi = new TH1F("mMidPTPhi", "#phi hist at 50<p_{T}<150 MeV; #phi (rad); counts", 100, -3.5, 3.5);

    int n_events = 10000000;
    double m_pi = 0.139;
    for (int i = 0; i < n_events; i++) {
        TLorentzVector *lv0 = random_LV(0, 0.2, -6, 6, -3.141592, 3.141592, 2*m_pi, 1.5);
        TLorentzVector *lv1 = decay_LV( lv0, m_pi);

        TLorentzVector rho, pi1, pi2, gauspi1, gauspi2, lvRecon, exagpi1, exagpi2, lvExag;
        rho.SetPtEtaPhiM( lv0->Pt(), lv0->Eta(), lv0->Phi(), lv0->M() );
        pi1.SetPxPyPzE( lv1->Px(), lv1->Py(), lv1->Pz(), lv1->E() );
        pi2 = rho - pi1;
        double noise1 = rng.Gaus( 0, 0.01*pi1.Pt() );
        double noise2 = rng.Gaus( 0, 0.01*pi2.Pt() );
        gauspi1.SetPtEtaPhiM( pi1.Pt() + noise1, pi1.Eta(), pi1.Phi(), pi1.M() );
        gauspi2.SetPtEtaPhiM( pi2.Pt() + noise2, pi2.Eta(), pi2.Phi(), pi2.M() );
        exagpi1.SetPtEtaPhiM( pi1.Pt() + rng.Gaus( 0, 0.2*pi1.Pt() ), pi1.Eta(), pi1.Phi(), pi1.M() );
        exagpi2.SetPtEtaPhiM( pi2.Pt() + rng.Gaus( 0, 0.2*pi2.Pt() ), pi2.Eta(), pi2.Phi(), pi2.M() );
        lvRecon = gauspi1 + gauspi2;
        lvExag = exagpi1 + exagpi2;
        
        mGaus->Fill( pi1.Pt(), noise1 );
        mPairPT->Fill( (pi1 + pi2).Pt() );
        mNoisyPairPT->Fill( lvRecon.Pt() );
        
        mPairPhi->Fill( calc_Phi(pi1, pi2) );
        mNoisyPairPhi->Fill( calc_Phi(gauspi1, gauspi2) );

        mRhoM->Fill(rho.M());
        mRhoPT->Fill(rho.Pt());

        if ( gauspi1.Pt() > 0.1 && gauspi2.Pt() > 0.1 ) {
            double pairPhi = calc_Phi( gauspi1, gauspi2 );
            double reconPT = lvRecon.Pt();
            mReconstructedM->Fill(lvRecon.M()); 
            mCutPairPT->Fill(reconPT);
            mCutPairPhi->Fill( pairPhi );
            mPhivsPT->Fill( pairPhi, reconPT);
            mCos2phivsPT->Fill( 2*cos(2* pairPhi), reconPT);

            mPolarPhivsPT->Fill( reconPT*cos(pairPhi), reconPT*sin(pairPhi));
            mExagPolarPhivsPT->Fill( lvExag.Pt()*cos( calc_Phi( exagpi1, exagpi2)), lvExag.Pt()*sin(calc_Phi( exagpi1, exagpi2)));
            if ( reconPT < 0.01 ) { mLowPTPhi->Fill(pairPhi); }
            if ( reconPT > 0.05 && reconPT < 0.15 ) { mMidPTPhi->Fill(pairPhi); }
        }
    }

TH1F* mMassRCbyMC = (TH1F*)mReconstructedM->Clone("mMassMCbyRC");
TH1F* mPTRCbyMC = (TH1F*)mPairPT->Clone("mPTMCbyRC");

auto *mToyPTCos2phiMoments = mCos2phivsPT->ProfileY("mToyPTCos2phiMoments", 1, -1);

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
mNoisyPairPT->SetLineColor(kBlue);
mCutPairPT->SetLineColor(kOrange);
mNoisyPairPT->Draw();
mPairPT->Draw("SAME");
mCutPairPT->Draw("SAME");
mPairPT->SetMinimum(0);
auto legend = new TLegend(0.70,0.1,0.9,0.3);
legend->SetHeader("Legend","C"); // option "C" allows to center the header
legend->AddEntry(mPairPT,"No noise, no cut");
legend->AddEntry(mNoisyPairPT,"With noise, no cut");
legend->AddEntry(mCutPairPT,"With cut and noise");
legend->Draw();
gPad->Print( "plot_mToyPairPT.pdf" );
gPad->Print( "plot_mToyPairPT.png" );

makeCanvas();
mMassRCbyMC->Divide(mRhoM);
mMassRCbyMC->SetLineColor(kBlack);
mMassRCbyMC->SetTitle("#pi^{+}+#pi^{-}/#rho^{0}; Mass (GeV); ratio of counts");
mMassRCbyMC->Draw();
gPad->Print( "plot_mMassRCbyMC.pdf" );
gPad->Print( "plot_mMassRCbyMC.png" );

makeCanvas();
mPTRCbyMC->Divide(mRhoPT);
mPTRCbyMC->SetLineColor(kBlack);
mPTRCbyMC->SetTitle("#pi^{+}+#pi^{-}/#rho^{0}; P_{T} (GeV); ratio of counts");
mPTRCbyMC->SetMinimum(0);
mPTRCbyMC->Draw();
gPad->Print( "plot_mPTRCbyMC.pdf" );
gPad->Print( "plot_mPTRCbyMC.png" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mNoisyPairPhi->SetLineColor(kBlue);
mCutPairPhi->SetLineColor(kOrange);
gStyle->SetOptFit();
mNoisyPairPhi->SetMinimum(16000);
mNoisyPairPhi->Draw();
mPairPhi->Draw("Same");
mCutPairPhi->Draw("Same");
auto legend2 = new TLegend(0.7,0.1,0.9,0.3);
legend2->SetHeader("Legend","C"); // option "C" allows to center the header
legend2->AddEntry(mPairPhi,"No noise, no cut");
legend2->AddEntry(mNoisyPairPhi,"With noise, no cut");
legend2->AddEntry(mCutPairPhi,"With cut and noise");
legend2->Draw();
gPad->Print( "plot_mToyPairPhi.pdf" );
gPad->Print( "plot_mToyPairPhi.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsPT->Draw("colz");
gPad->Print( "plot_mToyPhivsPT.pdf" );
gPad->Print( "plot_mToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mPolarPhivsPT->Draw("colz");
gPad->Print( "plot_mToyPolarPhivsPT.pdf" );
gPad->Print( "plot_mToyPolarPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mExagPolarPhivsPT->Draw("colz");
gPad->Print( "plot_mToyExagPolarPhivsPT.pdf" );
gPad->Print( "plot_mToyExagPolarPhivsPT.png" );

makeCanvas();
mCos2phivsPT->Draw("colz");
gPad->SetLogz();
gPad->Print( "plot_mToyCos2phivsPT.pdf" );
gPad->Print( "plot_mToyCos2phivsPT.png" );

makeCanvas();
mToyPTCos2phiMoments->SetLineColor(kBlack);
mToyPTCos2phiMoments->SetTitle("Strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mToyPTCos2phiMoments->Draw();
gPad->Print( "plot_mToyPTCos2phiMoments.pdf" );
gPad->Print( "plot_mToyPTCos2phiMoments.png" );

makeCanvas();
mLowPTPhi->SetLineColor(kBlack);
mLowPTPhi->Draw();
gPad->Print( "plot_mToyLowPTPhi.pdf" );
gPad->Print( "plot_mToyLowPTPhi.png" );

makeCanvas();
mMidPTPhi->SetLineColor(kBlack);
mMidPTPhi->Draw();
gPad->Print( "plot_mToyMidPTPhi.pdf" );
gPad->Print( "plot_mToyMidPTPhi.png" );

fo->Write();
}
