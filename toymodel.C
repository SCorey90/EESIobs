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
TRandom ptr_rng(global_rng());

TLorentzVector *random_LV( double ptmin, double ptmax, double etamin, double etamax, double phimin, double phimax, double Mmin, double Mmax) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = rng.Uniform(ptmin, ptmax);
    double eta = rng.Uniform(etamin, etamax);
    double phi = rng.Uniform(phimin, phimax);
    double mass = rng.Uniform(Mmin, Mmax);
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
}

//TLorentzVector *hist_sample_LV( TH1F * PtHist, double etamin, double etamax, double phimin, double phimax, double Mmin, double Mmax) {
//    TLorentzVector *lv0 = new TLorentzVector;
//    double pt = PtHist->GetRandom();
//    double eta = rng.Uniform(etamin, etamax);
//    double phi = rng.Uniform(phimin, phimax);
//    double mass = rng.Uniform(Mmin, Mmax);
//    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
//    return lv0;
//}

TLorentzVector *hist_sample_LV( TH1F * PtHist, TH1F * EtaHist, TH1F * LVPhiHist, TH1F * MassHist) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = PtHist->GetRandom();
    double eta = EtaHist->GetRandom();
    double phi = LVPhiHist->GetRandom();
    double mass = MassHist->GetRandom();
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
    lv1->Boost((lv0->BoostVector()));
    return lv1;
}

TLorentzVector *asym_decay( TLorentzVector lv1, double mass1, double mass2, double p_decay ) {
    double decay_check = rng.Uniform(0, 1);
    TLorentzVector *new_lv = new TLorentzVector;
    if ( decay_check >= p_decay ) { new_lv->SetPxPyPzE( lv1.Px(), lv1.Py(), lv1.Pz(), lv1.E() ); }
    if ( decay_check < p_decay ) {
        double daughter1_E = ((lv1.M())*(lv1.M()) + (mass1*mass1) - (mass2*mass2)) / (2*lv1.M());
        double daughter1_absP = sqrt( (daughter1_E*daughter1_E) - (mass1*mass1) );
        double phi = rng.Uniform(-3.141592, 3.141592);
        double theta = acos(rng.Uniform(-1, 1));
        new_lv->SetPxPyPzE( daughter1_absP*sin(theta)*cos(phi) , daughter1_absP*sin(theta)*sin(phi) , daughter1_absP*cos(theta) , daughter1_E );
        new_lv->Boost((lv1.BoostVector()));
    }
    return new_lv;
}

double calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    lv1.Boost(-lvPlus.BoostVector());
    lv2.Boost(-lvPlus.BoostVector());
    TLorentzVector lvMinus = lv1 - lv2;
    double Px = lvPlus.Px();
    double Py = lvPlus.Py();
    double Qx = lvMinus.Px();
    double Qy = lvMinus.Py();
    double PcrossQ = (Px*Qy) - (Py*Qx);
    double cosphi = (Px*Qx + Py*Qy) / (lvPlus.Pt()*lvMinus.Pt());
    double PairPhi = acos(cosphi);
    if ( PcrossQ > 0 ){ 
        return PairPhi - 3.141592;
    } else { 
        return 3.141592 - PairPhi;
    }
}
    

void toymodel() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * dataHists = TFile::Open("EESIplots.root");

    dataHists->ls();
    auto * mDataMass = (TH1F*)dataHists->Get("mMass");
    auto * mDataPairPT = (TH1F*)dataHists->Get("mPperp");
    auto * mDataEta = (TH1F*)dataHists->Get("mEta");
    auto * mDataLVPhi = (TH1F*)dataHists->Get("mLVPhi");

    TFile *fo = new TFile( "toymodelplots.root", "RECREATE" );

    auto *mRhoM = new TH1F("mRhoM", "#rho^{0} mass; mass (GeV); counts", 1000, 0, 1.5);
    auto *mReconstructedM = new TH1F("mReconstructedM", "Sum of #pi^{+}#pi^{-} mass; mass (GeV); counts", 500, 0, 2);
    auto *mRhoPT = new TH1F("mRhoPT", "#rho^{0} P_{T}; #rho^{0} P_{T} (GeV/c)`; counts", 100, 0, 1.5);
    auto *mPairPT = new TH1F("mPairPT", "P_{T}(#pi^{+}+#pi^{-}); P_{T} (Rad); counts", 500, 0, 1.5);
    auto *mPairPTwMu = new TH1F("mPairPTwMu", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV/c); counts", 500, 0, 1.5);
    auto *mNoisyPairPT = new TH1F("mNoisyPairPT", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV/c); counts", 500, 0, 1.5);
    auto *mCutPairPT = new TH1F("mCutPairPT", "#pi^{+}#pi^{-} pair P_{T}; pair P_{T} (GeV/c); counts", 500, 0, 1.5);

    auto * mPi1PT = new TH1F("mPi1PT", "model #pi^{#pm} P_{T} distribution;P_{T} (GeV/c);# events", 500, 0, 1.5);
    auto * mPi1PTwMu = new TH1F("mPi1PTwMu", "model #pi^{#pm} P_{T} distribution;P_{T} (GeV/c);# events", 500, 0, 1.5);
    auto * mNoisyPi1PT = new TH1F("mNoisyPi1PT", "model #pi^{#pm} P_{T} distribution;P_{T} (GeV/c);# events", 500, 0, 1.5);
    auto * mCutPi1PT = new TH1F("mCutPi1PT", "model #pi^{#pm} P_{T} distribution;P_{T} (GeV/c);# events", 500, 0, 1.5);

    auto * mPairQT = new TH1F("mPairQT", "P_{T}((#pi^{+} - #pi^{-}); P_{T} (GeV/c); counts", 100, 0, 1.5);    
    auto * mPairQTwMu = new TH1F("mPairQTwMu", "P_{T}((#pi^{+} - #pi^{-}); P_{T} (GeV/c); counts", 100, 0, 1.5);    
    auto * mNoisyPairQT = new TH1F("mNoisyPairQT", "P_{T}((#pi^{+} - #pi^{-}); P_{T} (GeV/c); counts", 100, 0, 1.5);    
    auto * mCutPairQT = new TH1F("mCutPairQT", "P_{T}((#pi^{+} - #pi^{-}); P_{T} (GeV/c); counts", 100, 0, 1.5);    

    auto * mPairPhi = new TH1F("mPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);
    auto * mPairPhiwMu = new TH1F("mPairPhiwMu", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);
    auto * mNoisyPairPhi = new TH1F("mNoisyPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);
    auto * mCutPairPhi = new TH1F("mCutPairPhi", "model #pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);

    auto * mPhivsPT = new TH2F("mPhivsPT", "model #pi^{#pm} #phi distribution vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV/c); Counts", 100, -3.14, 3.14, 100, 0, 0.25);
    auto * mPolarPhivsPT = new TH2F("mPolarPhivsPT", "model #pi^{#pm} #phi distribution vs. parent P_{T}; P_{T}*cos(#phi) (GeV); P_{T}*sin(#phi); Counts", 100, -0.025, 0.025, 100, -0.025, 0.025);
    auto * mExagPolarPhivsPT = new TH2F("mExagPolarPhivsPT", "model #pi^{#pm} #phi distribution with 20% P_{T} smear vs. parent P_{T}; P_{T}*cos(#phi) (GeV/c); P_{T}*sin(#phi) (GeV/c); Counts", 100, -0.25, 0.25, 100, -0.25, 0.25);
    
    auto * mCos2flatPhi = new TH1F("mCos2flatPhi", "Cos2#phi from flat #phi distribution; #phi (rad); counts", 100, -1, 1);
    auto * mCos2flatphivsPT = new TH2F("mCos2flatvsPT", "cos2#phi from flat #phi distribution vs P_{T}; 2cos2#phi; Parent P_{T} (GeV/c); Counts", 100, -2, 2, 100, 0, 0.5);
    auto * mCos2decayphivsPT = new TH2F("mCos2decayphivsPT", "cos2#phi from flat #phi distribution w/ possibility of #mu decay vs P_{T}; 2cos2#phi; Parent P_{T} (GeV/c); Counts", 100, -2, 2, 100, 0, 0.5);
    auto * mCos2cutflatphivsPT = new TH2F("mCos2cutflatvsPT", "cos2#phi from flat #phi distribution vs P_{T} w/ acceptance cut; 2cos2#phi; Parent P_{T} (GeV/c); Counts", 100, -2, 2, 100, 0, 0.5);
    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "model cos2#phi distribution vs P_{T}; 2cos2#phi; Parent P_{T} (GeV/c); Counts", 100, -2, 2, 100, 0, 0.5);
    auto * mCos4phivsPT = new TH2F("mCos4phivsPT", "model cos4#phi distribution vs P_{T}; 4cos4#phi; Parent P_{T} (GeV/c); Counts", 100, -2, 2, 100, 0, 0.5);

    auto * mCos2flatphivsMass = new TH2F("mCos2flatvsMass", "cos2#phi from flat #phi distribution vs. mass; 2cos2#phi; Parent Mass (GeV); Counts", 100, -2, 2, 100, 0.2, 2);
    auto * mCos2cutflatphivsMass = new TH2F("mCos2cutflatvsMass", "cos2#phi from flat #phi distribution vs. mass w/ acceptance cut; 2cos2#phi; Parent Mass (GeV); Counts", 100, -2, 2, 100, 0.2, 2);
    auto * mPhivsMass = new TH2F("mPhivsMass", "model #pi^{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 100, 0, 2);
    auto * mCos2phivsMass = new TH2F("mCos2phivsMass", "model cos2#phi distribution vs mass; 2cos2#phi; Parent Mass (GeV); Counts", 100, -2, 2, 100, 0.2, 2);
    auto * mCos4phivsMass = new TH2F("mCos4phivsMass", "model cos4#phi distribution vs mass; 4cos4#phi; Parent Mass (GeV); Counts", 100, -2, 2, 100, 0.2, 2);

    auto * mGaus = new TH2F("mGaus", "Gaussian Testing; P_{T}; Gaus(0, 0.1*P_{T})", 100, 0, 0.2, 100, -0.2, 0.2);
    auto * mLowPTPhi = new TH1F("mLowPTPhi", "#phi hist at p_{T}<60 MeV/c; #phi (rad); counts", 100, -3.5, 3.5); 
    auto * mMidPTPhi = new TH1F("mMidPTPhi", "#phi hist at p_{T}>900 MeV/c; #phi (rad); counts", 100, -3.5, 3.5);
    auto * mHighMidPTPhi = new TH1F("mMidPTPhi", "#phi hist at p_{T}>900 MeV/c; #phi (rad); counts", 100, -3.5, 3.5);

    int n_events = 10000000;
    double m_pi = 0.139;
    double m_mu = 0.105;
    double m_nu_mu = 1.2 * pow(10, -10);
    for (int i = 0; i < n_events; i++) {
        //random parent 4-vector and simulate decay
        TLorentzVector *lv0 = hist_sample_LV(mDataPairPT, mDataEta, mDataLVPhi, mDataMass);
        TLorentzVector *lv1 = decay_LV( lv0, m_pi);

        TLorentzVector rho, pi1, pi2, gauspi1, gauspi2, lvRecon, exagpi1, exagpi2, lvExag, pi_or_mu1, pi_or_mu2;
        rho.SetPtEtaPhiM( lv0->Pt(), lv0->Eta(), lv0->Phi(), lv0->M() );
        pi1.SetPxPyPzE( lv1->Px(), lv1->Py(), lv1->Pz(), lv1->E() );
        pi2 = rho - pi1;

        //add possibility of pi->mu+nu
        TLorentzVector *decay_lv1 = asym_decay( pi1, m_mu, m_nu_mu, 0.01 );
        TLorentzVector *decay_lv2 = asym_decay( pi2, m_mu, m_nu_mu, 0.01 );
        pi_or_mu1.SetPtEtaPhiM( decay_lv1->Pt(), decay_lv1->Eta(), decay_lv1->Phi(), decay_lv1->M() );
        pi_or_mu2.SetPtEtaPhiM( decay_lv2->Pt(), decay_lv2->Eta(), decay_lv2->Phi(), decay_lv2->M() );

        //add gaussian noise to pi P_T
        double noise1 = rng.Gaus( 0, 0.04*pi1.Pt() );
        double noise2 = rng.Gaus( 0, 0.04*pi2.Pt() );
        gauspi1.SetPtEtaPhiM( pi1.Pt() + noise1, pi1.Eta(), pi1.Phi(), pi1.M() );
        gauspi2.SetPtEtaPhiM( pi2.Pt() + noise2, pi2.Eta(), pi2.Phi(), pi2.M() );
        pi_or_mu1.SetPtEtaPhiM( pi_or_mu1.Pt() + rng.Gaus( 0, 0.04*pi_or_mu1.Pt() ), pi_or_mu1.Eta(), pi_or_mu1.Phi(), pi_or_mu1.M() );
        pi_or_mu2.SetPtEtaPhiM( pi_or_mu2.Pt() + rng.Gaus( 0, 0.04*pi_or_mu2.Pt() ), pi_or_mu2.Eta(), pi_or_mu2.Phi(), pi_or_mu2.M() );
        lvRecon = gauspi1 + gauspi2;
        lvExag = exagpi1 + exagpi2;
        
        //fill histograms without pT cut 
        mGaus->Fill( pi1.Pt(), noise1 );
        mPairPT->Fill( (pi1 + pi2).Pt() );
        mPairPTwMu->Fill( (pi_or_mu1 + pi_or_mu2).Pt() );
        mNoisyPairPT->Fill( lvRecon.Pt() );
        
        mPi1PT->Fill( pi1.Pt() );
        mPi1PTwMu->Fill( pi_or_mu1.Pt() );
        mNoisyPi1PT->Fill( gauspi1.Pt() );

        mPairQT->Fill( (pi1 - pi2).Pt() );
        mPairQTwMu->Fill( (pi_or_mu1-pi_or_mu2).Pt() );
        mNoisyPairQT->Fill( (gauspi1 - gauspi2).Pt() );

        mPairPhi->Fill( calc_Phi(pi1, pi2) );
        mCos2flatPhi->Fill( cos(2*calc_Phi(pi1, pi2)) );
        mCos2flatphivsPT->Fill( 2*cos(2*calc_Phi(pi1, pi2)), (pi1+pi2).Pt());

        //mCos2decayphivsPT->Fill( 2*cos(2*calc_Phi(pi_or_mu1, pi_or_mu2)) , (pi_or_mu1+pi_or_mu2).Pt());

        mCos2flatphivsMass->Fill( 2*cos(2*calc_Phi(pi1, pi2)), (pi1+pi2).M());        

        if ( pi1.Pt() > 0.1 && pi2.Pt() > 0.1 ) {
             mCos2cutflatphivsPT->Fill( 2*cos(2*calc_Phi(pi1, pi2)), (pi1+pi2).Pt());
             mCos2cutflatphivsMass->Fill( 2*cos(2*calc_Phi(pi1, pi2)), (pi1+pi2).M());
        }
        
        mPairPhiwMu->Fill( calc_Phi(pi_or_mu1, pi_or_mu2) );
        mNoisyPairPhi->Fill( calc_Phi(gauspi1, gauspi2) );

        mRhoM->Fill(rho.M());
        mRhoPT->Fill(rho.Pt());

        //"" with pT cut
        if ( gauspi1.Pt() > 0.1 && gauspi2.Pt() > 0.1 ) {
            double pairPhi = calc_Phi( gauspi1, gauspi2 );
            double reconPT = lvRecon.Pt();
            mReconstructedM->Fill(lvRecon.M()); 
            mCutPi1PT->Fill( gauspi1.Pt() );
            mCutPairPT->Fill(reconPT);
            mCutPairQT->Fill( (gauspi1 - gauspi2).Pt() );
            
            //if (lvRecon.M() > 0.65 && lvRecon.M() < 0.9){
                if ( reconPT < 0.06 ) { mCutPairPhi->Fill( pairPhi ); }
                //mCutPairPhi->Fill( pairPhi );
                mPhivsPT->Fill( pairPhi, reconPT);
                mCos2phivsPT->Fill( 2*cos(2* pairPhi), reconPT);
                mCos4phivsPT->Fill( 2*cos(4* pairPhi), reconPT);
            //}

            //if ( reconPT < 0.06 ){
                mPhivsMass->Fill( pairPhi, lvRecon.M());
                mCos2phivsMass->Fill( 2*cos(2* pairPhi), lvRecon.M());
                mCos4phivsMass->Fill( 2*cos(4* pairPhi), lvRecon.M());

               mPolarPhivsPT->Fill( reconPT*cos(pairPhi), reconPT*sin(pairPhi));
            //}
            mExagPolarPhivsPT->Fill( lvExag.Pt()*cos( calc_Phi( exagpi1, exagpi2)), lvExag.Pt()*sin(calc_Phi( exagpi1, exagpi2)));
            if ( reconPT < 0.06 ) { mLowPTPhi->Fill(pairPhi); }
            if ( reconPT > 0.9 ) { mMidPTPhi->Fill(pairPhi); }
        }
        if ( pi_or_mu1.Pt() > 0.1 && pi_or_mu2.Pt() > 0.1 ) { mCos2decayphivsPT->Fill( 2*cos(2*calc_Phi(pi_or_mu1, pi_or_mu2)) , (pi_or_mu1+pi_or_mu2).Pt()); }
    }

TH1F* mMassRCbyMC = (TH1F*)mReconstructedM->Clone("mMassMCbyRC");
TH1F* mPTRCbyMC = (TH1F*)mPairPT->Clone("mPTMCbyRC");

auto *mToyPTCos2flatphiMoments = mCos2flatphivsPT->ProfileY("mToyPTCos2flatphiMoments", 1, -1);

auto *mToyPTCos2decayphiMoments = mCos2decayphivsPT->ProfileY("mToyPTCos2decayphiMoments", 1, -1);

auto *mToyPTCos2cutflatphiMoments = mCos2cutflatphivsPT->ProfileY("mToyPTCos2cutflatphiMoments", 1, -1);

auto *mToyPTCos2phiMoments = mCos2phivsPT->ProfileY("mToyPTCos2phiMoments", 1, -1);
auto *mToyPTCos4phiMoments = mCos4phivsPT->ProfileY("mToyPTCos4phiMoments", 1, -1);

auto *mToyMassCos2flatphiMoments = mCos2flatphivsMass->ProfileY("mToyMassCos2flatphiMoments", 1, -1);

auto *mToyMassCos2cutflatphiMoments = mCos2cutflatphivsMass->ProfileY("mToyMassCos2cutflatphiMoments", 1, -1);

auto *mToyMassCos2phiMoments = mCos2phivsMass->ProfileY("mToyMassCos2phiMoments", 1, -1);
auto *mToyMassCos4phiMoments = mCos4phivsMass->ProfileY("mToyMassCos4phiMoments", 1, -1);

fo -> cd();

makeCanvas();
mRhoM->SetLineColor(kBlack);
mRhoM->Draw();
gPad->Print( "plots/model/Mass/plot_mToyRhoM.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyRhoM.png" );

makeCanvas();
mReconstructedM->SetLineColor(kBlack);
mReconstructedM->Draw();
gPad->Print( "plots/model/Mass/plot_mToyReconstructedM.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyReconstructedM.png" );

makeCanvas();
mRhoPT->SetLineColor(kBlack);
mRhoPT->Draw();
gPad->Print( "plots/model/PT/plot_mToyRhoPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyRhoPT.png" );

makeCanvas();
mPairPT->SetLineColor(kBlack);
mPairPTwMu->SetLineColor(kRed);
mNoisyPairPT->SetLineColor(kBlue);
mCutPairPT->SetLineColor(kOrange);
mNoisyPairPT->Draw();
mPairPT->Draw("SAME");
mPairPTwMu->Draw("SAME");
mCutPairPT->Draw("SAME");
mNoisyPairPT->SetMinimum(0);;
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->SetHeader("Legend","C"); // option "C" allows to center the header
legend->AddEntry(mPairPT,"No #pi->#mu decay, no noise, no cut");
legend->AddEntry(mPairPTwMu,"With #pi->#mu decay, no noise, no cut");
legend->AddEntry(mNoisyPairPT,"With #pi->#mu decay, with noise, no cut");
legend->AddEntry(mCutPairPT,"With #pi->#mu decay, with noise, with cut");
legend->Draw();
gPad->Print( "plots/model/PT/plot_mToyPairPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPairPT.png" );

makeCanvas();
mPi1PT->SetLineColor(kBlack);
mPi1PTwMu->SetLineColor(kRed);
mNoisyPi1PT->SetLineColor(kBlue);
mCutPi1PT->SetLineColor(kOrange);
mPi1PT->Draw();
mNoisyPi1PT->Draw("SAME");
mPi1PTwMu->Draw("SAME");
mCutPi1PT->Draw("SAME");
mPi1PT->SetMinimum(0);
mPi1PT->SetMaximum(70000);
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->SetHeader("Legend","C"); // option "C" allows to center the header
legend3->AddEntry(mPi1PT,"No #pi->#mu decay, no noise, no cut");
legend3->AddEntry(mPi1PTwMu,"With #pi->#mu decay, no noise, no cut");
legend3->AddEntry(mNoisyPi1PT,"With #pi->#mu decay, with noise, no cut");
legend3->AddEntry(mCutPi1PT,"With #pi->#mu decay, with noise, with cut");
legend3->Draw();
gPad->Print( "plots/model/PT/plot_mToyPi1PT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPi1PT.png" );

makeCanvas();
mPairQT->SetLineColor(kBlack);
mPairQTwMu->SetLineColor(kRed);
mNoisyPairQT->SetLineColor(kBlue);
mCutPairQT->SetLineColor(kOrange);
mPairQT->Draw();
mNoisyPairQT->Draw("SAME");
mPairQTwMu->Draw("SAME");
mCutPairQT->Draw("SAME");
mPairQT->SetMinimum(0);
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->SetHeader("Legend","C"); // option "C" allows to center the header
legend4->AddEntry(mPairQT,"No #pi->#mu decay, no noise, no cut");
legend4->AddEntry(mPairQTwMu,"With #pi->#mu decay, no noise, no cut");
legend4->AddEntry(mNoisyPairQT,"With #pi->#mu decay, with noise, no cut");
legend4->AddEntry(mCutPairQT,"With #pi->#mu decay, with noise, with cut");
legend4->Draw();
gPad->Print( "plots/model/PT/plot_mToyPairQT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPairQT.png" );

makeCanvas();
mMassRCbyMC->Divide(mRhoM);
mMassRCbyMC->SetLineColor(kBlack);
mMassRCbyMC->SetTitle("#pi^{+}+#pi^{-}/#rho^{0}; Mass (GeV); ratio of counts");
mMassRCbyMC->Draw();
gPad->Print( "plots/model/Mass/plot_mMassRCbyMC.pdf" );
gPad->Print( "plots/model/Mass/plot_mMassRCbyMC.png" );

makeCanvas();
mPTRCbyMC->Divide(mRhoPT);
mPTRCbyMC->SetLineColor(kBlack);
mPTRCbyMC->SetTitle("#pi^{+}+#pi^{-}/#rho^{0}; P_{T} (GeV); ratio of counts");
mPTRCbyMC->SetMinimum(0);
mPTRCbyMC->Draw();
gPad->Print( "plots/model/PT/plot_mPTRCbyMC.pdf" );
gPad->Print( "plots/model/PT/plot_mPTRCbyMC.png" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mPairPhiwMu->SetLineColor(kRed);
mNoisyPairPhi->SetLineColor(kBlue);
mCutPairPhi->SetLineColor(kOrange);
gStyle->SetOptFit();
mNoisyPairPhi->SetMinimum(10000);
mNoisyPairPhi->Draw();
mPairPhi->Draw("Same");
mPairPhiwMu->Draw("Same");
mCutPairPhi->Draw("Same");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->SetHeader("Legend","C"); // option "C" allows to center the header
legend2->AddEntry(mPairPhi,"No #pi->#mu decay, no noise, no cut");
legend2->AddEntry(mPairPhiwMu,"With #pi->#mu decay, no noise, no cut");
legend2->AddEntry(mNoisyPairPhi,"With #pi->#mu decay, with noise, no cut");
legend2->AddEntry(mCutPairPhi,"With #pi->#mu decay, with noise, with cut");
legend2->Draw();
gPad->Print( "plots/model/plot_mToyPairPhi.pdf" );
gPad->Print( "plots/model/plot_mToyPairPhi.png" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mPairPhi->Draw();
gPad->Print( "plots/model/plot_mToyPurePhi.pdf" );
gPad->Print( "plots/model/plot_mToyPurePhi.png" );

makeCanvas();
mCos2flatPhi->SetLineColor(kBlack);
mCos2flatPhi->Draw();
gPad->Print( "plots/model/plot_mToyPureCos2Phi.pdf" );
gPad->Print( "plots/model/plot_mToyPureCos2Phi.png" );

makeCanvas();
mCos2flatphivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyCos2flatphivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyCos2flatphivsPT.png" );

makeCanvas();
mToyPTCos2flatphiMoments->SetLineColor(kBlack);
mToyPTCos2flatphiMoments->SetTitle("No detector effects strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mToyPTCos2flatphiMoments->Draw();
gPad->Print( "plots/model/PT/plot_mToyPTCos2flatphiMoments.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPTCos2flatphiMoments.png" );

makeCanvas();
mCos2cutflatphivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyCos2cutflatphivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyCos2cutflatphivsPT.png" );

makeCanvas();
mToyPTCos2cutflatphiMoments->SetLineColor(kBlack);
mToyPTCos2cutflatphiMoments->SetTitle("P_{T} acceptance cut strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mToyPTCos2cutflatphiMoments->Draw();
gPad->Print( "plots/model/PT/plot_mToyPTCos2cutflatphiMoments.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPTCos2cutflatphiMoments.png" );

makeCanvas();
mCutPairPhi->SetLineColor(kBlack);
mCutPairPhi->Draw();
gPad->Print( "plots/model/plot_mToyCutPhi.pdf" );
gPad->Print( "plots/model/plot_mToyCutPhi.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyPhivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mPolarPhivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyPolarPhivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPolarPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mExagPolarPhivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyExagPolarPhivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyExagPolarPhivsPT.png" );

makeCanvas();
mCos2phivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyCos2phivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyCos2phivsPT.png" );

makeCanvas();
mCos4phivsPT->Draw("colz");
gPad->Print( "plots/model/PT/plot_mToyCos4phivsPT.pdf" );
gPad->Print( "plots/model/PT/plot_mToyCos4phivsPT.png" );

makeCanvas();
mToyPTCos2flatphiMoments->SetLineColor(kBlack);
mToyPTCos2decayphiMoments->SetLineColor(kMagenta);
mToyPTCos2cutflatphiMoments->SetLineColor(kOrange);
mToyPTCos2phiMoments->SetLineColor(kBlue);
mToyPTCos2phiMoments->SetTitle("Strength of cos(2#phi) signal vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mToyPTCos2phiMoments->SetMinimum(-1);
mToyPTCos2phiMoments->SetMaximum(0.25);
mToyPTCos2phiMoments->Draw();
mToyPTCos2flatphiMoments->Draw("SAME");
mToyPTCos2decayphiMoments->Draw("SAME");
mToyPTCos2cutflatphiMoments->Draw("SAME");
auto legend5 = new TLegend(0.65,0.1,0.95,0.4);
legend5->SetHeader("Legend","C"); // option "C" allows to center the header
legend5->AddEntry(mToyPTCos2flatphiMoments,"No effects (completely flat #phi hist)");
legend5->AddEntry(mToyPTCos2decayphiMoments,"#pi->#mu+#nu with P_{T} smear and detector acceptance");
legend5->AddEntry(mToyPTCos2cutflatphiMoments,"Detector acceptance");
legend5->AddEntry(mToyPTCos2phiMoments,"P_{T} smear, detector acceptance");
legend5->Draw();
gPad->Print( "plots/model/PT/plot_mToyPTCos2phiMoments.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPTCos2phiMoments.png" );

makeCanvas();
mToyPTCos4phiMoments->SetLineColor(kOrange);
mToyPTCos4phiMoments->SetTitle("Strength of cos(4#phi) signal vs. P_{T}; P_{T} (GeV); 4<cos(4#phi)>");
mToyPTCos4phiMoments->Draw();
gPad->Print( "plots/model/PT/plot_mToyPTCos4phiMoments.pdf" );
gPad->Print( "plots/model/PT/plot_mToyPTCos4phiMoments.png" );

makeCanvas();
mLowPTPhi->SetLineColor(kBlack);
mLowPTPhi->Draw();
gPad->Print( "plots/model/PT/plot_mToyLowPTPhi.pdf" );
gPad->Print( "plots/model/PT/plot_mToyLowPTPhi.png" );

makeCanvas();
mMidPTPhi->SetLineColor(kBlack);
mMidPTPhi->Draw();
gPad->Print( "plots/model/PT/plot_mToyMidPTPhi.pdf" );
gPad->Print( "plots/model/PT/plot_mToyMidPTPhi.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsMass->Draw("colz");
gPad->Print( "plots/model/Mass/plot_mToyPhivsMass.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyPhivsMass.png" );

makeCanvas();
mCos2phivsMass->Draw("colz");
gPad->Print( "plots/model/Mass/plot_mToyCos2phivsMass.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyCos2phivsMass.png" );

makeCanvas();
mToyMassCos2flatphiMoments->SetLineColor(kBlack);
//mToyMassCos2decayphiMoments->SetLineColor(kPink);
mToyMassCos2cutflatphiMoments->SetLineColor(kOrange);
mToyMassCos2phiMoments->SetLineColor(kBlue);
mToyMassCos2phiMoments->SetTitle("Strength of cos(2#phi) signal vs. Mass; Mass (GeV); 2<cos(2#phi)>");
mToyMassCos2phiMoments->Draw();
mToyMassCos2flatphiMoments->Draw("SAME");
//mToyMassCos2decayphiMoments->Draw("SAME");
mToyMassCos2cutflatphiMoments->Draw("SAME");
auto legend6 = new TLegend(0.65,0.1,0.95,0.4);
legend6->SetHeader("Legend","C"); // option "C" allows to center the header
legend6->AddEntry(mToyMassCos2flatphiMoments,"No effects (completely flat #phi hist)");
//legend6->AddEntry(mToyMassCos2decayphiMoments,"#pi->#mu+#nu decay only)");
legend6->AddEntry(mToyMassCos2cutflatphiMoments,"Detector acceptance");
legend6->AddEntry(mToyMassCos2phiMoments,"P_{T} smear, detector acceptance");
legend6->Draw();
gPad->Print( "plots/model/Mass/plot_mToyMassCos2phiMoments.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyMassCos2phiMoments.png" );

makeCanvas();
mCos4phivsMass->Draw("colz");
gPad->Print( "plots/model/Mass/plot_mToyCos4phivsMass.pdf" );
gPad->Print( "plots/model/Mass/plot_mToyCos4phivsMass.png" );

fo->Write();
}
