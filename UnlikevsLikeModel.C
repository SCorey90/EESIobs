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

TLorentzVector *hist_sample_PT( TH1F * PtHist, double etamin, double etamax, double phimin, double phimax, double mass) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = PtHist->GetRandom();
    double eta = rng.Uniform(etamin, etamax);
    double phi = rng.Uniform(phimin, phimax);
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
}

TLorentzVector *hist_sample_Eta( double ptmin, double ptmax, TH1F * EtaHist, double phimin, double phimax, double mass ) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = rng.Uniform(ptmin, ptmax);
    double eta = EtaHist->GetRandom();
    double phi = rng.Uniform(phimin, phimin);
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
}

TLorentzVector *hist_sample_Phi( double ptmin, double ptmax, double etamin, double etamax, TH1F * LVPhiHist, double mass) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = rng.Uniform(ptmin, ptmax);
    double eta = rng.Uniform(etamin, etamax);
    double phi = LVPhiHist->GetRandom();
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
}

TLorentzVector *hist_sample_all( TH1F * PtHist, TH1F * EtaHist, TH1F * LVPhiHist, TH1F * MassHist) {
    TLorentzVector *lv0 = new TLorentzVector;
    double pt = PtHist->GetRandom();
    double eta = EtaHist->GetRandom();
    double phi = LVPhiHist->GetRandom();
    double mass = MassHist->GetRandom();
    lv0->SetPtEtaPhiM( pt, eta, phi, mass);
    return lv0;
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

void UnlikevsLikeModel() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * unlikeHists = TFile::Open("EESIplots.root");
    TFile * likeHists = TFile::Open("MixedEventplots.root");
    
    unlikeHists->ls();
    auto * mUnlikePiPT = (TH1F*)unlikeHists->Get("mPiPT");
    auto * mUnlikePiEta = (TH1F*)unlikeHists->Get("mPiEta");
    auto * mUnlikePiAngle = (TH1F*)unlikeHists->Get("mPiAngle");

    likeHists->ls();
    auto * mlikePiPT = (TH1F*)likeHists->Get("mPiDataLikePT");
    auto * mlikePiEta = (TH1F*)likeHists->Get("mPiDataLikeEta");
    auto * mlikePiAngle = (TH1F*)likeHists->Get("mPiDataLikeAngle");

    TFile *fo = new TFile( "UnlikevsLikeModelplots.root", "RECREATE" );

    int n_events = 1000000;
    double m_pi = 0.139;
    double PI = 3.1415926535;

    auto * mflatPT = new TH1F("mflatPT", "Random Daughters (uniform)->Parent P_{T}; P_{T} (GeV/c); counts", 500, 0, 4);
    auto * mflatMass = new TH1F("mflatMass", "Random Daughters (uniform)->Parent Mass; Mass (GeV); counts", 500, 0, 5);
    auto * mflatPhi = new TH1F("mflatPhi", "Random Daughters (uniform)->Parent #phi; #phi (rad); counts", 500, -PI, PI);

    auto * mlikePThistPT = new TH1F("mlikePThistPT", "Random Daughters (like sign P_{T} hist)->Parent P_{T}; P_{T} (GeV/c); counts", 500, 0, 4);
    auto * mlikePThistMass = new TH1F("mlikePThistMass", "Random Daughters (like sign P_{T} hist)->Parent Mass; Mass (GeV); counts", 500, 0, 5);
    auto * mlikePThistPhi = new TH1F("mlikePThistPhi", "Random Daughters (like sign P_{T} hist)->Parent #phi; #phi (rad); counts", 500, -PI, PI);

    auto * munlikePThistPT = new TH1F("munlikePThistPT", "Random Daughters (unlike sign P_{T} hist)->Parent P_{T}; P_{T} (GeV/c); counts", 500, 0, 4);
    auto * munlikePThistMass = new TH1F("munlikePThistMass", "Random Daughters (unlike sign P_{T} hist)->Parent Mass; Mass (GeV); counts", 500, 0, 5);
    auto * munlikePThistPhi = new TH1F("munlikePThistPhi", "Random Daughters (unlike sign P_{T} hist)->Parent #phi; #phi (rad); counts", 500, -PI, PI);

    auto * munlikePThistCos2phivsPT = new TH2F("munlikePThistCos2phivsPT", "Random Daughters (unlike sign P_{T} hist)->Parent cos2#phi vs P_{T}; 2cos2#phi; P_{T}; counts", 500, -2, 2, 500, 0, 2);

    for (int i = 0; i < n_events; i++) {
        TLorentzVector flatLV1, flatLV2, flatLVSum, likePThistLV1, likePThistLV2, likePThistLVSum, unlikePThistLV1, unlikePThistLV2, unlikePThistLVSum;

        TLorentzVector * ptr_flatLV1 = random_LV(0.05, 2, -1.5, 1.5, -PI, PI, m_pi, m_pi);
        TLorentzVector * ptr_flatLV2 = random_LV(0.05, 2, -1.5, 1.5, -PI, PI, m_pi, m_pi);
        flatLV1.SetPtEtaPhiM( ptr_flatLV1->Pt(), ptr_flatLV1->Eta(), ptr_flatLV1->Phi(), m_pi );
        flatLV2.SetPtEtaPhiM( ptr_flatLV2->Pt(), ptr_flatLV2->Eta(), ptr_flatLV2->Phi(), m_pi );
        flatLVSum = flatLV1 + flatLV2;

        TLorentzVector *ptr_likePThistLV1 = hist_sample_PT( mlikePiPT, -1.5, 1.5, -PI, PI, m_pi);
        TLorentzVector *ptr_likePThistLV2 = hist_sample_PT( mlikePiPT, -1.5, 1.5, -PI, PI, m_pi);
        likePThistLV1.SetPtEtaPhiM( ptr_likePThistLV1->Pt(), ptr_likePThistLV1->Eta(), ptr_likePThistLV1->Phi(), m_pi );
        likePThistLV2.SetPtEtaPhiM( ptr_likePThistLV2->Pt(), ptr_likePThistLV2->Eta(), ptr_likePThistLV2->Phi(), m_pi );
        likePThistLVSum = likePThistLV1 + likePThistLV2;

        TLorentzVector *ptr_unlikePThistLV1 = hist_sample_PT( mUnlikePiPT, -1.5, 1.5, -PI, PI, m_pi);
        TLorentzVector *ptr_unlikePThistLV2 = hist_sample_PT( mUnlikePiPT, -1.5, 1.5, -PI, PI, m_pi);
        unlikePThistLV1.SetPtEtaPhiM( ptr_unlikePThistLV1->Pt(), ptr_unlikePThistLV1->Eta(), ptr_unlikePThistLV1->Phi(), m_pi );
        unlikePThistLV2.SetPtEtaPhiM( ptr_unlikePThistLV2->Pt(), ptr_unlikePThistLV2->Eta(), ptr_unlikePThistLV2->Phi(), m_pi );
        unlikePThistLVSum = unlikePThistLV1 + unlikePThistLV2;

        mflatPT->Fill( flatLVSum.Pt() );
        mflatMass->Fill( flatLVSum.M() );
        mflatPhi->Fill(calc_Phi(flatLV1, flatLV2));

        mlikePThistPT->Fill( likePThistLVSum.Pt() );
        mlikePThistMass->Fill( likePThistLVSum.M() );
        mlikePThistPhi->Fill(calc_Phi(likePThistLV1, likePThistLV2));
  
        munlikePThistPT->Fill( unlikePThistLVSum.Pt() );
        munlikePThistMass->Fill( unlikePThistLVSum.M() );
        munlikePThistPhi->Fill(calc_Phi(unlikePThistLV1, unlikePThistLV2));

        munlikePThistCos2phivsPT->Fill( 2*cos(2*calc_Phi(unlikePThistLV1, unlikePThistLV2)), unlikePThistLVSum.Pt());

    }

auto * munlikePThistCos2phiMomentsvsPT = munlikePThistCos2phivsPT->ProfileY("munlikePThistCos2phiMomentsvsPT", -1, 1);

fo->cd();

makeCanvas();
mflatPT->SetLineColor(kBlack);
mflatPT->Draw();
gPad->Print("plots/UnlikevsLikeModel/PT/plot_mflatPT.png");
gPad->Print("plots/UnlikevsLikeModel/PT/plot_mflatPT.pdf");

makeCanvas();
mflatMass->SetLineColor(kBlack);
mflatMass->Draw();
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_mflatMass.png");
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_mflatMass.pdf");

makeCanvas();
mflatPhi->SetLineColor(kBlack);
mflatPhi->Draw();
gPad->Print("plots/UnlikevsLikeModel/plot_mflatPhi.png");
gPad->Print("plots/UnlikevsLikeModel/plot_mflatPhi.pdf");

makeCanvas();
mlikePThistPT->SetLineColor(kBlack);
mlikePThistPT->Draw();
gPad->Print("plots/UnlikevsLikeModel/PT/plot_mlikePThistPT.png");
gPad->Print("plots/UnlikevsLikeModel/PT/plot_mlikePThistPT.pdf");

makeCanvas();
mlikePThistMass->SetLineColor(kBlack);
mlikePThistMass->Draw();
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_mlikePThistMass.png");
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_mlikePThistMass.pdf");

makeCanvas();
mlikePThistPhi->SetLineColor(kBlack);
mlikePThistPhi->Draw();
gPad->Print("plots/UnlikevsLikeModel/plot_mlikePThistPhi.png");
gPad->Print("plots/UnlikevsLikeModel/plot_mlikePThistPhi.pdf");

makeCanvas();
munlikePThistPT->SetLineColor(kBlack);
munlikePThistPT->Draw();
gPad->Print("plots/UnlikevsLikeModel/PT/plot_munlikePThistPT.png");
gPad->Print("plots/UnlikevsLikeModel/PT/plot_munlikePThistPT.pdf");

makeCanvas();
munlikePThistMass->SetLineColor(kBlack);
munlikePThistMass->Draw();
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_munlikePThistMass.png");
gPad->Print("plots/UnlikevsLikeModel/Mass/plot_munlikePThistMass.pdf");

makeCanvas();
munlikePThistPhi->SetLineColor(kBlack);
munlikePThistPhi->Draw();
gPad->Print("plots/UnlikevsLikeModel/plot_munlikePThistPhi.png");
gPad->Print("plots/UnlikevsLikeModel/plot_munlikePThistPhi.pdf");

makeCanvas();
munlikePThistCos2phiMomentsvsPT->SetLineColor(kBlack);
munlikePThistCos2phiMomentsvsPT->Draw();
gPad->Print("plots/UnlikevsLikeModel/PT/plot_munlikePThistCos2phiMomentsvsPT.png");
gPad->Print("plots/UnlikevsLikeModel/PT/plot_munlikePThistCos2phiMomentsvsPT.pdf");

fo->Write();

}
