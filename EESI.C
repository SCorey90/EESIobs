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


double* histMoments( TH2F* hist , int n) {
    int nbinsy = hist->GetNbinsY();
    int nbinsx = hist->GetNbinsX();
    double phivals[nbinsx];
    double weights[nbinsx];
    double *cosnphi_moments = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
	for (int j = 1; j < nbinsx+1; j++) {
            phivals[j-1] = onedhist->GetXaxis()->GetBinCenter(j);
            weights[j-1] = onedhist->GetBinContent(j);
        }
        cosnphi_moments[i-1] = 0;
        for (int k = 0; k < nbinsx; k++){
            if ( onedhist->GetEntries() > 0) {
                cosnphi_moments[i-1] += weights[k] * n * cos(n*phivals[k]) / (onedhist->GetEntries());
            } else { cosnphi_moments[i-1] += 0; }
        }
    }
    
    return cosnphi_moments;
} 

double* moment_error( TH2F* hist, int n ) {
    int nbinsx = hist->GetNbinsX();
    int nbinsy = hist->GetNbinsY();
    double *moment_errors = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        auto* onedhist = hist->ProjectionX("1dhist", i, i);
        auto* momenthist = new TH1F("momenthist", "ncosnphi moments", nbinsx, -1*n, 1*n);
        for (int j = 1; j < nbinsx+1; j++) {
            for (int k = 0; k < onedhist->GetBinContent(j); k++) {
                momenthist->Fill( n * cos( n * (onedhist->GetXaxis()->GetBinCenter(j))));
            }
        }
        moment_errors[i-1] = momenthist->GetMeanError();
    }
    return moment_errors;
}       

double* histYbins( TH2F* hist) {
    int nbinsy = hist->GetNbinsY();
    double *ybinvals = new double[nbinsy];
    for (int i = 1; i < nbinsy+1; i++) {
        ybinvals[i-1] = hist->GetYaxis()->GetBinCenter(i);
    }
    return ybinvals;
}

void EESI() {
    TH1F("h1", "ntuple", 100, -4, 4);

    TFile * fo = new TFile( "EESIplots.root", "RECREATE" );

    auto * mMass = new TH1F("mMass", "Parent (#rho^{0})  Mass; Mass (GeV); Counts", 500, 0, 2);
    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1);
    auto * mZDCTotal = new TH1F("mZDCtotal", "ZDC East+West readout; ZDC readout; counts", 500, 0, 1200);
    auto * mPairPhi = new TH1F("mPairPhi", "#pi_{#pm} #phi distribution;#phi (rad);# events", 500, -3.13, 3.13);
    auto * mPxvsPy = new TH2F("mPxvsPy", "#rho^{0} 2D momentum dist; P_{x} (GeV); P_{y} (GeV); Counts", 200, -0.1, 0.1, 200, -0.1, 0.1);
    auto * mPhivsMass = new TH2F("mPhivsMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 50, 0.3, 1.35);
    auto * mPhivsPT = new TH2F("mPhivsPT", "#pi^{#pm} #phi distribution vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV); Counts", 100, -3.14, 3.14, 100, 0, 0.25);
    auto * mPhivsRapidity = new TH2F("mPhivsRapidity", "#pi^{#pm} #phi distribution vs. Rapidity; #phi (rad); Rapidity (GeV); Counts", 100, -3.14, 3.14, 100, -2, 2);
    auto * mPhivsLowMass = new TH2F("mPhivsLowMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 100, 0.3, 0.55);
    auto * mPhivsZDC = new TH2F("mPhivsZDC", "#pi^{#pm} #phi distribution vs. ZDC readout; #phi (rad); (ZDC East^{2} + ZDC West^{2})^{1/2}; Counts", 100, -3.14, 3.14, 50, 0, 1200);
    auto * mPhivsEastZDC = new TH2F("mPhivsEastZDC", "#pi^{#pm} #phi distribution vs. ZDC East readout; #phi (rad); ZDC East counts; Counts", 100, -3.14, 3.14, 50, 0, 700);
    auto * mPhivsWestZDC = new TH2F("mPhivsWestZDC", "#pi^{#pm} #phi distribution vs. ZDC West readout; #phi (rad); ZDC West counts; Counts", 100, -3.14, 3.14, 50, 0, 700);

    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -1, 1, 100, 0, 0.25);
    auto * mCos4phivsPT = new TH2F("mCos4phivsPT", "cos4#phi distribution vs P_{T}", 100, -1, 1, 100, 0, 0.25);

    auto * mCos2phivsMass = new TH2F("mCos2phivsMass", "cos2#phi distribution vs Mass", 100, -1, 1, 100, 0.3, 1.35);
    auto * mCos4phivsMass = new TH2F("mCos4phivsMass", "cos4#phi distribution vs Mass", 100, -1, 1, 100, 0.3, 1.35);

    auto * mPhiFit = new TF1("mPhiFit", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x)", -3.14, 3.14);
    auto * mLinearFit = new TF1("mLinearFit", "[0] + [1]*x", -5, 5);
    auto * mQuadraticFit = new TF1("mQuadraticFit", "[0] + [1]*x + [2]*x*x", -5, 5);
    //Open pairDST
    TFile *myFile = TFile::Open("/Users/samcorey/code/data/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        double EastZDC = pair->mZDCEast;
        double WestZDC = pair->mZDCWest;
        double TotalZDC = sqrt((EastZDC*EastZDC + WestZDC*WestZDC));
        Float_t mRapidity = pair->mRapidity;
        
        Float_t mMassVal = pair->mMass; 
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
        lv = lv1+lv2;
        double absPperp = lv.Pt();
        double PairPhi = calc_Phi( lv1, lv2);

        if ( chipipi <10 && dca1 <1 && dca2 <1 ){
		mPperp->Fill( absPperp );
		mMass->Fill( lv.M() );
                mZDCTotal->Fill( TotalZDC );
            if ( lv.M() > 0.65 && lv.M() <0.9){
                if ( absPperp < 0.06){

                    mPairPhi->Fill ( PairPhi );
                    mPhivsRapidity->Fill ( PairPhi, mRapidity);
                    mPhivsZDC->Fill ( PairPhi, TotalZDC );
                    mPhivsEastZDC->Fill ( PairPhi, EastZDC );
                    mPhivsWestZDC->Fill ( PairPhi, WestZDC );
                }
                mPxvsPy->Fill( absPperp*cos(PairPhi), absPperp*sin(PairPhi));
                mPhivsPT->Fill ( PairPhi, absPperp);
                mCos2phivsPT->Fill( 2*cos(2 *(PairPhi)), absPperp );
                mCos4phivsPT->Fill( 4*cos(4 *(PairPhi)), absPperp );
            }
            if ( lv.M() > 0.2 && lv.M() <1.5 && absPperp < 0.06) { 
                mPhivsMass->Fill ( PairPhi, lv.M());
                mCos2phivsMass->Fill( 2*cos(2 *(PairPhi)), lv.M() );
                mCos4phivsMass->Fill( 4*cos(4 *(PairPhi)), lv.M() ); 
            }
            if ( lv.M() > 0.2 && lv.M() <0.55 && absPperp < 0.06) { mPhivsLowMass->Fill ( PairPhi, lv.M()); }
        }
    }
//mPairPhi->Fit(mPhiFit);

auto *mPTmomentsplot = new TGraphErrors(mPhivsPT->GetNbinsY(), histYbins(mPhivsPT), histMoments(mPhivsPT, 2), 0, moment_error(mPhivsPT, 2));
auto *mMassmomentsplot = new TGraphErrors(mPhivsMass->GetNbinsY(), histYbins(mPhivsMass), histMoments(mPhivsMass, 2), 0, moment_error(mPhivsMass, 2));
auto *mLowMassmomentsplot = new TGraphErrors(mPhivsLowMass->GetNbinsY(), histYbins(mPhivsLowMass), histMoments(mPhivsLowMass, 2), 0, moment_error(mPhivsLowMass, 2));
auto *mRapiditymomentsplot = new TGraphErrors(mPhivsRapidity->GetNbinsY(), histYbins(mPhivsRapidity), histMoments(mPhivsRapidity, 2), 0, moment_error(mPhivsRapidity, 2));
auto *mZDCmomentsplot = new TGraphErrors(mPhivsZDC->GetNbinsY(), histYbins(mPhivsZDC), histMoments(mPhivsZDC, 2), 0, moment_error(mPhivsZDC, 2));
auto *mPTcos4phimoments = new TGraphErrors(mPhivsPT->GetNbinsY(), histYbins(mPhivsPT), histMoments(mPhivsPT, 4), 0, moment_error(mPhivsPT, 4));
auto *mMasscos4phimoments = new TGraphErrors(mPhivsMass->GetNbinsY(), histYbins(mPhivsMass), histMoments(mPhivsMass, 4), 0, moment_error(mPhivsMass, 4));
auto *mRapiditycos4phimoments = new TGraphErrors(mPhivsRapidity->GetNbinsY(), histYbins(mPhivsRapidity), histMoments(mPhivsRapidity, 4), 0, moment_error(mPhivsRapidity, 4));

auto *mv2PTcos2phimoments = mCos2phivsPT->ProfileY("mv2PTcos2phimoments", 1, -1);
auto *mv2PTcos4phimoments = mCos4phivsPT->ProfileY("mv2PTcos4phimoments", 1, -1);

TGraph *mRapiditymomentsQuad = (TGraph*)mRapiditymomentsplot->Clone("mRapiditymomentsQuad");
mRapiditymomentsplot->Fit(mLinearFit,"", "", -0.5, 0.5);
mRapiditymomentsQuad->Fit(mQuadraticFit,"", "", -1.6, 1.6);

fo -> cd();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plots/plot_mMass.pdf" );

makeCanvas();
mPperp->SetLineColor(kBlack);
mPperp->Draw();
gPad->Print( "plots/plot_mPperp.pdf" );

makeCanvas();
mZDCTotal->SetLineColor(kBlack);
mZDCTotal->Draw();
gPad->Print( "plots/plot_mZDCTotal.pdf" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
//gStyle->SetOptFit();
mPairPhi->Draw();
gPad->Print( "plots/plot_mPairPhi.pdf" );
gPad->Print( "plots/plot_mPairPhi.png" );

makeCanvas();
mPhiFit->SetLineColor(kBlack);
mPhiFit->Draw();
gPad->Print( "plots/plot_mPhiFit.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPxvsPy->Draw("colz");
gPad->Print( "plots/plot_mPxvsPy.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsMass->Draw("colz");
gPad->Print( "plots/plot_mPhivsMass.pdf" );
gPad->Print( "plots/plot_mPhivsMass.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsLowMass->Draw("colz");
gPad->Print( "plots/plot_mPhivsLowMass.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsPT->Draw("colz");
gPad->Print( "plots/plot_mPhivsPT.pdf" );
gPad->Print( "plots/plot_mPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsRapidity->Draw("colz");
gPad->Print( "plots/plot_mPhivsRapidity.pdf" );
gPad->Print( "plots/plot_mPhivsRapidity.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsEastZDC->Draw("colz");
gPad->Print( "plots/plot_mPhivsEastZDC.pdf" );
gPad->Print( "plots/plot_mPhivsEastZDC.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsWestZDC->Draw("colz");
gPad->Print( "plots/plot_mPhivsWestZDC.pdf" );
gPad->Print( "plots/plot_mPhivsWestZDC.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsZDC->Draw("colz");
gPad->Print( "plots/plot_mPhivsZDC.pdf" );
gPad->Print( "plots/plot_mPhivsZDC.png" );

makeCanvas();
mPTmomentsplot->SetTitle("cos(2#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mPTmomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mPTmomentsplot.pdf" );
gPad->Print( "plots/plot_mPTmomentsplot.png" );

makeCanvas();
mMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
mMassmomentsplot->SetMinimum(0);
mMassmomentsplot->SetMaximum(1);
mMassmomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mMassmomentsplot.pdf" );
gPad->Print( "plots/plot_mMassmomentsplot.png" );

makeCanvas();
mRapiditymomentsplot->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mRapiditymomentsplot.pdf" );
gPad->Print( "plots/plot_mRapiditymomentsplot.png" );
mRapiditymomentsplot->SetMinimum(0);
mRapiditymomentsplot->SetMaximum(0.5);
gStyle->SetOptFit(1111);
mRapiditymomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mZoomedRapiditymomentsplot.pdf" );
gPad->Print( "plots/plot_mZoomedRapiditymomentsplot.png" );

makeCanvas();
mRapiditymomentsQuad->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsQuad->Draw("AC*");
gPad->Print( "plots/plot_mRapiditymomentsQuad.pdf" );
gPad->Print( "plots/plot_mRapiditymomentsQuad.png" );
mRapiditymomentsQuad->SetMinimum(0);
mRapiditymomentsQuad->SetMaximum(0.5);
gStyle->SetOptFit(1111);
mRapiditymomentsQuad->Draw("AC*");
gPad->Print( "plots/plot_mZoomedRapiditymomentsQuad.pdf" );
gPad->Print( "plots/plot_mZoomedRapiditymomentsQuad.png" );

makeCanvas();
mZDCmomentsplot->SetTitle("cos(2#phi) moments vs. ZDC; (ZDC East^{2} + ZDC West^{2})^{1/2}; 2<cos(2#phi)>");
mZDCmomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mZDCmomentsplot.pdf" );
gPad->Print( "plots/plot_mZDCmomentsplot.png" );

makeCanvas();
mPTcos4phimoments->SetTitle("cos(4#phi) moments vs. P_{T}; P_{T} (GeV); 4<cos(4#phi)>");
mPTcos4phimoments->Draw("AC*");
gPad->Print( "plots/plot_mPTcos4phimoments.pdf" );

makeCanvas();
mMasscos4phimoments->SetTitle("cos(4#phi) moments vs. Mass; Mass (GeV); 4<cos(4#phi)>");
mMasscos4phimoments->Draw("AC*");
gPad->Print( "plots/plot_mMasscos4phimoments.pdf" );

makeCanvas();
mRapiditycos4phimoments->SetTitle("cos(4#phi) moments vs. Rapidity; Rapidity; 4<cos(4#phi)>");
mRapiditycos4phimoments->Draw("AC*");
gPad->Print( "plots/plot_mRapiditycos4phimoments.pdf" );

makeCanvas();
mLowMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
//mLowMassmomentsplot->SetMinimum(0);
//mLowMassmomentsplot->SetMaximum(1);
mLowMassmomentsplot->Draw("AC*");
gPad->Print( "plots/plot_mLowMassmomentsplot.pdf" );
gPad->Print( "plots/plot_mLowMassmomentsplot.png" );

makeCanvas();
mv2PTcos2phimoments->SetLineColor(kBlack);
mv2PTcos2phimoments->Draw();
gPad->Print( "plots/plot_mv2PTcos2phimoments.pdf" );
gPad->Print( "plots/plot_mv2PTcos2phimoments.png" );

makeCanvas();
mv2PTcos4phimoments->SetLineColor(kBlack);
mv2PTcos4phimoments->Draw();
gPad->Print( "plots/plot_mv2PTcos4phimoments.pdf" );
gPad->Print( "plots/plot_mv2PTcos4phimoments.png" );

makeCanvas();
mCos2phivsPT->Draw("colz");
gPad->Print( "plots/plot_mCos2phivsPT.pdf" );

makeCanvas();
mCos4phivsPT->Draw("colz");
gPad->Print( "plots/plot_mCos4phivsPT.pdf" );

makeCanvas();
mCos2phivsMass->Draw("colz");
gPad->Print( "plots/plot_mCos2phivsMass.pdf" );

makeCanvas();
mCos4phivsMass->Draw("colz");
gPad->Print( "plots/plot_mCos4phivsMass.pdf" );

fo->Write();
}
