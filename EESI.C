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

double lab_calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
    //lv1.Boost(-lvPlus.BoostVector());
    //lv2.Boost(-lvPlus.BoostVector());
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
                cosnphi_moments[i-1] += weights[k] * 2 * cos(n*phivals[k]) / (onedhist->GetEntries());
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
        auto* momenthist = new TH1F("momenthist", "2cosnphi moments", nbinsx, -1*n, 1*n);
        for (int j = 1; j < nbinsx+1; j++) {
            for (int k = 0; k < onedhist->GetBinContent(j); k++) {
                momenthist->Fill( 2 * cos( n * (onedhist->GetXaxis()->GetBinCenter(j))));
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
    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mEta = new TH1F("mEta", "Parent (#rho^{0}) Eta; Eta; counts", 500, -6, 6);
    auto * mLVPhi = new TH1F("mLVPhi", "Parent (#rho^{0}) LV.Phi(); #phi (rad); counts", 500, -3.15, 3.15);
    auto * mZDCEast = new TH1F("mZDCEast", "ZDC East; ZDC East; counts", 200, 0, 1300);
    auto * mZDCWest = new TH1F("mZDCWest", "ZDC West; ZDC West; counts", 200, 0, 1300);

    auto * mZDCTotal = new TH1F("mZDCtotal", "ZDC (East^{2} + West^{2})^{1/2}; ZDC (East^{2} + West^{2})^{1/2}; counts", 1000, 0, 1900);
    auto * mPairPhi = new TH1F("mPairPhi", "#pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);
    auto * mLabPairPhi = new TH1F("mLabPairPhi", "#pi_{#pm} #phi distribution;#phi (rad);# events", 100, -3.13, 3.13);
    auto * mPxvsPy = new TH2F("mPxvsPy", "#rho^{0} 2D momentum dist; P_{x} (GeV/c); P_{y} (GeV/c); Counts", 200, -0.1, 0.1, 200, -0.1, 0.1);
    auto * mPhivsMass = new TH2F("mPhivsMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 50, 0.3, 1.35);
    auto * mPhivsPT = new TH2F("mPhivsPT", "#pi^{#pm} #phi distribution vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV/c); Counts", 100, -3.14, 3.14, 100, 0, 0.25);
    auto * mPhivsRapidity = new TH2F("mPhivsRapidity", "#pi^{#pm} #phi distribution vs. Rapidity; #phi (rad); Rapidity (GeV); Counts", 100, -3.14, 3.14, 100, -2, 2);
    auto * mPhivsLowMass = new TH2F("mPhivsLowMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 100, 0.3, 0.55);
    
    auto * mPhivsZDC = new TH2F("mPhivsZDC", "#pi^{#pm} #phi distribution vs. ZDC readout; #phi (rad); (ZDC East^{2} + ZDC West^{2})^{1/2}; Counts", 100, -3.14, 3.14, 500, 0, 1900);
    auto * mPhivsEastZDC = new TH2F("mPhivsEastZDC", "#pi^{#pm} #phi distribution vs. ZDC East readout; #phi (rad); ZDC East counts; Counts", 100, -3.14, 3.14, 50, 0, 700);
    auto * mPhivsWestZDC = new TH2F("mPhivsWestZDC", "#pi^{#pm} #phi distribution vs. ZDC West readout; #phi (rad); ZDC West counts; Counts", 100, -3.14, 3.14, 50, 0, 700);
    auto * mPhi1n1n = new TH1F("mPhi1n1n", "#phi distribution of 1n1n ZDC peak; #phi (rad); counts", 20, -3.14, 3.14);
    auto * mPhi2n1n = new TH1F("mPhi2n1n", "#phi distribution of 2n1n ZDC peaks; #phi (rad); counts", 20, -3.14, 3.14);
    auto * mPhi2n2n = new TH1F("mPhi2n2n", "#phi distribution of 2n2n ZDC peak; #phi (rad); counts", 20, -3.14, 3.14);

    auto * mCos2phivsEastZDC = new TH2F("mCos2phivsEastZDC", "cos2#phi distribution vs East ZDC; 2cos2#phi; East ZDC readout; counts", 100, -2, 2, 50, 0, 1300);
    auto * mCos2phivsWestZDC = new TH2F("mCos2phivsWestZDC", "cos2#phi distribution vs West ZDC; 2cos2#phi; West ZDC readout; counts", 100, -2, 2, 50, 0, 1300);
    auto * mCos2phivsTotalZDC = new TH2F("mCos2phivsTotalZDC", "cos2#phi distribution vs ZDC (East^{2} + West^{2})^{1/2}; 2cos2#phi; ZDC (East^{2} + West^{2})^{1/2}; counts", 100, -2, 2, 50, 0, 1700);

    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.5);
    auto * mLabCos2phivsPT = new TH2F("mLabCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.5);
    auto * mCos4phivsPT = new TH2F("mCos4phivsPT", "cos4#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.25);

    auto * mLowZDCcos2phivsPT = new TH2F("mLowZDCcos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.25);
    auto * mHighZDCcos2phivsPT = new TH2F("mHighZDCcos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.25);

    auto * mCos2phivsMass = new TH2F("mCos2phivsMass", "cos2#phi distribution vs Mass", 100, -2, 2, 100, 0.3, 1.35);
    auto * mCos4phivsMass = new TH2F("mCos4phivsMass", "cos4#phi distribution vs Mass", 100, -2, 2, 100, 0.3, 1.35);

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
                mEta->Fill( lv.Eta() );
                mLVPhi->Fill( lv.Phi() );
                mZDCEast->Fill( EastZDC );
                mZDCWest->Fill( WestZDC );
                mZDCTotal->Fill( TotalZDC );
            if ( lv.M() > 0.65 && lv.M() <0.9){
                if ( absPperp < 0.06){

                    mPairPhi->Fill ( PairPhi );
                    mLabPairPhi->Fill ( lab_calc_Phi(lv1, lv2) );
                    mPhivsRapidity->Fill ( PairPhi, mRapidity);
                    mPhivsZDC->Fill ( PairPhi, TotalZDC );
                    mPhivsEastZDC->Fill ( PairPhi, EastZDC );
                    mPhivsWestZDC->Fill ( PairPhi, WestZDC );

                    mCos2phivsEastZDC->Fill ( 2*cos(2*PairPhi), EastZDC );
                    mCos2phivsWestZDC->Fill ( 2*cos(2*PairPhi), WestZDC );
                    mCos2phivsTotalZDC->Fill ( 2*cos(2*PairPhi), TotalZDC );
                }
                if ( WestZDC > 30 && WestZDC < 70 && EastZDC > 30 && EastZDC < 70 ){mPhi1n1n->Fill(PairPhi);}
                if ( WestZDC > 30 && WestZDC < 70 && EastZDC > 100 && EastZDC < 135 ){mPhi2n1n->Fill(PairPhi);}
                if ( WestZDC > 100 && WestZDC < 135 && EastZDC > 30 && EastZDC < 70 ){mPhi2n1n->Fill(PairPhi);}
                if ( WestZDC > 100 && WestZDC < 135 && EastZDC > 100 && EastZDC < 130 ){mPhi2n2n->Fill(PairPhi);}
                if ( TotalZDC < 300 ) { mLowZDCcos2phivsPT->Fill( 2*cos(2*PairPhi), absPperp); }
                if ( TotalZDC > 600 ) { mHighZDCcos2phivsPT->Fill( 2*cos(2*PairPhi), absPperp); }
                mPxvsPy->Fill( absPperp*cos(PairPhi), absPperp*sin(PairPhi));
                mPhivsPT->Fill ( PairPhi, absPperp);
                mCos2phivsPT->Fill( 2*cos(2 *(PairPhi)), absPperp );
                mLabCos2phivsPT->Fill( 2*cos(2 *(lab_calc_Phi(lv1, lv2))), absPperp );
                mCos4phivsPT->Fill( 4*cos(4 *(PairPhi)), absPperp );
            }
            if ( lv.M() > 0.2 && lv.M() <1.5 && absPperp < 0.06) { 
                mPhivsMass->Fill ( PairPhi, lv.M());
                mCos2phivsMass->Fill( 2*cos(2 *(PairPhi)), lv.M() );
                mCos4phivsMass->Fill( 2*cos(4 *(PairPhi)), lv.M() ); 
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

//mPhi1n1n->Scale(2*3.1415/(mPhi1n1n->Integral("width")));
//mPhi2n1n->Scale(2*3.1415/(mPhi2n1n->Integral("width")));
//mPhi2n2n->Scale(2*3.1415/(mPhi2n2n->Integral("width")));

auto *mv2PTcos2phimoments = mCos2phivsPT->ProfileY("mv2PTcos2phimoments", 1, -1);
auto *mlabPTcos2phimoments = mLabCos2phivsPT->ProfileY("mlabPTcos2phimoments", 1, -1);
auto *mv2PTcos4phimoments = mCos4phivsPT->ProfileY("mv2PTcos4phimoments", 1, -1);

auto *mEastZDCcos2phimoments = mCos2phivsEastZDC->ProfileY("mEastZDCcos2phimoments", 1, -1);
auto *mWestZDCcos2phimoments = mCos2phivsWestZDC->ProfileY("mWestZDCcos2phimoments", 1, -1);
auto *mTotalZDCcos2phimoments = mCos2phivsTotalZDC->ProfileY("mTotalZDCcos2phimoments", 1, -1);
auto *mLowZDCPTcos2phimoments = mLowZDCcos2phivsPT->ProfileY("mLowZDCPTcos2phimoments", 1, -1);
auto *mHighZDCPTcos2phimoments = mHighZDCcos2phivsPT->ProfileY("mHighZDCPTcos2phimoments", 1, -1);

TGraph *mRapiditymomentsQuad = (TGraph*)mRapiditymomentsplot->Clone("mRapiditymomentsQuad");
mRapiditymomentsplot->Fit(mLinearFit,"", "", -0.5, 0.5);
mRapiditymomentsQuad->Fit(mQuadraticFit,"", "", -1.6, 1.6);

mTotalZDCcos2phimoments->Fit("pol0", "","", 600, 1600);

fo -> cd();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plots/data/Mass/plot_mMass.pdf" );

makeCanvas();
mPperp->SetLineColor(kBlack);
mPperp->Draw();
gPad->Print( "plots/data/PT/plot_mPperp.pdf" );

makeCanvas();
mZDCEast->SetLineColor(kBlack);
mZDCEast->Draw();
gPad->Print( "plots/data/ZDC/plot_mZDCEast.pdf" );
gPad->Print( "plots/data/ZDC/plot_mZDCEast.png" );

makeCanvas();
mZDCWest->SetLineColor(kBlack);
mZDCWest->Draw();
gPad->Print( "plots/data/ZDC/plot_mZDCWest.pdf" );
gPad->Print( "plots/data/ZDC/plot_mZDCWest.png" );

makeCanvas();
mZDCTotal->SetLineColor(kBlack);
mZDCTotal->Draw();
gPad->Print( "plots/data/ZDC/plot_mZDCTotal.pdf" );
gPad->Print( "plots/data/ZDC/plot_mZDCTotal.png" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
//gStyle->SetOptFit();
mPairPhi->Draw();
gPad->Print( "plots/data/plot_mPairPhi.pdf" );
gPad->Print( "plots/data/plot_mPairPhi.png" );

makeCanvas();
mPhiFit->SetLineColor(kBlack);
mPhiFit->Draw();
gPad->Print( "plots/data/plot_mPhiFit.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPxvsPy->Draw("colz");
gPad->Print( "plots/data/PT/plot_mPxvsPy.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsMass->Draw("colz");
gPad->Print( "plots/data/Mass/plot_mPhivsMass.pdf" );
gPad->Print( "plots/data/Mass/plot_mPhivsMass.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsLowMass->Draw("colz");
gPad->Print( "plots/data/Mass/plot_mPhivsLowMass.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsPT->Draw("colz");
gPad->Print( "plots/data/PT/plot_mPhivsPT.pdf" );
gPad->Print( "plots/data/PT/plot_mPhivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsRapidity->Draw("colz");
gPad->Print( "plots/data/Rapidity/plot_mPhivsRapidity.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mPhivsRapidity.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsEastZDC->Draw("colz");
gPad->Print( "plots/data/ZDC/plot_mPhivsEastZDC.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhivsEastZDC.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsWestZDC->Draw("colz");
gPad->Print( "plots/data/ZDC/plot_mPhivsWestZDC.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhivsWestZDC.png" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsZDC->Draw("colz");
gPad->Print( "plots/data/ZDC/plot_mPhivsZDC.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhivsZDC.png" );

makeCanvas();
mPhi1n1n->SetLineColor(kBlack);
mPhi1n1n->Draw("e1 P*");
gPad->Print( "plots/data/ZDC/plot_mPhi1n1n.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhi1n1n.png" );

makeCanvas();
mPhi2n1n->SetLineColor(kBlack);
mPhi2n1n->Draw("e1 P*");
gPad->Print( "plots/data/ZDC/plot_mPhi2n1n.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhi2n1n.png" );

makeCanvas();
mPhi2n2n->SetLineColor(kBlack);
mPhi2n2n->Draw("e1 P*");
gPad->Print( "plots/data/ZDC/plot_mPhi2n2n.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhi2n2n.png" );

//makeCanvas();
//mPhi1n1n->SetLineColor(kBlack);
//mPhi1n1n->SetTitle("#phi at various ZDC peaks; #phi (rad); counts");
//mPhi2n1n->SetLineColor(kOrange);
//mPhi2n2n->SetLineColor(kBlue);
//mPhi1n1n->Draw();
//mPhi2n1n->Draw("SAME");
//mPhi2n2n->Draw("SAME");
//auto legend = new TLegend(0.65,0.1,0.95,0.4);
//legend->AddEntry(mPhi1n1n,"1n1n peaks");
//legend->AddEntry(mPhi2n1n,"2n1n peaks");
//legend->AddEntry(mPhi2n2n,"2n2n peak");
//legend->Draw();
//gPad->Print( "plots/data/ZDC/plot_mPhiNeutronPeaks.pdf" );
//gPad->Print( "plots/data/ZDC/plot_mPhiNeutronPeaks.png" );

makeCanvas();
mPTmomentsplot->SetTitle("cos(2#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mPTmomentsplot->Draw("AC*");
gPad->Print( "plots/data/PT/plot_mPTmomentsplot.pdf" );
gPad->Print( "plots/data/PT/plot_mPTmomentsplot.png" );

makeCanvas();
mMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
mMassmomentsplot->SetMinimum(0);
mMassmomentsplot->SetMaximum(1);
mMassmomentsplot->Draw("AC*");
gPad->Print( "plots/data/Mass/plot_mMassmomentsplot.pdf" );
gPad->Print( "plots/data/Mass/plot_mMassmomentsplot.png" );

makeCanvas();
mRapiditymomentsplot->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsplot->Draw("AC*");
gPad->Print( "plots/data/Rapidity/plot_mRapiditymomentsplot.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mRapiditymomentsplot.png" );
mRapiditymomentsplot->SetMinimum(0);
mRapiditymomentsplot->SetMaximum(0.5);
gStyle->SetOptFit(1111);
mRapiditymomentsplot->Draw("AC*");
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsplot.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsplot.png" );

makeCanvas();
mRapiditymomentsQuad->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsQuad->Draw("AC*");
gPad->Print( "plots/plot_mRapiditymomentsQuad.pdf" );
gPad->Print( "plots/plot_mRapiditymomentsQuad.png" );
mRapiditymomentsQuad->SetMinimum(0);
mRapiditymomentsQuad->SetMaximum(0.5);
gStyle->SetOptFit(1111);
mRapiditymomentsQuad->Draw("AC*");
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsQuad.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsQuad.png" );

makeCanvas();
mZDCmomentsplot->SetTitle("cos(2#phi) moments vs. ZDC; (ZDC East^{2} + ZDC West^{2})^{1/2}; 2<cos(2#phi)>");
mZDCmomentsplot->Draw("AC*");
gPad->Print( "plots/data/ZDC/plot_mZDCmomentsplot.pdf" );
gPad->Print( "plots/data/ZDC/plot_mZDCmomentsplot.png" );

makeCanvas();
mPTcos4phimoments->SetTitle("cos(4#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(4#phi)>");
mPTcos4phimoments->Draw("AC*");
gPad->Print( "plots/data/PT/plot_mPTcos4phimoments.pdf" );

makeCanvas();
mMasscos4phimoments->SetTitle("cos(4#phi) moments vs. Mass; Mass (GeV); 2<cos(4#phi)>");
mMasscos4phimoments->Draw("AC*");
gPad->Print( "plots/data/Mass/plot_mMasscos4phimoments.pdf" );

makeCanvas();
mRapiditycos4phimoments->SetTitle("cos(4#phi) moments vs. Rapidity; Rapidity; 2<cos(4#phi)>");
mRapiditycos4phimoments->Draw("AC*");
gPad->Print( "plots/data/Rapidity/plot_mRapiditycos4phimoments.pdf" );

makeCanvas();
mLowMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
//mLowMassmomentsplot->SetMinimum(0);
//mLowMassmomentsplot->SetMaximum(1);
mLowMassmomentsplot->Draw("AC*");
gPad->Print( "plots/data/Mass/plot_mLowMassmomentsplot.pdf" );
gPad->Print( "plots/data/Mass/plot_mLowMassmomentsplot.png" );

makeCanvas();
mv2PTcos2phimoments->SetLineColor(kBlack);
mv2PTcos2phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mv2PTcos2phimoments.pdf" );
gPad->Print( "plots/data/PT/plot_mv2PTcos2phimoments.png" );

makeCanvas();
mv2PTcos4phimoments->SetLineColor(kBlack);
mv2PTcos4phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mv2PTcos4phimoments.pdf" );
gPad->Print( "plots/data/PT/plot_mv2PTcos4phimoments.png" );

makeCanvas();
mCos2phivsPT->Draw("colz");
gPad->Print( "plots/data/PT/plot_mCos2phivsPT.pdf" );

makeCanvas();
mCos4phivsPT->Draw("colz");
gPad->Print( "plots/data/PT/plot_mCos4phivsPT.pdf" );

makeCanvas();
mCos2phivsMass->Draw("colz");
gPad->Print( "plots/data/Mass/plot_mCos2phivsMass.pdf" );

makeCanvas();
mCos4phivsMass->Draw("colz");
gPad->Print( "plots/data/Mass/plot_mCos4phivsMass.pdf" );

makeCanvas();
mEastZDCcos2phimoments->SetLineColor(kBlack);
mWestZDCcos2phimoments->SetLineColor(kMagenta);
mEastZDCcos2phimoments->SetTitle("cos(2#phi) moments vs. ZDC readout; East readout; 2<cos2#phi>");
mEastZDCcos2phimoments->Draw();
mWestZDCcos2phimoments->Draw("SAME");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->AddEntry(mEastZDCcos2phimoments,"East ZDC");
legend2->AddEntry(mWestZDCcos2phimoments,"West ZDC");
legend2->Draw();
gPad->Print( "plots/data/ZDC/plot_mEWZDCcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mEWZDCcos2phimoments.png" );

makeCanvas();
mTotalZDCcos2phimoments->SetLineColor(kBlack);
mTotalZDCcos2phimoments->SetTitle("cos(2#phi) moments vs. ZDC (East^{2} + West^{2})^{1/2}; ZDC (East^{2} + West^{2})^{1/2}; 2<cos2#phi>");
mTotalZDCcos2phimoments->SetMinimum(0.15);
mTotalZDCcos2phimoments->SetMaximum(0.4);
mTotalZDCcos2phimoments->Draw();
gPad->Print( "plots/data/ZDC/plot_mTotalZDCcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mTotalZDCcos2phimoments.png" );

makeCanvas();
mLowZDCPTcos2phimoments->SetLineColor(kBlack);
mLowZDCPTcos2phimoments->SetTitle("cos(2#phi) moments vs. P_{T} at ZDC (East^{2} + West^{2})^{1/2} < 300; P_{T} (GeV/c); 2<cos2#phi>");
mLowZDCPTcos2phimoments->Draw();
gPad->Print( "plots/data/ZDC/plot_mLowZDCPTcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mLowZDCPTcos2phimoments.png" );

makeCanvas();
mHighZDCPTcos2phimoments->SetLineColor(kBlack);
mHighZDCPTcos2phimoments->SetTitle("cos(2#phi) moments vs. P_{T} at ZDC (East^{2} + West^{2})^{1/2} > 600; P_{T} (GeV/c); 2<cos2#phi>");
mHighZDCPTcos2phimoments->Draw();
gPad->Print( "plots/data/ZDC/plot_mHighZDCPTcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mHighZDCPTcos2phimoments.png" );

makeCanvas();
mLowZDCPTcos2phimoments->SetLineColor(kBlack);
mHighZDCPTcos2phimoments->SetLineColor(kMagenta);
mLowZDCPTcos2phimoments->SetTitle("cos(2#phi) moments vs. P_{T}; P_{T} (GeV/c); 2<cos2#phi>");
mLowZDCPTcos2phimoments->Draw();
mHighZDCPTcos2phimoments->Draw("SAME");
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->AddEntry(mLowZDCPTcos2phimoments,"ZDC (East^{2} + West^{2})^{1/2} < 300");
legend3->AddEntry(mHighZDCPTcos2phimoments,"ZDC (East^{2} + West^{2})^{1/2} > 600");
legend3->Draw();
gPad->Print( "plots/data/ZDC/plot_mDiffZDCPTcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mDiffZDCPTcos2phimoments.png" );

makeCanvas();
mv2PTcos2phimoments->SetLineColor(kBlack);
mlabPTcos2phimoments->SetLineColor(kGreen);
mv2PTcos2phimoments->SetTitle("Strength of cos2#phi modulation; P_{T} (GeV/c); 2<cos2#phi>");
mv2PTcos2phimoments->Draw();
mlabPTcos2phimoments->Draw("SAME");
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->AddEntry(mv2PTcos2phimoments,"Rest frame 2<cos2#phi>");
legend4->AddEntry(mlabPTcos2phimoments,"Lab frame 2<cos2#phi>");
legend4->Draw();
gPad->Print( "plots/data/PT/plot_mLabVRestCos2phivsPT.pdf" );
gPad->Print( "plots/data/PT/plot_mLabVRestCos2phivsPT.png" );

fo->Write();
}
