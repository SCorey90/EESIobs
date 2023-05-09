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
    auto * mPairPhi = new TH1F("mPairPhi", "#pi_{#pm} #phi distribution;#phi (rad);# events", 500, -3.13, 3.13);
    auto * mPxvsPy = new TH2F("mPxvsPy", "#rho^{0} 2D momentum dist; P_{x} (GeV); P_{y} (GeV); Counts", 200, -0.1, 0.1, 200, -0.1, 0.1);
    auto *mPhivsMass = new TH2F("mPhivsMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 50, 0.3, 1.35);
    auto *mPhivsPT = new TH2F("mPhivsPT", "#pi^{#pm} #phi distribution vs. parent P_{T}; #phi (rad); Parent P_{T} (GeV); Counts", 100, -3.14, 3.14, 100, 0, 0.25);
    auto *mPhivsRapidity = new TH2F("mPhivsRapidity", "#pi^{#pm} #phi distribution vs. Rapidity; #phi (rad); Rapidity (GeV); Counts", 100, -3.14, 3.14, 100, -2, 2);
    auto *mPhivsLowMass = new TH2F("mPhivsLowMass", "#pi_{#pm} #phi distribution vs. parent mass; #phi (rad); Parent Mass (GeV); Counts", 100, -3.14, 3.14, 100, 0.3, 0.55);

    auto * mPhiFit = new TF1("mPhiFit", "[0] + [1]*cos(x) + [2]*cos(2*x) + [3]*cos(3*x) + [4]*cos(4*x)", -3.14, 3.14);
    //Open pairDST
    TFile *myFile = TFile::Open("/Users/samcorey/code/data/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        Float_t mRapidity = pair->mRapidity;
        
        Float_t mMassVal = pair->mMass; 
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
        lv = lv1+lv2;
        lvn = lv1 - lv2;
        double Px = lv.Px();
        double Py = lv.Py();
        double Qx = lvn.Px();
        double Qy = lvn.Py();
        double absPperp = pow((Px*Px)+(Py*Py), 0.5);
        double absQperp = pow((Qx*Qx)+(Qy*Qy), 0.5);
        double PcrossQ = (Px*Qy) - (Py*Qx);
        double PdotQ = (Px*Qx) + (Py*Qy);
        double cosphi = (Px*Qx + Py*Qy) / (absPperp*absQperp);
        double PairPhi = acos(cosphi);

        if ( chipipi <10 && dca1 <1 && dca2 <1 ){
		mPperp->Fill( absPperp );
		mMass->Fill( lv.M() );
            if ( lv.M() > 0.65 && lv.M() <0.9){
                if ( PcrossQ > 0 ){
                    mPxvsPy->Fill( absPperp*cos(PairPhi), absPperp*sin(PairPhi));
                    mPhivsPT->Fill ( PairPhi - 3.1415, absPperp);
                }
		if ( PcrossQ < 0 ){
                    mPxvsPy->Fill( absPperp*cos(PairPhi), -absPperp*sin(PairPhi)); 
                    mPhivsPT->Fill ( 3.1415 - PairPhi, absPperp); 
                }
                if ( PcrossQ < 0 && absPperp < 0.06){ 
                    mPairPhi->Fill ( PairPhi - 3.1415); 
                    mPhivsRapidity->Fill ( PairPhi - 3.1415, mRapidity);
                }
                if ( PcrossQ > 0 && absPperp < 0.06){
                    mPairPhi->Fill ( 3.1415 - PairPhi );  
                    mPhivsRapidity->Fill ( 3.1415 - PairPhi, mRapidity); 
                }
            }
            if ( lv.M() > 0.2 && lv.M() <1.5){
                if ( PcrossQ < 0 && absPperp < 0.06){
                    mPhivsMass->Fill ( PairPhi - 3.1415, lv.M());
                }
                if ( PcrossQ > 0 && absPperp < 0.06){
                    mPhivsMass->Fill ( 3.1415 - PairPhi, lv.M());
                }
            }
            if ( lv.M() > 0.2 && lv.M() <0.55){
                if ( PcrossQ < 0 && absPperp < 0.06){
                    mPhivsLowMass->Fill ( PairPhi - 3.1415, lv.M());
                }
                if ( PcrossQ > 0 && absPperp < 0.06){
                    mPhivsLowMass->Fill ( 3.1415 - PairPhi, lv.M());
                }
            }

        }
    }
mPairPhi->Fit(mPhiFit);

auto *mPTmomentsplot = new TGraph(mPhivsPT->GetNbinsY(), histYbins(mPhivsPT), histMoments(mPhivsPT, 2));
auto *mMassmomentsplot = new TGraph(mPhivsMass->GetNbinsY(), histYbins(mPhivsMass), histMoments(mPhivsMass, 2));
auto *mLowMassmomentsplot = new TGraph(mPhivsLowMass->GetNbinsY(), histYbins(mPhivsLowMass), histMoments(mPhivsLowMass, 2));
auto *mRapiditymomentsplot = new TGraph(mPhivsRapidity->GetNbinsY(), histYbins(mPhivsRapidity), histMoments(mPhivsRapidity, 2));
auto *mPTcos4phimoments = new TGraph(mPhivsPT->GetNbinsY(), histYbins(mPhivsPT), histMoments(mPhivsPT, 4));
auto *mMasscos4phimoments = new TGraph(mPhivsMass->GetNbinsY(), histYbins(mPhivsMass), histMoments(mPhivsMass, 4));
auto *mRapiditycos4phimoments = new TGraph(mPhivsRapidity->GetNbinsY(), histYbins(mPhivsRapidity), histMoments(mPhivsRapidity, 4));

fo -> cd();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plot_mMass.pdf" );

makeCanvas();
mPperp->SetLineColor(kBlack);
mPperp->Draw();
gPad->Print( "plot_mPperp.pdf" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
gStyle->SetOptFit();
mPairPhi->Draw();
gPad->Print( "plot_mPairPhi.pdf" );

makeCanvas();
mPhiFit->SetLineColor(kBlack);
mPhiFit->Draw();
gPad->Print( "plot_mPhiFit.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPxvsPy->Draw("colz");
gPad->Print( "plot_mPxvsPy.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsMass->Draw("colz");
gPad->Print( "plot_mPhivsMass.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsLowMass->Draw("colz");
gPad->Print( "plot_mPhivsLowMass.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsPT->Draw("colz");
gPad->Print( "plot_mPhivsPT.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPhivsRapidity->Draw("colz");
gPad->Print( "plot_mPhivsRapidity.pdf" );

makeCanvas();
mPTmomentsplot->SetTitle("cos(2#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mPTmomentsplot->Draw("AC*");
gPad->Print( "plot_mPTmomentsplot.pdf" );

makeCanvas();
mMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
mMassmomentsplot->SetMinimum(0);
mMassmomentsplot->SetMaximum(1);
mMassmomentsplot->Draw("AC*");
gPad->Print( "plot_mMassmomentsplot.pdf" );

makeCanvas();
mRapiditymomentsplot->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsplot->Draw("AC*");
gPad->Print( "plot_mRapiditymomentsplot.pdf" );

makeCanvas();
mPTcos4phimoments->SetTitle("cos(4#phi) moments vs. P_{T}; P_{T} (GeV); 4<cos(4#phi)>");
mPTcos4phimoments->Draw("AC*");
gPad->Print( "plot_mPTcos4phimoments.pdf" );

makeCanvas();
mMasscos4phimoments->SetTitle("cos(4#phi) moments vs. Mass; Mass (GeV); 4<cos(4#phi)>");
mMasscos4phimoments->Draw("AC*");
gPad->Print( "plot_mMasscos4phimoments.pdf" );

makeCanvas();
mRapiditycos4phimoments->SetTitle("cos(4#phi) moments vs. Rapidity; Rapidity; 4<cos(4#phi)>");
mRapiditycos4phimoments->Draw("AC*");
gPad->Print( "plot_mRapiditycos4phimoments.pdf" );

makeCanvas();
mLowMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
//mLowMassmomentsplot->SetMinimum(0);
//mLowMassmomentsplot->SetMaximum(1);
mLowMassmomentsplot->Draw("AC*");
gPad->Print( "plot_mLowMassmomentsplot.pdf" );


fo->Write();
}
