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
    //can->SetRightMargin(0.01);
}

void EESI() {
    TH1F("h1", "ntuple", 100, -4, 4);

    TFile * fo = new TFile( "EESIplots.root", "RECREATE" );

    auto * mMass = new TH1F("mMass", "Parent (#rho^{0})  Mass; Mass (GeV); Counts", 500, 0.65, 0.9);
    auto * mPairPhi = new TH1F("mPairPhi", "#pi_{#pm} #phi distribution;#phi (rad);# events", 500, -6.28, 6.28);
    auto * mPxvsPy = new TH2F("mPxvsPy", "#rho^{0} 2D momentum dist; P_{x} (GeV); P_{y} (GeV); Counts", 200, -0.1, 0.1, 200, -0.1, 0.1);

    //Open pairDST
    TFile *myFile = TFile::Open("/Users/samcorey/code/data/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    TLorentzVector lv1, lv2, lv, lvn;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        
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
            if ( lv.M() > 0.65 && lv.M() <0.9){
                mMass->Fill( lv.M() );
                if ( PcrossQ > 0 ){mPxvsPy->Fill( absPperp*cos(PairPhi), absPperp*sin(PairPhi)); }
		if ( PcrossQ < 0 ){mPxvsPy->Fill( absPperp*cos(PairPhi), -absPperp*sin(PairPhi)); }
                if ( PcrossQ < 0 && absPperp < 0.06){ mPairPhi->Fill ( PairPhi - 3.1415); }
                if ( PcrossQ > 0 && absPperp < 0.06){ mPairPhi->Fill ( 3.1415 - PairPhi ); }
            }
        }
    }

fo -> cd();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plot_mMass.pdf" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mPairPhi->Draw();
gPad->Print( "plot_mPairPhi.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mPxvsPy->Draw("colz");
gPad->Print( "plot_mPxvsPy.pdf" );

fo->Write();
}
