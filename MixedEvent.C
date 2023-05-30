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

std::random_device global_rng;
TRandom3 rng(global_rng());

vector<TLorentzVector> update_buffer( vector<TLorentzVector> buffer, TLorentzVector current_ptcl ) {
    if ( buffer.size() < 10 ) { buffer.push_back(current_ptcl); }
    if ( buffer.size() >= 10 ) {
        buffer.resize(10);
        int rand_int = rng.Integer(10);
        buffer[rand_int] = current_ptcl;
    }
    return buffer;
}

void MixedEvent() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "MixedEventplots.root", "RECREATE" );
    
    //define histograms
    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mMass = new TH1F("mMass", "Parent (#rho^{0})  Mass; Mass (GeV); Counts", 500, 0, 2);
    
    TFile *myFile = TFile::Open("/Users/samcorey/code/data/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    vector<TLorentzVector> posBuffer;
    vector<TLorentzVector> negBuffer;

    vector<TLorentzVector> posPtcls;
    vector<TLorentzVector> negPtcls;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
            
        if ( chipipi <10 && dca1 <1 && dca2 <1 ){
            if (posBuffer.size() == 10 && negBuffer.size() == 10) {
                for (int i = 0; i < 10; i++ ){
                    posPtcls.push_back(lv1);
                    negPtcls.push_back(negBuffer[i]);
                }
                for (int i = 0; i < 10; i++ ){
                    negPtcls.push_back(lv2);
                    posPtcls.push_back(posBuffer[i]);
                }
            }
            posBuffer = update_buffer( posBuffer, lv1 );
            negBuffer = update_buffer( negBuffer, lv2 );
        }
    }
    for ( int i = 0; i < posPtcls.size(); i++ ) {
        TLorentzVector posLV = posPtcls[i];
        TLorentzVector negLV = negPtcls[i];
        mPperp->Fill( (posLV+negLV).Pt());
        mMass->Fill( (posLV+negLV).M());
    }
makeCanvas();
mPperp->SetLineColor(kBlack);
mPperp->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mPperp.png" );
gPad->Print( "plots/mixedevent/PT/plot_mPperp.pdf" );

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();
gPad->Print( "plots/mixedevent/Mass/plot_mMass.png" );
gPad->Print( "plots/mixedevent/Mass/plot_mMass.pdf" );
}
