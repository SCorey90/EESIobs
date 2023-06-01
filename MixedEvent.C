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
    if ( buffer.size() < 5 ) { buffer.push_back(current_ptcl); }
    if ( buffer.size() >= 5 ) {
        buffer.resize(5);
        int rand_int = rng.Integer(5);
        buffer[rand_int] = current_ptcl;
    }
    return buffer;
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

void MixedEvent() {
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "MixedEventplots.root", "RECREATE" );
    
    //define histograms
    auto * mPiDataLikePT = new TH1F("mPiDataLikePT", "#pi Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mPiDataLikeEta = new TH1F("mPiDataLikeEta", "#pi Rapidity; Rapidity; counts", 250, -2.5, 2.5);
    auto * mPiDataLikeAngle = new TH1F("mPiDataLikeAngle", "#pi Azimuthal angle; #phi (rad); counts", 100, -3.14, 3.14);
    auto * mPiDataLikeMass = new TH1F("mPiDataLikeMass", "Like Sign #pi Mass; Mass (GeV); counts", 500, 0, 1.5);

    auto * mDataLikePhi = new TH1F("mDataLikePhi", "Same sign data #phi distribution; #phi (rad); counts", 200, -3.13, 3.13);
    auto * mDataLikeCos2PhivsPT = new TH2F("mDataLikeCos2PhivsPT", "Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mMass = new TH1F("mMass", "Parent (#rho^{0})  Mass; Mass (GeV); Counts", 500, 0, 2);

    auto * mPairPhi = new TH1F("mPairPhi", "#pi^{+}#pi^{-} mixed event #phi distribution; #phi (rad); counts", 200, -3.13, 3.13);
    auto * mCos2PhivsPT = new TH2F("mCos2PhivsPT", "Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);
    auto * mLikePhi = new TH1F("mLikePhi", "#pi^{+}#pi^{-} like sign #phi distribution; #phi (rad); counts", 200, -3.13, 3.13);
    auto * mLikeCos2PhivsPT = new TH2F("mLikeCos2PhivsPT", "Cos2#phi signal vs. parent P_{T} (like sign); 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    TFile *myFile = TFile::Open("/Users/samcorey/code/data/pair_dst_Run12UU.root");
    TTreeReader myReader("PairDst", myFile);
    TTreeReaderValue<FemtoPair> pair(myReader, "Pairs");

    vector<TLorentzVector> posBuffer;
    vector<TLorentzVector> negBuffer;

    vector<TLorentzVector> posPtcls;
    vector<TLorentzVector> negPtcls;

    vector<TLorentzVector> likesign1;
    vector<TLorentzVector> likesign2;

    while (myReader.Next()) {
        double chipipi = pow( pair->d1_mNSigmaPion, 2) + pow( pair->d2_mNSigmaPion, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
            
        if ( chipipi <10 && dca1 <1 && dca2 <1 ){
            if (abs(pair->mChargeSum) == 2 ) {
                mPiDataLikePT->Fill(lv1.Pt());
                mPiDataLikeEta->Fill(lv1.Eta());
                mPiDataLikeAngle->Fill(lv1.Phi());
                mPiDataLikeMass->Fill(lv1.M());
                mDataLikePhi->Fill(calc_Phi(lv1, lv2));
                mDataLikeCos2PhivsPT->Fill( 2*cos(2*calc_Phi(lv1,lv2)), (lv1+lv2).Pt());
            }
            if (posBuffer.size() == 5 && negBuffer.size() == 5 && abs(pair->mChargeSum) == 2) {
                //opposite sign pairs
                for (int i = 0; i < 5; i++ ){
                    posPtcls.push_back(lv1);
                    negPtcls.push_back(negBuffer[i]);
                }
                for (int i = 0; i < 5; i++ ){
                    negPtcls.push_back(lv2);
                    posPtcls.push_back(posBuffer[i]);
                }
                //like sign pairs
                for (int i = 0; i < 5; i++ ){
                    likesign1.push_back(lv1);
                    likesign2.push_back(posBuffer[i]);
                }
                for (int i = 0; i < 5; i++ ){
                    likesign2.push_back(lv2);
                    likesign1.push_back(negBuffer[i]);
                }
            }
            posBuffer = update_buffer( posBuffer, lv1 );
            negBuffer = update_buffer( negBuffer, lv2 );
        }
    }
    for ( int i = 0; i < posPtcls.size(); i++ ) {
        TLorentzVector posLV = posPtcls[i];
        TLorentzVector negLV = negPtcls[i];
        TLorentzVector likeLV1 = likesign1[i];
        TLorentzVector likeLV2 = likesign2[i];
        
        //calc useful values
        double pairPT = (posLV + negLV).Pt();
        double pairMass = (posLV + negLV).M();
        double pairPhi = calc_Phi( posLV, negLV );
        double cos2phi = cos( 2 * pairPhi );

        //like sign pair values
        double likePT = (likeLV1 + likeLV2).Pt();
        double likePhi = calc_Phi( likeLV1, likeLV2 );
        double likecos2phi = cos( 2 * likePhi );

        //basic histograms
        mPperp->Fill( pairPT );
        mMass->Fill( pairMass );

        //phi histograms
        mPairPhi->Fill( pairPhi );
        mCos2PhivsPT->Fill( (2 * cos2phi), pairPT );
        mLikePhi->Fill( likePhi );
        mLikeCos2PhivsPT->Fill( (2 * likecos2phi), likePT );
        
    }

//Profile histograms
auto * mCos2PhivsPTmoments = mCos2PhivsPT->ProfileY("mCos2PhivsPTmoments", 1, -1);
auto * mLikeCos2PhivsPTmoments = mLikeCos2PhivsPT->ProfileY("mLikeCos2PhivsPTmoments", 1, -1);
auto * mDataLikeCos2PhivsPTmoments = mDataLikeCos2PhivsPT->ProfileY("mDataLikeCos2PhivsPTmoments", 1, -1);

//Draw histograms
makeCanvas();
mPiDataLikePT->SetLineColor(kBlack);
mPiDataLikePT->Draw();

makeCanvas();
mPiDataLikeEta->SetLineColor(kBlack);
mPiDataLikeEta->Draw();

makeCanvas();
mPiDataLikeAngle->SetLineColor(kBlack);
mPiDataLikeAngle->Draw();

makeCanvas();
mPiDataLikeMass->SetLineColor(kBlack);
mPiDataLikeMass->Draw();

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

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mPairPhi->Draw();
gPad->Print( "plots/mixedevent/plot_mPairPhi.png" );
gPad->Print( "plots/mixedevent/plot_mPairPhi.pdf" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
mLikePhi->SetLineColor(kGreen);
mPairPhi->Draw();
mLikePhi->Draw("SAME");
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->SetHeader("Legend","C"); // option "C" allows to center the header
legend->AddEntry(mPairPhi,"Opposite sign pairs");
legend->AddEntry(mLikePhi,"Like sign pairs");
legend->Draw();
gPad->Print( "plots/mixedevent/plot_mOppositevsLikePhi.png" );
gPad->Print( "plots/mixedevent/plot_mOppositevsLikePhi.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mCos2PhivsPT->Draw("colz");
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPT.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPT.pdf" );

makeCanvas();
mCos2PhivsPTmoments->SetLineColor(kBlack);
mCos2PhivsPTmoments->SetTitle( "Cos2#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos2#phi>" );
mCos2PhivsPTmoments->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPTmoments.pdf" );

makeCanvas();
mCos2PhivsPTmoments->SetLineColor(kBlack);
mDataLikeCos2PhivsPTmoments->SetLineColor(kGreen);
mCos2PhivsPTmoments->SetTitle( "Cos2#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos2#phi>" );
mCos2PhivsPTmoments->Draw();
mDataLikeCos2PhivsPTmoments->Draw("Same");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->SetHeader("Legend","C"); // option "C" allows to center the header
legend2->AddEntry(mCos2PhivsPTmoments,"Mixed event");
legend2->AddEntry(mDataLikeCos2PhivsPTmoments,"Like sign pairs from data");
legend2->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mOppositevsLikeCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mOppositevsLikeCos2PhivsPTmoments.pdf" );

fo->Write();
}
