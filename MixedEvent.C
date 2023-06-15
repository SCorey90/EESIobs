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

vector<TLorentzVector> update_buffer( vector<TLorentzVector> buffer, TLorentzVector current_ptcl, int buffer_size ) {
    if ( buffer.size() < buffer_size ) { buffer.push_back(current_ptcl); }
    if ( buffer.size() >= buffer_size ) {
        buffer.resize(buffer_size);
        int rand_int = rng.Integer(buffer_size);
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

    double PI = 3.1415926535;
 
    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "MixedEventplots.root", "RECREATE" );
    
    //define histograms
    auto * mPiDataLikePT = new TH1F("mPiDataLikePT", "#pi Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mPiDataLikeEta = new TH1F("mPiDataLikeEta", "#pi Rapidity; Rapidity; counts", 250, -2.5, 2.5);
    auto * mPiDataLikeAngle = new TH1F("mPiDataLikeAngle", "#pi Azimuthal angle; #phi (rad); counts", 100, -3.14, 3.14);
    auto * mPiDataLikeMass = new TH1F("mPiDataLikeMass", "Like Sign #pi Mass; Mass (GeV); counts", 500, 0, 1.5);
    auto * mDataLikePiPTvsRhoPT = new TH2F("mDataLikePiPTvsRhoPT", "Like sign pairs #pi vs. #rho^{0} P_{T}; #pi P_T; #rho P_{T}, counts", 500, 0, 1.5, 500, 0, 1.5);

    auto * mDataLikePT = new TH1F("mDataLikePT", "Same sign data P_{T} distribution; P_{T} (GeV/c); counts", 500, 0, 1.5);
    auto * mDataLikePhi = new TH1F("mDataLikePhi", "Same sign data #phi distribution; #phi (rad); counts", 200, -3.13, 3.13);
    auto * mDataLikeCos2PhivsPT = new TH2F("mDataLikeCos2PhivsPT", "Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 250, 0, 0.5);
    
    auto * mDataLikeCos2PhivsMass = new TH2F("mDataLikeCos2PhivsMass", "Cos2#phi signal vs. parent Mass; 2cos2#phi;  Mass (GeV); counts", 100, -2, 2, 250, 0.25, 2);

    auto * mDataLikeCos2PhivsRapidity = new TH2F("mDataLikeCos2PhivsRapidity", "Cos2#phi signal vs. parent Rapidity; 2cos2#phi;  Rapidity; counts", 100, -2, 2, 250, -2, 2);

    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum; Transverse Momentum (GeV); counts", 500, 0, 1.5);
    auto * mMass = new TH1F("mMass", "Parent (#rho^{0})  Mass; Mass (GeV); Counts", 500, 0, 2);
    auto * mRapidity = new TH1F("mRapidity", "Parent (#rho^{0})  Rapidity; Rapidity; Counts", 500, -6, 6);

    auto * mPairPhi = new TH1F("mPairPhi", "#pi^{+}#pi^{-} mixed event #phi distribution; #phi (rad); counts", 200, -3.13, 3.13);
    auto * mCos1PhivsPT = new TH2F("mCos1PhivsPT", "Cos1#phi signal vs. parent P_{T}; 2cos1#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 250, 0, 0.5);
    auto * mCos2PhivsPT = new TH2F("mCos2PhivsPT", "Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 250, 0, 0.5);
    auto * mCos3PhivsPT = new TH2F("mCos3PhivsPT", "Cos3#phi signal vs. parent P_{T}; 2cos3#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 250, 0, 0.5);
    auto * mCos4PhivsPT = new TH2F("mCos4PhivsPT", "Cos4#phi signal vs. parent P_{T}; 2cos4#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 250, 0, 0.5);
    auto * mCos2PhivsMass = new TH2F("mCos2PhivsMass", "Cos2#phi signal vs. parent Mass; 2cos2#phi;  Mass (GeV); counts", 100, -2, 2, 250, 0.25, 2);
    auto * mCos2PhivsRapidity = new TH2F("mCos2PhivsRapidity", "Cos2#phi signal vs. parent Rapidity; 2cos2#phi;  Rapidity; counts", 100, -2, 2, 250, -2, 2);
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

        int buffer_size = 10;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
            
        if ( chipipi <10 && dca1 <1 && dca2 <1 ){
            if (abs(pair->mChargeSum) == 0 ) {
                mDataLikePT->Fill((lv1+lv2).Pt());
                mPiDataLikePT->Fill(lv1.Pt());
                mPiDataLikeEta->Fill(lv1.Eta());
                mPiDataLikeAngle->Fill(lv1.Phi());
                mPiDataLikeMass->Fill(lv1.M());
                mDataLikePiPTvsRhoPT->Fill( lv1.Pt(), (lv1 + lv2).Pt() );
                if ( (lv1+lv2).M() < 0.5 && (lv1+lv2).M() > 0.4) {
                    mDataLikePhi->Fill(calc_Phi(lv1, lv2));
                    mDataLikeCos2PhivsPT->Fill( 2*cos(2*calc_Phi(lv1,lv2)), (lv1+lv2).Pt());
                }
                if ( (lv1+lv2).Pt() < 0.6 ) {
                    mDataLikeCos2PhivsMass->Fill( 2*cos(2*calc_Phi(lv1,lv2)), (lv1+lv2).M());
                }
                if ( (lv1+lv2).Pt() < 0.6 && (lv1+lv2).M() < 0.5 && (lv1+lv2).M() > 0.4  ) {
                    mDataLikeCos2PhivsRapidity->Fill( 2*cos(2*calc_Phi(lv1,lv2)), (lv1+lv2).Rapidity());
                }
            }
            if (posBuffer.size() == buffer_size && negBuffer.size() == buffer_size && abs(pair->mChargeSum) == 0 /*&& (lv1+lv2).Pt() > 0.125*/ ) {
                TLorentzVector parentLV = lv1 + lv2;
                //opposite sign pairs
                for (int i = 0; i < buffer_size; i++ ){
                    TLorentzVector Sum1 = lv1 + negBuffer[i];
                    lv1.Boost(-Sum1.BoostVector());
                    negBuffer[i].Boost(-Sum1.BoostVector());
                    //if ( parentLV.Pt() > Sum1.Pt()  ) {
                    //if ( parentLV.Pt() <= Sum1.Pt()  ) {
                    //if ( abs(parentLV.Pt()-Sum1.Pt()) < 0.05*parentLV.Pt()  ) {
                        lv1.Boost(Sum1.BoostVector());
                        negBuffer[i].Boost(Sum1.BoostVector()); 
                        posPtcls.push_back(lv1);
                        negPtcls.push_back(negBuffer[i]);
                    //}
                }
                for (int i = 0; i < buffer_size; i++ ){
                    TLorentzVector Sum2 = lv2 + posBuffer[i];
                    lv2.Boost(-Sum2.BoostVector());
                    posBuffer[i].Boost(-Sum2.BoostVector());
                    //if ( parentLV.Pt() > Sum2.Pt() ) {
                    //if ( parentLV.Pt() <= Sum2.Pt()  ) {
                    //if ( abs(parentLV.Pt()-Sum2.Pt()) < 0.05*parentLV.Pt()  ) {
                        lv2.Boost(Sum2.BoostVector());
                        posBuffer[i].Boost(Sum2.BoostVector());
                        negPtcls.push_back(lv2);
                        posPtcls.push_back(posBuffer[i]);
                    //}
                }
                //like sign pairs
                for (int i = 0; i < buffer_size; i++ ){
                    //if ( abs( lv1.Pt() - posBuffer[i].Pt() ) < 0.05*lv1.Pt() ) {
                        likesign1.push_back(lv1);
                        likesign2.push_back(posBuffer[i]);
                    //}
                }
                for (int i = 0; i < buffer_size; i++ ){
                    //if ( abs( lv2.Pt() - negBuffer[i].Pt() ) < 0.05*lv2.Pt() ) {
                        likesign2.push_back(lv2);
                        likesign1.push_back(negBuffer[i]);
                    //}
                }
            }
            posBuffer = update_buffer( posBuffer, lv1, buffer_size );
            negBuffer = update_buffer( negBuffer, lv2, buffer_size );
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
        double pairRapidity = (posLV + negLV).Rapidity();
        double pairPhi = calc_Phi( posLV, negLV );
        double cos2phi = cos( 2 * pairPhi );

        //like sign pair values
        double likePT = (likeLV1 + likeLV2).Pt();
        double likePhi = calc_Phi( likeLV1, likeLV2 );
        double likecos2phi = cos( 2 * likePhi );

        //basic histograms
        mPperp->Fill( pairPT );
        mMass->Fill( pairMass );
        mRapidity->Fill( (posLV + negLV).Rapidity() );

        //phi histograms
        if ( pairMass > 0.4 && pairMass < 0.5 ) {
            mPairPhi->Fill( pairPhi );
            mCos1PhivsPT->Fill( 2 * cos( pairPhi ), pairPT );
            mCos2PhivsPT->Fill( (2 * cos2phi), pairPT );
            mCos3PhivsPT->Fill( 2 * cos( 3*pairPhi ), pairPT );
            mCos4PhivsPT->Fill( 2 * cos( 4*pairPhi ), pairPT );
            mLikePhi->Fill( likePhi );
            mLikeCos2PhivsPT->Fill( (2 * likecos2phi), likePT );
        }
        if ( pairPT < 0.6 ) {
            mCos2PhivsMass->Fill( (2 * cos2phi), pairMass );
        }
        if ( pairPT < 0.6 && pairMass > 0.4 && pairMass < 0.5 ) {
            mCos2PhivsRapidity->Fill( (2*cos2phi), pairRapidity );
        }
    }

//Profile histograms
TProfile * mCos1PhivsPTmoments = mCos1PhivsPT->ProfileY("mCos1PhivsPTmoments", 1, -1);
TProfile * mCos2PhivsPTmoments = mCos2PhivsPT->ProfileY("mCos2PhivsPTmoments", 1, -1);
TProfile * mCos3PhivsPTmoments = mCos3PhivsPT->ProfileY("mCos3PhivsPTmoments", 1, -1);
TProfile * mCos4PhivsPTmoments = mCos4PhivsPT->ProfileY("mCos4PhivsPTmoments", 1, -1);
auto * mLikeCos2PhivsPTmoments = mLikeCos2PhivsPT->ProfileY("mLikeCos2PhivsPTmoments", 1, -1);
TProfile * mDataLikeCos2PhivsPTmoments = mDataLikeCos2PhivsPT->ProfileY("mDataLikeCos2PhivsPTmoments", 1, -1);

auto * mCos2PhivsMassmoments = mCos2PhivsMass->ProfileY("mCos2PhivsMassmoments", 1, -1);
auto * mDataLikeCos2PhivsMassmoments = mDataLikeCos2PhivsMass->ProfileY("mDataLikeCos2PhivsMassmoments", 1, -1);

auto * mCos2PhivsRapiditymoments = mCos2PhivsRapidity->ProfileY("mCos2PhivsRapiditymoments", 1, -1);
auto * mDataLikeCos2PhivsRapiditymoments = mDataLikeCos2PhivsRapidity->ProfileY("mDataLikeCos2PhivsRapiditymoments", 1, -1);

//Trying corrections
mCos2PhivsPT->Scale(2*PI*0.5/mPairPhi->Integral("width"));
mDataLikeCos2PhivsPT->Scale(2*PI*0.5/mDataLikePhi->Integral("width"));
mDataLikeCos2PhivsPT->Add(mDataLikeCos2PhivsPT, mCos2PhivsPT, 1, -1);

auto * mCorrectedCos2PhivsPTmoments = new TH1F("mCorrectedCos2PhivsPTmoments", "Corrected  (0.4 < M_{#rho^{0}} < 0.5 GeV); P_{T} (GeV/c); 2<cos2#phi>", 250, 0, 0.5);
for (int i = 0; i < (mCorrectedCos2PhivsPTmoments->GetNbinsX()) - 1; i++) { 
    double gamma_2 = mDataLikeCos2PhivsPTmoments->GetBinContent( i+1 );
    double omega_2 = mCos2PhivsPTmoments->GetBinContent( i+1 );
    double delta_gamma_2 =  mDataLikeCos2PhivsPTmoments->GetBinError( i+1 );
    double delta_omega_2 =  mCos2PhivsPTmoments->GetBinError( i+1 );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * alpha_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedCos2PhivsPTmoments->SetBinContent( i+1, alpha_2 );
    mCorrectedCos2PhivsPTmoments->SetBinError( i+1, delta_alpha_2);
}

auto * mCorrectedCos2PhivsMassmoments = new TH1F("mCorrectedCos2PhivsMassmoments", " Corrected (P_{T}(#rho^{0}) < 0.6 GeV/c); Pair Mass (GeV); 2<cos2#phi>", 250, 0.25, 2);
for (int i = 0; i < (mCorrectedCos2PhivsMassmoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mDataLikeCos2PhivsMassmoments->GetBinContent( i+1 );
    double omega_2 = mCos2PhivsMassmoments->GetBinContent( i+1 );
    double delta_gamma_2 =  mDataLikeCos2PhivsMassmoments->GetBinError( i+1 );
    double delta_omega_2 =  mCos2PhivsMassmoments->GetBinError( i+1 );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * alpha_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedCos2PhivsMassmoments->SetBinContent( i+1, alpha_2 );
    mCorrectedCos2PhivsMassmoments->SetBinError( i+1, delta_alpha_2);
}

auto * mCorrectedCos2PhivsRapiditymoments = new TH1F("mCorrectedCos2PhivsRapiditymoments", "Corrected  (P_{T}(#rho^{0}) < 0.6 GeV/c, 0.4 < M_{#rho^{0}} < 0.5 GeV); Pair Rapidity; 2<cos2#phi>", 250, -2, 2);
for (int i = 0; i < (mCorrectedCos2PhivsRapiditymoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mDataLikeCos2PhivsRapiditymoments->GetBinContent( i+1 );
    double omega_2 = mCos2PhivsRapiditymoments->GetBinContent( i+1 );
    double delta_gamma_2 =  mDataLikeCos2PhivsRapiditymoments->GetBinError( i+1 );
    double delta_omega_2 =  mCos2PhivsRapiditymoments->GetBinError( i+1 );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (4 - 2*(omega_2 * omega_2)) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * alpha_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedCos2PhivsRapiditymoments->SetBinContent( i+1, alpha_2 );
    mCorrectedCos2PhivsRapiditymoments->SetBinError( i+1, delta_alpha_2);
}


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
gStyle->SetPalette(1);
mDataLikePiPTvsRhoPT->Draw("colz");
gPad->Print( "plots/data/PT/plot_mLikePiPTvsRhoPt.pdf" );
gPad->Print( "plots/data/PT/plot_mLikePiPTvsRhoPt.png" );

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
mCos1PhivsPTmoments->SetLineColor(kBlack);
mCos1PhivsPTmoments->SetTitle( "Cos1#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos1#phi>" );
mCos1PhivsPTmoments->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mCos1PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos1PhivsPTmoments.pdf" );

makeCanvas();
mCos2PhivsPTmoments->SetLineColor(kBlack);
mCos2PhivsPTmoments->SetTitle( "Cos2#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos2#phi>" );
mCos2PhivsPTmoments->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos2PhivsPTmoments.pdf" );

makeCanvas();
mCos3PhivsPTmoments->SetLineColor(kBlack);
mCos3PhivsPTmoments->SetTitle( "Cos3#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos3#phi>" );
mCos3PhivsPTmoments->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mCos3PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos3PhivsPTmoments.pdf" );

makeCanvas();
mCos4PhivsPTmoments->SetLineColor(kBlack);
mCos4PhivsPTmoments->SetTitle( "Cos4#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos4#phi>" );
mCos4PhivsPTmoments->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mCos4PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCos4PhivsPTmoments.pdf" );

makeCanvas();
mCos2PhivsPTmoments->SetLineColor(kBlack);
mDataLikeCos2PhivsPTmoments->SetLineColor(kYellow);
mCorrectedCos2PhivsPTmoments->SetLineColor(kBlue);
mCos2PhivsPTmoments->SetTitle( "0.4 < M_{#rho^{0}} < 0.5 GeV; P_{T} (GeV/c); 2<cos2#phi>" );
mCos2PhivsPTmoments->SetMaximum(1);
mCos2PhivsPTmoments->Draw();
mDataLikeCos2PhivsPTmoments->Draw("SAME");
mCorrectedCos2PhivsPTmoments->Draw("SAME");
auto legend2 = new TLegend(0.65,0.1,0.95,0.4);
legend2->SetHeader("Legend","C"); // option "C" allows to center the header
legend2->AddEntry(mCos2PhivsPTmoments,"Mixed event unlike sign source");
legend2->AddEntry(mDataLikeCos2PhivsPTmoments,"Unike sign pairs from data");
legend2->AddEntry(mCorrectedCos2PhivsPTmoments,"Baseline subtracted signal");
legend2->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mOppositevsLikeCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mOppositevsLikeCos2PhivsPTmoments.pdf" );

makeCanvas();
mCos2PhivsMassmoments->SetLineColor(kBlack);
mDataLikeCos2PhivsMassmoments->SetLineColor(kGreen);
mCorrectedCos2PhivsMassmoments->SetLineColor(kBlue);
mCos2PhivsMassmoments->SetTitle( "P_{T}(#rho^{0}) < 0.6 GeV/c; Pair Mass (GeV); 2<cos2#phi>" );
mCos2PhivsMassmoments->SetMaximum(1);
mCos2PhivsMassmoments->Draw();
mDataLikeCos2PhivsMassmoments->Draw("Same");
mCorrectedCos2PhivsMassmoments->Draw("Same");
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->SetHeader("Legend","C"); // option "C" allows to center the header
legend3->AddEntry(mCos2PhivsMassmoments,"Mixed event unlike sign source");
legend3->AddEntry(mDataLikeCos2PhivsMassmoments,"Unlike sign pairs from data");
legend3->AddEntry(mCorrectedCos2PhivsMassmoments,"Baseline subtracted signal");
legend3->Draw();
gPad->Print( "plots/mixedevent/Mass/plot_mDataLikevsMixedCos2PhivsMassmoments.png" );
gPad->Print( "plots/mixedevent/Mass/plot_mDataLikevsMixedCos2PhivsMassmoments.pdf" );

makeCanvas();
mCos2PhivsRapiditymoments->SetLineColor(kBlack);
mDataLikeCos2PhivsRapiditymoments->SetLineColor(kGreen);
mCorrectedCos2PhivsRapiditymoments->SetLineColor(kBlue);
mCos2PhivsRapiditymoments->SetTitle( "P_{T}(#rho^{0}) < 0.6 GeV/c, 0.4 < M_{#rho^{0}} < 0.5 GeV; Pair Rapidity; 2<cos2#phi>" );
mCos2PhivsRapiditymoments->SetMaximum(1);
mCos2PhivsRapiditymoments->Draw();
mDataLikeCos2PhivsRapiditymoments->Draw("Same");
mCorrectedCos2PhivsRapiditymoments->Draw("Same");
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->SetHeader("Legend","C"); // option "C" allows to center the header
legend4->AddEntry(mCos2PhivsRapiditymoments,"Mixed event unlike sign source");
legend4->AddEntry(mDataLikeCos2PhivsRapiditymoments,"Unike sign pairs from data");
legend4->AddEntry(mCorrectedCos2PhivsRapiditymoments,"Baseline subtracted signal");
legend4->Draw();
gPad->Print( "plots/mixedevent/Mass/plot_mDataLikevsMixedCos2PhivsRapiditymoments.png" );
gPad->Print( "plots/mixedevent/Mass/plot_mDataLikevsMixedCos2PhivsRapiditymoments.pdf" );

/*makeCanvas();
mCorrectedCos2PhivsPTmoments->SetLineColor(kBlack);
mCorrectedCos2PhivsPTmoments->SetMaximum(1);
mCorrectedCos2PhivsPTmoments->SetMinimum(-1.5);
mCorrectedCos2PhivsPTmoments->Draw("E1");
gPad->Print( "plots/mixedevent/PT/plot_mCorrectedCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mCorrectedCos2PhivsPTmoments.pdf" );

makeCanvas();
mCorrectedCos2PhivsMassmoments->SetLineColor(kBlack);
mCorrectedCos2PhivsMassmoments->SetMaximum(1);
mCorrectedCos2PhivsMassmoments->SetMinimum(-1.5);
mCorrectedCos2PhivsMassmoments->Draw("E1");
gPad->Print( "plots/mixedevent/mass/plot_mCorrectedCos2PhivsMassmoments.png" );
gPad->Print( "plots/mixedevent/mass/plot_mCorrectedCos2PhivsMassmoments.pdf" );

makeCanvas();
mCorrectedCos2PhivsRapiditymoments->SetLineColor(kBlack);
mCorrectedCos2PhivsRapiditymoments->SetMaximum(1);
mCorrectedCos2PhivsRapiditymoments->SetMinimum(-1.5);
mCorrectedCos2PhivsRapiditymoments->Draw("E1");
gPad->Print( "plots/mixedevent/Rapidity/plot_mCorrectedCos2PhivsRapiditymoments.png" );
gPad->Print( "plots/mixedevent/Rapidity/plot_mCorrectedCos2PhivsRapiditymoments.pdf" );
*/
fo->Write();
}
