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

vector<TLorentzVector> update_buffer( vector<TLorentzVector> buffer, TLorentzVector current_ptcl, int buffer_size, int index ) {
    if ( buffer.size() < buffer_size ) { buffer.push_back(current_ptcl); }
    if ( buffer.size() >= buffer_size ) {
        buffer.resize(buffer_size);
        buffer[index] = current_ptcl;
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

double ddTOF ( double mass, double p1, double tof1, double l1,  double p2, double tof2, double l2 ) {

    double c = 3.0e1; //in cm/ns
    double me2 = pow(mass,2);
    double dTOF = tof2 - tof1 ;
    double tof1exp = (l1 / c) * sqrt( 1 + ( (mass*mass)/(p1*p1) ) );
    double tof2exp = (l2 / c) * sqrt( 1 + ( (mass*mass)/(p2*p2) ) );
    double dTOFexp = tof2exp - tof1exp ;
    double ddTOF = dTOF - dTOFexp ;

    return ddTOF;
}

void MixedEventv2() {
    
    double PI = 3.14159265;

    TH1F("h1", "ntuple", 100, -4, 4);
    TFile * fo = new TFile( "MixedEventplotsv2.root", "RECREATE" );

    //data histograms
    auto * mDataCos2phivsPT = new TH2F("mDataCos2phivsPT", "Data Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    //mixed histograms
    auto * mMixedCos2phivsPT = new TH2F("mMixedCos2phivsPT", "Mixed Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    //data analysis and filling buffers
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
        double eZDC = pair->mZDCEast;
        double wZDC = pair->mZDCWest;
        double vertex = pair->mVertexZ;

        int buffer_size = 10;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );

        double ddTOFval = ddTOF( 0.135, lv1.P(), pair->d1_mTof, pair->d1_mLength, lv2.P(), pair->d2_mTof, pair->d2_mLength );

        if ( chipipi <10 && dca1 <1 && dca2 <1 && abs(pair->mChargeSum)==0 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.3  ){
            
            //fill histograms

            double dataPairPhi = calc_Phi(lv1, lv2);
            if ( (lv1+lv2).M() > 0.65 && (lv1+lv2).M() < 0.9 ) {
                mDataCos2phivsPT->Fill( 2*cos(2*dataPairPhi), (lv1+lv2).Pt() );
            }

            //add to mixed particle list
            if ( posBuffer.size() == buffer_size && negBuffer.size() == buffer_size ) {
                for ( int i = 0; i < buffer_size; i++ ) {
                    posPtcls.push_back(lv1);
                    negPtcls.push_back(negBuffer[i]);
                }
                for ( int i = 0; i < buffer_size; i++ ) {
                    posPtcls.push_back(posBuffer[i]);
                    negPtcls.push_back(lv2);
                }
            }

            //update buffers
            int buff_index = rng.Integer(buffer_size - 1 );
            posBuffer = update_buffer( posBuffer, lv1, buffer_size, buff_index );
            negBuffer = update_buffer( negBuffer, lv2, buffer_size, buff_index );
        }
    }
    cout << posBuffer.size();
    //analyze mixed pairs
    for ( int i = 0; i < posPtcls.size(); i++ ) {
        TLorentzVector mixed_parent = posPtcls[i] + negPtcls[i] ;
        double mixedPairPhi = calc_Phi( posPtcls[i], negPtcls[i] );

        //fill histograms
        if ( mixed_parent.M() > 0.65 && mixed_parent.M() < 0.9 ) {
            mMixedCos2phivsPT->Fill( 2*cos(2*mixedPairPhi), mixed_parent.Pt() );
        }
    }

//profile histograms

TProfile * mDataCos2phivsPTmoments = mDataCos2phivsPT->ProfileY("mDataCos2phivsPTmoments", 1, -1);

TProfile * mMixedCos2phivsPTmoments = mMixedCos2phivsPT->ProfileY("mMixedCos2phivsPTmoments", 1, -1);

//data corrections

auto * mCorrectedCos2phivsPTmoments = new TH1F("mCorrectedCos2phivsPTmoments", "Corrected  (0.65 < M_{#rho^{0}} < 0.9 GeV); P_{T} (GeV/c); 2<cos2#phi>", 100, 0, 0.5);
for (int i = 0; i < (mCorrectedCos2phivsPTmoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mDataCos2phivsPTmoments->GetBinContent( i+1 )/2;
    double omega_2 = mMixedCos2phivsPTmoments->GetBinContent( i+1 )/2;
    double delta_gamma_2 =  mDataCos2phivsPTmoments->GetBinError( i+1 )/2;
    double delta_omega_2 =  mMixedCos2phivsPTmoments->GetBinError( i+1 )/2;
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * gamma_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedCos2phivsPTmoments->SetBinContent( i+1, 2*alpha_2 );
    mCorrectedCos2phivsPTmoments->SetBinError( i+1, delta_alpha_2);
}

//draw histograms

makeCanvas();
mDataCos2phivsPTmoments->SetLineColor(kBlack);
mDataCos2phivsPTmoments->SetTitle( "Cos2#phi signal strength vs P_{T} (data); P_{T} (GeV/c); 2<cos2#phi>" );
mDataCos2phivsPTmoments->Draw();

makeCanvas();
mMixedCos2phivsPTmoments->SetLineColor(kBlack);
mMixedCos2phivsPTmoments->SetTitle( "Cos2#phi signal strength vs P_{T} (mixed event); P_{T} (GeV/c); 2<cos2#phi>" );
mMixedCos2phivsPTmoments->Draw();

makeCanvas();
mMixedCos2phivsPTmoments->SetLineColor(kBlack);
mDataCos2phivsPTmoments->SetLineColor(kGreen-2);
mCorrectedCos2phivsPTmoments->SetLineColor(kMagenta);
mMixedCos2phivsPTmoments->SetTitle( "0.65 < M_{#rho^{0}} < 0.9 GeV; P_{T} (GeV/c); 2<cos2#phi>" );
mMixedCos2phivsPTmoments->SetMaximum(1);
mMixedCos2phivsPTmoments->Draw();
mDataCos2phivsPTmoments->Draw("SAME");
mCorrectedCos2phivsPTmoments->Draw("SAME");
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->SetHeader("Legend","C"); // option "C" allows to center the header
legend->AddEntry(mMixedCos2phivsPTmoments,"Mixed event unlike sign source");
legend->AddEntry(mDataCos2phivsPTmoments,"Unike sign pairs from data");
legend->AddEntry(mCorrectedCos2phivsPTmoments,"Baseline subtracted signal");
legend->Draw();
gPad->Print( "plots/mixedevent/PT/plot_mMixedCorrCos2PhivsPTmoments.png" );
gPad->Print( "plots/mixedevent/PT/plot_mMixedCorrCos2PhivsPTmoments.pdf" );

fo->Write();

}
