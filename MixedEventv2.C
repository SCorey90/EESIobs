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
    //auto * mDataCos2phivsPT = new TH2F("mDataCos2phivsPT", "Data Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    auto * mDataCos1phivsPTvsMass = new TH3F("mDataCos1phivsPTvsMass", "Data Cos1#phi signal vs. parent P_{T}; 2cos1#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mDataCos2phivsPTvsMass = new TH3F("mDataCos2phivsPTvsMass", "Data Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mDataCos3phivsPTvsMass = new TH3F("mDataCos3phivsPTvsMass", "Data Cos3#phi signal vs. parent P_{T}; 2cos3#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mDataCos4phivsPTvsMass = new TH3F("mDataCos4phivsPTvsMass", "Data Cos4#phi signal vs. parent P_{T}; 2cos4#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );


    //mixed histograms
    //auto * mMixedCos2phivsPT = new TH2F("mMixedCos2phivsPT", "Mixed Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); counts", 100, -2, 2, 100, 0, 0.5);

    auto * mMixedCos1phivsPTvsMass = new TH3F("mMixedCos1phivsPTvsMass", "Mixed Cos1#phi signal vs. parent P_{T}; 2cos1#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mMixedCos2phivsPTvsMass = new TH3F("mMixedCos2phivsPTvsMass", "Mixed Cos2#phi signal vs. parent P_{T}; 2cos2#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mMixedCos3phivsPTvsMass = new TH3F("mMixedCos3phivsPTvsMass", "Mixed Cos3#phi signal vs. parent P_{T}; 2cos3#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );
    auto * mMixedCos4phivsPTvsMass = new TH3F("mMixedCos4phivsPTvsMass", "Mixed Cos4#phi signal vs. parent P_{T}; 2cos4#phi;  P_{T} (GeV/c); Mass (GeV); counts", 100, -2, 2, 100, 0, 0.5, 100, 0.25, 2 );

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

        int buffer_size = 200;

        TLorentzVector lv1, lv2;
        lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
        lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );

        double ddTOFval = ddTOF( 0.135, lv1.P(), pair->d1_mTof, pair->d1_mLength, lv2.P(), pair->d2_mTof, pair->d2_mLength );

        if ( chipipi <10 && dca1 <1 && dca2 <1 && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.3  ){
            
            //fill histograms
            if ( abs(pair->mChargeSum)==0 && abs((lv1+lv2).Rapidity()) < 1) {
                TLorentzVector lv = lv1 + lv2 ;
                double dataPairPhi = calc_Phi(lv1, lv2);
                mDataCos1phivsPTvsMass->Fill( 2*cos(1*dataPairPhi), lv.Pt(), lv.M() );
                mDataCos2phivsPTvsMass->Fill( 2*cos(2*dataPairPhi), lv.Pt(), lv.M() );
                mDataCos3phivsPTvsMass->Fill( 2*cos(3*dataPairPhi), lv.Pt(), lv.M() );
                mDataCos4phivsPTvsMass->Fill( 2*cos(4*dataPairPhi), lv.Pt(), lv.M() );
            }

            //add to mixed particle list
            if ( abs(pair->mChargeSum)==2 ) {
                for ( TLorentzVector mixed_lv1 : posBuffer ) {
                    TLorentzVector mixed_lv = mixed_lv1 + lv2 ;
                    double mixedPairPhi = calc_Phi(mixed_lv1, lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 ) {
                        //if (mixed_lv.M() > 0.65 && mixed_lv.M() < 0.66 ) { mMixedCos2phivsPT->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt() ); }

                        mMixedCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }    
                }
                for ( TLorentzVector mixed_lv2 : negBuffer ) {
                    TLorentzVector mixed_lv = lv1 + mixed_lv2 ;
                    double mixedPairPhi = calc_Phi(lv1, mixed_lv2);
                    if ( abs(mixed_lv.Rapidity()) < 1 ) {
                        //if (mixed_lv.M() > 0.65 && mixed_lv.M() < 0.66 ) { mMixedCos2phivsPT->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt() ); }
 
                        mMixedCos1phivsPTvsMass->Fill( 2*cos(1*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos2phivsPTvsMass->Fill( 2*cos(2*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos3phivsPTvsMass->Fill( 2*cos(3*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                        mMixedCos4phivsPTvsMass->Fill( 2*cos(4*mixedPairPhi), mixed_lv.Pt(), mixed_lv.M() );
                    }
                }
            }

            //update buffers
            int buff_index = rng.Integer(buffer_size - 1 );
            posBuffer = update_buffer( posBuffer, lv1, buffer_size, buff_index );
            negBuffer = update_buffer( negBuffer, lv2, buffer_size, buff_index );
        }
    }
/*
auto * mDataTest = mDataCos2phivsPT->ProfileY("mDataTest", 1, -1);
auto * mMixedTest = mMixedCos2phivsPT->ProfileY("mMixedTest", 1, -1);

makeCanvas();
mDataTest->SetLineColor(kBlack);
mMixedTest->SetLineColor(kBlue);
mDataTest->Draw();
mMixedTest->Draw("SAME");
*/

fo->Write();

}
