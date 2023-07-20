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

int ican2 = 0;
void makeCanvas()  {
    TCanvas * can = new TCanvas( TString::Format( "can%d", ican2++ ), "", 900, 600);
    //can->SetTopMargin(0.04);
    can->SetRightMargin(0.35);
}

std::random_device global_rng;
TRandom3 rng(global_rng());

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
    return PairPhi;
}

double lab_calc_Phi( TLorentzVector lv1, TLorentzVector lv2) {
    TLorentzVector lvPlus = lv1 + lv2;
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

void EESI() {
    TH1F("h1", "ntuple", 100, -4, 4);

    TFile * fo = new TFile( "EESIplots.root", "RECREATE" );

    auto * mChiChi = new TH2F("mChiChi", " ; #Chi_{#pi#pi} ; #Chi_{ee} ; counts", 100, 0, 10, 100, 0, 200);
    auto * mddTOF = new TH1F("mddTOF", " #Delta#DeltaTOF ; #Delta#DeltaTOF; counts", 500, -1, 1 );

    auto * mDCA1 = new TH1F("mDCA1", "DCA1", 500, 0, 10);
    auto * mDCA1v2 = new TH2F("mDCA1v2", "DCA1v2", 500, 0, 3, 500, 0, 3);

    auto * mPiPT = new TH1F("mPiPT", "#pi P_{T}; P_{T} (GeV/c); counts", 500, 0, 1.5);
    auto * mPiEta = new TH1F("mPiEta", "#pi Rapidity; P_{T} (GeV/c); counts", 250, -2.5, 2.5);
    auto * mPiAngle = new TH1F("mPiAngle", "#pi Azimuthal Angle; #phi (GeV/c); counts", 100, -3.14, 3.14);
    auto * mPiMass = new TH1F("mPiMass", "#pi Mass; Mass (GeV); counts", 500, 0, 1.5);
    auto * mUnlikePiPTvsRhoPT = new TH2F("mUnlikePiPTvsRhoPT", "Unlike sign pairs #pi vs. #rho^{0} P_{T}; #pi P_T; #rho P_{T}, counts", 500, 0, 1.5, 500, 0, 1.5);

    auto * mMass = new TH1F("mMass", "M_{#pi#pi} (U+U); Mass (GeV); Counts", 500, 0, 2);
    auto * mPperp = new TH1F("mPperp", "Parent (#rho^{0}) Transverse Momentum (U+U); Transverse Momentum (GeV); counts", 500, 0, 0.3);
    auto * mEta = new TH1F("mEta", "Parent (#rho^{0}) Eta; Eta; counts", 500, -6, 6);
    auto * mLVPhi = new TH1F("mLVPhi", "Parent (#rho^{0}) LV.Phi(); #phi (rad); counts", 500, -3.15, 3.15);
    auto * mZDCEast = new TH1F("mZDCEast", "ZDC East; ZDC East; counts", 200, 0, 1300);
    auto * mZDCWest = new TH1F("mZDCWest", "ZDC West; ZDC West; counts", 200, 0, 1300);

    auto * mZDCEastvsWest = new TH2F("mZDCEastvsWest", "ZDC East vs. West; ZDC East; ZDC West counts", 200, 0, 1300, 200, 0, 1300);
    auto * mZDCTotal = new TH1F("mZDCtotal", "ZDC (East^{2} + West^{2})^{1/2}; ZDC (East^{2} + West^{2})^{1/2}; counts", 1000, 0, 1900);
    
    auto * mPairPhi = new TH1F("mPairPhi", "0.65 < M_{#pi#pi} < 0.9 GeV, P_{T} < 0.06 MeV/c (U+U) ;#phi (rad);norm. counts", 100, -3.13, 3.13);
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

    auto * mCos2phivsPT = new TH2F("mCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 25, 0, 0.25);
    auto * mCos1phivsPT = new TH2F("mCos1phivsPT", "cos1#phi distribution vs P_{T}", 100, -2, 2, 10, 0, 0.25);
    auto * mCos3phivsPT = new TH2F("mCos3phivsPT", "cos3#phi distribution vs P_{T}", 100, -2, 2, 10, 0, 0.25);
    auto * mLabCos2phivsPT = new TH2F("mLabCos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 25, 0, 0.25);
    auto * mCos4phivsPT = new TH2F("mCos4phivsPT", "cos4#phi distribution vs P_{T}", 100, -2, 2, 10, 0, 0.25);

    auto * mLikeCos2phivsPT = new TH2F("mLikeCos2phivsPT", "cos2#phi distribution vs P_{T} (Like sign)", 100, -2, 2, 5, 0, 0.25);
    auto * mLikeCos1phivsPT = new TH2F("mLikeCos1phivsPT", "cos1#phi distribution vs P_{T}", 100, -2, 2, 8, 0, 0.25);
    auto * mLikeCos3phivsPT = new TH2F("mLikeCos3phivsPT", "cos3#phi distribution vs P_{T}", 100, -2, 2, 8, 0, 0.25);
    auto * mLikeCos4phivsPT = new TH2F("mLikeCos4phivsPT", "cos4#phi distribution vs P_{T}", 100, -2, 2, 8, 0, 0.25);

    auto * mLowZDCcos2phivsPT = new TH2F("mLowZDCcos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.25);
    auto * mHighZDCcos2phivsPT = new TH2F("mHighZDCcos2phivsPT", "cos2#phi distribution vs P_{T}", 100, -2, 2, 100, 0, 0.25);

    auto * mCos2phivsMass = new TH2F("mCos2phivsMass", "cos2#phi distribution vs Mass", 100, -2, 2, 25, 0.25, 1);
    auto * mCos4phivsMass = new TH2F("mCos4phivsMass", "cos4#phi distribution vs Mass", 100, -2, 2, 100, 0.3, 1.35);

    auto * mLikeCos2phivsMass = new TH2F("mLikeCos2phivsMass", "cos2#phi distribution vs Mass (like sign)", 100, -2, 2, 5, 0.25, 1);

    auto * mCos2phivsRapidity = new TH2F("mCos2phivsRapidity", "cos2#phi distribution vs Rapidity", 100, -2, 2, 100, -2, 2);
    auto * mCos4phivsRapidity = new TH2F("mCos4phivsRapidity", "cos4#phi distribution vs Rapidity", 100, -2, 2, 100, -2, 2);

    auto * mLikeCos2phivsRapidity = new TH2F("mLikeCos2phivsRapidity", "cos2#phi distribution vs Rapidity (like sign)", 100, -2, 2, 5, -2, 2);

    auto * mCorrectedPTcos2phimoments = new TH1F("mCorrectedPTcos2phimoments", "Corrected  (0.65 < M_{#pi#pi} < 0.9 GeV); P_{T} (GeV/c); 2<cos2#phi>", 25, 0, 0.25);

    auto * mCos2phivsPTvsMass = new TH3F("mCos2phivsPTvsMass", ";pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>; counts", 25, 0.25, 0.7, 25, 0, 0.5, 100, -2, 2);

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
        double chiee = pow( pair->d1_mNSigmaElectron, 2) + pow( pair->d2_mNSigmaElectron, 2);
        double dca1 = pair->d1_mDCA;
        double dca2 = pair->d2_mDCA;
        double EastZDC = pair->mZDCEast;
        double WestZDC = pair->mZDCWest;
        double TotalZDC = sqrt((EastZDC*EastZDC + WestZDC*WestZDC));
        Float_t mRapidity = pair->mRapidity;

        if ( chipipi <10 ){
            mDCA1->Fill(dca1);
            mDCA1v2->Fill(dca1, dca2);
        }       

        double ddTOFval = ddTOF( 0.135, lv1.P(), pair->d1_mTof, pair->d1_mLength, lv2.P(), pair->d2_mTof, pair->d2_mLength );
 
        Float_t mMassVal = pair->mMass;
        double rand_val = 0; //rng.Uniform(0, 1);
        if ( rand_val < 0.5 ) {
            lv1.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
            lv2.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
        } else {
            lv2.SetPtEtaPhiM( pair->d1_mPt, pair->d1_mEta, pair->d1_mPhi, 0.135 );
            lv1.SetPtEtaPhiM( pair->d2_mPt, pair->d2_mEta, pair->d2_mPhi, 0.135 );
        }
        lv = lv1+lv2;
        double absPperp = lv.Pt();
        double PairPhi = calc_Phi( lv1, lv2);
        //like sign
        if ( chipipi <10 && dca1 <1 && dca2 <1 && abs(pair->mChargeSum)==2 /* && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.3 */  ){
            if ( lv.M() > 0.65 && lv.M() <0.9 ){
                mLikeCos1phivsPT->Fill( 2*cos(1*calc_Phi(lv1, lv2)), lv.Pt() );
                mLikeCos2phivsPT->Fill( 2*cos(2*calc_Phi(lv1, lv2)), lv.Pt() );
                mLikeCos3phivsPT->Fill( 2*cos(3*calc_Phi(lv1, lv2)), lv.Pt() );
                mLikeCos4phivsPT->Fill( 2*cos(4*calc_Phi(lv1, lv2)), lv.Pt() );
            }
            if ( lv.Pt() < 0.06 ) {
                mLikeCos2phivsMass->Fill( 2*cos(2*calc_Phi(lv1, lv2)), lv.M() );
            }
            if ( lv.M() > 0.65 && lv.M() <0.9 && lv.Pt() < 0.06 ) {
                mLikeCos2phivsRapidity->Fill( 2*cos(2*calc_Phi(lv1, lv2)), mRapidity );
            }
        }
        
        if ( pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 ) {
            mddTOF->Fill( ddTOFval );
        }
        if ( abs(ddTOFval) < 0.3  && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 ) {
            mChiChi->Fill( chipipi, chiee );
        }
 
        //signal
        if ( chipipi <10 && dca1 <1 && dca2 <1 && pair->mChargeSum==0  /* && pair->d1_mMatchFlag !=0 && pair->d2_mMatchFlag !=0 && abs(ddTOFval) < 0.3 */ ){

            //mChiChi->Fill( chipipi, chiee );
      
	    mPperp->Fill( lv.Pt() );
            mMass->Fill( lv.M() );
            mEta->Fill( lv.Eta() );
            mLVPhi->Fill( lv.Phi() );
            mPiPT->Fill( lv1.Pt() );
            mPiEta->Fill( lv1.Eta() );
            mPiAngle->Fill( lv1.Phi() );
            mPiMass->Fill( lv1.M() );
            mZDCEast->Fill( EastZDC );
                    mZDCWest->Fill( WestZDC );
            mZDCTotal->Fill( TotalZDC );
            mUnlikePiPTvsRhoPT->Fill( lv1.Pt(), lv.Pt() );
            mCos2phivsPTvsMass->Fill( lv.M(), lv.Pt(), 2*cos( PairPhi ));

            if ( lv.M() > 0.65 && lv.M() <0.9 ){

                if ( absPperp < 0.06){

                    mPairPhi->Fill ( PairPhi );
                    mLabPairPhi->Fill ( lab_calc_Phi(lv1, lv2) );
                    mPhivsRapidity->Fill ( PairPhi, mRapidity);
                    mPhivsZDC->Fill ( PairPhi, TotalZDC );
                    mPhivsEastZDC->Fill ( PairPhi, EastZDC );
                    mPhivsWestZDC->Fill ( PairPhi, WestZDC );

                    mCos2phivsRapidity->Fill( 2*cos(2*PairPhi), lv.Rapidity() );
                    mCos4phivsRapidity->Fill( 2*cos(4*PairPhi), lv.Rapidity() );

                    mCos2phivsEastZDC->Fill ( 2*cos(2*PairPhi), EastZDC );
                    mCos2phivsWestZDC->Fill ( 2*cos(2*PairPhi), WestZDC );
                    mCos2phivsTotalZDC->Fill ( 2*cos(2*PairPhi), TotalZDC );

                }

                if ( WestZDC > 30 && WestZDC < 70 && EastZDC > 30 && EastZDC < 70 ){mPhi1n1n->Fill(PairPhi);}
                if ( WestZDC > 30 && WestZDC < 70 && EastZDC > 100 && EastZDC < 135 ){mPhi2n1n->Fill(PairPhi);}
                if ( WestZDC > 100 && WestZDC < 135 && EastZDC > 30 && EastZDC < 70 ){mPhi2n1n->Fill(PairPhi);}
                if ( WestZDC > 100 && WestZDC < 135 && EastZDC > 100 && EastZDC < 130 ){mPhi2n2n->Fill(PairPhi);}
                if ( TotalZDC < 150 ) { mLowZDCcos2phivsPT->Fill( 2*cos(2*PairPhi), absPperp); }
                if ( TotalZDC > 150 ) { mHighZDCcos2phivsPT->Fill( 2*cos(2*PairPhi), absPperp); }

                if ( TotalZDC < 150 ) { mZDCEastvsWest->Fill(EastZDC, WestZDC); }

                mPxvsPy->Fill( absPperp*cos(PairPhi), absPperp*sin(PairPhi));
                mPhivsPT->Fill ( PairPhi, absPperp);
                mCos1phivsPT->Fill( 2*cos(PairPhi), absPperp );
                mCos2phivsPT->Fill( 2*cos(2 *PairPhi), absPperp );
                mLabCos2phivsPT->Fill( 2*cos(2 *(lab_calc_Phi(lv1, lv2))), absPperp );
                mCos3phivsPT->Fill( 2*cos(3 * PairPhi), absPperp);
                mCos4phivsPT->Fill( 2*cos(4 * PairPhi), absPperp );
            }

            if ( lv.M() > 0.2 && lv.M() <1.5 && absPperp < 0.06) { 

                mPhivsMass->Fill ( PairPhi, lv.M());
                mCos2phivsMass->Fill( 2*cos(2 *(PairPhi)), lv.M() );
                mCos4phivsMass->Fill( 2*cos(4 *(PairPhi)), lv.M() ); 

            }
        }
    }
mPairPhi->Scale(2*3.1415/(mPairPhi->Integral("width")));
mPairPhi->Fit(mPhiFit);

mPhi1n1n->Scale(2*3.1415/(mPhi1n1n->Integral("width")));
mPhi2n1n->Scale(2*3.1415/(mPhi2n1n->Integral("width")));
mPhi2n2n->Scale(2*3.1415/(mPhi2n2n->Integral("width")));

auto * mPTcos1phimoments = mCos1phivsPT->ProfileY("mPTcos1phimoments", 1, -1);
auto * mPTcos2phimoments = mCos2phivsPT->ProfileY("mPTcos2phimoments", 1, -1);
auto * mlabPTcos2phimoments = mLabCos2phivsPT->ProfileY("mlabPTcos2phimoments", 1, -1);
auto * mPTcos3phimoments = mCos3phivsPT->ProfileY("mPTcos3phimoments", 1, -1);
auto * mPTcos4phimoments = mCos4phivsPT->ProfileY("mPTcos4phimoments", 1, -1);

auto * mLikePTcos1phimoments = mLikeCos1phivsPT->ProfileY("mLikePTcos1phimoments", 1, -1);
auto * mLikePTcos2phimoments = mLikeCos2phivsPT->ProfileY("mLikePTcos2phimoments", 1, -1);
auto * mLikePTcos3phimoments = mLikeCos3phivsPT->ProfileY("mLikePTcos3phimoments", 1, -1);
auto * mLikePTcos4phimoments = mLikeCos4phivsPT->ProfileY("mLikePTcos4phimoments", 1, -1);

auto * mMasscos2phimoments = mCos2phivsMass->ProfileY("mMasscos2phimoments", 1, -1);
auto * mMasscos4phimoments = mCos4phivsMass->ProfileY("mMasscos4phimoments", 1, -1);

auto * mLikeMasscos2phimoments = mLikeCos2phivsMass->ProfileY("mLikeMasscos2phimoments", 1, -1);

auto * mRapiditycos2phimoments = mCos2phivsRapidity->ProfileY("mRapiditycos2phimoments", 1, -1);
auto * mRapiditycos4phimoments = mCos4phivsRapidity->ProfileY("mRapiditycos4phimoments", 1, -1);

auto * mLikeRapiditycos2phimoments = mLikeCos2phivsRapidity->ProfileY("mLikeRapiditycos2phimoments", 1, -1);

auto * mEastZDCcos2phimoments = mCos2phivsEastZDC->ProfileY("mEastZDCcos2phimoments", 1, -1);
auto * mWestZDCcos2phimoments = mCos2phivsWestZDC->ProfileY("mWestZDCcos2phimoments", 1, -1);
auto * mTotalZDCcos2phimoments = mCos2phivsTotalZDC->ProfileY("mTotalZDCcos2phimoments", 1, -1);
auto * mLowZDCPTcos2phimoments = mLowZDCcos2phivsPT->ProfileY("mLowZDCPTcos2phimoments", 1, -1);
auto * mHighZDCPTcos2phimoments = mHighZDCcos2phivsPT->ProfileY("mHighZDCPTcos2phimoments", 1, -1);

auto * mPTvsMasscos2phimoments = mCos2phivsPTvsMass->Project3DProfile("xy");

mRapiditycos2phimoments->Fit("pol0","", "", -0.8, 0.8);

mTotalZDCcos2phimoments->Fit("pol0", "","", 600, 1600);

//auto * mCorrectedPTcos2phimoments = new TH1F("mCorrectedPTcos2phimoments", "Corrected  (0.65 < M_{#pi#pi} < 0.9 GeV); P_{T} (GeV/c); 2<cos2#phi>", 25, 0, 0.25);
for (int i = 0; i < (mCorrectedPTcos2phimoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mPTcos2phimoments->GetBinContent( i+1 )/2;
    double omega_2 = mLikePTcos2phimoments->GetBinContent( mLikePTcos2phimoments->FindBin( mPTcos2phimoments->GetBinCenter( i+1 ) ) )/2;
    double delta_gamma_2 =  mPTcos2phimoments->GetBinError( i+1 );
    double delta_omega_2 =  mLikePTcos2phimoments->GetBinError( mLikePTcos2phimoments->FindBin( mPTcos2phimoments->GetBinCenter( i+1 ) ) );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * gamma_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedPTcos2phimoments->SetBinContent( i+1, 2*alpha_2 );
    mCorrectedPTcos2phimoments->SetBinError( i+1, delta_alpha_2);
}

auto * mCorrectedMasscos2phimoments = new TH1F("mCorrectedMasscos2phimoments", "Corrected  (P_{T} < 0.06 GeV/c); Mass (GeV); 2<cos2#phi>", 25, 0.25, 1);
for (int i = 0; i < (mCorrectedMasscos2phimoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mMasscos2phimoments->GetBinContent( i+1 )/2;
    double omega_2 = mLikeMasscos2phimoments->GetBinContent( mLikeMasscos2phimoments->FindBin( mMasscos2phimoments->GetBinCenter( i+1 ) ) )/2;
    double delta_gamma_2 =  mMasscos2phimoments->GetBinError( i+1 );
    double delta_omega_2 =  mLikeMasscos2phimoments->GetBinError( mLikeMasscos2phimoments->FindBin( mMasscos2phimoments->GetBinCenter( i+1 ) ) );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * gamma_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedMasscos2phimoments->SetBinContent( i+1, 2*alpha_2 );
    mCorrectedMasscos2phimoments->SetBinError( i+1, delta_alpha_2);
}

auto * mCorrectedRapiditycos2phimoments = new TH1F("mCorrectedRapiditycos2phimoments", "Corrected  (0.6 < M_{#pi#pi} < 0.8 GeV, P_{T} < 0.06 MeV); Rapidity; 2<cos2#phi>", 25, -2, 2);
for (int i = 0; i < (mCorrectedRapiditycos2phimoments->GetNbinsX()) - 1; i++) {
    double gamma_2 = mRapiditycos2phimoments->GetBinContent( i+1 )/2;
    double omega_2 = mLikeRapiditycos2phimoments->GetBinContent( mLikeRapiditycos2phimoments->FindBin( mRapiditycos2phimoments->GetBinCenter( i+1 ) ) )/2;
    double delta_gamma_2 =  mRapiditycos2phimoments->GetBinError( i+1 );
    double delta_omega_2 =  mLikeRapiditycos2phimoments->GetBinError( mLikeRapiditycos2phimoments->FindBin( mRapiditycos2phimoments->GetBinCenter( i+1 ) ) );
    double dadw = (2*(gamma_2 * gamma_2) - 4) / pow( (omega_2*gamma_2 - 2), 2);
    double dadg = (-2*(omega_2 * omega_2) + 4) / pow( (omega_2*gamma_2 - 2), 2);

    double alpha_2 = ( 2 * (omega_2 - gamma_2) ) / ( (omega_2 * gamma_2) -2 );
    double delta_alpha_2 = sqrt( (dadw*dadw*delta_omega_2*delta_omega_2) + (dadg*dadg*delta_gamma_2*delta_gamma_2) );

    mCorrectedRapiditycos2phimoments->SetBinContent( i+1, 2*alpha_2 );
    mCorrectedRapiditycos2phimoments->SetBinError( i+1, delta_alpha_2);
}


//make plots
fo -> cd();

makeCanvas();
gStyle->SetPalette(1);
mChiChi->Draw("colz");
gPad->SetLogz();
gPad->Print( "plots/plot_mChiChi.pdf" );
gPad->Print( "plots/plot_mChiChi.png" );

makeCanvas();
mddTOF->SetLineColor(kBlack);
gPad->SetLogy();
mddTOF->Draw();
gPad->Print( "plots/data/plot_mddTOF.pdf" );
gPad->Print( "plots/data/plot_mddTOF.png" );

makeCanvas();
mDCA1->SetLineColor(kBlack);
mDCA1->Draw();
gPad->Print( "plots/data/plot_mDCA1.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mDCA1v2->Draw("colz");
gPad->Print( "plots/data/plot_mDCA1v2.pdf" );

makeCanvas();
gStyle->SetPalette(1);
mZDCEastvsWest->Draw("colz");
gPad->Print( "plots/data/ZDC/plot_mZDCEastvsWest.pdf" );
gPad->Print( "plots/data/ZDC/plot_mZDCEastvsWest.png" );

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->SetMaximum(22500);
mMass->Draw();
TLine *l1 = new TLine(0.65,0,0.65,22500);
l1->SetLineColor(kBlue);
l1->SetLineWidth(2);
l1->Draw("Same");
TLine *l2 = new TLine(0.9,0,0.9,22500);
l2->SetLineColor(kRed-2);
l2->SetLineWidth(2);
l2->Draw("Same");
auto masslegend = new TLegend(0.65,0.1,0.95,0.4);
masslegend->SetHeader("Legend","C"); // option "C" allows to center the header
masslegend->AddEntry(l1,"Lower mass cut (650 MeV)");
masslegend->AddEntry(l2,"Upper mass cut (900 MeV)");
masslegend->Draw();
gPad->Print( "plots/data/Mass/plot_mMass.pdf" );

makeCanvas();
mPperp->SetLineColor(kBlack);
mPperp->SetMaximum(10000);
mPperp->Draw();
TLine *l3 = new TLine(0.06,0,0.06,10000);
l3->SetLineColor(kRed-2);
l3->SetLineWidth(2);
l3->Draw("Same");
auto ptlegend = new TLegend(0.65,0.1,0.95,0.4);
ptlegend->SetHeader("Legend","C"); // option "C" allows to center the header
ptlegend->AddEntry(l3,"Upper P_{T} cut (60 MeV)");
ptlegend->Draw();
gPad->Print( "plots/data/PT/plot_mPperp.pdf" );

makeCanvas();
mPairPhi->SetLineColor(kBlack);
gStyle->SetOptFit();
mPairPhi->Draw();
gPad->Print( "plots/data/plot_mPairPhi.pdf" );
gPad->Print( "plots/data/plot_mPairPhi.png" );

/*
makeCanvas();
mEta->SetLineColor(kBlack);
mLVPhi->Draw();

makeCanvas();
mMass->SetLineColor(kBlack);
mMass->Draw();

makeCanvas();
gStyle->SetPalette(1);
mUnlikePiPTvsRhoPT->Draw("colz");
gPad->Print( "plots/data/PT/plot_mUnlikePiPTvsRhoPT.pdf" );
gPad->Print( "plots/data/PT/plot_mUnlikePiPTvsRhoPT.png" );

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
gStyle->SetOptFit();
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

makeCanvas();
mPhi1n1n->SetLineColor(kBlack);
mPhi1n1n->SetTitle("#phi at various ZDC peaks; #phi (rad); counts");
mPhi2n1n->SetLineColor(kOrange);
mPhi2n2n->SetLineColor(kBlue);
mPhi1n1n->Draw();
mPhi2n1n->Draw("SAME");
mPhi2n2n->Draw("SAME");
auto legend = new TLegend(0.65,0.1,0.95,0.4);
legend->AddEntry(mPhi1n1n,"1n1n peaks");
legend->AddEntry(mPhi2n1n,"2n1n peaks");
legend->AddEntry(mPhi2n2n,"2n2n peak");
legend->Draw();
gPad->Print( "plots/data/ZDC/plot_mPhiNeutronPeaks.pdf" );
gPad->Print( "plots/data/ZDC/plot_mPhiNeutronPeaks.png" );

makeCanvas();
mPTcos2phimoments->SetTitle("cos(2#phi) moments vs. P_{T}; P_{T} (GeV); 2<cos(2#phi)>");
mPTcos2phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mPTmomentsplot.pdf" );
gPad->Print( "plots/data/PT/plot_mPTmomentsplot.png" );

makeCanvas();
mMassmomentsplot->SetTitle("cos(2#phi) moments vs. Mass; Mass (GeV); 2<cos(2#phi)>");
mMassmomentsplot->SetMinimum(-0.5);
mMassmomentsplot->SetMaximum(1.0);
mMassmomentsplot->Draw();
gPad->Print( "plots/data/Mass/plot_mMassmomentsplot.pdf" );
gPad->Print( "plots/data/Mass/plot_mMassmomentsplot.png" );

makeCanvas();
mRapiditymomentsplot->SetTitle("cos(2#phi) moments vs. Rapidity; Rapidity; 2<cos(2#phi)>");
mRapiditymomentsplot->Draw();
gPad->Print( "plots/data/Rapidity/plot_mRapiditymomentsplot.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mRapiditymomentsplot.png" );
mRapiditymomentsplot->SetMinimum(0);
mRapiditymomentsplot->SetMaximum(0.5);
gStyle->SetOptFit(1111);
mRapiditymomentsplot->Draw();
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsplot.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mZoomedRapiditymomentsplot.png" );

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
*/
makeCanvas();
mLowZDCPTcos2phimoments->SetLineColor(kBlack);
mHighZDCPTcos2phimoments->SetLineColor(kMagenta);
mLowZDCPTcos2phimoments->SetTitle("cos(2#phi) moments vs. Pair P_{T} (U+U); Pair P_{T} (GeV/c); 2<cos2#phi>");
mLowZDCPTcos2phimoments->Draw();
mHighZDCPTcos2phimoments->Draw("SAME");
auto legend3 = new TLegend(0.65,0.1,0.95,0.4);
legend3->AddEntry(mLowZDCPTcos2phimoments,"ZDC (East^{2} + West^{2})^{1/2} < 150");
legend3->AddEntry(mHighZDCPTcos2phimoments,"ZDC (East^{2} + West^{2})^{1/2} > 150");
legend3->Draw();
gPad->Print( "plots/data/ZDC/plot_mDiffZDCPTcos2phimoments.pdf" );
gPad->Print( "plots/data/ZDC/plot_mDiffZDCPTcos2phimoments.png" );

makeCanvas();
mPTcos2phimoments->SetLineColor(kBlack);
mlabPTcos2phimoments->SetLineColor(kMagenta);
mPTcos2phimoments->SetTitle("Strength of cos2#phi modulation; P_{T} (GeV/c); 2<cos2#phi>");
mPTcos2phimoments->SetMinimum(-0.5);
mPTcos2phimoments->Draw();
mlabPTcos2phimoments->Draw("SAME");
auto legend4 = new TLegend(0.65,0.1,0.95,0.4);
legend4->AddEntry(mPTcos2phimoments,"Rest frame 2<cos2#phi>");
legend4->AddEntry(mlabPTcos2phimoments,"Lab frame 2<cos2#phi>");
legend4->Draw();
gPad->Print( "plots/data/PT/plot_mLabVRestCos2phivsPT.pdf" );
gPad->Print( "plots/data/PT/plot_mLabVRestCos2phivsPT.png" );

makeCanvas();
gStyle->SetPalette(1);
mPTvsMasscos2phimoments->SetTitle("; pair P_{T} (GeV/c); pair Mass (GeV); 2<cos2#phi>");
mPTvsMasscos2phimoments->Draw("colz");
gPad->Print( "plots/data/plot_mPTvsMasscos2phimoments.pdf" );
gPad->Print( "plots/data/plot_mPTvsMasscos2phimoments.png" );

makeCanvas();
mPTcos1phimoments->SetTitle("cos(1#phi) moments vs. Pair P_{T} (U+U); Pair P_{T} (GeV/c); 2<cos(1#phi)>");
mPTcos1phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mPTcos1phimoments.pdf" );

makeCanvas();
mPTcos2phimoments->SetTitle("cos(2#phi) moments vs. Pair P_{T} (U+U); Pair P_{T} (GeV/c); 2<cos(2#phi)>");
mPTcos2phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mPTcos2phimoments.pdf" );

makeCanvas();
mPTcos3phimoments->SetTitle("cos(3#phi) moments vs. Pair P_{T} (U+U); Pair P_{T} (GeV/c); 2<cos(3#phi)>");
mPTcos3phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mPTcos3phimoments.pdf" );

makeCanvas();
mPTcos4phimoments->SetTitle("cos(4#phi) moments vs. Pair P_{T} (U+U); Pair P_{T} (GeV/c); 2<cos(4#phi)>");
mPTcos4phimoments->Draw();
gPad->Print( "plots/data/PT/plot_mPTcos4phimoments.pdf" );

makeCanvas();
mMasscos2phimoments->SetTitle("cos(2#phi) moments vs. Mass (U+U); Mass (GeV); 2<cos(2#phi)>");
mMasscos2phimoments->SetMinimum(0);
mMasscos2phimoments->SetMaximum(1.2);
mMasscos2phimoments->Draw();
gPad->Print( "plots/data/Mass/plot_mMasscos2phimoments.pdf" );
gPad->Print( "plots/data/Mass/plot_mMasscos2phimoments.png" );

makeCanvas();
mRapiditycos2phimoments->SetTitle("cos(2#phi) moments vs. Rapidity (U+U); Rapidity; 2<cos(2#phi)>");
mRapiditycos2phimoments->SetMinimum(-0.5);
mRapiditycos2phimoments->SetMaximum(0.5);
gStyle->SetOptFit();
mRapiditycos2phimoments->Draw();
gPad->Print( "plots/data/Rapidity/plot_mRapiditycos2phimoments.pdf" );
gPad->Print( "plots/data/Rapidity/plot_mRapiditycos2phimoments.png" );

makeCanvas();
mLikePTcos1phimoments->SetLineColor(kBlack);
mPTcos1phimoments->SetLineColor(kGreen+2);
//mCorrectedPTcos1phimoments->SetLineColor(kMagenta);
mLikePTcos1phimoments->SetTitle( "2<cos1#phi> (0.65 < M_{#pi#pi} < 0.9 GeV) (U+U); Pair P_{T} (GeV/c); 2<cos#phi>" );
mLikePTcos1phimoments->SetLineWidth(2);
mPTcos1phimoments->SetLineWidth(2);
//mCorrectedPTcos1phimoments->SetLineWidth(2);
mLikePTcos1phimoments->Draw();
mPTcos1phimoments->Draw("SAME");
//mCorrectedPTcos1phimoments->Draw("SAME");
auto legend6 = new TLegend(0.65,0.1,0.95,0.4);
legend6->SetHeader("Legend","C"); // option "C" allows to center the header
legend6->AddEntry(mLikePTcos1phimoments,"Like sign pairs");
legend6->AddEntry(mPTcos1phimoments,"Unike sign pairs");
//legend6->AddEntry(mCorrectedPTcos1phimoments,"Baseline subtracted signal");
legend6->Draw();
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos1PhivsPTmoments.png" );
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos1PhivsPTmoments.pdf" );

makeCanvas();
mLikePTcos2phimoments->SetLineColor(kBlack);
mPTcos2phimoments->SetLineColor(kGreen+2);
mCorrectedPTcos2phimoments->SetLineColor(kMagenta);
mLikePTcos2phimoments->SetTitle( "2<cos2#phi> (0.65 < M_{#pi#pi} < 0.9 GeV) (U+U); Pair P_{T} (GeV/c); 2<cos2#phi>" );
mLikePTcos2phimoments->SetMaximum(0.5);
mLikePTcos2phimoments->SetLineWidth(2);
mPTcos2phimoments->SetLineWidth(2);
mCorrectedPTcos2phimoments->SetLineWidth(2);
mLikePTcos2phimoments->Draw();
mPTcos2phimoments->Draw("SAME");
mCorrectedPTcos2phimoments->Draw("SAME");
auto legend5 = new TLegend(0.65,0.1,0.95,0.4);
legend5->SetHeader("Legend","C"); // option "C" allows to center the header
legend5->AddEntry(mLikePTcos2phimoments,"Like sign pairs");
legend5->AddEntry(mPTcos2phimoments,"Unike sign pairs");
legend5->AddEntry(mCorrectedPTcos2phimoments,"Baseline subtracted signal");
legend5->Draw();
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos2PhivsPTmoments.png" );
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos2PhivsPTmoments.pdf" );

makeCanvas();
mLikePTcos3phimoments->SetLineColor(kBlack);
mPTcos3phimoments->SetLineColor(kGreen+2);
//mCorrectedPTcos3phimoments->SetLineColor(kMagenta);
mLikePTcos3phimoments->SetTitle( "2<cos3#phi> (0.65 < M_{#pi#pi} < 0.9 GeV) (U+U); Pair P_{T} (GeV/c); 2<cos3#phi>" );
mLikePTcos3phimoments->SetLineWidth(2);
mPTcos3phimoments->SetLineWidth(2);
//mCorrectedPTcos1phimoments->SetLineWidth(2);
mLikePTcos3phimoments->Draw();
mPTcos3phimoments->Draw("SAME");
//mCorrectedPTcos1phimoments->Draw("SAME");
auto legend7 = new TLegend(0.65,0.1,0.95,0.4);
legend7->SetHeader("Legend","C"); // option "C" allows to center the header
legend7->AddEntry(mLikePTcos3phimoments,"Like sign pairs");
legend7->AddEntry(mPTcos3phimoments,"Unike sign pairs");
//legend7->AddEntry(mCorrectedPTcos1phimoments,"Baseline subtracted signal");
legend7->Draw();
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos3PhivsPTmoments.png" );
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos3PhivsPTmoments.pdf" );

makeCanvas();
mLikePTcos4phimoments->SetLineColor(kBlack);
mPTcos4phimoments->SetLineColor(kGreen+2);
//mCorrectedPTcos4phimoments->SetLineColor(kMagenta);
mLikePTcos4phimoments->SetTitle( "2<cos4#phi> (0.65 < M_{#pi#pi} < 0.9 GeV) (U+U); Pair P_{T} (GeV/c); 2<cos4#phi>" );
mLikePTcos4phimoments->SetMaximum(0.15);
mLikePTcos4phimoments->SetMinimum(-0.15);
mLikePTcos4phimoments->SetLineWidth(2);
mPTcos4phimoments->SetLineWidth(2);
//mCorrectedPTcos4phimoments->SetLineWidth(2);
mLikePTcos4phimoments->Draw();
mPTcos4phimoments->Draw("SAME");
//mCorrectedPTcos1phimoments->Draw("SAME");
auto legend8 = new TLegend(0.65,0.1,0.95,0.4);
legend8->SetHeader("Legend","C"); // option "C" allows to center the header
legend8->AddEntry(mLikePTcos4phimoments,"Like sign pairs");
legend8->AddEntry(mPTcos4phimoments,"Unike sign pairs");
//legend8->AddEntry(mCorrectedPTcos4phimoments,"Baseline subtracted signal");
legend8->Draw();
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos4PhivsPTmoments.png" );
gPad->Print( "plots/data/PT/plot_mLikeCorrectionCos4PhivsPTmoments.pdf" );

makeCanvas();
mLikeMasscos2phimoments->SetLineColor(kBlack);
mMasscos2phimoments->SetLineColor(kGreen+2);
mCorrectedMasscos2phimoments->SetLineColor(kMagenta);
mLikeMasscos2phimoments->SetTitle( "Pair P_{T} < 0.06 GeV/c (U+U); M_{#pi#pi} (GeV); 2<cos2#phi>" );
mLikeMasscos2phimoments->SetMaximum(1);
mLikeMasscos2phimoments->SetLineWidth(2);
mMasscos2phimoments->SetLineWidth(2);
mCorrectedMasscos2phimoments->SetLineWidth(2);
mLikeMasscos2phimoments->Draw();
mMasscos2phimoments->Draw("SAME");
mCorrectedMasscos2phimoments->Draw("SAME");
auto legend9 = new TLegend(0.65,0.1,0.95,0.4);
legend9->SetHeader("Legend","C"); // option "C" allows to center the header
legend9->AddEntry(mLikeMasscos2phimoments,"Like sign pairs");
legend9->AddEntry(mMasscos2phimoments,"Unike sign pairs");
legend9->AddEntry(mCorrectedMasscos2phimoments,"Baseline subtracted signal");
legend9->Draw();
gPad->Print( "plots/data/Mass/plot_mLikeCorrectionCos2PhivsMassmoments.png" );
gPad->Print( "plots/data/Mass/plot_mLikeCorrectionCos2PhivsMassmoments.pdf" );

makeCanvas();
mLikeRapiditycos2phimoments->SetLineColor(kBlack);
mRapiditycos2phimoments->SetLineColor(kGreen+2);
//mCorrectedRapiditycos2phimoments->SetLineColor(kMagenta);
mLikeRapiditycos2phimoments->SetTitle( "0.65 < M_{#pi#pi} < 0.9 GeV, Pair P_{T} < 0.06 GeV/c (U+U); Rapidity; 2<cos2#phi>" );
mLikeRapiditycos2phimoments->SetMaximum(0.5);
mLikeRapiditycos2phimoments->SetLineWidth(2);
mRapiditycos2phimoments->SetLineWidth(2);
//mCorrectedRapiditycos2phimoments->SetLineWidth(2);
gStyle->SetOptFit();
mLikeRapiditycos2phimoments->Draw();
mRapiditycos2phimoments->Draw("SAME");
//mCorrectedRapiditycos2phimoments->Draw("SAME");
auto legend10 = new TLegend(0.65,0.1,0.95,0.4);
legend10->SetHeader("Legend","C"); // option "C" allows to center the header
legend10->AddEntry(mLikeRapiditycos2phimoments,"Like sign pairs");
legend10->AddEntry(mRapiditycos2phimoments,"Unike sign pairs");
//legend10->AddEntry(mCorrectedRapiditycos2phimoments,"Baseline subtracted signal");
legend10->Draw();
gPad->Print( "plots/data/Rapidity/plot_mLikeCorrectionCos2PhivsRapiditymoments.png" );
gPad->Print( "plots/data/Rapidity/plot_mLikeCorrectionCos2PhivsRapiditymoments.pdf" );

fo->Write();
}
