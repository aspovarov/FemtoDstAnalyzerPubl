/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading STAR FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \authors: Grigory Nigmatkulov, Povarov Alexey, Demanov Alexandr
 * \date May 29, 2018
 */

// This is needed for calling standalone classes (not needed on RACF)
#define _VANILLA_ROOT_

// C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream> 

// ROOT headers
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"

// FemtoDst headers

#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoDstReader.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoDst.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoEvent.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoTrack.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoV0.h"
#include "/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/StFemtoXi.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/libStFemtoDst)
#endif

const Int_t n = 8;    // eta-gap by TPC + BBC and ZDC

// Forward declarations
// Check event and track
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZ, Double_t VtxR, Double_t delta_VtxY);
Bool_t isGoodTrack(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_EVENT);
Bool_t isGoodTrackFlow(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_FLOW);

// Calculate Q-vector by TPC
Bool_t CalculateTPCQVec( TVector2 Q2vec[][n], TVector2 Q3vec[][n], StFemtoDst *dst,Int_t ngap, Double_t DCA_EVENT);

// Calculate Q-vector by ZDC
Float_t ZDCSMD(Int_t eastwest, Int_t verthori, Int_t strip, StFemtoEvent *event );
Float_t ZDCSMD_GetPosition( Int_t eastwest, Int_t verthori, Int_t strip);
TVector2 CalculateZDCQVec(StFemtoEvent *event, Int_t ew);

// Calculate Q-vector by BBC
TVector2 CalculateBBCQVec(StFemtoEvent *event, Int_t ew, Float_t harm);
Float_t GetBBCTilePhi(const Int_t e_w, const Int_t iTile);

// Recentering for TPC, BBC, ZDC
void Recentering( TVector2 Q1vec[][2], TVector2 Q2vec[][n], TVector2 Q3vec[][n], Int_t RunID, Int_t cent,TProfile2D *tpQx1[][2],
TProfile2D *tpQy1[][2],TProfile2D *tpQx2[][n],TProfile2D *tpQy2[][n],TProfile2D *tpQx3[][n],TProfile2D *tpQy3[][n] );

const Float_t electron_mass = 0.0005485799;
const Float_t pion_mass = 0.13957061;
const Float_t kaon_mass = 0.493677;
const Float_t proton_mass = 0.9382720813;

const Float_t electron_mass_sqr = 0.000000301;
const Float_t pion_mass_sqr = 0.019479955;
const Float_t kaon_mass_sqr = 0.24371698;
const Float_t proton_mass_sqr = 0.880354499;

Double_t mZDCSMDCenterex = 4.72466;
Double_t mZDCSMDCenterey = 5.53629;
Double_t mZDCSMDCenterwx = 4.39604;
Double_t mZDCSMDCenterwy = 5.19968;



// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.FemtoDst.root files
//_________________
void FemtoDstAnalyzer(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                           const Char_t *oFileName = "oTest.root",
                           const Char_t *mode = "QAmode",
                           const Char_t *energy = "39GeV") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  
  //QA//0 - QA mode
  //raw//1 - raw mode
  //rec//2 - recentering mode
  //flow//3 - flattening and flow mode
  //Int_t mode = 1;

  //false - RunQA OFF
  //true - RunQA ON
  Bool_t useRunQA = false;

  Double_t VtxZ;
  Double_t VtxR;
  Double_t DCA_EVENT;
  Double_t DCA_FLOW;
  Double_t DCA_SYS;
  Double_t delta_VtxY;
  Double_t nSigm;
  
  Double_t sys_VtxZ = 45.0;
  Int_t sys_nHits = 18;
  Double_t sys_DCA = 1.5;

  Int_t runIdBins;// = 20000;
  Int_t runIdRange[2];// = { 11095000, 11115000 };

  if(strncmp(energy, "39GeV",5)==0){
    runIdRange[0]=11095000;
    runIdRange[1]=11115000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=40.0;
    VtxR=2.0;
    delta_VtxY=0.0;
    DCA_EVENT=2.0;
    DCA_FLOW=1.0;

    useRunQA = true;
    DCA_SYS=2.0;
    Double_t sys_VtxZ=35.0;
    Int_t sys_nHits=18;
    Double_t sys_DCA=1.5;
    nSigm = 3.0;

  }
  if(strncmp(energy, "27GeV",5)==0){
    runIdRange[0]=12171000;
    runIdRange[1]=12180000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=70.0;
    VtxR=2.0;
    delta_VtxY=0.0;
    DCA_EVENT=2.0;
    DCA_FLOW=1.0;
    
    DCA_SYS=2.0;
    Double_t sys_VtxZ = 45.0;
    Int_t sys_nHits = 18;
    Double_t sys_DCA = 1.5;
    useRunQA = true;
    nSigm = 1.5;

  }
  if(strncmp(energy, "19GeV",5)==0){
    runIdRange[0]=12110000;
    runIdRange[1]=12123000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=70.0;
    VtxR=2.0;
    delta_VtxY=0.0;
    DCA_EVENT=2.0;
    DCA_FLOW=1.0;
    
    DCA_SYS=2.0;
    Double_t sys_VtxZ = 45.0;
    Int_t sys_nHits = 18;
    Double_t sys_DCA = 1.5;
    useRunQA = true;
    nSigm = 3.0;

  }
  if(strncmp(energy, "14GeV",5)==0){
    runIdRange[0]=15045000;
    runIdRange[1]=15075000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=70.0;
    VtxR=1.0;
    delta_VtxY=0.8847;
    DCA_EVENT=2.0;
    DCA_FLOW=3.0;
    
    DCA_SYS=3.0;
    Double_t sys_VtxZ = 65.0;
    Int_t sys_nHits = 18;
    Double_t sys_DCA = 2.5;
    useRunQA = true;
    nSigm = 3.0;

  }
  if(strncmp(energy, "11GeV",5)==0){
    runIdRange[0]=11145000;
    runIdRange[1]=11165000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=50.0;
    VtxR=2.0;
    delta_VtxY=0.0;
    DCA_EVENT=2.0;
    DCA_FLOW=1.0;
    
    DCA_SYS=2.0;
    Double_t sys_VtxZ = 45.0;
    Int_t sys_nHits = 18;
    Double_t sys_DCA = 1.5;
    useRunQA = true;
    nSigm = 3.0;

  }
  if(strncmp(energy, "7GeV",4)==0){
    runIdRange[0]=11110000;
    runIdRange[1]=11150000;
    runIdBins=runIdRange[1]-runIdRange[0];
    VtxZ=50.0;
    VtxR=2.0;
    delta_VtxY=0.0;
    DCA_EVENT=2.0;
    DCA_FLOW=1.0;
    
    DCA_SYS=2.0;
    Double_t sys_VtxZ = 45.0;
    Int_t sys_nHits = 18;
    Double_t sys_DCA = 1.5;
  }

  std::cout<< DCA_SYS<<"\t"<<DCA_FLOW<<"\t"<<DCA_EVENT<<"\t"<< sys_VtxZ<<"\t"<<nSigm<<std::endl;

  Int_t eventCounter = 0;
  Int_t hundredIter = 0;

  Bool_t badEvent = true; 
  Int_t c = 9;    // cent9()
  Int_t cent, RunID;
  Double_t w = 0.0;
  Double_t Phi, pt, q, SqM;

  Double_t etagap[] ={0.05, 0.075, 0.1, 0.15, 0.2, 0.5};

  const Char_t *direction[] = {"east","west"};
  const Char_t *detector[] = {"TPC","BBC","ZDC"};
  const Char_t *ngap[] = {"Eta01","Eta15","Eta02","Eta03","Eta04","Eta1","BBC","ZDC"};
  const Char_t *gap[] = {"gap 0.1","gap 0.15","gap 0.2","gap 0.3","gap 0.4","gap 1.0","BBC","ZDC"};
  const Char_t *sign[] = {"Pos","Neg"};
  const Char_t *particles[] = {"pion","kaon","proton","hadrons"};
  const Char_t *partLateX[] = {"#pi^{+}","pi^{-}","K^{+}","K^{-}","p","#bar{p}"};
  const Char_t *systematic[] = {"VtxZ", "nHits", "DCA"};
  const Char_t *systematic_text[] = {"45", "18", "1.5"};

  gSystem->Load("./libStFemtoDst.so");
  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("./libStFemtoDst.so");
  #endif

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  //Long64_t events2read = femtoReader->chain()->GetEntries();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInt_tree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInt_tree: "  << eventsInt_tree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  //Reaction plane and flow

  //histogram for Q-vectors and event planes with eta-gap and without error
  TH1D *h_Qx1[2][2][c], *h_Qy1[2][2][c]; 
  TH1D *h_Qx2[2][n][c], *h_Qy2[2][n][c], *h_Qx3[2][n][c], *h_Qy3[2][n][c];
  TH1D *h_Psi1[2][2][c],*h_Psi2[2][n][c], *h_Psi3[2][n][c];
  //histogram for check resolution
  TH1D *h_sinPsi2westTPCE[n][c], *h_cosPsi2westTPCE[n][c], *h_sinPsi3westTPCE[n][c], *h_cosPsi3westTPCE[n][c];
    
    // loop by direction
    for(Int_t l = 0; l < 2; l++) {
      //loop by cent                                                                                                  
      for(Int_t j = 0; j < c; j++) {
        //loop by eta-gap + BBC and ZDC
        for(Int_t i = 0; i < n; i++) {

          //historam of Q-vectors for east eta-gap 
          h_Qx2[l][i][j] = new TH1D(Form("h_Qx2%s%scent%i",direction[l],ngap[i],j),
              Form("Q_{x} for #psi_{2} %s %s cent %i;Q_{x}",direction[l],gap[i],j),300,-1.5,1.5);
          h_Qy2[l][i][j] = new TH1D(Form("h_Qy2%s%scent%i",direction[l],ngap[i],j),
              Form("Q_{y} for #psi_{2} %s %s cent %i;Q_{y}",direction[l],gap[i],j),300,-1.5,1.5);
          h_Qx3[l][i][j] = new TH1D(Form("h_Qx3%s%scent%i",direction[l],ngap[i],j),
              Form("Q_{x} for #psi_{3} %s %s cent %i;Q_{x}",direction[l],gap[i],j),300,-1.5,1.5);
          h_Qy3[l][i][j] = new TH1D(Form("h_Qy3%s%scent%i",direction[l],ngap[i],j),
              Form("Q_{y} for #psi_{3} %s %s cent %i;Q_{y}",direction[l],gap[i],j),300,-1.5,1.5);
                
          //historam of event planes for east eta-gap 
          h_Psi2[l][i][j] = new TH1D(Form("h_Psi2%s%scent%i",direction[l],ngap[i],j),
              Form("#psi_{2} %s %s cent %i;#psi_{2}",direction[l],gap[i],j),100,-0.05,3.2); 
          h_Psi3[l][i][j] = new TH1D(Form("h_Psi3%s%scent%i",direction[l],ngap[i],j),
              Form("#psi_{3} %s %s cent %i;#psi_{3}",direction[l],gap[i],j),100,-0.05,2.15);
        }// loop by n

        for(Int_t det = 0; det < 2; det++) {
          h_Qx1[l][det][j] = new TH1D(Form("h_Qx1%s%scent%i",direction[l],ngap[det+6],j),
            Form("Q_{x} for #psi_{1} %s %s cent %i;Q_{x}",direction[l],gap[det+6],j),300,-1.5,1.5);
          h_Qy1[l][det][j] = new TH1D(Form("h_Qy1%s%scent%i",direction[l],ngap[det+6],j),
            Form("Q_{y} for #psi_{1} %s %s cent %i;Q_{y}",direction[l],gap[det+6],j),300,-1.5,1.5);

          h_Psi1[l][det][j] = new TH1D(Form("h_Psi1%s%scent%i",direction[l],ngap[det+6],j),
              Form("#psi_{1} %s %s cent %i;#psi_{1}",direction[l],gap[det+6],j),200,-0.1,6.35);
        }

      }// loop by cent
    }// loop by direstion
  

  
  // RunQA
  std::vector<Int_t> BadRuns;
  if(useRunQA == true) { 
    Int_t buff = 0;
    BadRuns.reserve(1000);
    ifstream BadRunList(Form("/mnt/pool/1/aspovarov/basov/test/Bad_Run_%s.txt", energy));
    while( !BadRunList.eof() ) {
        BadRunList >> buff;
        BadRuns.push_back( buff );
    }
    std::cout << "BadRunLists yep" << std::endl;
  }//RunQA
  
  
  //QA mode /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( strncmp(mode, "QA",2)==0 ){
    TFile *oFile = new TFile(oFileName, "RECREATE");


    // Histogramming
    // Event
    TH1D *hRefMult = new TH1D("hRefMult", "Reference multiplicity;RefMult;Entries",
                              600, -0.5, 599.5);
    TH2D *hVtxXvsY = new TH2D("hVtxXvsY", "hVtxXvsY;x (cm);y (cm)",
                              200,-10.,10.,200,-10.,10.);
    TH1D *hVtxZ = new TH1D("hVtxZ","hVtxZ;z (cm); Entries",
                           140, -70., 70.);
    TH1D *hRefMult2 = new TH1D("hRefMult2","Reference multiplicity in |#eta|<1;RefMult2;Entries",
                               600, -0.5, 599.5);
    TH1D *hGRefMult = new TH1D("hGRefMult","Reference multiplicity of global tracks;gRefMult;Entries",
                               800, -0.5, 799.5);
    TH1D *hNumberOfPrimaries = new TH1D("hNumberOfPrimaries","Number of primary tracks;Number of primary tracks;Entries",
                                        600, -0.5, 599.5);
    TH1D *hNumberOfGlobals = new TH1D("hNumberOfGlobals","Number of global tracks;Number of global tracks;Entries",
                                      600, -0.5, 599.5);
    TH1D *hCent9 = new TH1D("hCent9","Centralitity;Cent9;Entries",
                            13, -1.5, 11.5);
    TH1D *hCent16 = new TH1D("hCent16","Centralitity;Cent16;Entries",
                            19, -1.5, 17.5);
    TH1D *hBTofHit = new TH1D("hBTofHit","Number of hits in TOF;bTofTrayMult;Entries",
                              600, -0.5, 599.5);
    TH1D *hBTofMatched = new TH1D("hBTofMatched","Number of TOF-matched tracks;bTofMatched;Entries",
                                  400, -0.5, 399.5);
    TH1D *hBemcMatched = new TH1D("hBemcMatched","Number of BEMC-matched tracks;bEmcMatched;Entries",
                                  400, -0.5, 399.5);
    TH1D *hRanking = new TH1D("hRanking","Primary vertex ranking;Primary vertex ranking;Entries",
                              21, -10.5, 10.5);
    TH2D *hVpdVzDiffVsVz = new TH2D("hVpdVzDiffVsVz","v_{z}(TPC) - v_{z}(VPD) vs. v_{z}(TPC);v_{z}(TPC);v_{z}(TPC) - v_{z}(VPD)",
                                    280, -70., 70., 80, -20., 20.);
    TH2D *hBTofTrayMultVsRefMult = new TH2D("hBTofTrayMultVsRefMult","TOF tray multiplicity vs. refMult;refMult;bTofTrayMult",
                                            600, -0.5, 599.5, 600, -0.5, 599.5);
    TH2D *hBTofMatchedVsRefMult = new TH2D("hBTofMatchedVsRefMult","TOF-matched tracks vs. refMult;refMult;TOF-matched",
                                            600, -0.5, 599.5, 400, -0.5, 399.5);
    TH1D *hTransSphericity = new TH1D("hTransSphericity","Transverse sphericity;Sphericity;Entries",
                                      10, 0., 1.);
    TH1D *hTransSphericity2 = new TH1D("hTransSphericity2","Transverse sphericity in |#eta|<1;Sphericity;Entries",
                                       10, 0., 1.);
    TH1D *hNumberOfVertices = new TH1D("hNumberOfVertices","Number of primary vertices;Number of primary vertices;Entries",
                                       15, -0.5, 14.5);
    TProfile *hEventProfile[8];
    hEventProfile[0] = new TProfile("hEventProfile_0","Profile of refMult;Run ID;<refMult>",
                                    runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[1] = new TProfile("hEventProfile_1","Profile of TOF tray multiplicity;Run ID;<bTofTrayMultiplicity>",
                                    runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[2] = new TProfile("hEventProfile_2","Profile of TOF-matched tracks;Run ID;<bTofMatched>",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[3] = new TProfile("hEventProfile_3","Profile of number of primary tracks;Run ID;<nPrimTracks>",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[4] = new TProfile("hEventProfile_4","Profile of number of global tracks;Run ID;<nGlobTracks>",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[5] = new TProfile("hEventProfile_5","Profile of ZDC ADC;Run ID; <ADC_{ZDC}>",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[6] = new TProfile("hEventProfile_6","Profile of BBC ADC;Run ID; <ADC_{BBC}>",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    hEventProfile[7] = new TProfile("hEventProfile_7","Profile of primary vertex Z position;Run ID; <VtxZ> (cm)",
                                              runIdBins, runIdRange[0], runIdRange[1] );
    for(Int_t iHist=0; iHist<8; iHist++) {
      Int_t mColor = 1;
      hEventProfile[iHist]->SetMarkerStyle(20);    // filled circle
      hEventProfile[iHist]->SetMarkerColor(mColor);
      hEventProfile[iHist]->SetMarkerSize(1.1);
      hEventProfile[iHist]->SetLineWidth(2);
      hEventProfile[iHist]->SetLineColor(mColor);  // black
    }

    // Track
    TH1D *hGlobalPtot = new TH1D("hGlobalPtot","Global track momentum;p (GeV/c);Entries",
                                 200, 0., 2. );
    TH1D *hPrimaryPtot = new TH1D("hPrimaryPtot","Primary track momentum;p (GeV/c);Entries",
                                  200, 0., 2. );
    TH1D *hGlobalPt = new TH1D("hGlobalPt","Global track transverse momentum;p_{T} (GeV/c)",
                                200, 0., 2.);
    TH1D *hPrimaryPt = new TH1D("hPrimaryPt","Primary track transverse momentum;p_{T} (GeV/c)",
                                200, 0., 2.);
    TH1D *hNHits = new TH1D("hNHits","Number of hits;nHits;Entries", 80, -0.5, 79.5);
    TH1D *hNHitsRatio = new TH1D("hNHitsRatio","nHitsFit to nHitsPoss ratio;nHitsFit/nHitsRatio;Entries",
                                 10, 0., 1. );
    TH1D *hChi2 = new TH1D("hChi2","#chi^{2} of the track;#chi^{2};Entries",
                           200, 0., 20.);
    TH1D *hDca = new TH1D("hDca","DCA to primary vertex;DCA (cm);Entries",
                          100, 0., 10.);
    TH2D *hDcaVsPt = new TH2D("hDcaVsPt","charge*p_{T} vs. DCA;charge * p_{T} (GeV/c);DCA (cm)",
                              840, -2.1, 2.1, 100, 0., 10.);
    TH1D *hPhi = new TH1D("hPhi","Azimuthal angle distribution;#phi;Entries",
                          640, -3.2, 3.2 );
    TH1D *hEta = new TH1D("hEta","Track pseudorapidity;#eta;Entries", 220, -1.1, 1.1 );
    TH1D *hEtaG = new TH1D("hEtaG","Track pseudorapidity of global track;#eta;Entires", 220, -1., 1. );
    TH2D *hPtVsEta = new TH2D("hPtVsEta","p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)",
                              220, -1.1, 1.1, 80, 0.05, 2.05);
    TH2D *hPrimaryPhiVsPt[2];
    for(Int_t i=0; i<2; i++) {
      hPrimaryPhiVsPt[i] = new TH2D(Form("hPrimaryPhiVsPt_%d",i),
           Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
           300, 0., 3., 630, -3.15, 3.15 );
    }
    TH1D* hDedx = new TH1D("hDedx","dE/dx;dE/dx (keV/cm);Entries",
                           125, 0., 12.5);
    TH2D *hDedxVsPt = new TH2D("hDedxVsPt", "dE/dx vs. charge*p_{T};charge * p_{T} (GeV/c);dE/dx (keV/cm)",
                               840, -2.1, 2.1, 600, 0., 12.);
    TH2D *hNSigmaPionVsPt = new TH2D("hNSigmaPionVsPt","n#sigma(#pi) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(#pi)",
                                     840, -2.1, 2.1, 200, -10., 10.);
    TH2D *hNSigmaElectronVsPt = new TH2D("hNSigmaElectronVsPt","n#sigma(e) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(e)",
                                         840, -2.1, 2.1, 200, -10., 10.);
    TH2D *hNSigmaKaonVsPt = new TH2D("hNSigmaKaonVsPt","n#sigma(K) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(K)",
                                     840, -2.1, 2.1, 200, -10., 10.);
    TH2D *hNSigmaProtonVsPt = new TH2D("hNSigmaProtonVsPt","n#sigma(p) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(p)",
                                       840, -2.1, 2.1, 200, -10., 10.);

    TH2D *hDedxVsPtPID[4];
    for ( Int_t i=0; i<4; i++ ) {
      TString name = "hDedxVsPtPID_";
      name += i;
      TString title = "dE/dx vs. charge*p_{T} ";
      switch (i) {
        case 0: title += "|n#sigma(e)| #leq 2;"; break;
        case 1: title += "|n#sigma(#pi)| #leq 2;"; break;
        case 2: title += "|n#sigma(K)| #leq 2;"; break;
        case 3: title += "|n#sigma(p)| #leq 2;"; break;
        default: title += "unknown PID;";
      }
      title += "charge*p_{T} (GeV/c);dE/dx (keV/cm)";
      hDedxVsPtPID[i] = new TH2D(name.Data(), title.Data(),
                                 840, -2.1, 2.1, 600, 0., 12.);
    }

    // TofPidTrait
    TH1D *hTofBeta = new TH1D("hTofBeta","BTofPidTraits #beta;#beta",
                              2000, 0., 2.);
    TH2D *hInvBetaVsPt = new TH2D("hInvBetaVsPt","1/#beta vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta",
                                  840, -2.1, 2.1, 200, 0.8, 2.8);
    TH1D *hMassSqr = new TH1D("hMassSqr","m^{2};m^{2} (GeV/c^{2})^{2};dN/dm^{2} (entries)",
                              520, -0.1, 5.1 );
    TH2D *hMassSqrVsPt = new TH2D("hMassSqrVsPt","m^{2} vs. charge*p_{T};charge * p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",
                                  840, -2.1, 2.1, 200, -0.2, 1.8);
    TH2D *hDedxVsMassSqr[2];
    hDedxVsMassSqr[0] = new TH2D("hDedxVsMassSqr_0","dE/dx vs. mass^{2} charge>0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
               440, -0.4, 1.8, 250, 0., 12.5 );
    hDedxVsMassSqr[1] = new TH2D("hDedxVsMassSqr_1","dE/dx vs. mass^{2} charge<0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
               440, -0.4, 1.8, 250, 0., 12.5 );
    TH2D *hInvBetaDiffElectronVsPt = new TH2D("hInvBetaDiffElectronVsPt","1/#beta - 1/#beta(electron) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(e)",
                                              840, -2.1, 2.1, 200, -0.1, 0.1);
    TH2D *hInvBetaDiffPionVsPt = new TH2D("hInvBetaDiffPionVsPt","1/#beta - 1/#beta(pion) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(#pi)",
                                              840, -2.1, 2.1, 200, -0.1, 0.1);
    TH2D *hInvBetaDiffKaonVsPt = new TH2D("hInvBetaDiffKaonVsPt","1/#beta - 1/#beta(kaon) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(K)",
                                              840, -2.1, 2.1, 200, -0.1, 0.1);
    TH2D *hInvBetaDiffProtonVsPt = new TH2D("hInvBetaDiffProtonVsPt","1/#beta - 1/#beta(p) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(p)",
                                            840, -2.1, 2.1, 200, -0.1, 0.1);
    TProfile *hTrackProfile[6];
    hTrackProfile[0] = new TProfile("hTrackProfile_0","Profile of track #phi;Run ID;<#phi>",
                                             runIdBins, runIdRange[0], runIdRange[1] );
    hTrackProfile[1] = new TProfile("hTrackProfile_1","Profile of track p_{T};Run ID;<p_{T}>",
                                             runIdBins, runIdRange[0], runIdRange[1] );
    hTrackProfile[2] = new TProfile("hTrackProfile_2","Profile of track nHits;Run ID;<nHits>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
    hTrackProfile[3] = new TProfile("hTrackProfile_3","Profile of track DCA;Run ID;<DCA>",
                                            runIdBins, runIdRange[0], runIdRange[1] );
    hTrackProfile[4] = new TProfile("hTrackProfile_4","Profile of track #beta;Run ID;<#beta>",
                                               runIdBins, runIdRange[0], runIdRange[1] );
    hTrackProfile[5] = new TProfile("hTrackProfile_5","Profile of track dE/dx;Run ID;<dE/dx> (keV/cm)",
                                              runIdBins, runIdRange[0], runIdRange[1] );

    for(Int_t iTrk=0; iTrk<6; iTrk++) {
      Int_t mColor = 1;
      hTrackProfile[iTrk]->SetMarkerStyle(20);    // filled circle
      hTrackProfile[iTrk]->SetMarkerColor(mColor);
      hTrackProfile[iTrk]->SetMarkerSize(1.1);
      hTrackProfile[iTrk]->SetLineWidth(2);
      hTrackProfile[iTrk]->SetLineColor(mColor);  // black
    }
  
    // Loop over events
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

      eventCounter++;
      if( eventCounter >= 100000 ) {
        eventCounter = 0;
        hundredIter++;
        std::cout << "Working on event #[" << (hundredIter * 100000)
          	      << "/" << events2read << "]" << std::endl;
      }


      Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
      if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
        break;
      }

      // Retrieve femtoDst
      StFemtoDst *dst = femtoReader->femtoDst();

      // Retrieve event information
      StFemtoEvent *event = dst->event();
      if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
        break;
      }

      TVector3 pVtx = event->primaryVertex();

      // Simple event cut
      if ( !isGoodEvent( event, VtxZ, VtxR, delta_VtxY ) ) continue;
      
      //RunQA
      if( useRunQA == true && std::find(BadRuns.begin(), BadRuns.end(), event -> runId()) != BadRuns.end() ) continue; 
      

      // Fill event histograms
      hRefMult->Fill( event->refMult() );
      hVtxXvsY->Fill( event->primaryVertex().X(), event->primaryVertex().Y() );
      hVtxZ->Fill( event->primaryVertex().Z() );
      hRefMult2->Fill( event->refMult2() );
      hGRefMult->Fill( event->gRefMult() );
      hNumberOfPrimaries->Fill( event->numberOfPrimaryTracks() );
      hNumberOfGlobals->Fill( event->numberOfGlobalTracks() );
      hCent9->Fill( event->cent9() );
      hCent16->Fill( event->cent16() );
      hBTofHit->Fill( event->numberOfBTofHit() );
      hBTofMatched->Fill( event->numberOfTofMatched() );
      //hBemcMatched->Fill( event->numberOfBEMCMatched() );
      hRanking->Fill( event->ranking() );
      hVpdVzDiffVsVz->Fill( event->primaryVertex().Z(),
                            event->primaryVertex().Z() - event->vpdVz() );
      hBTofTrayMultVsRefMult->Fill( event->refMult(),
                                    event->numberOfBTofHit() );
      hBTofMatchedVsRefMult->Fill( event->refMult(),
                                  event->numberOfTofMatched() );
      hTransSphericity->Fill( event->transverseSphericity() );
      hTransSphericity2->Fill( event->transverseSphericity2() );
      hNumberOfVertices->Fill( event->numberOfPrimaryVertices() );
      hEventProfile[0]->Fill( event->runId(), event->gRefMult() );
      hEventProfile[1]->Fill( event->runId(), event->numberOfBTofHit() );
      hEventProfile[2]->Fill( event->runId(), event->numberOfTofMatched() );
      hEventProfile[3]->Fill( event->runId(), event->numberOfPrimaryTracks() );
      hEventProfile[4]->Fill( event->runId(), event->numberOfGlobalTracks() );
      hEventProfile[5]->Fill( event->runId(), event->zdcSumAdcEast() + event->zdcSumAdcWest() );
      Float_t bbcAdcSum = 0;
      for( Int_t iTile=0; iTile<24; iTile++ ) {
        bbcAdcSum += event->bbcAdcEast(iTile);
        bbcAdcSum += event->bbcAdcWest(iTile);
      }
      hEventProfile[6]->Fill( event->runId(), bbcAdcSum );
      hEventProfile[7]->Fill( event->runId(), pVtx.Z() );

      // Track analysis
      Int_t nTracks = dst->numberOfTracks();

      // Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

        // Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst->track(iTrk);

        if (!femtoTrack) continue;
        //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;

        if ( !femtoTrack->isPrimary() ) continue;

        // Simple single-track cut
        if( !isGoodTrack(event, femtoTrack, DCA_EVENT) ) continue;

        // Fill track histograms
        hGlobalPtot->Fill( femtoTrack->gMom().Mag() );
        hPrimaryPtot->Fill( femtoTrack->pMom().Mag() );
        hGlobalPt->Fill( femtoTrack->gMom().Pt() );
        hPrimaryPt->Fill( femtoTrack->pMom().Pt() );
        hNHits->Fill( femtoTrack->nHits() );
        hNHitsRatio->Fill( (Float_t)femtoTrack->nHitsFit()/femtoTrack->nHitsPoss() );
        hChi2->Fill( femtoTrack->chi2() );
        hDca->Fill( femtoTrack->gDCA( pVtx.X(), pVtx.Y(), pVtx.Z() ) );
        hDcaVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                        femtoTrack->gDCA( pVtx.X(), pVtx.Y(), pVtx.Z() ) );
        hPtVsEta->Fill( femtoTrack->pMom().Eta(), femtoTrack->pMom().Pt() );
        hPhi->Fill( femtoTrack->pMom().Phi() );
        hEta->Fill( femtoTrack->pMom().Eta() );
        hEtaG->Fill( femtoTrack->gMom().Eta() );
        hDedx->Fill( femtoTrack->dEdx() * 1e6 );
        hDedxVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                         femtoTrack->dEdx() * 1e6 );

        // If electron has passed PID nsigma cut
        if ( TMath::Abs( femtoTrack->nSigmaElectron() ) <= 2. ) {
          hDedxVsPtPID[0]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                 femtoTrack->dEdx() * 1e6 );
        }

        // If pion has passed PID nsigma cut
        if ( TMath::Abs( femtoTrack->nSigmaPion() ) <= 2. ) {
          hDedxVsPtPID[1]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                 femtoTrack->dEdx() * 1e6 );
        }

        // If kaon has passed PID nsigma cut
        if ( TMath::Abs( femtoTrack->nSigmaKaon() ) <= 2. ) {
          hDedxVsPtPID[2]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                 femtoTrack->dEdx() * 1e6 );
        }

        // If proton has passed PID nsigma cut
        if ( TMath::Abs( femtoTrack->nSigmaProton() ) <= 2. ) {
          hDedxVsPtPID[3]->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                 femtoTrack->dEdx() * 1e6 );
        }

        if( femtoTrack->charge() > 0 ) {
        	hPrimaryPhiVsPt[0]->Fill( femtoTrack->pMom().Pt(),
                                    femtoTrack->pMom().Phi() );
        }
        else {
  	       hPrimaryPhiVsPt[1]->Fill( femtoTrack->pMom().Pt(),
  				                           femtoTrack->pMom().Phi() );
        }
        hNSigmaElectronVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                   femtoTrack->nSigmaElectron() );
        hNSigmaPionVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->nSigmaPion() );
        hNSigmaKaonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                               femtoTrack->nSigmaKaon() );
        hNSigmaProtonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                 femtoTrack->nSigmaProton() );
        hTrackProfile[0]->Fill( event->runId(), femtoTrack->pMom().Phi() );
        hTrackProfile[1]->Fill( event->runId(), femtoTrack->pMom().Pt() );
        hTrackProfile[2]->Fill( event->runId(), femtoTrack->nHits() );
        hTrackProfile[3]->Fill( event->runId(),
                                femtoTrack->gDCA(pVtx.X(), pVtx.Y(), pVtx.Z() ) );
        hTrackProfile[5]->Fill( event->runId(),
                                femtoTrack->dEdx() * 1e6 );

        // Check if track has TOF signal
        if ( femtoTrack->isTofTrack() ) {
  	      hTofBeta->Fill( femtoTrack->beta() );
          hInvBetaVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                              femtoTrack->invBeta() );
          hMassSqr->Fill( femtoTrack->massSqr() );
          hMassSqrVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                              femtoTrack->massSqr() );
          hInvBetaDiffElectronVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                          femtoTrack->invBeta() - TMath::Sqrt(electron_mass_sqr +
                                          femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
          hInvBetaDiffPionVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                      femtoTrack->invBeta() - TMath::Sqrt(pion_mass_sqr +
                                          femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
          hInvBetaDiffKaonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                      femtoTrack->invBeta() - TMath::Sqrt(kaon_mass_sqr +
                                          femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
          hInvBetaDiffProtonVsPt->Fill( femtoTrack->charge() * femtoTrack->pMom().Pt(),
                                        femtoTrack->invBeta() - TMath::Sqrt(proton_mass_sqr +
                                          femtoTrack->pMom().Mag2() )/ femtoTrack->pMom().Mag() );
          if ( femtoTrack->charge() > 0 ) {
            hDedxVsMassSqr[0]->Fill( femtoTrack->massSqr(), femtoTrack->dEdx() * 1e6 );
          }
          else {
            hDedxVsMassSqr[1]->Fill( femtoTrack->massSqr(), femtoTrack->dEdx() * 1e6 );
          }
          hTrackProfile[4]->Fill( event->runId(), femtoTrack->beta() );
        } //if( isTofTrack() )

      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
    oFile->Write();
    oFile->Close();
  }//QA mode
  
  //raw mode /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( strncmp(mode, "raw",3)==0 ){ 
    TFile *oFile = new TFile(oFileName, "RECREATE");

    //TProfile2D for recentering
    TProfile2D *tp_Qx1[2][2], *tp_Qy1[2][2];
    TProfile2D *tp_Qx2[2][n], *tp_Qx3[2][n], *tp_Qy2[2][n], *tp_Qy3[2][n];

    // loop by direction 
    for(Int_t l = 0; l < 2; l++) {
      // loop by eta-gap of TPC + BBC and ZDC
      for(Int_t i = 0; i < n; i++) {
      
        //TProfile2D for recentering 
        tp_Qx2[l][i] = new TProfile2D(Form("tp_Qx2%s%s",direction[l],ngap[i]),
        Form("<Q_{x}> for #psi_{2} %s %s;RunID;cent",direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

        tp_Qy2[l][i] = new TProfile2D(Form("tp_Qy2%s%s",direction[l],ngap[i]),
        Form("<Q_{y}> for #psi_{2} %s %s;RunID;cent",direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

        tp_Qx3[l][i] = new TProfile2D(Form("tp_Qx3%s%s",direction[l],ngap[i]),
        Form("<Q_{x}> for #psi_{3} %s %s;RunID;cent",direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

        tp_Qy3[l][i] = new TProfile2D(Form("tp_Qy3%s%s",direction[l],ngap[i]),
        Form("<Q_{y}> for #psi_{3} %s %s;RunID;cent",direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

      }// loop by eta-gap of TPC + BBC and ZDC

      for(Int_t det = 0; det < 2; det++) {

        tp_Qx1[l][det] = new TProfile2D(Form("tp_Qx1%s%s",direction[l],ngap[det+6]),
        Form("<Q_{x}> for #psi_{1} %s %s;RunID;cent",direction[l],gap[det+6]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

        tp_Qy1[l][det] = new TProfile2D(Form("tp_Qy1%s%s",direction[l],ngap[det+6]),
        Form("<Q_{y}> for #psi_{1} %s %s;RunID;cent",direction[l],gap[det+6]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);
      }
    }// loop by direction

    //loop by event
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

      if ( iEvent % 10000 == 0) {
        std::cout << "Working on event #[" << (iEvent+1)
     	      	    << "/" << events2read << "]" << std::endl;
      }

      if (iEvent == events2read-1) {
        std::cout << "Working on event #[" << (events2read)
     	 	          << "/" << events2read << "]" << std::endl;
      }

    	Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    	if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      	break;
    	}

    	// Retrieve femtoDst
    	StFemtoDst *dst = femtoReader->femtoDst();

      // Retrieve event information
    	StFemtoEvent *event = dst->event();
    	if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      	break;
      }

    	// Simple event cut
    	if( !isGoodEvent( event, VtxZ, VtxR, delta_VtxY ) ) continue;

    	// Track analysis
    	Int_t nTracks = dst->numberOfTracks();
            
      //RunQA
      if( useRunQA == true && std::find(BadRuns.begin(), BadRuns.end(), event -> runId()) != BadRuns.end() ) continue; 

      TVector2 Q1vec[2][2];
      TVector2 Q2vec[2][8], Q3vec[2][8];

      for(Int_t l = 0; l < 2; l++) { 
        for(Int_t det = 0; det < 2; det++) {
          Q1vec[l][det].Set(0.,0.);
        }
        for(Int_t i = 0; i < n-1; i++) {
          Q2vec[l][i].Set(0.,0.);
          Q3vec[l][i].Set(0.,0.);
        }
      }

      for(Int_t dir = 0; dir < 2; dir++) {
        Q1vec[dir][0] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[dir][1] = CalculateZDCQVec(event, dir);
      }
      for(Int_t dir = 0; dir < 2; dir++) {
        Q2vec[dir][6] = CalculateBBCQVec(event, dir, 2.0);
        Q3vec[dir][6] = CalculateBBCQVec(event, dir, 3.0);
      }
      badEvent = CalculateTPCQVec(Q2vec, Q3vec, dst, n-2, DCA_EVENT);
            
    	cent = event -> cent9();
    	RunID = event -> runId();

      for(Int_t l = 0; l < 2; l ++) {
        for(Int_t i = 0; i < n-1; i++) {  

          // Fill only TPC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i < 6) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

        	  h_Psi2[l][i][cent] -> Fill( 1.0/2.0 * Q2vec[l][i].Phi()  );
        	  h_Psi3[l][i][cent] -> Fill( 1.0/3.0 * Q3vec[l][i].Phi()  ); 

   		      tp_Qx2[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vec[l][i].X() );
   		      tp_Qy2[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vec[l][i].Y() );
   		      tp_Qx3[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vec[l][i].X() );
   		      tp_Qy3[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vec[l][i].Y() );
          }
          // Fill BBc and ZDC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i > 5) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

            h_Psi2[l][i][cent] -> Fill( 1.0/2.0 * Q2vec[l][i].Phi()  );
            h_Psi3[l][i][cent] -> Fill( 1.0/3.0 * Q3vec[l][i].Phi()  ); 

            tp_Qx2[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vec[l][i].X() );
            tp_Qy2[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vec[l][i].Y() );
            tp_Qx3[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vec[l][i].X() );
            tp_Qy3[l][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vec[l][i].Y() );
          }	      
        } //for(Int_t i = 0; i < n-1; i++)  

        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[l][det].Mod() != 0. ) {
            h_Qx1[l][det][cent] -> Fill( Q1vec[l][det].X() );
            h_Qy1[l][det][cent] -> Fill( Q1vec[l][det].Y() );

            h_Psi1[l][det][cent] -> Fill( Q1vec[l][det].Phi() - TMath::Pi() );

            tp_Qx1[l][det] -> Fill( (Double_t)RunID, (Double_t)cent, Q1vec[l][det].X() );
            tp_Qy1[l][det] -> Fill( (Double_t)RunID, (Double_t)cent, Q1vec[l][det].Y() );
          }

        }

      } //for(Int_t l = 0; l < 2; l ++) 
    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    for(Int_t l = 0; l < 2; l++) {
      for(Int_t det = 0; det < 2; det++) {
        for(Int_t j = 0; j < c; j++) { 
          h_Qx1[l][det][j] -> Write();
          h_Qy1[l][det][j] -> Write();
          h_Psi1[l][det][j] -> Write();
        }
        tp_Qx1[l][det] -> Write();
        tp_Qy1[l][det] -> Write();
      }
      for(Int_t i = 0; i < n; i++) {
        for(Int_t j = 0; j < c; j++ ) {

          h_Qx2[l][i][j] -> Write();
          h_Qy2[l][i][j] -> Write();
          h_Qx3[l][i][j] -> Write();
          h_Qy3[l][i][j] -> Write();

          h_Psi2[l][i][j] -> Write();
          h_Psi3[l][i][j] -> Write();
        }
        tp_Qx2[l][i] -> Write();
        tp_Qy2[l][i] -> Write();
        tp_Qx3[l][i] -> Write();
        tp_Qy3[l][i] -> Write();
      }
    }
  	oFile->Close();
  }//raw mode 
      
	// Recentering mode ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( strncmp(mode, "rec",3)==0 ){ 

    TProfile2D *tp_sinPsi1[2][2][20], *tp_cosPsi1[2][2][20];
    TProfile2D *tp_sinPsi2[2][n][4], *tp_cosPsi2[2][n][4], *tp_sinPsi3[2][n][4], *tp_cosPsi3[2][n][4];
    for(Int_t l = 0; l < 2; l++) {
      for(Int_t i = 0; i < n; i++) {
        for(Int_t j = 0; j < 4; j++) {
          //TProfile2D for flattening 
          tp_sinPsi2[l][i][j] = new TProfile2D(Form("tp_%isinPsi2%s%s",j+1,direction[l],ngap[i]),
          Form("<sin(%i*2#psi_{2})> %s %s",j+1,direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

          tp_cosPsi2[l][i][j] = new TProfile2D(Form("tp_%icosPsi2%s%s",j+1,direction[l],ngap[i]),
          Form("<cos(%i*2#psi_{2})> %s %s",j+1,direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

          tp_sinPsi3[l][i][j] = new TProfile2D(Form("tp_%isinPsi3%s%s",j+1,direction[l],ngap[i]),
          Form("<sin(%i*3#psi_{3})> %s %s",j+1,direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

          tp_cosPsi3[l][i][j] = new TProfile2D(Form("tp_%icosPsi3%s%s",j+1,direction[l],ngap[i]),
          Form("<cos(%i*3#psi_{3})> %s %s",j+1,direction[l],gap[i]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

        }
      }
      for(Int_t det = 0; det < 2; det++) {
        for(Int_t j = 0; j < 20; j++) {
          //TProfile2D for flattening 
          tp_sinPsi1[l][det][j] = new TProfile2D(Form("tp_%isinPsi1%s%s",j+1,direction[l],ngap[det+6]),
          Form("<sin(%i*#psi_{1})> %s %s",j+1,direction[l],gap[det+6]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

          tp_cosPsi1[l][det][j] = new TProfile2D(Form("tp_%icosPsi1%s%s",j+1,direction[l],ngap[det+6]),
          Form("<cos(%i*#psi_{1})> %s %s",j+1,direction[l],gap[det+6]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);
        }
      }
    }

    TFile *f1 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/NoRe_%s.root",energy),"READ");
    TProfile2D *tpQx1[2][2], *tpQy1[2][2];
    TProfile2D *tpQx2[2][n], *tpQy2[2][n], *tpQx3[2][n], *tpQy3[2][n];
    for(Int_t l = 0; l < 2; l++) {
      for(Int_t i = 0; i < n; i++) {
        tpQx2[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx2%s%s",direction[l],ngap[i]) );
        tpQx3[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx3%s%s",direction[l],ngap[i]) );
        tpQy2[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy2%s%s",direction[l],ngap[i]) );
        tpQy3[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy3%s%s",direction[l],ngap[i]) );
      }
      for(Int_t det = 0; det < 2; det++) {
        tpQx1[l][det] = (TProfile2D*) f1 -> Get( Form("tp_Qx1%s%s",direction[l],ngap[det+6]) );
        tpQy1[l][det] = (TProfile2D*) f1 -> Get( Form("tp_Qy1%s%s",direction[l],ngap[det+6]) );
      }
    }

    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
      if ( iEvent % 10000 == 0) {
        std::cout << "Working on event #[" << (iEvent+1)
     	      		  << "/" << events2read << "]" << std::endl;
	    }
	    if (iEvent == events2read-1) {
        std::cout << "Working on event #[" << (events2read)
     	            << "/" << events2read << "]" << std::endl;
	    }

    	Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    	if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
        break;
    	}

    	// Retrieve femtoDst
    	StFemtoDst *dst = femtoReader->femtoDst();

    	// Retrieve event information
    	StFemtoEvent *event = dst->event();
    	if( !event ) {
      	std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      	break;
    	}

    	// Simple event cut
    	if ( !isGoodEvent( event , VtxZ, VtxR, delta_VtxY) ) continue;

      // Track analysis
    	Int_t nTracks = dst->numberOfTracks();

      //RunQA
      if( useRunQA == true && std::find(BadRuns.begin(), BadRuns.end(), event -> runId()) != BadRuns.end() ) continue; 

      TVector2 Q1vec[2][2];
      TVector2 Q2vec[2][8], Q3vec[2][8];

      for(Int_t l = 0; l < 2; l++) { 
        for(Int_t det = 0; det < 2; det++) {
          Q1vec[l][det].Set(0.,0.);
        }
        for(Int_t i = 0; i < n-1; i++) {
          Q2vec[l][i].Set(0.,0.);
          Q3vec[l][i].Set(0.,0.);
        }
      }

      for(Int_t dir = 0; dir < 2; dir++) {
        Q1vec[dir][0] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[dir][1] = CalculateZDCQVec(event, dir);
      }

      for(Int_t l = 0; l < 2; l++) {
        Q2vec[l][6] = CalculateBBCQVec(event, l, 2.0);
        Q3vec[l][6] = CalculateBBCQVec(event, l, 3.0);
      }
      badEvent = CalculateTPCQVec(Q2vec, Q3vec, dst, n-2, DCA_EVENT);
             
      cent = (Double_t)(event -> cent9());
      RunID = event -> runId();

      Recentering(Q1vec,Q2vec,Q3vec,RunID,cent,tpQx1,tpQy1,tpQx2,tpQy2,tpQx3,tpQy3);   

      for(Int_t l = 0; l < 2; l ++) { 
        for(Int_t i = 0; i < n-1; i++) {  

          // Fill only TPC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i < 6) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

            h_Psi2[l][i][cent] -> Fill( 1.0/2.0 * Q2vec[l][i].Phi()  );
            h_Psi3[l][i][cent] -> Fill( 1.0/3.0 * Q3vec[l][i].Phi()  );

            for(Int_t j =0; j < 4; j++) {
              tp_sinPsi2[l][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q2vec[l][i].Phi() ) );
              tp_cosPsi2[l][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q2vec[l][i].Phi() ) );
              tp_sinPsi3[l][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q3vec[l][i].Phi() ) ); 
              tp_cosPsi3[l][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q3vec[l][i].Phi() ) );
            }
          }
          // Fill BBc and ZDC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i > 5) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

            h_Psi2[l][i][cent] -> Fill( 1.0/2.0 * Q2vec[l][i].Phi()  );
            h_Psi3[l][i][cent] -> Fill( 1.0/3.0 * Q3vec[l][i].Phi()  );

            for(Int_t j =0; j < 4; j++) { 
              tp_sinPsi2[l][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q2vec[l][i].Phi() ) );
              tp_cosPsi2[l][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q2vec[l][i].Phi() ) ); 
              tp_sinPsi3[l][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q3vec[l][i].Phi() ) );
              tp_cosPsi3[l][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q3vec[l][i].Phi() ) );
            }
          }           
        } //for(Int_t i = 0; i < n-1; i++)  

        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[l][det].Mod() != 0. ) {
            h_Qx1[l][det][cent] -> Fill( Q1vec[l][det].X() );
            h_Qy1[l][det][cent] -> Fill( Q1vec[l][det].Y() );

            h_Psi1[l][det][cent] -> Fill( Q1vec[l][det].Phi() );

            for(Int_t j = 0; j < 20; j++) {
              tp_sinPsi1[l][det][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q1vec[l][det].Phi() ) );
              tp_cosPsi1[l][det][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q1vec[l][det].Phi() ) ); 
            }
          }
        }

      } //for(Int_t l = 0; l < 2; l ++)
    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    f1 -> Close();
    TFile *oFile = new TFile(oFileName, "RECREATE");
    for(Int_t l = 0; l < 2; l++) {
      for(Int_t det = 0; det < 2; det++) {
        for(Int_t j = 0; j < c; j++) { 
          h_Qx1[l][det][j] -> Write();
          h_Qy1[l][det][j] -> Write();
          h_Psi1[l][det][j] -> Write();
        }
        for(Int_t j = 0; j < 20; j++) { 
          tp_sinPsi1[l][det][j] -> Write();
          tp_cosPsi1[l][det][j] -> Write();
        }
      }
      for(Int_t i = 0; i < n; i++) {
        for(Int_t j = 0; j < c; j++ ) {

          h_Qx2[l][i][j] -> Write();
          h_Qy2[l][i][j] -> Write();
          h_Qx3[l][i][j] -> Write();
          h_Qy3[l][i][j] -> Write();

          h_Psi2[l][i][j] -> Write();
          h_Psi3[l][i][j] -> Write();
        }
        for(Int_t j = 0; j < 4; j++) {
          tp_sinPsi2[l][i][j] -> Write();
          tp_cosPsi2[l][i][j] -> Write();
          tp_sinPsi3[l][i][j] -> Write();
          tp_cosPsi3[l][i][j] -> Write();
        }
      }
    }
    oFile->Close();
  }// Recentering mode 

  // Flattening and flow mode ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  if( strncmp(mode, "flow",4)==0 ){ 
  
    Double_t SqMdown[] = {-0.15, 0.2, 0.74};
    Double_t SqMup[] = {0.1, 0.32, 1.20};
    Double_t nSigma[3];
    

    //TProfile for resolution for Psi1
    TProfile *tp_SqRes1[2];
    for(Int_t det = 0; det < 2; det++) {
      tp_SqRes1[det] = new TProfile(Form("tp_SqRes1%s",ngap[det+6]),Form("Resolution^{2} for v_{1} %s",gap[det+6]),9,0,9);
    }

    //TProfile for resolution for Psi2 and Psi3 
    TProfile *tp_SqRes2[n], *tp_SqRes3[n];
    for(Int_t i = 0; i < n; i++) {
      tp_SqRes2[i] = new TProfile(Form("tp_SqRes2%s",ngap[i]),Form("Resolution^{2} for v_{2} %s",gap[i]),9,0,9);
      tp_SqRes3[i] = new TProfile(Form("tp_SqRes3%s",ngap[i]),Form("Resolution^{2} for v_{3} %s",gap[i]),9,0,9);
    }

    Int_t n_sys=3;
    // 0 - VtZ
    // 1 - nHits
    // 2 - DCA

    TH1D *h_pt_pion[2];
    for(Int_t i=0; i<2; i++){
      h_pt_pion[i] = new TH1D(Form("h_pt_pion%s",sign[i]), 
                              Form("Pt mean Pion %s",sign[i]), 60,0.2,3.2);
    }
    
    TProfile *tp_v1[2][9];
    TProfile *tp_v2cent[3][n-2], *tp_v3cent[3][n-2];
    TProfile2D *tp_v2[3][n-2], *tp_v3[3][n-2];
    TProfile2D *tp_v2_sys[3][n-2][3], *tp_v3_sys[3][n-2][3];
    TProfile2D *tp_v2PID[2][3][3][n-2], *tp_v3PID[2][3][3][n-2];
    TProfile2D *tp_v2PID_sys[2][3][3][n-2][3], *tp_v3PID_sys[2][3][3][n-2][3];
    //TProfile for mean pt for PID and hadrons
    TProfile *tp_meanPt_PID[2][3][3][n-2];
    TProfile *tp_meanPt_PID_sys[2][3][3][n-2][3];
    TProfile *tp_meanPt_hadrons[3][n-2];
    TProfile *tp_meanPt_hadrons_sys[3][n-2][3];

    for(Int_t det = 0; det < 2; det++) {
      for(Int_t j = 0; j < c; j++) {
        tp_v1[det][j] = new TProfile( Form( "tp_v1%scent%i",detector[det+1],j), 
                                 Form( "v_{1} of pseudorapidity cent %i by %s",j,detector[det+1]),20,-1,1);
      }
    }

    for(Int_t i = 0; i < n-2; i++) {
      for(Int_t j = 0; j < 3; j++) {
        tp_v2cent[j][i] = new TProfile(Form("tp_v2cent%s%s",detector[j],ngap[i]),
                              Form("v_{2} of cent %s by %s;cent;v_{2}",gap[i],detector[j]),9,0,9);
        tp_v3cent[j][i] = new TProfile(Form("tp_v3cent%s%s",detector[j],ngap[i]),
                              Form("v_{3} of cent %s by %s;cent;v_{3}",gap[i],detector[j]),9,0,9);
        tp_v2[j][i] = new TProfile2D(Form("tp_v2_hadrons%s%s",detector[j],ngap[i]),
                          Form("v_{2} of p_{t} and cent %s by %s;p_{t} [GeV/c];cent",gap[i],detector[j]),30,0.2,3.2,9,0,9);
        tp_v3[j][i] = new TProfile2D(Form("tp_v3_hadrons%s%s",detector[j],ngap[i]),
                          Form("v_{3} of p_{t} and cent %s by %s;p_{t} [GeV/c];cent",gap[i],detector[j]),30,0.2,3.2,9,0,9);
        tp_meanPt_hadrons[j][i] = new TProfile(Form("tp_meanPt_hadrons%s%s",detector[j],ngap[i]),
                          Form("Mean p_{t} for bins v_{2} #Delta#eta-gap=%s %s; bin; p_{t} [GeV/c]",gap[i],detector[j]), 30, 0.2, 3.2);

        if(i < 3 && j==0){
          for(Int_t nsys = 0; nsys < 3; nsys++){
            
            tp_v2_sys[j][i][nsys] = new TProfile2D(Form("tp_v2_hadrons%s%s%s",detector[j],ngap[i], systematic[nsys]),
                              Form("v_{2} of p_{t} and cent %s by %s Systematic:%s;p_{t} [GeV/c];cent",gap[i],detector[j], systematic[nsys]),30,0.2,3.2,9,0,9);
            tp_v3_sys[j][i][nsys] = new TProfile2D(Form("tp_v3_hadrons%s%s%s",detector[j],ngap[i],systematic[nsys]),
                              Form("v_{3} of p_{t} and cent %s by %s Systematic:%s;p_{t} [GeV/c];cent",gap[i],detector[j],systematic[nsys]),30,0.2,3.2,9,0,9);
            tp_meanPt_hadrons_sys[j][i][nsys] = new TProfile(Form("tp_meanPt_hadrons%s%s%s",detector[j],ngap[i],systematic[nsys]),
                          Form("Mean p_{t} for bins v_{2} #Delta#eta-gap=%s %s Systematic:%s; bin; p_{t} [GeV/c]",gap[i],detector[j],systematic[nsys]), 30, 0.2, 3.2);
          }
        }

        Int_t part = 0;
        for(Int_t k = 0; k < 3; k++) {
          for(Int_t l = 0; l < 2; l++) {
            tp_v2PID[l][k][j][i] = new TProfile2D( Form("tp_v2_%s%s%s%s",particles[k],sign[l],detector[j],ngap[i]), 
              Form("v_{2} for %s of p_{t} and cent %s by %s;p_{t} [GeV/c];cent",partLateX[part],gap[i],detector[j]),30,0.2,3.2,9,0,9);
            tp_v3PID[l][k][j][i] = new TProfile2D( Form("tp_v3_%s%s%s%s",particles[k],sign[l],detector[j],ngap[i]), 
              Form("v_{3} for %s of p_{t} and cent %s by %s;p_{t} [GeV/c];cent",partLateX[part],gap[i],detector[j]),30,0.2,3.2,9,0,9);
            tp_meanPt_PID[l][k][j][i] = new TProfile(Form("tp_meanPt_%s%s%s%s",particles[k],sign[l],detector[j],ngap[i]),
              Form("Mean p_{t} for %s of p_{t} and cent %s by %s;p_{t} [GeV/c];cent",particles[k],sign[l],detector[j],ngap[i]), 30, 0.2, 3.2);
            
            if(i < 3 && j==0){
              for(Int_t nsys = 0; nsys < 3; nsys++){
                tp_v2PID_sys[l][k][j][i][nsys] = new TProfile2D( Form("tp_v2_%s%s%s%s%s",particles[k],sign[l],detector[j],ngap[i], systematic[nsys]), 
                  Form("v_{2} for %s of p_{t} and cent %s by %s Systematic: %s;p_{t} [GeV/c];cent",partLateX[part],gap[i],detector[j], systematic[nsys]),30,0.2,3.2,9,0,9);
                tp_v3PID_sys[l][k][j][i][nsys] = new TProfile2D( Form("tp_v3_%s%s%s%s%s",particles[k],sign[l],detector[j],ngap[i], systematic[nsys]), 
                  Form("v_{3} for %s of p_{t} and cent %s by %s Systematic: %s;p_{t} [GeV/c];cent",partLateX[part],gap[i],detector[j], systematic[nsys]),30,0.2,3.2,9,0,9);
                tp_meanPt_PID_sys[l][k][j][i][nsys] = new TProfile(Form("tp_meanPt_%s%s%s%s%s",particles[k],sign[l],detector[j],ngap[i], systematic[nsys]),
                  Form("Mean p_{t} for %s of p_{t} and cent %s by %s Systematic: %s;p_{t} [GeV/c];cent",particles[k],sign[l],detector[j],ngap[i],systematic[nsys]), 30, 0.2, 3.2);
              }
            }
            part++;
          }
        }
      } //for(Int_t i = 0; i < n; i++)
    } //for(Int_t j = 0; j < 3; j++)

    TProfile *tp_meanPt = new TProfile("tp_meanPt","Mean p_{t} for bins v_{2}; bin; p_{t} [GeV/c]", 14, 0, 14);
    
    TFile *f1 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/NoRe_%s.root",energy),"READ");
    TProfile2D *tpQx1[2][2], *tpQy1[2][2];
    TProfile2D *tpQx2[2][n], *tpQy2[2][n], *tpQx3[2][n], *tpQy3[2][n];
    for(Int_t l = 0; l < 2; l++) {   
      for(Int_t i = 0; i < n; i++) {
        tpQx2[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx2%s%s",direction[l],ngap[i]) );
        tpQx3[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx3%s%s",direction[l],ngap[i]) );
        tpQy2[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy2%s%s",direction[l],ngap[i]) );
        tpQy3[l][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy3%s%s",direction[l],ngap[i]) );
      }
      for(Int_t det = 0; det < 2; det++) {
        tpQx1[l][det] = (TProfile2D*) f1 -> Get( Form("tp_Qx1%s%s",direction[l],ngap[det+6]) );
        tpQy1[l][det] = (TProfile2D*) f1 -> Get( Form("tp_Qy1%s%s",direction[l],ngap[det+6]) );
      }
    }

    TFile *f2 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/Re_%s.root",energy),"READ");
    TProfile2D *tpsinPsi1[2][2][20], *tpcosPsi1[2][2][20];
    TProfile2D *tpsinPsi2[2][n][4], *tpcosPsi2[2][n][4], *tpsinPsi3[2][n][4], *tpcosPsi3[2][n][4];
    for(Int_t l = 0; l < 2; l++) {
      for(Int_t i = 0; i < n; i++) {
        for(Int_t j = 0; j < 4; j++) {
          tpsinPsi2[l][i][j] = (TProfile2D*) f2 -> Get( Form("tp_%isinPsi2%s%s",j+1,direction[l],ngap[i]) );
          tpcosPsi2[l][i][j] = (TProfile2D*) f2 -> Get( Form("tp_%icosPsi2%s%s",j+1,direction[l],ngap[i]) );
          tpsinPsi3[l][i][j] = (TProfile2D*) f2 -> Get( Form("tp_%isinPsi3%s%s",j+1,direction[l],ngap[i]) );
          tpcosPsi3[l][i][j] = (TProfile2D*) f2 -> Get( Form("tp_%icosPsi3%s%s",j+1,direction[l],ngap[i]) );
        }
      }
      for(Int_t det = 0; det < 2; det++) {
        for(Int_t j = 0; j < 20; j++) { 
          tpsinPsi1[l][det][j] = (TProfile2D*) f2 -> Get( Form("tp_%isinPsi1%s%s",j+1,direction[l],ngap[det+6]) );
          tpcosPsi1[l][det][j] = (TProfile2D*) f2 -> Get( Form("tp_%icosPsi1%s%s",j+1,direction[l],ngap[det+6]) );
        }
      }
    }
 
    for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
      if ( iEvent % 10000 == 0) {
        std::cout << "Working on event #[" << (iEvent+1)
     	      		  << "/" << events2read << "]" << std::endl;
      }
      if (iEvent == events2read-1) {
        std::cout << "Working on event #[" << (events2read)
     	            << "/" << events2read << "]" << std::endl;
      }

    	Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
      if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
        break;
    	}

    	// Retrieve femtoDst
    	StFemtoDst *dst = femtoReader->femtoDst();

    	// Retrieve event information
    	StFemtoEvent *event = dst->event();
    	if( !event ) {
      	std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      	break;
    	}

    	// Simple event cut
      if ( !isGoodEvent( event, VtxZ, VtxR, delta_VtxY ) ) continue;

      TVector3 pVtx = event->primaryVertex();

      // Track analysis
      Int_t nTracks = dst->numberOfTracks();

      //RunQA
      if(useRunQA == true) {
        if( std::find(BadRuns.begin(), BadRuns.end(), event -> runId()) != BadRuns.end() ) continue; 
      }

      Double_t sinPsi1 = 0., cosPsi1 = 0., sinPsi2 = 0., cosPsi2 = 0., sinPsi3 = 0., cosPsi3 = 0.;
      Double_t dPsi1 = 0., dPsi2 = 0., dPsi3 = 0.;
      Double_t Psi1[2][2], Psi2[2][n], Psi3[2][n];
      TVector2 Q1vec[2][2];
      TVector2 Q2vec[2][8], Q3vec[2][8];

      for(Int_t l = 0; l < 2; l++) { 
        for(Int_t det = 0; det < 2; det++) {
          Q1vec[l][det].Set(0.,0.);
        }
        for(Int_t i = 0; i < n-1; i++) {
          Q2vec[l][i].Set(0.,0.);
          Q3vec[l][i].Set(0.,0.);
        }
      }

      for(Int_t dir = 0; dir < 2; dir++) {
        Q1vec[dir][0] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[dir][1] = CalculateZDCQVec(event, dir);
      }

      for(Int_t l = 0; l < 2; l++) {
        Q2vec[l][6] = CalculateBBCQVec(event, l, 2.0);
        Q3vec[l][6] = CalculateBBCQVec(event, l, 3.0);
      }
      badEvent = CalculateTPCQVec(Q2vec, Q3vec, dst, n-2, DCA_EVENT);

      RunID = event -> runId();
      cent = event -> cent9();

      Recentering(Q1vec,Q2vec,Q3vec,RunID,cent,tpQx1,tpQy1,tpQx2,tpQy2,tpQx3,tpQy3);

      for(Int_t l = 0; l < 2; l ++) { 
        for(Int_t i = 0; i < n-1; i++) {  

          // Fill only TPC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i < 6) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

            Psi2[l][i] = 1.0/2.0 * Q2vec[l][i].Phi();
            Psi3[l][i] = 1.0/3.0 * Q3vec[l][i].Phi();
          }
          // Fill BBc and ZDC
          if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i > 5) {
            h_Qx2[l][i][cent] -> Fill( Q2vec[l][i].X() );
            h_Qy2[l][i][cent] -> Fill( Q2vec[l][i].Y() );
            h_Qx3[l][i][cent] -> Fill( Q3vec[l][i].X() );
            h_Qy3[l][i][cent] -> Fill( Q3vec[l][i].Y() );

            Psi2[l][i] = 1.0/2.0 * Q2vec[l][i].Phi();
            Psi3[l][i] = 1.0/3.0 * Q3vec[l][i].Phi();
          }       
        } //for(Int_t i = 0; i < n-1; i++)
        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[l][det].Mod() != 0. ) {
            h_Qx1[l][det][cent] -> Fill( Q1vec[l][det].X() );
            h_Qy1[l][det][cent] -> Fill( Q1vec[l][det].Y() );

            Psi1[l][det] = Q1vec[l][det].Phi();
          }
        }  
      
      } //for(Int_t l = 0; l < 2; l ++)   

      // Flattening stage 
      for(Int_t l = 0; l < 2; l++) {
        
          for(Int_t det = 0; det < 2; det++) {
            if( Q1vec[l][det].Mod() != 0. ) {
            for(Int_t j = 0; j < 20; j++) { 
               
              sinPsi1 = tpsinPsi1[l][det][j]->GetBinContent( tpsinPsi1[l][det][j]->FindBin(RunID, cent) );
              cosPsi1 = tpcosPsi1[l][det][j]->GetBinContent( tpcosPsi1[l][det][j]->FindBin(RunID, cent) );
              //std::cout << "sinPsi = "  << sinPsi1 << std::endl;
              //std::cout << "cosPsi = "  << cosPsi1 << std::endl;
   
              dPsi1 += -2.0*( sinPsi1 * TMath::Cos( (Double_t)(j+1)*Psi1[l][det] ) )/( (Double_t)(j+1) ) 
                       +2.0*( cosPsi1 * TMath::Sin( (Double_t)(j+1)*Psi1[l][det] ) )/( (Double_t)(j+1) );
            }
            //std::cout << "dPsi1 = "  << dPsi1 << std::endl; 
            Psi1[l][det] += dPsi1;
            dPsi1 = 0.;
          } 
        }//for(Int_t det = 0; det < 2; det++)

        
          for(Int_t i = 0; i < n; i++) {
            if( Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. ) {
            for(Int_t k = 0; k < 4; k++) {

              sinPsi2 = tpsinPsi2[l][i][k]->GetBinContent( tpsinPsi2[l][i][k] -> FindBin( RunID, cent ) );
              cosPsi2 = tpcosPsi2[l][i][k]->GetBinContent( tpcosPsi2[l][i][k] -> FindBin( RunID, cent ) );
              sinPsi3 = tpsinPsi3[l][i][k]->GetBinContent( tpsinPsi3[l][i][k] -> FindBin( RunID, cent ) );
              cosPsi3 = tpcosPsi3[l][i][k]->GetBinContent( tpcosPsi3[l][i][k] -> FindBin( RunID, cent ) );

              dPsi2 += -2.0*( sinPsi2 * TMath::Cos( (Double_t)(k+1)*2.0*Psi2[l][i] ) )/( 2.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi2 * TMath::Sin( (Double_t)(k+1)*2.0*Psi2[l][i] ) )/( 2.0*(Double_t)(k+1) );

              dPsi3 += -2.0*( sinPsi3 * TMath::Cos( (Double_t)(k+1)*3.0*Psi3[l][i] ) )/( 3.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi3 * TMath::Sin( (Double_t)(k+1)*3.0*Psi3[l][i] ) )/( 3.0*(Double_t)(k+1) );

            } //for(Int_t k = 0; k < 4; k++)
            Psi2[l][i] += dPsi2;
            Psi3[l][i] += dPsi3;
            dPsi2 = dPsi3 = 0.;
          } 
        }//for(Int_t i = 0; i < n; i++)

      } //for(Int_t l = 0; l < 2; l++)

      for(Int_t l = 0; l < 2; l ++) {  
        for(Int_t i = 0; i < n-1; i++) {  
         

          // Fill only TPC
          if(Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i < 6) {
            h_Psi2[l][i][cent] -> Fill( Psi2[l][i]  );
            h_Psi3[l][i][cent] -> Fill( Psi3[l][i]  );
          }
          // Fill BBc and ZDC
          if(Q2vec[l][i].Mod() != 0. && Q3vec[l][i].Mod() != 0. && i > 5) {
            h_Psi2[l][i][cent] -> Fill( Psi2[l][i]  );
            h_Psi3[l][i][cent] -> Fill( Psi3[l][i]  );
          }
        } //for(Int_t i = 0; i < n-1; i++)

        for(Int_t det = 0; det < 2; det++) {
          if(Q1vec[l][det].Mod() != 0.) {
            h_Psi1[l][det][cent] -> Fill( Psi1[l][det] );
          }
        }

      } //for(Int_t l = 0; l < 2; l ++)  

      for(Int_t det = 0; det < 2; det++) {
        if(Q1vec[0][det].Mod() != 0. && Q1vec[1][det].Mod() != 0.) {
          tp_SqRes1[det] -> Fill(cent, TMath::Cos( (Psi1[1][det] - Psi1[0][det]) ) );
        }
      }  

      for(Int_t i = 0; i < n-1; i++) {
        if(Q2vec[0][i].Mod() != 0. && Q3vec[0][i].Mod() != 0. && Q2vec[1][i].Mod() != 0. && Q3vec[1][i].Mod() != 0.) {
          tp_SqRes2[i] -> Fill(cent, TMath::Cos( 2*(Psi2[1][i] - Psi2[0][i]) ) );
          tp_SqRes3[i] -> Fill(cent, TMath::Cos( 3*(Psi3[1][i] - Psi3[0][i]) ) );
        }
      }

      //Track loop for Flow 
      for(Int_t iTrk = 0; iTrk < nTracks; iTrk++) {

        //Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst -> track(iTrk);

        // Simple single-track cut
        if( !isGoodTrackFlow(event, femtoTrack, DCA_SYS) ) continue;

        if(femtoTrack -> pt() < 2.0) {
          w = femtoTrack -> pt();
        }
        else w = 2.0;         

        Phi = femtoTrack -> phi();
        pt = femtoTrack -> pt();
        q = femtoTrack -> charge();
        Int_t q_charge;
        if(femtoTrack -> charge() > 0){
          q_charge=0;
        }
        if(femtoTrack -> charge() < 0){
          q_charge=1;
        }
        
        SqM = femtoTrack -> massSqr();

        nSigma[0] = femtoTrack -> nSigmaPion();
        nSigma[1] = femtoTrack -> nSigmaKaon();
        nSigma[2] = femtoTrack -> nSigmaProton();

        if( pt >= 0.2 && pt < 0.4) tp_meanPt -> Fill(0., pt);
        if( pt >= 0.4 && pt < 0.6) tp_meanPt -> Fill(1., pt);
        if( pt >= 0.6 && pt < 0.8) tp_meanPt -> Fill(2., pt);
        if( pt >= 0.8 && pt < 1.0) tp_meanPt -> Fill(3., pt);
        if( pt >= 1.0 && pt < 1.2) tp_meanPt -> Fill(4., pt);
        if( pt >= 1.2 && pt < 1.4) tp_meanPt -> Fill(5., pt);
        if( pt >= 1.4 && pt < 1.6) tp_meanPt -> Fill(6., pt);
        if( pt >= 1.6 && pt < 1.8) tp_meanPt -> Fill(7., pt);
        if( pt >= 1.8 && pt < 2.0) tp_meanPt -> Fill(8., pt);
        if( pt >= 2.0 && pt < 2.2) tp_meanPt -> Fill(9., pt);
        if( pt >= 2.2 && pt < 2.4) tp_meanPt -> Fill(10., pt);
        if( pt >= 2.4 && pt < 2.6) tp_meanPt -> Fill(11., pt);
        if( pt >= 2.6 && pt < 2.8) tp_meanPt -> Fill(12., pt);
        if( pt >= 2.8 && pt < 3.0) tp_meanPt -> Fill(13., pt);

        Double_t v1 = 0.,v2 = 0., v3 = 0.;
 
          

            if( femtoTrack -> eta() < 0 ) {
              v1 = TMath::Cos( Phi - Psi1[0][0] );
              tp_v1[0][cent] -> Fill( femtoTrack -> eta(), v1);
              v1 = TMath::Cos( Phi - Psi1[0][1] );
              tp_v1[1][cent] -> Fill( femtoTrack -> eta(), v1);
            }
            if( femtoTrack -> eta() > 0 ) {
              v1 = TMath::Cos( Phi - Psi1[1][0] );
              tp_v1[0][cent] -> Fill( femtoTrack -> eta(), v1);
              v1 = TMath::Cos( Phi - Psi1[1][1] );
              tp_v1[1][cent] -> Fill( femtoTrack -> eta(), v1);
              
            }

        for(Int_t det = 0; det < 3; det++) {
          for(Int_t eta = 0; eta < n-2; eta++) {

            if( det == 0 ) {
              if( femtoTrack -> eta() < -etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[1][eta]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[1][eta]) );

                if(eta<3){
                  if( TMath::Abs( pVtx.Z() ) < sys_VtxZ){
                    tp_v2_sys[det][eta][0]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][0]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][0] -> Fill(pt,pt);
                  }
                  if( femtoTrack -> nHits() > sys_nHits){
                    tp_v2_sys[det][eta][1]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][1]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][1] -> Fill(pt,pt);
                  }
                  if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA){
                    tp_v2_sys[det][eta][2]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][2]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][2] -> Fill(pt,pt);
                  }  
                }

                if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
                  tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                  tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                  tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                  tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                  tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);
                }
                

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < nSigm && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
                      tp_v2PID[q_charge][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                      tp_v3PID[q_charge][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                      tp_meanPt_PID[q_charge][part][det][eta] -> Fill(pt,pt);
                    }
            
                    if(eta<3){
                      if( TMath::Abs( pVtx.Z() ) < sys_VtxZ){
                        tp_v2PID_sys[q_charge][part][det][eta][0] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][0] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][0] -> Fill(pt,pt);
                      }
                      if( femtoTrack -> nHits() > sys_nHits){
                        tp_v2PID_sys[q_charge][part][det][eta][1] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][1] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][1] -> Fill(pt,pt);
                      }
                      if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA){
                        tp_v2PID_sys[q_charge][part][det][eta][2] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][2] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][2] -> Fill(pt,pt);
                      }  
                    }
                  }
                }

              }
              if( femtoTrack -> eta() > etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[0][eta]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[0][eta]) );

                if(eta<3){
                  if( femtoTrack -> nHits() > sys_nHits){
                    tp_v2_sys[det][eta][0]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][0]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][0] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( pVtx.Z() ) < sys_VtxZ){
                    tp_v2_sys[det][eta][1]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][1]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][1] -> Fill(pt,pt);
                  }
                  if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA){
                    tp_v2_sys[det][eta][2]->Fill(pt, (Double_t)cent, v2, w );
                    tp_v3_sys[det][eta][2]->Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_hadrons_sys[det][eta][2] -> Fill(pt,pt);
                  }  
                }

                if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
                  tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                  tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                  tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                  tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                }
                
                tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < nSigm && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
                      tp_v2PID[q_charge][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                      tp_v3PID[q_charge][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                      tp_meanPt_PID[q_charge][part][det][eta] -> Fill(pt,pt);
                    }
            
                    if(eta<3){
                      if( TMath::Abs( pVtx.Z() ) < sys_VtxZ){
                        tp_v2PID_sys[q_charge][part][det][eta][0] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][0] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][0] -> Fill(pt,pt);
                      }
                      if( femtoTrack -> nHits() > sys_nHits){
                        tp_v2PID_sys[q_charge][part][det][eta][1] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][1] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][1] -> Fill(pt,pt);
                      }
                      if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA){
                        tp_v2PID_sys[q_charge][part][det][eta][2] -> Fill(pt, (Double_t)cent, v2, w );
                        tp_v3PID_sys[q_charge][part][det][eta][2] -> Fill(pt, (Double_t)cent, v3, w );
                        tp_meanPt_PID_sys[q_charge][part][det][eta][2] -> Fill(pt,pt);
                      }  
                    }
                  }
                }
              }
            } //if( det== 0)

            if( det==1 ) {
              if(femtoTrack -> eta() < -etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[1][6]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[1][6]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < nSigm && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < nSigm && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);
                  }
                }
              }

              if(femtoTrack -> eta() > etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[0][6]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[0][6]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < nSigm && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < nSigm && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);
                  }
                }
              }

            } //if( det==1 )
          } //for(Int_t eta = 0; eta < n-2; eta++)
        } //for(Int_t det = 0; det < 3; det++)

        /* no systematic
        for(Int_t det = 0; det < 3; det++) {
          for(Int_t eta = 0; eta < n-2; eta++) {

            if( det == 0 ) {
              if( femtoTrack -> eta() < -etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[1][eta]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[1][eta]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < 3 && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < 3 && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);

                  }
                }

              }
              if( femtoTrack -> eta() > etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[0][eta]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[0][eta]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < 3 && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < 3 && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);

                  }
                }
              }
            } //if( det== 0)

            if( det==1 ) {
              if(femtoTrack -> eta() < -etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[1][6]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[1][6]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                tp_meanPt_hadrons[det][eta] -> Fill(pt,pt);

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < 3 && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < 3 && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);
                  }
                }
              }

              if(femtoTrack -> eta() > etagap[eta] ) {
                v2 = TMath::Cos( 2.0*(Phi - Psi2[0][6]) );
                v3 = TMath::Cos( 3.0*(Phi - Psi3[0][6]) );

                tp_v2cent[det][eta] -> Fill( (Double_t)cent, v2, w );
                tp_v3cent[det][eta] -> Fill( (Double_t)cent, v3, w );
                tp_v2[det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                tp_v3[det][eta] -> Fill(pt, (Double_t)cent, v3, w );

                for(Int_t part = 0; part < 3; part++) {
                  if( TMath::Abs( nSigma[part] ) < 3 && q > 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[0][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[0][part][det][eta] -> Fill(pt,pt);
                  }
                  if( TMath::Abs( nSigma[part] ) < 3 && q < 0 && SqM > SqMdown[part] && SqM < SqMup[part]) {
                    tp_v2PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v2, w );
                    tp_v3PID[1][part][det][eta] -> Fill(pt, (Double_t)cent, v3, w );
                    tp_meanPt_PID[1][part][det][eta] -> Fill(pt,pt);
                  }
                }
              }

            } //if( det==1 )
          } //for(Int_t eta = 0; eta < n-2; eta++)
        } //for(Int_t det = 0; det < 3; det++)
        */  //no systematic

      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

    f1 -> Close();
    f2 -> Close();
    TFile *oFile = new TFile(oFileName, "RECREATE");
    tp_meanPt -> Write();
    for(Int_t det = 0; det < 2; det++) {
      for(Int_t l = 0; l < 2; l++) {
        for(Int_t j = 0; j < c; j++) { 
          h_Qx1[l][det][j] -> Write();
          h_Qy1[l][det][j] -> Write();
          h_Psi1[l][det][j] -> Write();
        }
      }
      tp_SqRes1[det] -> Write();
      for(Int_t j = 0; j < c; j++) {
        tp_v1[det][j] -> Write();
      }
    }
    for(Int_t i = 0; i < n; i++) {
      for(Int_t l = 0; l < 2; l++) {
        for(Int_t j = 0; j < c; j++ ) {

          h_Qx2[l][i][j] -> Write();
          h_Qy2[l][i][j] -> Write();
          h_Qx3[l][i][j] -> Write();
          h_Qy3[l][i][j] -> Write();

          h_Psi2[l][i][j] -> Write();
          h_Psi3[l][i][j] -> Write();
        }
      }
      tp_SqRes2[i] -> Write();
      tp_SqRes3[i] -> Write();
    }
        
    for(Int_t det = 0; det < 3; det++) {
      for(Int_t eta = 0; eta < n-2; eta++) {
        tp_v2cent[det][eta] -> Write();
        tp_v3cent[det][eta] -> Write();
        tp_v2[det][eta] -> Write();
        tp_v3[det][eta] -> Write();
        tp_meanPt_hadrons[det][eta] -> Write();

        if(eta < 3 && det==0){
          for(Int_t nsys = 0; nsys < 3; nsys++){
            tp_v2_sys[det][eta][nsys] -> Write();
            tp_v3_sys[det][eta][nsys] -> Write();
            tp_meanPt_hadrons_sys[det][eta][nsys] -> Write();
          }
        }

        for(Int_t part = 0; part < 3; part++) {
          for(Int_t l = 0; l < 2; l++) {
            tp_v2PID[l][part][det][eta] -> Write();
            tp_v3PID[l][part][det][eta] -> Write();
            tp_meanPt_PID[l][part][det][eta] -> Write();
            if(eta < 3 && det==0){
              for(Int_t nsys = 0; nsys < 3; nsys++){
                tp_v2PID_sys[l][part][det][eta][nsys]-> Write();
                tp_v3PID_sys[l][part][det][eta][nsys] -> Write();
                tp_meanPt_PID_sys[l][part][det][eta][nsys]-> Write();
              }
            }
          }
        }
      }
    }
    for(Int_t i=0;i<2;i++){
      h_pt_pion[i]->Write();
    }
    oFile->Close();
  }//flattening and flow mode 
		        
	  //////////////////////////////////////////////////////baaaaaaaad Ruuuuuuuuns!!!!!
    //BadRuns.clear();
    
  	std::cout << "BadRun clean!" << std::endl;
     // oFile->Write();
      //oFile->Close();
      
  	femtoReader->Finish();
  	std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	    	      << std::endl;
}


//_________________
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZ, Double_t VtxR, Double_t delta_VtxY) {
  Bool_t check = true; 
  TVector3 pVtx = event->primaryVertex();

  // Reject vertices that are far from the central membrane along the beam
  if( TMath::Abs( pVtx.Z() ) > VtxZ ) check = false;
  //if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + 0.8847, 2) ) > 1 ) check = false; 14.5
  if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + delta_VtxY, 2) ) > VtxR ) check = false;
  
  return check;
}

//_________________
Bool_t isGoodTrack(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_EVENT) {
  Bool_t check = true;
  TVector3 pVtx = event->primaryVertex();

  if ( !femtoTrack ) check = false;
  // Must be a primary track
  if ( !femtoTrack->isPrimary() ) check = false;
  if ( ( femtoTrack -> dEdx() ) == 0 ) check = false;
  // Simple single-track cut
  if( femtoTrack -> gMom().Mag() < 0.1 || femtoTrack -> gDCA(pVtx).Mag() > DCA_EVENT ) check = false;    
  if( TMath::Abs( femtoTrack -> eta() ) > 1.0 || femtoTrack -> nHits() < 15 || 
                  femtoTrack -> pt() < 0.2) check = false;
  if(femtoTrack -> pt() > 2.0) check = false; 
  if(  ( (Double_t)femtoTrack -> nHits() )/( (Double_t)femtoTrack -> nHitsPoss() )  < 0.52 ) check = false;

  return check; 
}

//_________________
Bool_t isGoodTrackFlow(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_FLOW) {
  Bool_t check = true;
  TVector3 pVtx = event->primaryVertex();

  if ( !femtoTrack ) check = false;
  // Must be a primary track
  if ( !femtoTrack->isPrimary() ) check = false;
  if ( ( femtoTrack -> dEdx() ) == 0 ) check = false;
  // Simple single-track cut
  if( femtoTrack -> gMom().Mag() < 0.1 || femtoTrack -> gDCA(pVtx).Mag() > DCA_FLOW ) check = false;   //1.5 
  if( TMath::Abs( femtoTrack -> eta() ) > 1.0 || femtoTrack -> nHits() < 15 || 
                  femtoTrack -> pt() < 0.2) check = false;
  if(femtoTrack -> pt() > 3.2) check = false; 
  if(  ( (Double_t)femtoTrack -> nHits() )/( (Double_t)femtoTrack -> nHitsPoss() )  < 0.52 ) check = false;

  return check; 
}

//_________________
void Recentering( TVector2 Q1vec[][2], TVector2 Q2vec[][n], TVector2 Q3vec[][n], Int_t RunID, Int_t cent,TProfile2D *tpQx1[][2], 
TProfile2D *tpQy1[][2], TProfile2D *tpQx2[][n], TProfile2D *tpQy2[][n],TProfile2D *tpQx3[][n], TProfile2D *tpQy3[][n] ) {
    
  TVector2 Q1mean, Q2mean, Q3mean;

  for(Int_t l = 0; l < 2; l++) {
    for(Int_t i = 0; i < n; i++) {
      Q2mean.Set(0.,0.);
      Q3mean.Set(0.,0.);

      Q2mean.Set( tpQx2[l][i] -> GetBinContent( tpQx2[l][i] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy2[l][i] -> GetBinContent( tpQy2[l][i] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      Q3mean.Set( tpQx3[l][i] -> GetBinContent( tpQx3[l][i] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy3[l][i] -> GetBinContent( tpQy3[l][i] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      //std::cout << "Q2meanX = " << Q2mean.X() << std::endl;
      //std::cout << "Q2meanX = " << Q2mean.Y() << std::endl;
      if( Q2vec[l][i].Mod() != 0.) { 
        Q2vec[l][i] -= Q2mean;
      }
      if( Q3vec[l][i].Mod() != 0.) {
        Q3vec[l][i] -= Q3mean;
      }
      //std::cout << "Recentering" << std::endl; 
      //std::cout << "Qx2[" << l << "][" << i << "] = " <<Q2vec[l][i].X() << std::endl;
      //std::cout << "Qy2[" << l << "][" << i << "] = " <<Q2vec[l][i].Y() << std::endl;
    }
    for(Int_t det = 0; det < 2; det++) { 
      Q1mean.Set(0.,0.);
      Q1mean.Set( tpQx1[l][det] -> GetBinContent( tpQx1[l][det] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy1[l][det] -> GetBinContent( tpQy1[l][det] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      if( Q1vec[l][det].Mod() != 0.) { 
        Q1vec[l][det] -= Q1mean;
      }
    }
  }
  /*std::cout << "Recentering" << std::endl; 
      for(Int_t l = 0; l < 2; l++ ){
          std::cout << "Qx2[" << l << "][" << 0 << "] = " <<Q2vec[l][0].X() << std::endl;
          std::cout << "Qy2[" << l << "][" << 0 << "] = " <<Q2vec[l][0].Y() << std::endl;
        }*/ 
}

Bool_t CalculateTPCQVec( TVector2 Q2vec[][8], TVector2 Q3vec[][8], StFemtoDst *dst, Int_t ngap, Double_t DCA_EVENT) {
  // 0 - this is east 
  // 1 - this is west
  Bool_t bad = true; 
  Double_t gap[6] = {0.05, 0.075, 0.1, 0.15, 0.2, 0.5};
  Float_t w[2][ngap];

  // Retrieve event information
  StFemtoEvent *event = dst->event();

  for(Int_t i = 0; i < 2; i++) {
    for(Int_t j = 0; j < ngap; j++) { 
      w[i][j] = 0.;
    }
  }

  Int_t nTracks = dst->numberOfTracks();

  // track loop 
  for(Int_t iTrk = 0; iTrk < nTracks; iTrk++) {

    // Retrieve i-th femto track
    StFemtoTrack *femtoTrack = dst -> track(iTrk);

    if( !isGoodTrack(event, femtoTrack, DCA_EVENT) ) continue;
    bad = false; 

    for(Int_t i = 0; i < ngap; i++) {

      if( femtoTrack -> eta() < -gap[i] ) {
        Q2vec[0][i].Set( Q2vec[0][i].X() + femtoTrack -> pt() * cos(2.0 * femtoTrack -> phi() ),
                         Q2vec[0][i].Y() + femtoTrack -> pt() * sin(2.0 * femtoTrack -> phi() ) );
        Q3vec[0][i].Set( Q3vec[0][i].X() + femtoTrack -> pt() * cos(3.0 * femtoTrack -> phi() ),
                         Q3vec[0][i].Y() + femtoTrack -> pt() * sin(3.0 * femtoTrack -> phi() ) );
        w[0][i]++; 
      }
      
      if( femtoTrack -> eta() > gap[i] ) {
        Q2vec[1][i].Set( Q2vec[1][i].X() + femtoTrack -> pt() * cos(2.0 * femtoTrack -> phi()),
                         Q2vec[1][i].Y() + femtoTrack -> pt() * sin(2.0 * femtoTrack -> phi() ) );
        Q3vec[1][i].Set( Q3vec[1][i].X() + femtoTrack -> pt() * cos(3.0 * femtoTrack -> phi() ),
                         Q3vec[1][i].Y() + femtoTrack -> pt() * sin(3.0 * femtoTrack -> phi() ) );
        w[1][i]++;
      }
    }
  }// track loop end

  for(Int_t j = 0; j < 2; j++) {
    for(Int_t i = 0; i < ngap; i++) {
      if( w[j][i] != 0 ) {
        Q2vec[j][i].Set( Q2vec[j][i].X()/w[j][i], Q2vec[j][i].Y()/w[j][i]);
        Q3vec[j][i].Set( Q3vec[j][i].X()/w[j][i], Q3vec[j][i].Y()/w[j][i]);
      }
    }
  }
  return bad;
} //CalculateTPCQVec

//_________________
Float_t ZDCSMD(Int_t eastwest, Int_t verthori, Int_t strip, StFemtoEvent *event ) {
  Float_t val = 0;

  if (eastwest == 0 && verthori == 0) val = event->zdcSmdEastVertical(strip-1);
  if (eastwest == 0 && verthori == 1) val = event->zdcSmdEastHorizontal(strip-1);
  if (eastwest == 1 && verthori == 0) val = event->zdcSmdWestVertical(strip-1);
  if (eastwest == 1 && verthori == 1) val = event->zdcSmdWestHorizontal(strip-1);
  
  return val;
}

Float_t ZDCSMD_GetPosition( Int_t eastwest, Int_t verthori, Int_t strip) {
  Float_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
  Float_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1]-mZDCSMDCenterex;
  if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x[strip-1];
  if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterey;
  if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterwy;

  return 0;
}

TVector2 CalculateZDCQVec(StFemtoEvent *event, Int_t ew) {
  Int_t nstrip = 0;
  Float_t Qvec[4]={0};
  TVector2 Q(0,0);

  if (ew == 0 && event->zdcSumAdcEast() < 1) return Q;
  if (ew == 1 && event->zdcSumAdcWest() < 1) return Q;
  
  for (Int_t xy = 0; xy < 2; xy++) {

    if (xy == 0) nstrip = 8;
    else         nstrip = 9;

    for( Int_t is=1; is<nstrip; is++ ){

      Float_t smd = ZDCSMD(ew,xy,is,event);
      Qvec[xy] += ZDCSMD_GetPosition(ew, xy, is)*smd; Qvec[xy+2] += smd;

    }

    if (Qvec[xy+2] > 0) Qvec[xy] /= Qvec[xy+2];
    //      cout << Qvec[xy] << endl;
  }//xy

  Q.Set(Qvec[0],Qvec[1]);
  return Q;

}//CalculateZDCQVec

//_________________
TVector2 CalculateBBCQVec(StFemtoEvent *event, Int_t ew, Float_t harm) {
   Float_t Qx = 0., Qy = 0., Qwgt = 0.;
   TVector2 Q(0.,0.);
    int nstrip = 16;
                        for( int is=0; is<nstrip; is++ ){
                                Float_t adc = 0;
                                if (ew == 0 ) adc = event->bbcAdcEast(is);//egain[is];
                                if (ew == 1 ) adc = event->bbcAdcWest(is);//wgain[is];
                                if (is == 6 || is == 11) adc *= 2;
                                Qx += cos(harm*GetBBCTilePhi(ew, is))*adc;
                                Qy += sin(harm*GetBBCTilePhi(ew, is))*adc;
                                Qwgt += adc;
                        }
                        if( Qwgt>0 ) {Qx /= Qwgt; Qy /= Qwgt;}
                        else                    return Q;
                Q.Set(Qx,Qy);
                return Q;
}//CalculateBBCQVec

Float_t GetBBCTilePhi(const Int_t e_w, const Int_t iTile){
  Double_t pi = TMath::Pi();
  Double_t phi_div  = pi/6.0;
  Float_t bbc_phi = phi_div;

  switch(iTile) {
    case 0: bbc_phi=3*phi_div;
    break;
    case 1: bbc_phi=phi_div;
    break;
    case 2: bbc_phi=-1*phi_div;
    break;
    case 3: bbc_phi=-3*phi_div;
    break;
    case 4: bbc_phi=-5*phi_div;
    break;
    case 5: bbc_phi=5*phi_div;
    break;
    case 6: bbc_phi= 3*phi_div;
    break;
    case 7: bbc_phi=3*phi_div;
    break;
    case 8: bbc_phi=phi_div;
    break;
    case 9: bbc_phi=0.;
    break;
    case 10: bbc_phi=-phi_div;
    break;
    case 11: bbc_phi=-3*phi_div;
    break;
    case 12: bbc_phi=-3*phi_div;
    break;
    case 13: bbc_phi=-5*phi_div;
    break;
    case 14: bbc_phi=pi;
    break;
    case 15: bbc_phi=5*phi_div;
    break;
  }

  if(e_w==0){
    if (bbc_phi > -0.001){ bbc_phi = pi-bbc_phi;}
    else {bbc_phi= -pi-bbc_phi;}
  }

if(bbc_phi<0.0) bbc_phi +=2*pi;
if(bbc_phi>2*pi) bbc_phi -=2*pi;

return bbc_phi;
}
