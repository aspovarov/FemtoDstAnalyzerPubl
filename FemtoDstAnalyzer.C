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

const Int_t n = 6;	// eta-gap by TPC 

// Forward declarations
// Check event and track
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZ, Double_t VtxR, Double_t delta_VtxY);
Bool_t isGoodTrack(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_EVENT);
Bool_t isGoodTrackFlow(StFemtoEvent *event, StFemtoTrack *femtoTrack, Double_t DCA_FLOW);

// Calculate Q-vector by TPC
void CalculateTPCQVec( TVector2 Q2vec[][n], TVector2 Q3vec[][n], StFemtoDst *dst,Int_t ngap, Double_t DCA_EVENT);

// Calculate Q-vector by ZDC
Float_t ZDCSMD(Int_t eastwest, Int_t verthori, Int_t strip, StFemtoEvent *event );
Float_t ZDCSMD_GetPosition( Int_t eastwest, Int_t verthori, Int_t strip);
TVector2 CalculateZDCQVec(StFemtoEvent *event, Int_t ew);

// Calculate Q-vector by BBC
TVector2 CalculateBBCQVec(StFemtoEvent *event, Int_t ew, Float_t harm);
Float_t GetBBCTilePhi(const Int_t e_w, const Int_t iTile);

// Recentering for TPC, BBC, ZDC
void Recentering( TVector2 Q1vec[][3], TVector2 Q2vec[][n], TVector2 Q3vec[][n], Int_t RunID, Int_t cent,TProfile2D *tpQx1[][3],
TProfile2D *tpQy1[][3],TProfile2D *tpQx2[][n],TProfile2D *tpQy2[][n],TProfile2D *tpQx3[][n],TProfile2D *tpQy3[][n] );

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
                      const Char_t *outFileName = "oTest.root",
                      const Char_t *mode = "QAmode",
                      const Char_t *energy = "39GeV") {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  gSystem->Load("/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/libStFemtoDst.so");
  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  gSystem->Load("/mnt/pool/rhic/1/nigmatkulov/soft/StFemtoEvent/libStFemtoDst.so");
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
  std::vector<Int_t> badRuns;
  
  Double_t sys_VtxZ = 45.0;
  Int_t sys_nHits = 18;
  Double_t sys_DCA = 1.5;

  Int_t runIdBins;
  Int_t runIdRange[2];

  if( strncmp(energy, "39GeV",5) == 0) {

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
    badRuns ={11099102, 11099103, 11099104, 11099106, 11099107, 
    		  11099125, 11100004, 11100005, 11100008, 11100010, 
    		  11100011, 11100016, 11100020, 11100071, 11101014, 
    		  11101104, 11102098, 11103009, 11103047, 11103065, 
    		  11105011, 11105029, 11106026, 11106027, 11106028, 
    		  11106029, 11106030, 11106040, 11106041, 11107008, 
    		  11107046, 11107083, 11108040, 11108053, 11108065, 
    		  11108075, 11109092, 11109102, 11109105, 11109104, 
    		  11110005, 11110041, 11110042, 11110086};

  }// if( strncmp(energy, "39GeV",5) == 0 )

  if( strncmp(energy, "27GeV",5) == 0 ){
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
    badRuns ={12172049, 12172056, 12173009, 12173018, 12173026, 
    		  12173053, 12173054, 12173055, 12173056, 12173057, 
    		  12173072, 12174077, 12174096, 12174109, 12175007, 
    		  12175030, 12175089, 12176046, 12176047, 12176067, 
    		  12176069, 12176104, 12178051, 12178093, 12179068, 
    		  12179083, 12179084, 12179085, 12179086};

  }// if( strncmp( energy, "27GeV",5) == 0 )

  if( strncmp(energy, "19GeV",5) == 0) {
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
    badRuns ={12113081, 12113084, 12114077, 12114079, 12114085, 
    	 	  12114088, 12114089, 12114091, 12114094, 12114095, 
    	 	  12114097, 12114098, 12114099, 12114100, 12114101, 
    	 	  12114102, 12114103, 12114104, 12114110, 12115025, 
    	 	  12115026, 12116015, 12116016, 12116063, 12116084, 
    	 	  12117009, 12118045, 12119008, 12119009, 12119011, 
    	 	  12119015, 12119016, 12119017, 12119019, 12119020, 
    	 	  12119021, 12119022, 12119023, 12119024, 12119025, 
    	 	  12119027, 12119028, 12119029, 12119030, 12119032, 
    	 	  12119035, 12119036, 12119039, 12120018, 12120073};

  }// if( strncmp(energy, "19GeV",5) == 0 )

  if( strncmp(energy, "14GeV",5) == 0 ) {
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
    badRuns ={15053027, 15053028, 15053029, 15053034, 15053035, 
			  15053048, 15053052, 15053053, 15053054, 15053055, 
			  15053056, 15053057, 15054019, 15054053, 15054054, 
			  15055018, 15055131, 15055133, 15055134, 15055135, 
			  15055136, 15055137, 15055138, 15055139, 15055140, 
			  15055141, 15056001, 15056002, 15056003, 15056004, 
			  15056005, 15056006, 15056007, 15056008, 15056009, 
			  15056113, 15056114, 15056116, 15056117, 15056124, 
			  15056125, 15057001, 15057003, 15057004, 15057006, 
			  15057007, 15057010, 15057011, 15057013, 15057014, 
			  15057018, 15057055, 15057059, 15058006, 15058011, 
			  15058018, 15060061, 15060069, 15061001, 15061002, 
			  15062006, 15065012, 15065014, 15066008, 15066013, 
			  15066017, 15068013, 15068014, 15068016, 15069034, 
			  15069036, 15070009, 15070010};

  }// if( strncmp(energy, "14GeV",5) == 0)

  if( strncmp(energy, "11GeV",5) == 0 ) {
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
    badRuns ={11148001, 11148008, 11148009, 11148010, 11148036, 
    		  11148055, 11149017, 11149018, 11149040, 11149043, 
    		  11150017, 11150029, 11151051, 11151057, 11153045, 
    		  11154026, 11154040, 11154059, 11156036, 11156043, 
    		  11156044, 11156045, 11157039};

  }// if( strncmp(energy, "11GeV",5) == 0 )

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

  }// if(strncmp(energy, "7GeV",4)==0)
 
  Int_t c = 9;    // cent9()
  Int_t cent, RunID;
  Double_t w = 0.0;
  Double_t Phi, pt, q, SqM;

  TVector2 Q1vec[2][3];
  TVector2 Q2vecTPC[3][n], Q3vecTPC[3][n];
  Double_t Psi1[2][3]; 
  Double_t Psi2TPC[3][n], Psi3TPC[3][n];

  Double_t etagap[] ={0.05, 0.075, 0.1, 0.15, 0.2, 0.5};
  Double_t SqMdown[] = {-0.15, 0.2, 0.74};
  Double_t SqMup[] = {0.1, 0.32, 1.20};
  Double_t nSigma[3];

	const Char_t *resol[] = {"ew","comb"};
  const Char_t *direction[] = {"east","west","comb"};
  const Char_t *detector[] = {"TPC","BBC","ZDC"};
  const Char_t *ngap[] = {"Eta01","Eta15","Eta02","Eta03","Eta04","Eta1"};
  const Char_t *gap[] = {"gap 0.1","gap 0.15","gap 0.2","gap 0.3","gap 0.4","gap 1.0"};
  const Char_t *sign[] = {"Pos","Neg"};
  const Char_t *particles[] = {"Pion","Kaon","Proton","hadrons"};
  const Char_t *partLateX[] = {"#pi^{+}","pi^{-}","K^{+}","K^{-}","p","#bar{p}"};
  const Char_t *systematic[] = {"VtxZ", "nHits", "DCA"};
  const Char_t *systematic_text[] = {"45", "18", "1.5"};

  TFile *outFile = new TFile(outFileName, "RECREATE");

  TFile *f1, *f2;
	TProfile2D *tpQx1[2][3], *tpQy1[2][3];
	TProfile2D *tpQx2[3][n], *tpQy2[3][n], *tpQx3[3][n], *tpQy3[3][n];
	TProfile2D *tpsinPsi1[2][3][20], *tpcosPsi1[2][3][20];
	TProfile2D *tpsinPsi2[3][n][4], *tpcosPsi2[3][n][4], *tpsinPsi3[3][n][4], *tpcosPsi3[3][n][4];

  //histogram for Q-vectors and event planes with eta-gap and without error
  TH1D *h_Qx1[2][3][c], *h_Qy1[2][3][c]; 
  TH1D *h_Qx2[3][n][c], *h_Qy2[3][n][c], *h_Qx3[3][n][c], *h_Qy3[3][n][c];
  TH1D *h_Psi1[2][3][c],*h_Psi2[3][n][c], *h_Psi3[3][n][c];
  //histogram for check resolution
  TH1D *h_sinPsi2westTPCE[n][c], *h_cosPsi2westTPCE[n][c], *h_sinPsi3westTPCE[n][c], *h_cosPsi3westTPCE[n][c];

    if( strncmp(mode, "QA",2) != 0 ) {
	    // loop by direction
	    for(Int_t l = 0; l < 3; l++) {
	      //loop by cent                                                                                                  
	      for(Int_t j = 0; j < c; j++) {
	        //loop by eta-gap + BBC and ZDC
	        for(Int_t i = 0; i < n; i++) {

	          //historam of Q-vectors for eta-gap 
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
	          h_Qx1[det][l][j] = new TH1D(Form("h_Qx1%s%scent%i",detector[det+1],direction[l],j),
	            Form("Q_{x} for #psi_{1} %s %s cent %i;Q_{x}",detector[det+1],direction[l],j),300,-1.5,1.5);
	          h_Qy1[det][l][j] = new TH1D(Form("h_Qy1%s%scent%i",detector[det+1],direction[l],j),
	            Form("Q_{y} for #psi_{1} %s %s cent %i;Q_{y}",detector[det+1],direction[l],j),300,-1.5,1.5);

	          h_Psi1[det][l][j] = new TH1D(Form("h_Psi1%s%scent%i",detector[det+1],direction[l],j),
	              Form("#psi_{1} %s %s cent %i;#psi_{1}",detector[det+1],direction[l],j),200,-0.1,6.35);
	        }

	      }// loop by cent
	    }// loop by direstion
	    
	}// if( strncmp(mode, "QA",2) != 0 )

  	//***************HISTOGRMS AND PROFILES FOR QA MODE***************//
  	//*****Event*****//
  	TH1D *hRefMult, *hVtxZ, *hRefMult2, *hGRefMult, *hNumberOfPrimaries, *hNumberOfGlobals;
  	TH1D *hCent9, *hCent16, *hBTofHit, *hBTofMatched, *hBemcMatched, *hRanking;
  	TH1D *hTransSphericity, *hTransSphericity2, *hNumberOfVertices;
  	TH2D *hVtxXvsY, *hVpdVzDiffVsVz, *hBTofTrayMultVsRefMult, *hBTofMatchedVsRefMult;
  	TProfile *hEventProfile[8];  
  	//*****Track*****// 
  	TH1D *hGlobalPtot, *hPrimaryPtot, *hGlobalPt, *hPrimaryPt, *hNHits, *hNHitsRatio;
  	TH1D *hChi2, *hDca, *hPhi, *hEta, *hEtaG, *hDedx;
  	TH2D *hDcaVsPt, *hPtVsEta, *hPrimaryPhiVsPt[2], *hDedxVsPt; 
  	TH2D *hNSigmaPionVsPt, *hNSigmaElectronVsPt, *hNSigmaKaonVsPt, *hNSigmaProtonVsPt;
  	TH2D *hDedxVsPtPID[4];
  	TProfile *hTrackProfile[6];
  	//*****TofPidTrait*****//
  	TH1D *hTofBeta, *hMassSqr;
  	TH2D *hInvBetaVsPt, *hMassSqrVsPt, *hDedxVsMassSqr[2]; 
  	TH2D *hInvBetaDiffElectronVsPt, *hInvBetaDiffPionVsPt, *hInvBetaDiffKaonVsPt, *hInvBetaDiffProtonVsPt;
  	
  	if( strncmp(mode, "QA",2) == 0 ) {
		//*****Event*****//
	   	hRefMult = new TH1D("hRefMult", "Reference multiplicity;RefMult;Entries",
	                              600, -0.5, 599.5);
	    hVtxXvsY = new TH2D("hVtxXvsY", "hVtxXvsY;x (cm);y (cm)",
	                              200,-10.,10.,200,-10.,10.);
	   	hVtxZ = new TH1D("hVtxZ","hVtxZ;z (cm); Entries",
	                           140, -70., 70.);
	   	hRefMult2 = new TH1D("hRefMult2","Reference multiplicity in |#eta|<1;RefMult2;Entries",
	                               600, -0.5, 599.5);
	   	hGRefMult = new TH1D("hGRefMult","Reference multiplicity of global tracks;gRefMult;Entries",
	                               800, -0.5, 799.5);
	   	hNumberOfPrimaries = new TH1D("hNumberOfPrimaries","Number of primary tracks;Number of primary tracks;Entries",
	                                        600, -0.5, 599.5);
	   	hNumberOfGlobals = new TH1D("hNumberOfGlobals","Number of global tracks;Number of global tracks;Entries",
	                                      600, -0.5, 599.5);
	   	hCent9 = new TH1D("hCent9","Centralitity;Cent9;Entries",
	                            13, -1.5, 11.5);
	   	hCent16 = new TH1D("hCent16","Centralitity;Cent16;Entries",
	                            19, -1.5, 17.5);
	   	hBTofHit = new TH1D("hBTofHit","Number of hits in TOF;bTofTrayMult;Entries",
	                              600, -0.5, 599.5);
	   	hBTofMatched = new TH1D("hBTofMatched","Number of TOF-matched tracks;bTofMatched;Entries",
	                                  400, -0.5, 399.5);
	   	hBemcMatched = new TH1D("hBemcMatched","Number of BEMC-matched tracks;bEmcMatched;Entries",
	                                  400, -0.5, 399.5);
	   	hRanking = new TH1D("hRanking","Primary vertex ranking;Primary vertex ranking;Entries",
	                              21, -10.5, 10.5);
	    hVpdVzDiffVsVz = new TH2D("hVpdVzDiffVsVz","v_{z}(TPC) - v_{z}(VPD) vs. v_{z}(TPC);v_{z}(TPC);v_{z}(TPC) - v_{z}(VPD)",
	                                    280, -70., 70., 80, -20., 20.);
	    hBTofTrayMultVsRefMult = new TH2D("hBTofTrayMultVsRefMult","TOF tray multiplicity vs. refMult;refMult;bTofTrayMult",
	                                            600, -0.5, 599.5, 600, -0.5, 599.5);
	    hBTofMatchedVsRefMult = new TH2D("hBTofMatchedVsRefMult","TOF-matched tracks vs. refMult;refMult;TOF-matched",
	                                            600, -0.5, 599.5, 400, -0.5, 399.5);
	   	hTransSphericity = new TH1D("hTransSphericity","Transverse sphericity;Sphericity;Entries",
	                                      10, 0., 1.);
	   	hTransSphericity2 = new TH1D("hTransSphericity2","Transverse sphericity in |#eta|<1;Sphericity;Entries",
	                                       10, 0., 1.);
	   	hNumberOfVertices = new TH1D("hNumberOfVertices","Number of primary vertices;Number of primary vertices;Entries",
	                                       15, -0.5, 14.5);
	    
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

	    //*****Track*****// 
	    hGlobalPtot = new TH1D("hGlobalPtot","Global track momentum;p (GeV/c);Entries",
	                                 200, 0., 2. );
	    hPrimaryPtot = new TH1D("hPrimaryPtot","Primary track momentum;p (GeV/c);Entries",
	                                  200, 0., 2. );
	    hGlobalPt = new TH1D("hGlobalPt","Global track transverse momentum;p_{T} (GeV/c)",
	                                200, 0., 2.);
	    hPrimaryPt = new TH1D("hPrimaryPt","Primary track transverse momentum;p_{T} (GeV/c)",
	                                200, 0., 2.);
	    hNHits = new TH1D("hNHits","Number of hits;nHits;Entries", 80, -0.5, 79.5);
	    hNHitsRatio = new TH1D("hNHitsRatio","nHitsFit to nHitsPoss ratio;nHitsFit/nHitsRatio;Entries",
	                                 10, 0., 1. );
	    hChi2 = new TH1D("hChi2","#chi^{2} of the track;#chi^{2};Entries",
	                           200, 0., 20.);
	    hDca = new TH1D("hDca","DCA to primary vertex;DCA (cm);Entries",
	                          100, 0., 10.);
	    hDcaVsPt = new TH2D("hDcaVsPt","charge*p_{T} vs. DCA;charge * p_{T} (GeV/c);DCA (cm)",
	                              840, -2.1, 2.1, 100, 0., 10.);
	    hPhi = new TH1D("hPhi","Azimuthal angle distribution;#phi;Entries",
	                          640, -3.2, 3.2 );
	    hEta = new TH1D("hEta","Track pseudorapidity;#eta;Entries", 220, -1.1, 1.1 );
	    hEtaG = new TH1D("hEtaG","Track pseudorapidity of global track;#eta;Entires", 220, -1., 1. );
	    hPtVsEta = new TH2D("hPtVsEta","p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)",
	                              220, -1.1, 1.1, 80, 0.05, 2.05);

	    for(Int_t i=0; i<2; i++) {
	      hPrimaryPhiVsPt[i] = new TH2D(Form("hPrimaryPhiVsPt_%d",i),
	           Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
	           300, 0., 3., 630, -3.15, 3.15 );
	    }

	    hDedx = new TH1D("hDedx","dE/dx;dE/dx (keV/cm);Entries",
	                           125, 0., 12.5);
	    hDedxVsPt = new TH2D("hDedxVsPt", "dE/dx vs. charge*p_{T};charge * p_{T} (GeV/c);dE/dx (keV/cm)",
	                               840, -2.1, 2.1, 600, 0., 12.);
	    hNSigmaPionVsPt = new TH2D("hNSigmaPionVsPt","n#sigma(#pi) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(#pi)",
	                                     840, -2.1, 2.1, 200, -10., 10.);
	    hNSigmaElectronVsPt = new TH2D("hNSigmaElectronVsPt","n#sigma(e) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(e)",
	                                         840, -2.1, 2.1, 200, -10., 10.);
	    hNSigmaKaonVsPt = new TH2D("hNSigmaKaonVsPt","n#sigma(K) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(K)",
	                                     840, -2.1, 2.1, 200, -10., 10.);
	    hNSigmaProtonVsPt = new TH2D("hNSigmaProtonVsPt","n#sigma(p) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(p)",
	                                       840, -2.1, 2.1, 200, -10., 10.);

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

	    //*****TofPidTrait*****//
	    hTofBeta = new TH1D("hTofBeta","BTofPidTraits #beta;#beta",
	                              2000, 0., 2.);
	    hInvBetaVsPt = new TH2D("hInvBetaVsPt","1/#beta vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta",
	                                  840, -2.1, 2.1, 200, 0.8, 2.8);
	    hMassSqr = new TH1D("hMassSqr","m^{2};m^{2} (GeV/c^{2})^{2};dN/dm^{2} (entries)",
	                              520, -0.1, 5.1 );
	    hMassSqrVsPt = new TH2D("hMassSqrVsPt","m^{2} vs. charge*p_{T};charge * p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",
	                                  840, -2.1, 2.1, 200, -0.2, 1.8);
	    hDedxVsMassSqr[0] = new TH2D("hDedxVsMassSqr_0","dE/dx vs. mass^{2} charge>0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
	               440, -0.4, 1.8, 250, 0., 12.5 );
	    hDedxVsMassSqr[1] = new TH2D("hDedxVsMassSqr_1","dE/dx vs. mass^{2} charge<0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)",
	               440, -0.4, 1.8, 250, 0., 12.5 );
	    hInvBetaDiffElectronVsPt = new TH2D("hInvBetaDiffElectronVsPt","1/#beta - 1/#beta(electron) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(e)",
	                                              840, -2.1, 2.1, 200, -0.1, 0.1);
	    hInvBetaDiffPionVsPt = new TH2D("hInvBetaDiffPionVsPt","1/#beta - 1/#beta(pion) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(#pi)",
	                                              840, -2.1, 2.1, 200, -0.1, 0.1);
	    hInvBetaDiffKaonVsPt = new TH2D("hInvBetaDiffKaonVsPt","1/#beta - 1/#beta(kaon) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(K)",
	                                              840, -2.1, 2.1, 200, -0.1, 0.1);
	    hInvBetaDiffProtonVsPt = new TH2D("hInvBetaDiffProtonVsPt","1/#beta - 1/#beta(p) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(p)",
	                                            840, -2.1, 2.1, 200, -0.1, 0.1);

  	}// if( strncmp(mode, "QA",2) == 0 )


  	//***************PROFILES FOR RAW MODE***************//
	TProfile2D *tp_Qx1[2][3], *tp_Qy1[2][3];
	TProfile2D *tp_Qx2[3][n], *tp_Qx3[3][n], *tp_Qy2[3][n], *tp_Qy3[3][n];

  	if( strncmp(mode, "raw",3) == 0 ) {

  		// loop by direction 
	    for(Int_t l = 0; l < 3; l++) {
	      // loop by eta-gap of TPC 
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

	      }// loop by eta-gap of TPC 

	      for(Int_t det = 0; det < 2; det++) {

	        tp_Qx1[det][l] = new TProfile2D(Form("tp_Qx1%s%s",detector[det+1],direction[l]),
	        Form("<Q_{x}> for #psi_{1} %s %s;RunID;cent",detector[det+1],direction[l]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

	        tp_Qy1[det][l] = new TProfile2D(Form("tp_Qy1%s%s",detector[det+1],direction[l]),
	        Form("<Q_{y}> for #psi_{1} %s %s;RunID;cent",detector[det+1],direction[l]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);
	      }
	    }// loop by direction

  	}// if( strncmp(mode, "raw",3) == 0 )


  	//***************PROFILES FOR RECENTERING MODE***************//
  	TProfile2D *tp_sinPsi1[2][3][20], *tp_cosPsi1[2][3][20];
		TProfile2D *tp_sinPsi2[3][n][4], *tp_cosPsi2[3][n][4], *tp_sinPsi3[3][n][4], *tp_cosPsi3[3][n][4];

  	if( strncmp(mode, "rec",3) == 0 ) { 

	    for(Int_t l = 0; l < 3; l++) {

	      for(Int_t i = 0; i < n; i++) {
	        for(Int_t j = 0; j < 4; j++) {
	          
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
	          
	          tp_sinPsi1[det][l][j] = new TProfile2D(Form("tp_%isinPsi1%s%s",j+1,detector[det+1],direction[l]),
	          Form("<sin(%i*#psi_{1})> %s %s",j+1,detector[det+1],direction[l]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);

	          tp_cosPsi1[det][l][j] = new TProfile2D(Form("tp_%icosPsi1%s%s",j+1,detector[det+1],direction[l]),
	          Form("<cos(%i*#psi_{1})> %s %s",j+1,detector[det+1],direction[l]),runIdBins ,runIdRange[0],runIdRange[1],9,0,9);
	        }
	      }

	    }// for(Int_t l = 0; l < 2; l++)

	    f1 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/NoRe_%s.root",energy),"READ");
	    for(Int_t dir = 0; dir < 3; dir++) {
	      for(Int_t i = 0; i < n; i++) {
	        tpQx2[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx2%s%s",direction[dir],ngap[i]) );
	        tpQx3[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx3%s%s",direction[dir],ngap[i]) );
	        tpQy2[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy2%s%s",direction[dir],ngap[i]) );
	        tpQy3[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy3%s%s",direction[dir],ngap[i]) );
	      }
	      for(Int_t det = 0; det < 2; det++) {
	        tpQx1[det][dir] = (TProfile2D*) f1 -> Get( Form("tp_Qx1%s%s",detector[det+1],direction[dir]) );
	        tpQy1[det][dir] = (TProfile2D*) f1 -> Get( Form("tp_Qy1%s%s",detector[det+1],direction[dir]) );
	      }
	    }

	}// if( strncmp(mode, "rec",3) == 0 )


	//***************PROFILES FOR FLATTENING AND FLOW MODE***************//
	Int_t n_sys=3; 
	// 0 - VtZ 
	// 1 - nHits
	// 2 - DCA
	//*****Resolution^2 for Psi1*****//
	TProfile *tp_SqRes1[2];
	//*****Resolution^2 for Psi2 and Psi3*****// 
	TProfile *tp_SqRes2[n], *tp_SqRes3[n];
	//*****Flows*****//
	TProfile *tp_v1ew[2][9], *tp_v1comb[2][9];
	TProfile2D *tp_v2hadrTPC[n][2], *tp_v3hadrTPC[n][2];
	TProfile2D *tp_v2hadrTPCsys[n][2][3], *tp_v3hadrTPCsys[n][2][3];

	TProfile2D *tp_v2pidTPC[n][2][2][3], *tp_v3pidTPC[n][2][2][3];

	TProfile2D *tp_meanPt_hadrons[n];
	TProfile2D *tp_meanPt_hadrons_sys[n][3];
	TProfile2D *tp_meanPt_PID[n][2][3];

	
	//*****Flows systematic*****//
	
	TProfile2D *tp_v2pidTPCsys[n][2][2][3][3]; 
	TProfile2D *tp_v3pidTPCsys[n][2][2][3][3];

	//*****Mean pt for PID and hadrons*****//
	
	TProfile2D *tp_meanPt_PID_sys[n][2][3][3];

	

	if( strncmp(mode, "flow",4) == 0 ) {

		for(Int_t det = 0; det < 2; det++) {
	  	tp_SqRes1[det] = new TProfile(Form("tp_SqRes1%s",detector[det+1]),
	  																Form("Resolution^{2} for v_{1} %s",detector[det+1]),9,0,9);
	  }

	  for(Int_t i = 0; i < n; i++) {
	  	
	    tp_SqRes2[i] = new TProfile(Form("tp_SqRes2%s",ngap[i]),
	    														Form("Resolution^{2} for v_{2} %s",gap[i]),9,0,9);
	    tp_SqRes3[i] = new TProfile(Form("tp_SqRes3%s",ngap[i]),
	    														Form("Resolution^{2} for v_{3} %s",gap[i]),9,0,9);
	  	
		}

	  for(Int_t det = 0; det < 2; det++) {
	  	for(Int_t j = 0; j < c; j++) {
	  		tp_v1ew[det][j] = new TProfile( Form( "tp_v1ew%scent%i",detector[det+1],j), 
	                                   		Form( "v_{1} of pseudorapidity cent %i ew by %s",j,detector[det+1]),20,-1,1);
	  		tp_v1comb[det][j] = new TProfile( Form( "tp_v1comb%scent%i",detector[det+1],j), 
	                                   			Form( "v_{1} of pseudorapidity cent %i comb by %s",j,detector[det+1]),20,-1,1);
	    }
	  }

	  for(Int_t i = 0; i < n; i++) {

	  	tp_meanPt_hadrons[i] = new TProfile2D(Form("tp_meanPt_hadrons%s",ngap[i]),
	    Form("Mean p_{t} for bins v_{2} %s; bin; p_{t} [GeV/c]",gap[i]),30,0.2,3.2,9,0,9);

	    for(Int_t nsys = 0; nsys < 3; nsys++) {
	    	tp_meanPt_hadrons_sys[i][nsys] = new TProfile2D(Form("tp_meanPt_hadrons_sys%s%s",ngap[i],systematic[nsys]),
	      Form("Mean p_{t} for bins v_{2} %s Systematic:%s; bin; p_{t} [GeV/c]",gap[i],systematic[nsys]),30,0.2,3.2,9,0,9);
	    }

	    Int_t p = 0; 
	    for(Int_t mark = 0; mark < 2; mark++) {
	    	for(Int_t part = 0; part < 3; part++) {

	    		tp_meanPt_PID[i][mark][part] = new TProfile2D(
	    		Form("tp_meanPt_PID%s%s%s",ngap[i],sign[mark],particles[part]),
	    		Form("Mean p_{t} for %s of p_{t} and cent %s;p_{t} [GeV/c];cent",partLateX[p],ngap[i]),
	    					30,0.2,3.2,9,0,9);

	    		for(Int_t nsys = 0; nsys < 3; nsys++) {

	    			tp_meanPt_PID_sys[i][mark][part][nsys] = new TProfile2D(
	    			Form("tp_meanPt_PID_sys%s%s%s%s",ngap[i],sign[mark],particles[part],systematic[nsys]),
	    			Form("Mean p_{t} for %s of p_{t} and cent %s %s;p_{t} [GeV/c];cent",partLateX[p],ngap[i],systematic[nsys]),
	    						30,0.2,3.2,9,0,9);

	    		}
	    		p++;
	    	}
	    }

	  	for(Int_t dir = 0; dir < 2; dir++) {

				tp_v2hadrTPC[i][dir] = new TProfile2D(Form("tp_v2hadrTPC%s%s",ngap[i],resol[dir]),
	                        Form("v_{2} of p_{t} and cent %s %s by TPC;p_{t} [GeV/c];cent",gap[i],resol[dir]),30,0.2,3.2,9,0,9);
	      tp_v3hadrTPC[i][dir] = new TProfile2D(Form("tp_v3hadrTPC%s%s",ngap[i],resol[dir]),
	                        Form("v_{3} of p_{t} and cent %s %s by TPC;p_{t} [GeV/c];cent",gap[i],resol[dir]),30,0.2,3.2,9,0,9);

	      for(Int_t nsys = 0; nsys < 3; nsys++) {

	        tp_v2hadrTPCsys[i][dir][nsys] = new TProfile2D(
	        Form("tp_v2hadrTPCsys%s%s%s",ngap[i],resol[dir],systematic[nsys]),
	        Form("v_{2} of p_{t} and cent %s %s by TPC Systematic: %s;p_{t} [GeV/c];cent",gap[i],resol[dir],systematic[nsys]),
	        			30,0.2,3.2,9,0,9);
	            
	        tp_v3hadrTPCsys[i][dir][nsys] = new TProfile2D(
	        Form("tp_v3hadrTPCsys%s%s%s",ngap[i],resol[dir],systematic[nsys]),
	        Form("v_{3} of p_{t} and cent %s %s by TPC Systematic: %s;p_{t} [GeV/c];cent",gap[i],resol[dir],systematic[nsys]),
	        			30,0.2,3.2,9,0,9);
	        
	      }

	      p = 0; 
	      for(Int_t mark = 0; mark < 2; mark++) {
	      	for(Int_t part = 0; part < 3; part++) {

	      		tp_v2pidTPC[i][dir][mark][part] = new TProfile2D(
	      		Form("tp_v2pidTPC%s%s%s%s",ngap[i],resol[dir],sign[mark],particles[part]), 
	          Form("v_{2} for %s of p_{t} and cent %s by TPC;p_{t} [GeV/c];cent",partLateX[p],gap[i]),
	          			30,0.2,3.2,9,0,9);

	          tp_v3pidTPC[i][dir][mark][part] = new TProfile2D(
	          Form("tp_v3pidTPC%s%s%s%s",ngap[i],resol[dir],sign[mark],particles[part]), 
	          Form("v_{3} for %s of p_{t} and cent %s by TPC;p_{t} [GeV/c];cent",partLateX[p],gap[i]),
	          			30,0.2,3.2,9,0,9);

	          for(Int_t nsys = 0; nsys < 3; nsys++) {

	            tp_v2pidTPCsys[i][dir][mark][part][nsys] = new TProfile2D( 
	            Form("tp_v2pidTPCsys%s%s%s%s%s",ngap[i],resol[dir],sign[mark],particles[part],systematic[nsys]), 
	            Form("v_{2} for %s of p_{t} and cent %s by TPC Systematic: %s;p_{t} [GeV/c];cent",partLateX[p],gap[i],systematic[nsys]),
	              		30,0.2,3.2,9,0,9);
	                
	            tp_v3pidTPCsys[i][dir][mark][part][nsys] = new TProfile2D( 
	            Form("tp_v3pidTPCsys%s%s%s%s%s",ngap[i],resol[dir],sign[mark],particles[part],systematic[nsys]), 
	            Form("v_{3} for %s of p_{t} and cent %s by TPC Systematic: %s;p_{t} [GeV/c];cent",partLateX[p],gap[i],systematic[nsys]),
	              		30,0.2,3.2,9,0,9);

	          }
	        	p++;
	      	} 
	      }

	    }// for(Int_t dir = 0; dir < 2; dir++)

	  }// for(Int_t i = 0; i < n; i++)

	  	f1 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/NoRe_%s.root",energy),"READ");
	    for(Int_t dir = 0; dir < 3; dir++) {
	      for(Int_t i = 0; i < n; i++) {
	        tpQx2[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx2%s%s",direction[dir],ngap[i]) );
	        tpQx3[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qx3%s%s",direction[dir],ngap[i]) );
	        tpQy2[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy2%s%s",direction[dir],ngap[i]) );
	        tpQy3[dir][i] = (TProfile2D*) f1 -> Get( Form("tp_Qy3%s%s",direction[dir],ngap[i]) );
	      }
	      for(Int_t det = 0; det < 2; det++) {
	        tpQx1[det][dir] = (TProfile2D*) f1 -> Get( Form("tp_Qx1%s%s",detector[det+1],direction[dir]) );
	        tpQy1[det][dir] = (TProfile2D*) f1 -> Get( Form("tp_Qy1%s%s",detector[det+1],direction[dir]) );
	      }
	    }

	    f2 = new TFile(Form("/mnt/pool/1/aspovarov/basov/test/OUT/Re_%s.root",energy),"READ");
	    for(Int_t l = 0; l < 3; l++) {
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
	          tpsinPsi1[det][l][j] = (TProfile2D*) f2 -> Get( Form("tp_%isinPsi1%s%s",j+1,detector[det+1],direction[l]) );
	          tpcosPsi1[det][l][j] = (TProfile2D*) f2 -> Get( Form("tp_%icosPsi1%s%s",j+1,detector[det+1],direction[l]) );
	        }
	      }
	    }

	}// if( strncmp(mode, "flow",4)==0 )


	/*////////////////////////////////////////////////////////////////////////////////////////*/
 /*________________________________START OF EVENT LOOP_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
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

    TVector3 pVtx = event->primaryVertex();

   	// Track analysis
    Int_t nTracks = dst->numberOfTracks();
            
    //RunQA
    if( useRunQA == true && std::find(badRuns.begin(), badRuns.end(), event -> runId()) != badRuns.end() ) continue;

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*______________________________________QA MODE___________________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
	if( strncmp(mode, "QA",2)==0 ){

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
	} //if( strncmp(mode, "QA",2)==0 )

  /*/////////////////////////////////////////////////////////////////////////////////////////*/
 /*______________________________________RAW MODE___________________________________________*/
/*/////////////////////////////////////////////////////////////////////////////////////////*/
    if( strncmp(mode, "raw",3) == 0 ) {

    	for(Int_t dir = 0; dir < 3; dir++) {
    		for(Int_t det = 0; det < 2; det++) {
    			Q1vec[det][dir].Set(0.,0.);
    		}
  			
  			for(Int_t i = 0; i < n; i++) {
  				Q2vecTPC[dir][i].Set(0.,0.);
  				Q3vecTPC[dir][i].Set(0.,0.);
  			}
  		}

      for(Int_t dir = 0; dir < 2; dir++) {
      	Q1vec[0][dir] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[1][dir] = CalculateZDCQVec(event, dir);
      }
      Q1vec[0][2] = Q1vec[0][1] - Q1vec[0][0];
      Q1vec[1][2] = Q1vec[1][1] - Q1vec[1][0];
      
      CalculateTPCQVec(Q2vecTPC, Q3vecTPC, dst, n, DCA_EVENT);
            
    	cent = event -> cent9();
    	RunID = event -> runId();

      for(Int_t dir = 0; dir < 3; dir ++) {
        for(Int_t i = 0; i < n; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            h_Qx2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].X() );
            h_Qy2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].Y() );
            h_Qx3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].X() );
            h_Qy3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].Y() );

        	  h_Psi2[dir][i][cent] -> Fill( 1.0/2.0 * Q2vecTPC[dir][i].Phi()  );
        	  h_Psi3[dir][i][cent] -> Fill( 1.0/3.0 * Q3vecTPC[dir][i].Phi()  ); 

   		      tp_Qx2[dir][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vecTPC[dir][i].X() );
   		      tp_Qy2[dir][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q2vecTPC[dir][i].Y() );
   		      tp_Qx3[dir][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vecTPC[dir][i].X() );
   		      tp_Qy3[dir][i] -> Fill( (Double_t)RunID, (Double_t)cent, Q3vecTPC[dir][i].Y() );
          }
        } //for(Int_t i = 0; i < n; i++)  

        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[det][dir].Mod() != 0. ) {
            h_Qx1[det][dir][cent] -> Fill( Q1vec[det][dir].X() );
            h_Qy1[det][dir][cent] -> Fill( Q1vec[det][dir].Y() );

            h_Psi1[det][dir][cent] -> Fill( Q1vec[det][dir].Phi() );

            tp_Qx1[det][dir] -> Fill( (Double_t)RunID, (Double_t)cent, Q1vec[det][dir].X() );
            tp_Qy1[det][dir] -> Fill( (Double_t)RunID, (Double_t)cent, Q1vec[det][dir].Y() );
          }
        }

      }// for(Int_t dir = 0; dir < 3; dir++)

		}// if( strncmp(mode, "raw",3) == 0 )


  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________________RECENTERING MODE_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
		if( strncmp(mode, "rec",3)==0 ) {

			

	    for(Int_t dir = 0; dir < 3; dir++) {
    		for(Int_t det = 0; det < 2; det++) {
    			Q1vec[det][dir].Set(0.,0.);
    		}
  			
  			for(Int_t i = 0; i < n; i++) {
  				Q2vecTPC[dir][i].Set(0.,0.);
  				Q3vecTPC[dir][i].Set(0.,0.);
  			}
  		}

      for(Int_t dir = 0; dir < 2; dir++) {
      	Q1vec[0][dir] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[1][dir] = CalculateZDCQVec(event, dir);
      }
      Q1vec[0][2] = Q1vec[0][1] - Q1vec[0][0];
      Q1vec[1][2] = Q1vec[1][1] - Q1vec[1][0];

      CalculateTPCQVec(Q2vecTPC, Q3vecTPC, dst, n, DCA_EVENT);
             
      cent = (Double_t)(event -> cent9());
      RunID = event -> runId();

      Recentering(Q1vec,Q2vecTPC,Q3vecTPC,RunID,cent,tpQx1,tpQy1,tpQx2,tpQy2,tpQx3,tpQy3);  

       for(Int_t dir = 0; dir < 3; dir ++) {
        for(Int_t i = 0; i < n; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            h_Qx2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].X() );
            h_Qy2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].Y() );
            h_Qx3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].X() );
            h_Qy3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].Y() );

        	  h_Psi2[dir][i][cent] -> Fill( 1.0/2.0 * Q2vecTPC[dir][i].Phi()  );
        	  h_Psi3[dir][i][cent] -> Fill( 1.0/3.0 * Q3vecTPC[dir][i].Phi()  ); 

        	  for(Int_t j =0; j < 4; j++) {
              tp_sinPsi2[dir][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q2vecTPC[dir][i].Phi() ) );
              tp_cosPsi2[dir][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q2vecTPC[dir][i].Phi() ) );
              tp_sinPsi3[dir][i][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q3vecTPC[dir][i].Phi() ) ); 
              tp_cosPsi3[dir][i][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q3vecTPC[dir][i].Phi() ) );
            }
   		      
          }
        } //for(Int_t i = 0; i < n; i++)  

        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[det][dir].Mod() != 0. ) {
            h_Qx1[det][dir][cent] -> Fill( Q1vec[det][dir].X() );
            h_Qy1[det][dir][cent] -> Fill( Q1vec[det][dir].Y() );

            h_Psi1[det][dir][cent] -> Fill( Q1vec[det][dir].Phi() );

            for(Int_t j = 0; j < 20; j++) {
              tp_sinPsi1[det][dir][j] -> Fill(RunID, cent, TMath::Sin( (Double_t)(j+1)*Q1vec[det][dir].Phi() ) );
              tp_cosPsi1[det][dir][j] -> Fill(RunID, cent, TMath::Cos( (Double_t)(j+1)*Q1vec[det][dir].Phi() ) ); 
            }
          }
        }

      }// for(Int_t dir = 0; dir < 3; dir++)

	} //if( strncmp(mode, "rec",3)==0 )

	/*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________FLATTENING AND FLOW MODE_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/
		if( strncmp(mode, "flow",4)==0 ){

			
	    
      Double_t sinPsi1 = 0., cosPsi1 = 0., sinPsi2 = 0., cosPsi2 = 0., sinPsi3 = 0., cosPsi3 = 0.;
      Double_t dPsi1 = 0., dPsi2 = 0., dPsi3 = 0.;

      for(Int_t dir = 0; dir < 3; dir++) {
    		for(Int_t det = 0; det < 2; det++) {
    			Q1vec[det][dir].Set(0.,0.);
    		}
  			
  			for(Int_t i = 0; i < n; i++) {
  				Q2vecTPC[dir][i].Set(0.,0.);
  				Q3vecTPC[dir][i].Set(0.,0.);
  			}
  		}

      for(Int_t dir = 0; dir < 2; dir++) {
      	Q1vec[0][dir] = CalculateBBCQVec(event, dir, 1.0);
        Q1vec[1][dir] = CalculateZDCQVec(event, dir);
      }
      Q1vec[0][2] = Q1vec[0][1] - Q1vec[0][0];
      Q1vec[1][2] = Q1vec[1][1] - Q1vec[1][0];

      CalculateTPCQVec(Q2vecTPC, Q3vecTPC, dst, n, DCA_EVENT);

      RunID = event -> runId();
      cent = event -> cent9();

      Recentering(Q1vec,Q2vecTPC,Q3vecTPC,RunID,cent,tpQx1,tpQy1,tpQx2,tpQy2,tpQx3,tpQy3);

      for(Int_t dir = 0; dir < 3; dir ++) {
        for(Int_t i = 0; i < n; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            h_Qx2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].X() );
            h_Qy2[dir][i][cent] -> Fill( Q2vecTPC[dir][i].Y() );
            h_Qx3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].X() );
            h_Qy3[dir][i][cent] -> Fill( Q3vecTPC[dir][i].Y() );

            Psi2TPC[dir][i] = 1.0/2.0 * Q2vecTPC[dir][i].Phi();
            Psi3TPC[dir][i] = 1.0/3.0 * Q3vecTPC[dir][i].Phi();

          }
        } //for(Int_t i = 0; i < n; i++)  

        for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[det][dir].Mod() != 0. ) {
            h_Qx1[det][dir][cent] -> Fill( Q1vec[det][dir].X() );
            h_Qy1[det][dir][cent] -> Fill( Q1vec[det][dir].Y() );

            Psi1[det][dir] = Q1vec[det][dir].Phi();

          }
        }

      }// for(Int_t dir = 0; dir < 3; dir++)

      // Flattening stage 
      for(Int_t dir = 0; dir < 3; dir++) {

      	for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[det][dir].Mod() != 0. ) {
            for(Int_t j = 0; j < 20; j++) { 
               
            	sinPsi1 = tpsinPsi1[det][dir][j]->GetBinContent( tpsinPsi1[det][dir][j]->FindBin(RunID, cent) );
              cosPsi1 = tpcosPsi1[det][dir][j]->GetBinContent( tpcosPsi1[det][dir][j]->FindBin(RunID, cent) );
            
              dPsi1 += -2.0*( sinPsi1 * TMath::Cos( (Double_t)(j+1)*Psi1[det][dir] ) )/( (Double_t)(j+1) ) 
                       +2.0*( cosPsi1 * TMath::Sin( (Double_t)(j+1)*Psi1[det][dir] ) )/( (Double_t)(j+1) );
            }
            Psi1[det][dir] += dPsi1;
            dPsi1 = 0.;
          } 
        }//for(Int_t det = 0; det < 2; det++)

        for(Int_t i = 0; i < n; i++) {
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            for(Int_t k = 0; k < 4; k++) {

              sinPsi2 = tpsinPsi2[dir][i][k]->GetBinContent( tpsinPsi2[dir][i][k] -> FindBin( RunID, cent ) );
              cosPsi2 = tpcosPsi2[dir][i][k]->GetBinContent( tpcosPsi2[dir][i][k] -> FindBin( RunID, cent ) );
              sinPsi3 = tpsinPsi3[dir][i][k]->GetBinContent( tpsinPsi3[dir][i][k] -> FindBin( RunID, cent ) );
              cosPsi3 = tpcosPsi3[dir][i][k]->GetBinContent( tpcosPsi3[dir][i][k] -> FindBin( RunID, cent ) );

              dPsi2 += -2.0*( sinPsi2 * TMath::Cos( (Double_t)(k+1)*2.0*Psi2TPC[dir][i] ) )/( 2.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi2 * TMath::Sin( (Double_t)(k+1)*2.0*Psi2TPC[dir][i] ) )/( 2.0*(Double_t)(k+1) );

              dPsi3 += -2.0*( sinPsi3 * TMath::Cos( (Double_t)(k+1)*3.0*Psi3TPC[dir][i] ) )/( 3.0*(Double_t)(k+1) ) 
                       +2.0*( cosPsi3 * TMath::Sin( (Double_t)(k+1)*3.0*Psi3TPC[dir][i] ) )/( 3.0*(Double_t)(k+1) );

            } //for(Int_t k = 0; k < 4; k++)
            Psi2TPC[dir][i] += dPsi2;
            Psi3TPC[dir][i] += dPsi3;
            dPsi2 = dPsi3 = 0.;
          } 
        }// for(Int_t i = 0; i < n; i++)

      }// for(Int_t dir = 0; dir < 3; dir++) {

      for(Int_t dir = 0; dir < 3; dir++) {  

      	for(Int_t det = 0; det < 2; det++) {
          if( Q1vec[det][dir].Mod() != 0. ) {
            h_Psi1[det][dir][cent] -> Fill( Psi1[det][dir] );
          }
        }

        for(Int_t i = 0; i < n; i++) {  
       
          if( Q2vecTPC[dir][i].Mod() != 0. && Q3vecTPC[dir][i].Mod() != 0.) {
            h_Psi2[dir][i][cent] -> Fill( Psi2TPC[dir][i]  );
            h_Psi3[dir][i][cent] -> Fill( Psi3TPC[dir][i]  );
          }
          
        }// for(Int_t i = 0; i < n; i++)

      }// for(Int_t dir = 0; dir < 3; dir++)  

      for(Int_t det = 0; det < 2; det++) {
        if(Q1vec[det][0].Mod() != 0. && Q1vec[det][1].Mod() != 0.) {
          tp_SqRes1[det] -> Fill(cent, TMath::Cos( (Psi1[det][1] - Psi1[det][0]) ) );
        }
      }  

      for(Int_t i = 0; i < n; i++) {
        if(Q2vecTPC[0][i].Mod() != 0. && Q3vecTPC[0][i].Mod() != 0. && Q2vecTPC[1][i].Mod() != 0. && Q3vecTPC[1][i].Mod() != 0.) {
          tp_SqRes2[i] -> Fill(cent, TMath::Cos( 2*(Psi2TPC[1][i] - Psi2TPC[0][i]) ) );
          tp_SqRes3[i] -> Fill(cent, TMath::Cos( 3*(Psi3TPC[1][i] - Psi3TPC[0][i]) ) );
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

        SqM = femtoTrack -> massSqr();
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

        nSigma[0] = femtoTrack -> nSigmaPion();
        nSigma[1] = femtoTrack -> nSigmaKaon();
        nSigma[2] = femtoTrack -> nSigmaProton();
        
        Double_t v1 = 0.,v2 = 0., v3 = 0.;

        for(Int_t det = 0; det < 2; det++) {

        	v1 = TMath::Cos( Phi - Psi1[det][2] );
          tp_v1comb[det][cent] -> Fill( femtoTrack -> eta(), v1);

	        if( femtoTrack -> eta() < 0 ) {
	        	v1 = TMath::Cos( Phi - Psi1[det][1] );
	          tp_v1ew[det][cent] -> Fill( femtoTrack -> eta(), v1);
	        }
	        if( femtoTrack -> eta() > 0 ) {
	          v1 = TMath::Cos( Phi - Psi1[det][0] );
	          tp_v1ew[det][cent] -> Fill( femtoTrack -> eta(), v1);
	        }
  
        }

        for(Int_t eta = 0; eta < n; eta++) {

        	if( femtoTrack -> eta() > abs( etagap[eta] ) ) {
          	v2 = TMath::Cos( 2.0*(Phi - Psi2TPC[2][eta]) );
          	v3 = TMath::Cos( 3.0*(Phi - Psi3TPC[2][eta]) );

          	tp_v2hadrTPC[eta][1] -> Fill(pt, (Double_t)cent, v2, w );
	        	tp_v3hadrTPC[eta][1] -> Fill(pt, (Double_t)cent, v3, w );
	        	tp_meanPt_hadrons[eta] -> Fill(pt,(Double_t)cent,pt);

	          if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
		        	tp_v2hadrTPCsys[eta][1][0] -> Fill(pt, (Double_t)cent, v2, w );
		        	tp_v3hadrTPCsys[eta][1][0] -> Fill(pt, (Double_t)cent, v3, w );
		        	tp_meanPt_hadrons_sys[eta][0] -> Fill(pt,(Double_t)cent,pt);
		        }
						if( femtoTrack -> nHits() > sys_nHits)  {
							tp_v2hadrTPCsys[eta][1][1] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][1][1] -> Fill(pt, (Double_t)cent, v3, w );
							tp_meanPt_hadrons_sys[eta][1] -> Fill(pt,(Double_t)cent,pt);
						}
						if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA) {
							tp_v2hadrTPCsys[eta][1][2] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][1][2] -> Fill(pt, (Double_t)cent, v3, w );
							tp_meanPt_hadrons_sys[eta][2] -> Fill(pt,(Double_t)cent,pt);
						}

						for(Int_t part = 0; part < 3; part++) {
							if( TMath::Abs( nSigma[part] ) < nSigm && SqM > SqMdown[part] && SqM < SqMup[part]) {
								if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
									tp_v2pidTPC[eta][1][q_charge][part] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPC[eta][1][q_charge][part] -> Fill(pt, (Double_t)cent, v3, w );
									tp_meanPt_PID[eta][q_charge][part] -> Fill(pt,(Double_t)cent,pt);
								}
								if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
									tp_v2pidTPCsys[eta][1][q_charge][part][0] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][1][q_charge][part][0] -> Fill(pt, (Double_t)cent, v3, w );
									tp_meanPt_PID_sys[eta][q_charge][part][0] -> Fill(pt,(Double_t)cent,pt);
								}
								if( femtoTrack -> nHits() > sys_nHits) {
									tp_v2pidTPCsys[eta][1][q_charge][part][1] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][1][q_charge][part][1] -> Fill(pt, (Double_t)cent, v3, w );
									tp_meanPt_PID_sys[eta][q_charge][part][1] -> Fill(pt,(Double_t)cent,pt);
								}
								if( femtoTrack -> gDCA(pVtx).Mag() < sys_DCA) {
									tp_v2pidTPCsys[eta][1][q_charge][part][2] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][1][q_charge][part][2] -> Fill(pt, (Double_t)cent, v3, w );
									tp_meanPt_PID_sys[eta][q_charge][part][2] -> Fill(pt,(Double_t)cent,pt);
								}
							}
						}// for(Int_t part = 0; part < 3; part++)
					}// if( abs( femtoTrack -> eta() ) > etagap[eta] 

					if( femtoTrack -> eta() < -etagap[eta] ) {

	        	v2 = TMath::Cos( 2.0*(Phi - Psi2TPC[1][eta]) );
          	v3 = TMath::Cos( 3.0*(Phi - Psi3TPC[1][eta]) );

	        	if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW) {
	        		tp_v2hadrTPC[eta][0] -> Fill(pt, (Double_t)cent, v2, w );
	        		tp_v3hadrTPC[eta][0] -> Fill(pt, (Double_t)cent, v3, w );
	        	}
	        	if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
	        		tp_v2hadrTPCsys[eta][0][0] -> Fill(pt, (Double_t)cent, v2, w );
	        		tp_v3hadrTPCsys[eta][0][0] -> Fill(pt, (Double_t)cent, v3, w );
	        	}
						if( femtoTrack -> nHits() > sys_nHits)  {
							tp_v2hadrTPCsys[eta][0][1] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][0][1] -> Fill(pt, (Double_t)cent, v3, w );
						}
						if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA) {
							tp_v2hadrTPCsys[eta][0][2] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][0][2] -> Fill(pt, (Double_t)cent, v3, w );
						}

						for(Int_t part = 0; part < 3; part++) {
							if( TMath::Abs( nSigma[part] ) < nSigm && SqM > SqMdown[part] && SqM < SqMup[part]) {
								if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
									tp_v2pidTPC[eta][0][q_charge][part] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPC[eta][0][q_charge][part] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
									tp_v2pidTPCsys[eta][0][q_charge][part][0] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][0] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( femtoTrack -> nHits() > sys_nHits) {
									tp_v2pidTPCsys[eta][0][q_charge][part][1] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][1] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA) {
									tp_v2pidTPCsys[eta][0][q_charge][part][2] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][2] -> Fill(pt, (Double_t)cent, v3, w );
								}
							}
						}// for(Int_t part = 0; part < 3; part++)

	        }// if( femtoTrack -> eta() < -etagap[eta] )

	        if( femtoTrack -> eta() > etagap[eta] ) {

	        	v2 = TMath::Cos( 2.0*(Phi - Psi2TPC[0][eta]) );
          	v3 = TMath::Cos( 3.0*(Phi - Psi3TPC[0][eta]) );

          	if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW) {
          		tp_v2hadrTPC[eta][0] -> Fill(pt, (Double_t)cent, v2, w );
	        		tp_v3hadrTPC[eta][0] -> Fill(pt, (Double_t)cent, v3, w );
	        	}
	        	if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
	        		tp_v2hadrTPCsys[eta][0][0] -> Fill(pt, (Double_t)cent, v2, w );
	        		tp_v3hadrTPCsys[eta][0][0] -> Fill(pt, (Double_t)cent, v3, w );
	        	}
						if( femtoTrack -> nHits() > sys_nHits) {
							tp_v2hadrTPCsys[eta][0][1] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][0][1] -> Fill(pt, (Double_t)cent, v3, w );
						}
						if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA) {
							tp_v2hadrTPCsys[eta][0][2] -> Fill(pt, (Double_t)cent, v2, w );
							tp_v3hadrTPCsys[eta][0][2] -> Fill(pt, (Double_t)cent, v3, w );
						}

						for(Int_t part = 0; part < 3; part++) {
							if( TMath::Abs( nSigma[part] ) < nSigm && SqM > SqMdown[part] && SqM < SqMup[part]) {
								if( femtoTrack->gDCA(pVtx).Mag() < DCA_FLOW){
									tp_v2pidTPC[eta][0][q_charge][part] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPC[eta][0][q_charge][part] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( TMath::Abs( pVtx.Z() ) < sys_VtxZ) {
									tp_v2pidTPCsys[eta][0][q_charge][part][0] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][0] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( femtoTrack -> nHits() > sys_nHits) {
									tp_v2pidTPCsys[eta][0][q_charge][part][1] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][1] -> Fill(pt, (Double_t)cent, v3, w );
								}
								if( femtoTrack->gDCA(pVtx).Mag() < sys_DCA) {
									tp_v2pidTPCsys[eta][0][q_charge][part][2] -> Fill(pt, (Double_t)cent, v2, w );
									tp_v3pidTPCsys[eta][0][q_charge][part][2] -> Fill(pt, (Double_t)cent, v3, w );
								}
							}
						}// for(Int_t part = 0; part < 3; part++)
	        		
	        }// if( femtoTrack -> eta() > etagap[eta] ) {

        }// for(Int_t eta = 0; eta < n; eta++)

      }// for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

		}// if(strncmp(mode, "flow", 4))

	}// for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->Write();
  outFile->Close();

  femtoReader->Finish();
  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << std::endl;

}// void FemtoDstAnalyzer()




	/*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

//********************CHECK EVENT ON GOOD********************//
Bool_t isGoodEvent(StFemtoEvent *event, Double_t VtxZ, Double_t VtxR, Double_t delta_VtxY) {
  Bool_t check = true; 
  TVector3 pVtx = event->primaryVertex();

  // Reject vertices that are far from the central membrane along the beam
  if( TMath::Abs( pVtx.Z() ) > VtxZ ) check = false;
  //if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + 0.8847, 2) ) > 1 ) check = false; 14.5
  if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + delta_VtxY, 2) ) > VtxR ) check = false;
  
  return check;
}// isGoodEvent(){}

//********************CHECK TRACK ON GOOD********************//
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
}// isGoodTrack(){}

//***************CHECK TRACK FOR FLOW ON GOOD***************//
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
}// isGoodTrackFlow(){}

//**********************RECENTERING**********************//
void Recentering( TVector2 Q1vec[][3], TVector2 Q2vec[][n], TVector2 Q3vec[][n], Int_t RunID, Int_t cent,TProfile2D *tpQx1[][3], 
TProfile2D *tpQy1[][3], TProfile2D *tpQx2[][n], TProfile2D *tpQy2[][n],TProfile2D *tpQx3[][n], TProfile2D *tpQy3[][n] ) {
    
  TVector2 Q1mean, Q2mean, Q3mean;

  for(Int_t dir = 0; dir < 3; dir++) {

    for(Int_t eta = 0; eta < n; eta++) {
      Q2mean.Set(0.,0.);
      Q3mean.Set(0.,0.);

      Q2mean.Set( tpQx2[dir][eta] -> GetBinContent( tpQx2[dir][eta] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy2[dir][eta] -> GetBinContent( tpQy2[dir][eta] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      Q3mean.Set( tpQx3[dir][eta] -> GetBinContent( tpQx3[dir][eta] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy3[dir][eta] -> GetBinContent( tpQy3[dir][eta] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      
      if( Q2vec[dir][eta].Mod() != 0.) { 
        Q2vec[dir][eta] -= Q2mean;
      }
      if( Q3vec[dir][eta].Mod() != 0.) {
        Q3vec[dir][eta] -= Q3mean;
      }
    }

    for(Int_t det = 0; det < 2; det++) { 
      Q1mean.Set(0.,0.);
      Q1mean.Set( tpQx1[det][dir] -> GetBinContent( tpQx1[det][dir] -> FindBin( (Double_t)RunID, (Double_t)cent )),
                  tpQy1[det][dir] -> GetBinContent( tpQy1[det][dir] -> FindBin( (Double_t)RunID, (Double_t)cent )) );
      if( Q1vec[det][dir].Mod() != 0.) { 
        Q1vec[det][dir] -= Q1mean;
      }
    }

  }// for(Int_t dir = 0; dir < 3; dir++)
  
}// void Recentering(){}

//**********************CALCULATE Q-VEC BY TPC**********************//
void CalculateTPCQVec( TVector2 Q2vec[][n], TVector2 Q3vec[][n], StFemtoDst *dst, Int_t ngap, Double_t DCA_EVENT) {
  // 0 - this is east 
  // 1 - this is west 
  Double_t gap[6] = {0.05, 0.075, 0.1, 0.15, 0.2, 0.5};
  Float_t w[2][ngap];

  // Retrieve event information
  StFemtoEvent *event = dst->event();

  for(Int_t dir = 0; dir < 2; dir++) {
    for(Int_t eta = 0; eta < ngap; eta++) { 
      w[dir][eta] = 0.;
    }
  }

  Int_t nTracks = dst->numberOfTracks();

  // track loop 
  for(Int_t iTrk = 0; iTrk < nTracks; iTrk++) {

    // Retrieve i-th femto track
    StFemtoTrack *femtoTrack = dst -> track(iTrk);

    if( !isGoodTrack(event, femtoTrack, DCA_EVENT) ) continue; 

    for(Int_t eta = 0; eta < ngap; eta++) {

      if( femtoTrack -> eta() < -gap[eta] ) {
        Q2vec[0][eta].Set( Q2vec[0][eta].X() + femtoTrack -> pt() * cos(2.0 * femtoTrack -> phi() ),
                         	 Q2vec[0][eta].Y() + femtoTrack -> pt() * sin(2.0 * femtoTrack -> phi() ) );
        Q3vec[0][eta].Set( Q3vec[0][eta].X() + femtoTrack -> pt() * cos(3.0 * femtoTrack -> phi() ),
                           Q3vec[0][eta].Y() + femtoTrack -> pt() * sin(3.0 * femtoTrack -> phi() ) );
        w[0][eta]++; 
      }
      
      if( femtoTrack -> eta() > gap[eta] ) {
        Q2vec[1][eta].Set( Q2vec[1][eta].X() + femtoTrack -> pt() * cos(2.0 * femtoTrack -> phi()),
                           Q2vec[1][eta].Y() + femtoTrack -> pt() * sin(2.0 * femtoTrack -> phi() ) );
        Q3vec[1][eta].Set( Q3vec[1][eta].X() + femtoTrack -> pt() * cos(3.0 * femtoTrack -> phi() ),
                           Q3vec[1][eta].Y() + femtoTrack -> pt() * sin(3.0 * femtoTrack -> phi() ) );
        w[1][eta]++;
      }
    }
  }// track loop end

	for(Int_t eta = 0; eta < ngap; eta++) {
  	for(Int_t dir = 0; dir < 2; dir++) {   
      if( w[dir][eta] != 0 ) {
        Q2vec[dir][eta].Set( Q2vec[dir][eta].X()/w[dir][eta], Q2vec[dir][eta].Y()/w[dir][eta]);
        Q3vec[dir][eta].Set( Q3vec[dir][eta].X()/w[dir][eta], Q3vec[dir][eta].Y()/w[dir][eta]);
      }
    }
		Q2vec[2][eta] = Q2vec[1][eta] - Q2vec[0][eta];
  	Q3vec[2][eta] = Q3vec[1][eta] - Q3vec[0][eta];
  }

}// void CalculateTPCQVec(){}

//********************RETURN SIGNAL FROM ZDC's STRIP**********************//
Float_t ZDCSMD(Int_t eastwest, Int_t verthori, Int_t strip, StFemtoEvent *event ) {
  Float_t val = 0;

  if (eastwest == 0 && verthori == 0) val = event->zdcSmdEastVertical(strip-1);
  if (eastwest == 0 && verthori == 1) val = event->zdcSmdEastHorizontal(strip-1);
  if (eastwest == 1 && verthori == 0) val = event->zdcSmdWestVertical(strip-1);
  if (eastwest == 1 && verthori == 1) val = event->zdcSmdWestHorizontal(strip-1);
  
  return val;
}// Float_t ZDCSMD(){}

//********************RETURN POSITION FOR ZDC's STRIP**********************//
Float_t ZDCSMD_GetPosition( Int_t eastwest, Int_t verthori, Int_t strip) {
  Float_t zdcsmd_x[7] = {0.5,2,3.5,5,6.5,8,9.5};
  Float_t zdcsmd_y[8] = {1.25,3.25,5.25,7.25,9.25,11.25,13.25,15.25};

  if(eastwest==0 && verthori==0) return zdcsmd_x[strip-1]-mZDCSMDCenterex;
  if(eastwest==1 && verthori==0) return mZDCSMDCenterwx-zdcsmd_x[strip-1];
  if(eastwest==0 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterey;
  if(eastwest==1 && verthori==1) return zdcsmd_y[strip-1]/sqrt(2.)-mZDCSMDCenterwy;

  return 0;
}// Float_t ZDCSMD_GetPosition(){}

//**********************CALCULATE Q-VEC BY ZDC**********************//
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

}// TVector2 CalculateZDCQVec(){}

//**********************CALCULATE Q-VEC BY BBC**********************//
TVector2 CalculateBBCQVec(StFemtoEvent *event, Int_t ew, Float_t harm) {
	Float_t Qx = 0., Qy = 0., Qwgt = 0.;
  TVector2 Q(0.,0.);
  Int_t nstrip = 16;
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
  else return Q;
  Q.Set(Qx,Qy);
  return Q;
}// TVector2 CalculateBBCQVec(){}

//********************RETURN PHI FOR BBC's TITLE**********************//
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
}// Float_t GetBBCTilePhi(){}