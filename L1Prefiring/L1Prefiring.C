#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <map>

#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile2D.h>


#include <time.h>       /* time_t, struct tm, difftime, time, mktime */

/* 
Written by David Petyt

Edited by Raman Khurana

11 December 2018
Added: information for the pre firing variables 
     : 1d histograms for each ieta ring, X axis is emulator index. for two set of selection, 
     : tp > 0 , tp > 0 and energy > 16 ADC 
     
12 December 2018: 
Added: configurable energy cut, to make it possible to run with different value of energy threshold. 

To be added: 


*/


using  namespace std;

bool is2017 = true; 


struct EcalAux
{
  int iMaskedTTeta[100];
  int iMaskedTTphi[100];
  int nMaskedCh;
  int runNb;
  int iMaskedChannelsFromStripsX[500];
  int iMaskedChannelsFromStripsY[500];
  int iMaskedChannelsFromStripsZ[500];
  int nMaskedChannelsFromStrips;
};


struct EcalTPGVariables
{
   
  // event variables
  UInt_t runNb ;
  ULong64_t evtNb ;
  UInt_t bxNb ;
  UInt_t lumiBlock ;
  ULong64_t orbitNb ;
  double timeStamp ; 
  UInt_t nbOfActiveTriggers ;
  UInt_t nbOfActiveTriggersBX ;
   
  int activeTriggers[128] ;
  int activeTriggersBX[128] ;
   
  int activeTechTriggers[64] ;
  UInt_t nbOfActiveTechTriggers ;
   
  // tower variables
  UInt_t nbOfTowers ; //max 4032 EB+EE
  int ieta[4032] ;
  int iphi[4032] ;
  int nbOfXtals[4032] ;
  int rawTPData[4032] ;
  int rawTPEmul1[4032] ;
  int rawTPEmul2[4032] ;
  int rawTPEmul3[4032] ;
  int rawTPEmul4[4032] ;
  int rawTPEmul5[4032] ;
  int rawTPEmulttFlag1[4032] ;
  int rawTPEmulttFlag2[4032] ;
  int rawTPEmulttFlag3[4032] ;
  int rawTPEmulttFlag4[4032] ;
  int rawTPEmulttFlag5[4032] ;
  int rawTPEmulsFGVB1[4032] ;
  int rawTPEmulsFGVB2[4032] ;
  int rawTPEmulsFGVB3[4032] ;
  int rawTPEmulsFGVB4[4032] ;
  int rawTPEmulsFGVB5[4032] ;
  
  float eRec[4032] ;
  int crystNb[4032] ;
  int spike[4032] ;
  int sevlv[4032];
  int ttFlag[4032];
  int trig_tower_adc[4032], trig_tower_sFGVB[4032]; 

   
  UInt_t nMaskedRCT ;
  int iMaskedRCTeta[100] ;
  int iMaskedRCTphi[100] ;
  int iMaskedRCTcrate[100] ;
   
  UInt_t nbOfL1IsoCands ;
  int L1IsoIeta[8] ;
  int L1IsoIphi[8] ;
  int L1IsoRank[8] ; 
  UInt_t nbOfL1NonisoCands ;
  int L1NonisoIeta[8] ;
  int L1NonisoIphi[8] ;
  int L1NonisoRank[8] ; 
   
  UInt_t nbOfL1preIsoCands ;
  int L1preIsoIeta[8] ;
  int L1preIsoIphi[8] ;
  int L1preIsoRank[8] ; 
  UInt_t nbOfL1preNonisoCands ;
  int L1preNonisoIeta[8] ;
  int L1preNonisoIphi[8] ;
  int L1preNonisoRank[8] ; 
   
  UInt_t nbOfL1postIsoCands ;
  int L1postIsoIeta[8] ;
  int L1postIsoIphi[8] ;
  int L1postIsoRank[8] ; 
  UInt_t nbOfL1postNonisoCands ;
  int L1postNonisoIeta[8] ;
  int L1postNonisoIphi[8] ;
  int L1postNonisoRank[8] ; 
   
} ;




bool isMasked(std::vector<float> etaV, std::vector<float> phiV, float ieta, float iphi){
  bool maskedstatus = false; 
  for (int i = 0 ; i < int( phiV.size()); i++) {
    
    int a = (int) etaV[i];
    int b = (int) phiV[i];
 
      if ((ieta==a) && (iphi==b) ) {
	maskedstatus = true;
	break;
      }
    
    
  }
  return maskedstatus; 
}





std::vector<float> MaskedCoordinate(int icol=1){
  std::vector<float> ietaV;
  ietaV.clear();
  string line;
  TString maskedfile;
  if (is2017)  maskedfile = "masked_TT_all_2017.txt";
  if (!is2017) maskedfile = "masked_TT_all_2018.txt";
  
  ifstream maskedTT (maskedfile);
  if (maskedTT.is_open()){
    while ( getline (maskedTT,line) ){
      
      std::stringstream ss(line);
      float a, b;
      ss >> a >> b;
      
      if (icol==1) ietaV.push_back(a);
      if (icol==2) ietaV.push_back(b);
      
    }
    
  }
  maskedTT.close();
  return ietaV;
}

void setBranchAddresses (TChain * chain, EcalTPGVariables & treeVars)
{
  chain->SetBranchAddress ("runNb",&treeVars.runNb) ; 
  chain->SetBranchAddress ("evtNb",&treeVars.evtNb) ; 
  chain->SetBranchAddress ("bxNb",&treeVars.bxNb) ; 
  chain->SetBranchAddress ("lumiBlock",&treeVars.lumiBlock) ;
  chain->SetBranchAddress ("orbitNb",&treeVars.orbitNb) ; 
  chain->SetBranchAddress ("timeStamp",&treeVars.timeStamp) ; 
  chain->SetBranchAddress ("nbOfActiveTriggers",&treeVars.nbOfActiveTriggers) ; 
  chain->SetBranchAddress ("activeTriggers",treeVars.activeTriggers) ; 
  chain->SetBranchAddress ("activeTriggersBX",treeVars.activeTriggersBX) ;
   
  chain->SetBranchAddress ("nbOfActiveTechTriggers",&treeVars.nbOfActiveTechTriggers) ; 
  chain->SetBranchAddress ("activeTechTriggers",treeVars.activeTechTriggers) ; 
   
  chain->SetBranchAddress ("nbOfTowers",&treeVars.nbOfTowers) ; 
  chain->SetBranchAddress ("ieta",treeVars.ieta) ; 
  chain->SetBranchAddress ("iphi",treeVars.iphi) ; 
  chain->SetBranchAddress ("nbOfXtals",treeVars.nbOfXtals) ; 
  chain->SetBranchAddress ("rawTPData",treeVars.rawTPData) ; 
  chain->SetBranchAddress ("rawTPEmul1",treeVars.rawTPEmul1) ; 
  chain->SetBranchAddress ("rawTPEmul2",treeVars.rawTPEmul2) ; 
  chain->SetBranchAddress ("rawTPEmul3",treeVars.rawTPEmul3) ; 
  chain->SetBranchAddress ("rawTPEmul4",treeVars.rawTPEmul4) ; 
  chain->SetBranchAddress ("rawTPEmul5",treeVars.rawTPEmul5) ; 
  chain->SetBranchAddress ("eRec",treeVars.eRec) ; 
  chain->SetBranchAddress ("crystNb",treeVars.crystNb) ;
  chain->SetBranchAddress ("spike",treeVars.spike) ;
  chain->SetBranchAddress ("sevlv", treeVars.sevlv);
  chain->SetBranchAddress ("ttFlag", treeVars.ttFlag);
  chain->SetBranchAddress ("trig_tower_adc",treeVars.trig_tower_adc) ; 
  chain->SetBranchAddress ("trig_tower_sFGVB",treeVars.trig_tower_sFGVB) ; 
  chain->SetBranchAddress ("rawTPEmulsFGVB1",treeVars.rawTPEmulsFGVB1) ; 
  chain->SetBranchAddress ("rawTPEmulsFGVB2",treeVars.rawTPEmulsFGVB2) ; 
  chain->SetBranchAddress ("rawTPEmulsFGVB3",treeVars.rawTPEmulsFGVB3) ; 
  chain->SetBranchAddress ("rawTPEmulsFGVB4",treeVars.rawTPEmulsFGVB4) ; 
  chain->SetBranchAddress ("rawTPEmulsFGVB5",treeVars.rawTPEmulsFGVB5) ; 
  chain->SetBranchAddress ("rawTPEmulttFlag1",treeVars.rawTPEmulttFlag1) ; 
  chain->SetBranchAddress ("rawTPEmulttFlag2",treeVars.rawTPEmulttFlag2) ; 
  chain->SetBranchAddress ("rawTPEmulttFlag3",treeVars.rawTPEmulttFlag3) ; 
  chain->SetBranchAddress ("rawTPEmulttFlag4",treeVars.rawTPEmulttFlag4) ; 
  chain->SetBranchAddress ("rawTPEmulttFlag5",treeVars.rawTPEmulttFlag5) ; 


   
  chain->SetBranchAddress ("nMaskedRCT",&treeVars.nMaskedRCT); //
  chain->SetBranchAddress ("iMaskedRCTeta", treeVars.iMaskedRCTeta);//
  chain->SetBranchAddress ("iMaskedRCTcrate", treeVars.iMaskedRCTcrate);//
  chain->SetBranchAddress ("iMaskedRCTphi", treeVars.iMaskedRCTphi);//
   
  chain->SetBranchAddress ("nbOfL1IsoCands",&treeVars.nbOfL1IsoCands); //
  chain->SetBranchAddress ("L1IsoIeta", treeVars.L1IsoIeta);//
  chain->SetBranchAddress ("L1IsoIphi", treeVars.L1IsoIphi);//
  chain->SetBranchAddress ("L1IsoRank", treeVars.L1IsoRank);//
   
  chain->SetBranchAddress ("nbOfL1NonisoCands",&treeVars.nbOfL1NonisoCands); //
  chain->SetBranchAddress ("L1NonisoIeta", treeVars.L1NonisoIeta);//
  chain->SetBranchAddress ("L1NonisoIphi", treeVars.L1NonisoIphi);//
  chain->SetBranchAddress ("L1NonisoRank", treeVars.L1NonisoRank);//
   
  chain->SetBranchAddress ("nbOfL1preIsoCands",&treeVars.nbOfL1preIsoCands); //
  chain->SetBranchAddress ("L1preIsoIeta", treeVars.L1preIsoIeta);//
  chain->SetBranchAddress ("L1preIsoIphi", treeVars.L1preIsoIphi);//
  chain->SetBranchAddress ("L1preIsoRank", treeVars.L1preIsoRank);//
   
  chain->SetBranchAddress ("nbOfL1preNonisoCands",&treeVars.nbOfL1preNonisoCands); //
  chain->SetBranchAddress ("L1preNonisoIeta", treeVars.L1preNonisoIeta);//
  chain->SetBranchAddress ("L1preNonisoIphi", treeVars.L1preNonisoIphi);//
  chain->SetBranchAddress ("L1preNonisoRank", treeVars.L1preNonisoRank);//
   
  chain->SetBranchAddress ("nbOfL1postIsoCands",&treeVars.nbOfL1postIsoCands); //
  chain->SetBranchAddress ("L1postIsoIeta", treeVars.L1postIsoIeta);//
  chain->SetBranchAddress ("L1postIsoIphi", treeVars.L1postIsoIphi);//
  chain->SetBranchAddress ("L1postIsoRank", treeVars.L1postIsoRank);//
   
  chain->SetBranchAddress ("nbOfL1postNonisoCands",&treeVars.nbOfL1postNonisoCands); //
  chain->SetBranchAddress ("L1postNonisoIeta", treeVars.L1postNonisoIeta);//
  chain->SetBranchAddress ("L1postNonisoIphi", treeVars.L1postNonisoIphi);//
  chain->SetBranchAddress ("L1postNonisoRank", treeVars.L1postNonisoRank);//
}



void setAuxBranchAddresses (TChain * chain, EcalAux & auxVars)
{
  chain->SetBranchAddress ("iMaskedTTeta",auxVars.iMaskedTTeta) ; 
  chain->SetBranchAddress ("iMaskedTTphi",auxVars.iMaskedTTphi) ;
  chain->SetBranchAddress ("nMaskedCh",&auxVars.nMaskedCh) ;
  chain->SetBranchAddress ("runNb",&auxVars.runNb) ;   
   
  chain->SetBranchAddress ("nMaskedChannelsFromStrips",&auxVars.nMaskedChannelsFromStrips) ;
  chain->SetBranchAddress ("iMaskedChannelsFromStripsX",auxVars.iMaskedChannelsFromStripsX) ;
  chain->SetBranchAddress ("iMaskedChannelsFromStripsY",auxVars.iMaskedChannelsFromStripsY) ;
  chain->SetBranchAddress ("iMaskedChannelsFromStripsZ",auxVars.iMaskedChannelsFromStripsZ) ;   
}


int getEt(UInt_t val) {return (val&0xff) ;}

UInt_t getFg(UInt_t val) {return ((val&0x100)!=0) ;}

UInt_t getTtf(UInt_t val) {return ((val>>9)&0x7) ;}


///////  Main program /////////

//void tpgreader()
void L1Prefiring(int threshold=16)
{  
  std::cout<<" threshold = "<<threshold << std::endl;

  time_t timer;
  double initial_t = time(&timer);
  
  TChain * chainAux = new TChain ("treeAux") ;
  TChain * chain = new TChain ("EcalTPGAnalysis") ;

  // input TPGTree file
  
  TString inputrootfile;
  if (is2017)  inputrootfile = "/eos/cms/store/user/pbarria/TPG/ECALTPGTree_ZeroBias2_Run2017F-Run306425_RAW-RECO.root"; // 2017 with prefiring 
  if (!is2017) inputrootfile = "/afs/cern.ch/user/s/sdutt/work/public/ECALTPGTree_subset.root" ;                          // 2018 without prefiring, full readout 
     
  
  chain->Add(inputrootfile);
  chainAux->Add(inputrootfile);
  


  // if set to true: print debug info on matched and unmatched towers
  
  bool debug=false;


  // define output histograms
  
  TH2F* ieta_vs_iphi_TP             = new TH2F("ieta_vs_iphi_TP","ieta_vs_iphi_TP", 57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_ETP            = new TH2F("ieta_vs_iphi_ETP","ieta_vs_iphi_ETP", 57,-28.5,28.5,72,0.5,72.5);
  
  //57,-28.5,28.5
  TH2F* ieta_vs_iphi_TP_ETP         = new TH2F("ieta_vs_iphi_TP_ETP","ieta_vs_iphi_TP_ETP",57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_TP_ETP_INTIME  =  new TH2F("ieta_vs_iphi_TP_ETP_INTIME","ieta_vs_iphi_TP_ETP_INTIME",57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_TP_ETP_OUTTIME = new TH2F("ieta_vs_iphi_TP_ETP_OUTTIME","ieta_vs_iphi_TP_ETP_OUTTIME",57,-28.5,28.5,72,0.5,72.5);

  
  TH2F* ieta_vs_iphi_TP0            = new TH2F("ieta_vs_iphi_TP0","ieta_vs_iphi_TP0",57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_TP0_ETP        = new TH2F("ieta_vs_iphi_TP0_ETP","ieta_vs_iphi_TP0_ETP",57,-28.5,28.5,72,0.5,72.5);


  TH2F* ieta_vs_idx_TP_ETP          = new TH2F("ieta_vs_idx_TP_ETP","ieta_vs_idx_TP_ETP",11,17.5,28.5, 5, -0.5, 4.5); 
  TH2F* ieta_vs_idx_TP0             = new TH2F("ieta_vs_idx_TP0","ieta_vs_idx_TP0", 11,17.5,28.5 , 5, -0.5, 4.5);
  TH2F* ieta_vs_idx_TP0_ETP         = new TH2F("ieta_vs_idx_TP0_ETP","ieta_vs_idx_TP0_ETP", 11,17.5,28.5 , 5, -0.5, 4.5);
  
  
  TH2F* ieta_vs_iphi_TP16           = new TH2F("ieta_vs_iphi_TP16","ieta_vs_iphi_TP16", 57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_ETP16_IDX2     = new TH2F("ieta_vs_iphi_ETP16_IDX2","ieta_vs_iphi_ETP16_IDX2",57,-28.5,28.5,72,0.5,72.5);
  TH2F* ieta_vs_iphi_ETP16_IDXAny   = new TH2F("ieta_vs_iphi_ETP16_IDXAny","ieta_vs_iphi_ETP16_IDXAny",57,-28.5,28.5,72,0.5,72.5);

  TH2F* idx_vs_ieta_TP16           = new TH2F("idx_vs_ieta_TP16","idx_vs_ieta_TP16", 11,17.5,28.5,5, -0.5, 4.5 );
  TH2F* idx_vs_ieta_ETP16_IDX2     = new TH2F("idx_vs_ieta_ETP16_IDX2","idx_vs_ieta_ETP16_IDX2",11,17.5,28.5, 5, -0.5, 4.5);
  TH2F* idx_vs_ieta_ETP16_IDXAny   = new TH2F("idx_vs_ieta_ETP16_IDXAny","idx_vs_ieta_ETP16_IDXAny",11,17.5,28.5,5, -0.5, 4.5);

  TH1F* idx_TP0                    = new TH1F("idx_TP0","idx_TP0_ETP",5, -0.5, 4.5);
  TH1F* idx_TP0_ETP                = new TH1F("idx_TP0_ETP","idx_TP0_ETP",5, -0.5, 4.5);
  
  TH1F* idx_TP0_vec[11];
  TH1F* idx_TP0_ETP_vec[11];
  
  for (int i=0; i<11; i++){
    TString postname;
    postname.Form("%d",i+18);
    std::cout<<" postname = "<<postname<<std::endl;
    idx_TP0_vec[i] = new TH1F("idx_TP0_vec_"+postname, "idx_TP0_vec_"+postname, 5, -0.5, 4.5);
    idx_TP0_ETP_vec[i] = new TH1F("idx_TP0_ETP_vec_"+postname, "idx_TP0_ETP_vec_"+postname, 5, -0.5, 4.5);
    
  }
  
  TH1F* ieta_TP_ETP_INTIME         = new TH1F("ieta_TP_ETP_INTIME","ieta_TP_ETP_INTIME",57,-28.5,28.5);
  TH1F* ieta_TP_ETP_OUTTIME        = new TH1F("ieta_TP_ETP_OUTTIME","ieta_TP_ETP_OUTTIME",57,-28.5,28.5);
  TH1F* ieta_TP0_ETP               = new TH1F("ieta_TP0_ETP","ieta_TP0_ETP",57,-28.5,28.5);
  TH1F* ieta_TP_ETP                = new TH1F("ieta_TP_ETP","ieta_TP_ETP",57,-28.5,28.5);
  TH1F* ieta_TP0                   = new TH1F("ieta_TP0","ieta_TP0",57,-28.5,28.5);

  
  /*
  
  TH1F *realTP=new TH1F("realTP","",256,0,256);
  TH1F *emulTP=new TH1F("emulTP","",256,0,256);

  TH1F *realTP_sev=new TH1F("realTP_sev","",256,0,256);
  TH1F *emulTP_sev=new TH1F("emulTP_sev","",256,0,256);

  TH2F *realvsemulTP=new TH2F("realVsEmulTP","",256,0,256,256,0,256);

  TH2F *etaphimismatch=new TH2F("etaphimismatch","",35,-17.5,17.5,72,0.5,72.5);

  TH2F *etaphimatch=new TH2F("etaphimatch","",35,-17.5,17.5,72,0.5,72.5);

  TH2F *realminusemulTPvsreal=new TH2F("realminusEmulTPvsreal","",32,0,256,21,-10.5,10.5);

  TH1F *realminusemulTP=new TH1F("realminusEmulTP","",21,-10.5,10.5);


  TH2F *realminusemulTPvsieta=new TH2F("realminusEmulTPvsieta","",57,-28.5,28.5,41,-20.5,20.5);

  TH2F *realminusemulTPratvsieta=new TH2F("realminusEmulTPratvsieta","",57,-28.5,28.5,41,-1.05,1.05);

  TH1F *realTPEE=new TH1F("realTPEE","",256,0,256);
  TH1F *emulTPEE=new TH1F("emulTPEE","",256,0,256);

  TH2F *realvsemulTPEE=new TH2F("realVsEmulTPEE","",256,0,256,256,0,256);

  TH2F *etaphimismatchEE=new TH2F("etaphimismatchEE","",57,-28.5,28.5,72,0.5,72.5);

  TH2F *etaphimatchEE=new TH2F("etaphimatchEE","",57,-28.5,28.5,72,0.5,72.5);

  TH2F *realminusemulTPvsrealEE=new TH2F("realminusEmulTPvsrealEE","",32,0,256,21,-10.5,10.5);

  TH1F *realminusemulTPEE=new TH1F("realminusEmulTPEE","",21,-10.5,10.5);
  */


  // set branch addresses to read trees
  
  EcalTPGVariables treeVars ;
  setBranchAddresses (chain, treeVars) ;
   
  EcalAux auxVars ;
  setAuxBranchAddresses (chainAux, auxVars) ;

  
  Int_t treeentries=chain->GetEntries();

  int count=0;
  int countsev=0;
  int countsev3=0;
  int countsev4=0;

  std::vector<float> etaV = MaskedCoordinate();
  std::vector<float> phiV = MaskedCoordinate();
  
  for (int entry = 0 ; entry < treeentries ; ++entry) {

    if (entry%5000==0) cout << entry << " / " << treeentries
			     << " events processed" << endl;


    //cout << entry << " / " << treeentries
    //	 << " events processed" << endl;
    
    chain->GetEntry (entry) ;

    UInt_t runNb = treeVars.runNb;
    ULong64_t evtNb = treeVars.evtNb;
    UInt_t bxNb = treeVars.bxNb;;
    UInt_t lumiBlock = treeVars.lumiBlock;
 


    for (UInt_t tower = 0 ; tower < treeVars.nbOfTowers ; tower++) {
      
      int ieta = treeVars.ieta[tower] ;
      int iphi = treeVars.iphi[tower] ;
      if (abs(ieta) <18) continue; 
      
      int tp = getEt(treeVars.rawTPData[tower]) ;
      int raw_spike = treeVars.trig_tower_sFGVB[tower];
      int ttflag=treeVars.ttFlag[tower];
      float erec_val=treeVars.eRec[tower];
      int severity_level=treeVars.sevlv[tower];
      int nxtals=treeVars.crystNb[tower];
      
      int emul[5] = {getEt(treeVars.rawTPEmul1[tower]),  
		     getEt(treeVars.rawTPEmul2[tower]),
		     getEt(treeVars.rawTPEmul3[tower]),
		     getEt(treeVars.rawTPEmul4[tower]),
		     getEt(treeVars.rawTPEmul5[tower])} ;
      int emulsFGVB[5] = {treeVars.rawTPEmulsFGVB1[tower],  
			  treeVars.rawTPEmulsFGVB2[tower],
			  treeVars.rawTPEmulsFGVB3[tower],
			  treeVars.rawTPEmulsFGVB4[tower],
			  treeVars.rawTPEmulsFGVB5[tower]} ;
      
      int emulTTF[5] = {treeVars.rawTPEmulttFlag1[tower],  
			    treeVars.rawTPEmulttFlag2[tower],
			treeVars.rawTPEmulttFlag3[tower],
			treeVars.rawTPEmulttFlag4[tower],
			treeVars.rawTPEmulttFlag5[tower]} ;
      
      
      
      int maxOfTPEmul = 0 ;
      int indexOfTPEmulMax = -1 ;
      int emulsfgvb=-1;
      int emulttf=-1;
      
      
      for (int i=0 ; i<5 ; i++) {
	if (emul[i]>maxOfTPEmul)
	  {
	    maxOfTPEmul = emul[i] ; 
	    indexOfTPEmulMax = i ;
	    emulsfgvb=emulsFGVB[i];
	    emulttf=emulTTF[i];
	  }
      }
      

      // get the masked towers.
      if (tp>0.0) ieta_vs_iphi_TP->Fill(ieta,iphi);
      if (maxOfTPEmul > 0.0) ieta_vs_iphi_ETP->Fill(ieta,iphi);
      // masked towers ends here. 
      
      

      // check if this is already tagged by COKE
      bool hasnoisyxtal  = false;
      if (is2017) hasnoisyxtal = false;
      if (!is2017) hasnoisyxtal  = (ttflag >3) ; 

      // check if the TT is already masked. 
      bool masked_ = isMasked(etaV, phiV, ieta, iphi);
      
      // skip the tower if this is msked or it is tagged by COKE as noisy. 
      if (masked_) continue;
      if (hasnoisyxtal) continue;
      
      
      
      // data TP has energy > 16 ADC
      if (tp>threshold){
	ieta_vs_iphi_TP16->Fill(ieta, iphi);
        idx_vs_ieta_TP16->Fill(ieta, indexOfTPEmulMax);

      }
      
      // emulated TP has energy > 16 ADC and it is in right bunch crossing 
      if(maxOfTPEmul>threshold && indexOfTPEmulMax==2){
	ieta_vs_iphi_ETP16_IDX2->Fill(ieta, iphi);
	idx_vs_ieta_ETP16_IDX2->Fill(ieta, indexOfTPEmulMax);
      }
      
      // emiulated Tp has energy > 16 ADC and it can be in any bunch crossing. 
      if(maxOfTPEmul>threshold){
	ieta_vs_iphi_ETP16_IDXAny->Fill(ieta, iphi);
	idx_vs_ieta_ETP16_IDXAny->Fill(ieta, indexOfTPEmulMax);
      }



      
      // denominator for condition 1 and 2 are same. The denominator for condition 3 is defined differently. 
      // I might have to change this  (ask David?) 
     

      // denominator seelction for 1 and 2
      if (tp>threshold && maxOfTPEmul>0){
	ieta_vs_iphi_TP_ETP->Fill(ieta, iphi);
	ieta_TP_ETP->Fill(ieta);

	
	// 
	ieta_vs_idx_TP_ETP->Fill(abs(ieta), indexOfTPEmulMax);
	
	
	
	//std::cout<<" denominator for the cond1 and cond2"<<std::endl;
	// condition 1: 
	// where both data and emulated TPs exist, and the emulated TP is in the correct bunch crossing
	if (tp>16 && maxOfTPEmul>0 && indexOfTPEmulMax==2){
	  ieta_vs_iphi_TP_ETP_INTIME->Fill(ieta, iphi);
	  ieta_TP_ETP_INTIME->Fill(ieta);
	  //std::cout<<" numerator for the cond1 and cond2"<<std::endl;
	}
	
	// condition 2: 
	// where both data and emulated TPs exist, but the emulated TP is in the wrong bunch crossing
	if (tp>threshold && maxOfTPEmul>0 && indexOfTPEmulMax!=2){
	  ieta_vs_iphi_TP_ETP_OUTTIME->Fill(ieta, iphi);
	  ieta_TP_ETP_OUTTIME->Fill(ieta);
	  //std::cout<<" numerator for the cond2"<<std::endl;
	}
      }// if (tp>16 && maxOfTPEmul>0){
      
      
      
      // denominator for condition 3: 
      if (tp==0  && !hasnoisyxtal ){
	//std::cout<<" this TT is masked "<<masked_<<" "<<ieta<<" " <<iphi<<" "<<std::endl;
	if (masked_) continue ;
	
	//
	idx_TP0->Fill(indexOfTPEmulMax);
	int ieta_index = abs(ieta) - 18;
	//std::cout<<" ieta, ieta_index "<< ieta<<"  "<<ieta_index<<std::endl;
	idx_TP0_vec[ieta_index]->Fill(indexOfTPEmulMax);
	
	ieta_vs_iphi_TP0->Fill(ieta, iphi);
	ieta_TP0->Fill(ieta) ;
	
	ieta_vs_idx_TP0->Fill(abs(ieta), indexOfTPEmulMax);
	// condition 3: 
	// where only an emulated TP exists  (here it will be necessary to remove towers that are masked)
	if ( tp==0 && maxOfTPEmul>threshold) {
	  
	  // 
	  ieta_vs_idx_TP0_ETP->Fill(abs(ieta), indexOfTPEmulMax);
	  
	  
	  //std::cout<<" TT flag = "<<ttflag<<std::endl;
	  
	  
	  ieta_vs_iphi_TP0_ETP->Fill(ieta, iphi);
	  idx_TP0_ETP->Fill(indexOfTPEmulMax);
	  
	  idx_TP0_ETP_vec[ieta_index]->Fill(indexOfTPEmulMax);
	  
	  ieta_TP0_ETP->Fill(ieta) ;
	  
	  //std::cout<<" numerator for the cond3"<<std::endl;

	}
      }//if (tp==0){
      
      
      
      // barrel: abs(ieta)  < 18
      // endcap: abs(ieta) >= 18 (10 TT on each side)
      
      
    }
    
  }
  
  TString outputfilename;
  TString cutval;
  cutval.Form("%d",threshold);
  if (is2017) outputfilename  = "PrefiringRateEE_2017data_cutval_"+cutval+".root";
  if (!is2017) outputfilename = "PrefiringRateEE_2018data_cutval_"+cutval+".root";
  
  TFile fout(outputfilename,"RECREATE");
  fout.cd();

  ieta_vs_iphi_TP->Write();
  ieta_vs_iphi_ETP->Write();
  
  ieta_vs_iphi_TP_ETP->Write();
  ieta_vs_iphi_TP_ETP_INTIME->Write();
  ieta_vs_iphi_TP_ETP_OUTTIME->Write();
  
  ieta_vs_iphi_TP0->Write();
  ieta_vs_iphi_TP0_ETP->Write();
  
  
  ieta_vs_idx_TP_ETP->Write();
  ieta_vs_idx_TP0->Write();
  ieta_vs_idx_TP0_ETP->Write();

  ieta_vs_iphi_TP16->Write();
  ieta_vs_iphi_ETP16_IDX2->Write();
  ieta_vs_iphi_ETP16_IDXAny->Write();

  idx_vs_ieta_TP16->Write();
  idx_vs_ieta_ETP16_IDX2->Write();
  idx_vs_ieta_ETP16_IDXAny->Write();
  

  idx_TP0->Write();
  idx_TP0_ETP->Write();


  //ieta_vs_idx_TP0->Write();
  //ieta_vs_idx_TP0_ETP->Write();
  
  // after writing the histograms, perform operations on them to extract fraction etc. 
  
  TH2F* frac_ieta_vs_iphi_TP_ETP_INTIME = ((TH2F*) fout.Get("ieta_vs_iphi_TP_ETP_INTIME"));
  frac_ieta_vs_iphi_TP_ETP_INTIME->SetName("frac_ieta_vs_iphi_TP_ETP_INTIME");//new TH2F("ieta_vs_iphi_TP_ETP_INTIME","ieta_vs_iphi_TP_ETP_INTIME",57,-28.5,28.5,72,0.5,72.5);
  frac_ieta_vs_iphi_TP_ETP_INTIME->Divide(ieta_vs_iphi_TP_ETP);
  frac_ieta_vs_iphi_TP_ETP_INTIME->Write();

  
  TH2F* frac_ieta_vs_iphi_TP_ETP_OUTTIME = (TH2F*) fout.Get("ieta_vs_iphi_TP_ETP_OUTTIME");
  frac_ieta_vs_iphi_TP_ETP_OUTTIME->SetName("frac_ieta_vs_iphi_TP_ETP_OUTTIME");
  frac_ieta_vs_iphi_TP_ETP_OUTTIME->Divide(ieta_vs_iphi_TP_ETP);
  frac_ieta_vs_iphi_TP_ETP_OUTTIME->Write();
  

  TH2F* frac_ieta_vs_iphi_TP0_ETP = (TH2F*) fout.Get("ieta_vs_iphi_TP0_ETP");
  frac_ieta_vs_iphi_TP0_ETP->SetName("frac_ieta_vs_iphi_TP0_ETP");
  frac_ieta_vs_iphi_TP0_ETP->Divide(ieta_vs_iphi_TP0);
  frac_ieta_vs_iphi_TP0_ETP->Write();
  

  TH1F* frac_idx_TP0_ETP = (TH1F*) fout.Get("idx_TP0_ETP");
  frac_idx_TP0_ETP->SetName("frac_idx_TP0_ETP");
  frac_idx_TP0_ETP->Divide(idx_TP0);
  frac_idx_TP0_ETP->Write();


  ieta_TP_ETP_INTIME->Write()  ;
  ieta_TP_ETP_OUTTIME->Write() ;
  ieta_TP0_ETP->Write()        ;
  ieta_TP_ETP->Write()         ;
  ieta_TP0->Write()            ;

  TH1F* frac_ieta_TP_ETP_INTIME = (TH1F*) fout.Get("ieta_TP_ETP_INTIME");
  frac_ieta_TP_ETP_INTIME->SetName("frac_ieta_TP_ETP_INTIME");
  frac_ieta_TP_ETP_INTIME->Divide(ieta_TP_ETP);
  frac_ieta_TP_ETP_INTIME->Write();


  TH1F* frac_ieta_TP_ETP_OUTTIME = (TH1F*) fout.Get("ieta_TP_ETP_OUTTIME");
  frac_ieta_TP_ETP_OUTTIME->SetName("frac_ieta_TP_ETP_OUTTIME");
  frac_ieta_TP_ETP_OUTTIME->Divide(ieta_TP_ETP);
  frac_ieta_TP_ETP_OUTTIME->Write();


  TH1F* frac_ieta_TP0_ETP = (TH1F*) fout.Get("ieta_TP0_ETP");
  frac_ieta_TP0_ETP->SetName("frac_ieta_TP0_ETP");
  frac_ieta_TP0_ETP->Divide(ieta_TP0);
  frac_ieta_TP0_ETP->Write();


  TH2F* frac_ieta_vs_iphi_TP16_DIV_ETP16_IDXAny = (TH2F*) fout.Get("ieta_vs_iphi_TP16");
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDXAny->SetName("frac_ieta_vs_iphi_TP16_DIV_ETP16_IDXAny");
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDXAny->Divide(ieta_vs_iphi_ETP16_IDXAny);
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDXAny->Write();

  TH2F* frac_ieta_vs_iphi_TP16_DIV_ETP16_IDX2 = (TH2F*) fout.Get("ieta_vs_iphi_TP16");
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDX2->SetName("frac_ieta_vs_iphi_TP16_DIV_ETP16_IDX2");
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDX2->Divide(ieta_vs_iphi_ETP16_IDX2);
  frac_ieta_vs_iphi_TP16_DIV_ETP16_IDX2->Write();


  TH2F* frac_idx_vs_ieta_TP16_DIV_ETP16_IDXAny = (TH2F*) fout.Get("idx_vs_ieta_TP16");
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDXAny->SetName("frac_idx_vs_ieta_TP16_DIV_ETP16_IDXAny");
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDXAny->Divide(idx_vs_ieta_ETP16_IDXAny);
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDXAny->Write();
  
  TH2F* frac_idx_vs_ieta_TP16_DIV_ETP16_IDX2 = (TH2F*) fout.Get("idx_vs_ieta_TP16");
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDX2->SetName("frac_idx_vs_ieta_TP16_DIV_ETP16_IDX2");
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDX2->Divide(idx_vs_ieta_ETP16_IDX2);
  frac_idx_vs_ieta_TP16_DIV_ETP16_IDX2->Write();
  
  
  TH2F* frac_ieta_vs_idx_TP0_ETP_DIV_TP0 = (TH2F*) fout.Get("ieta_vs_idx_TP0_ETP");
  frac_ieta_vs_idx_TP0_ETP_DIV_TP0->SetName("frac_ieta_vs_idx_TP0_ETP_DIV_TP0");
  frac_ieta_vs_idx_TP0_ETP_DIV_TP0->Divide(ieta_vs_idx_TP0);
  frac_ieta_vs_idx_TP0_ETP_DIV_TP0->Write();
  
  
  fout.Close();

  double final_t = time(&timer);
  
  double seconds = -difftime(initial_t,final_t);
  std::cout<<" time used to run the macro (in seconds):" << seconds<<std::endl; 
  
  
}

