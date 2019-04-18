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

22 January 2019: Resuming editing of the code following the suggestoion from David. 
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
   
  
  // new branches 

  UInt_t          nbOfL1preNonisoCandsm2;
  Int_t           L1preNonisoIetam2[9]; 
  Int_t           L1preNonisoIphim2[9]; 
  Int_t           L1preNonisoRankm2[9]; 
  UInt_t          nbOfL1preIsoCandsm2;
  Int_t           L1preIsoIetam2[12];   
  Int_t           L1preIsoIphim2[12];   
  Int_t           L1preIsoRankm2[12];   
  UInt_t          nbOfL1preNonisoCandsm1;
  Int_t           L1preNonisoIetam1[10];
  Int_t           L1preNonisoIphim1[10];
  Int_t           L1preNonisoRankm1[10];
  UInt_t          nbOfL1preIsoCandsm1;
  Int_t           L1preIsoIetam1[12];   
  Int_t           L1preIsoIphim1[12];   
  Int_t           L1preIsoRankm1[12];   
  UInt_t          nbOfL1preNonisoCandszero;
  Int_t           L1preNonisoIetazero[10];
  Int_t           L1preNonisoIphizero[10];
  Int_t           L1preNonisoRankzero[10];
  UInt_t          nbOfL1preIsoCandszero;
  Int_t           L1preIsoIetazero[12];   
  Int_t           L1preIsoIphizero[12];   
  Int_t           L1preIsoRankzero[12];   
  UInt_t          nbOfL1preNonisoCandsp1;
  Int_t           L1preNonisoIetap1[9];   
  Int_t           L1preNonisoIphip1[9];   
  Int_t           L1preNonisoRankp1[9];   
  UInt_t          nbOfL1preIsoCandsp1;
  Int_t           L1preIsoIetap1[12];   
  Int_t           L1preIsoIphip1[12];   
  Int_t           L1preIsoRankp1[12];   
  UInt_t          nbOfL1preNonisoCandsp2;
  Int_t           L1preNonisoIetap2[9]; 
  Int_t           L1preNonisoIphip2[9]; 
  Int_t           L1preNonisoRankp2[9]; 
  UInt_t          nbOfL1preIsoCandsp2;
  Int_t           L1preIsoIetap2[12];   
  Int_t           L1preIsoIphip2[12];   
  Int_t           L1preIsoRankp2[12];   
  
  
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

  
  // new variables 
  
  chain->SetBranchAddress("nbOfL1preNonisoCandsm2", &treeVars.nbOfL1preNonisoCandsm2);
  chain->SetBranchAddress("L1preNonisoIetam2", treeVars.L1preNonisoIetam2);
  chain->SetBranchAddress("L1preNonisoIphim2", treeVars.L1preNonisoIphim2);
  chain->SetBranchAddress("L1preNonisoRankm2", treeVars.L1preNonisoRankm2);
  chain->SetBranchAddress("nbOfL1preIsoCandsm2", &treeVars.nbOfL1preIsoCandsm2);
  chain->SetBranchAddress("L1preIsoIetam2", treeVars.L1preIsoIetam2);
  chain->SetBranchAddress("L1preIsoIphim2", treeVars.L1preIsoIphim2);
  chain->SetBranchAddress("L1preIsoRankm2", treeVars.L1preIsoRankm2);
  
  chain->SetBranchAddress("nbOfL1preNonisoCandsm1", &treeVars.nbOfL1preNonisoCandsm1);
  chain->SetBranchAddress("L1preNonisoIetam1", treeVars.L1preNonisoIetam1);
  chain->SetBranchAddress("L1preNonisoIphim1", treeVars.L1preNonisoIphim1);
  chain->SetBranchAddress("L1preNonisoRankm1", treeVars.L1preNonisoRankm1);
  chain->SetBranchAddress("nbOfL1preIsoCandsm1", &treeVars.nbOfL1preIsoCandsm1);
  chain->SetBranchAddress("L1preIsoIetam1", treeVars.L1preIsoIetam1);
  chain->SetBranchAddress("L1preIsoIphim1", treeVars.L1preIsoIphim1);
  chain->SetBranchAddress("L1preIsoRankm1", treeVars.L1preIsoRankm1);
  
  chain->SetBranchAddress("nbOfL1preNonisoCandszero", &treeVars.nbOfL1preNonisoCandszero);
  chain->SetBranchAddress("L1preNonisoIetazero", treeVars.L1preNonisoIetazero);
  chain->SetBranchAddress("L1preNonisoIphizero", treeVars.L1preNonisoIphizero);
  chain->SetBranchAddress("L1preNonisoRankzero", treeVars.L1preNonisoRankzero);
  chain->SetBranchAddress("nbOfL1preIsoCandszero", &treeVars.nbOfL1preIsoCandszero);
  chain->SetBranchAddress("L1preIsoIetazero", treeVars.L1preIsoIetazero);
  chain->SetBranchAddress("L1preIsoIphizero", treeVars.L1preIsoIphizero);
  chain->SetBranchAddress("L1preIsoRankzero", treeVars.L1preIsoRankzero);
  
  chain->SetBranchAddress("nbOfL1preNonisoCandsp1", &treeVars.nbOfL1preNonisoCandsp1);
  chain->SetBranchAddress("L1preNonisoIetap1", treeVars.L1preNonisoIetap1);
  chain->SetBranchAddress("L1preNonisoIphip1", treeVars.L1preNonisoIphip1);
  chain->SetBranchAddress("L1preNonisoRankp1", treeVars.L1preNonisoRankp1);
  chain->SetBranchAddress("nbOfL1preIsoCandsp1", &treeVars.nbOfL1preIsoCandsp1);
  chain->SetBranchAddress("L1preIsoIetap1", treeVars.L1preIsoIetap1);
  chain->SetBranchAddress("L1preIsoIphip1", treeVars.L1preIsoIphip1);
  chain->SetBranchAddress("L1preIsoRankp1", treeVars.L1preIsoRankp1);
  
  chain->SetBranchAddress("nbOfL1preNonisoCandsp2", &treeVars.nbOfL1preNonisoCandsp2);
  chain->SetBranchAddress("L1preNonisoIetap2", treeVars.L1preNonisoIetap2);
  chain->SetBranchAddress("L1preNonisoIphip2", treeVars.L1preNonisoIphip2);
  chain->SetBranchAddress("L1preNonisoRankp2", treeVars.L1preNonisoRankp2);
  chain->SetBranchAddress("nbOfL1preIsoCandsp2", &treeVars.nbOfL1preIsoCandsp2);
  chain->SetBranchAddress("L1preIsoIetap2", treeVars.L1preIsoIetap2);
  chain->SetBranchAddress("L1preIsoIphip2", treeVars.L1preIsoIphip2);
  chain->SetBranchAddress("L1preIsoRankp2", treeVars.L1preIsoRankp2);

  
  /*
  chain->SetBranchAddress("", & , ) ; 
  chain->SetBranchAddress("", , ) ;
  chain->SetBranchAddress("", , ) ;
  chain->SetBranchAddress("", , ) ;
  
  chain->SetBranchAddress("", , ) ;
  chain->SetBranchAddress("", , ) ;

  chain->SetBranchAddress("", & , ) ; 
  
  chain->SetBranchAddress("", , ) ;
  */
  
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
  TChain * chain ;
  
  if (is2017) chain =  new TChain ("tpAnalyzer/EcalTPGAnalysis") ;
  if (!is2017) chain =  new TChain ("EcalTPGAnalysis") ;

  // input TPGTree file
  
  TString inputrootfile;
  //if (is2017)  inputrootfile = "/eos/cms/store/user/pbarria/TPG/ECALTPGTree_ZeroBias2_Run2017F-Run306425_RAW-RECO.root"; // 2017 with prefiring 
  
  // rerun wioth fixed branches 
  if (is2017)  inputrootfile = "/eos/cms/store/user/khurana//Shalini/Run_2017/Merged_2017Data.root"; // 2017 with prefiring 
  
  //if (!is2017) inputrootfile = "/afs/cern.ch/user/s/sdutt/work/public/ECALTPGTree_subset.root" ;                          // 2018 without prefiring, full readout 
  //if (!is2017) inputrootfile = "/eos/cms/store/user/khurana/ECAL/TPGAnalysis/Merged_ECALTPGtree.root" ;                          // 2018 without prefiring, full readout 
   
  if (!is2017) inputrootfile = "/eos/cms/store/user/khurana//Shalini/Run_2018/Merged_2018.root";
  
  
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

  
  TH2F* ieta_vs_idx_L1Iso    = new TH2F("ieta_vs_idx_L1Iso","ieta_vs_idx_L1Iso",5, 0.5, 5.5, 140, -70.5,70.5);
  TH2F* ieta_vs_idx_L1NonIso = new TH2F("ieta_vs_idx_L1NonIso","ieta_vs_idx_L1NonIso",5, 0.5, 5.5, 140, -70.5,70.5);

  TH2F* ieta_vs_idx_L1Iso_binFix      = new TH2F("ieta_vs_idx_L1Iso_binFix","ieta_vs_idx_L1Iso_binFix",5, 0.5, 5.5, 140, -70.5,70.5);
  TH2F* ieta_vs_idx_L1NonIso_binFix   = new TH2F("ieta_vs_idx_L1NonIso_binFix","ieta_vs_idx_L1NonIso_binFix",5, 0.5, 5.5, 140, -70.5,70.5);
  

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
    //for (int entry = 15000 ; entry < 16000 ; ++entry) {
    
    chain->GetEntry (entry) ;
    
    
    if (entry%5000==0) cout << entry << " / " << treeentries
			    << " events processed" << endl;
    // for 2018: 4, 257, 336, 415, 811, 1151, 1230, 1309,  2045, 2124, 2302, 2596, 2939, 3117, 3097
    //     for 2017: 3, 46, 58, 70, 82, 121, 133, 145, 157, 172, 184, 196, 208, 223, 235, 247, 259, 274, 286, 298, 310, 349, 361, 373, 385, 400, 412, 424, 436, 451, 463, 475, 487, 502, 514, 526, 538, 577,
    // 589, 601, 613, 628, 640, 652, 664, 679, 691, 703, 715, 730, 742, 754, 766, 826, 838, 850, 862, 877, 889, 901, 913, 928, 940, 952, 964, 1003, 1015, 1027, 1039, 1054, 1066, 1078, 1090, 1105, 1117, 1129,
    // 1141, 1156, 1168, 1180, 1192, 1231, 1243, 1255, 1267, 1282, 1294, 1306, 1318, 1333, 1345, 1357, 1369, 1384, 1396, 1408, 1420, 1459, 1471, 1483, 1495, 1510, 1522, 1534, 1546, 1561, 1573, 1585, 1597, 
    //1612, 1624, 1636, 1648, 1720, 1732, 1744, 1756, 1771, 1783, 1795, 1807, 1822, 1834, 1846, 1858, 1897, 1909, 1921, 1933, 1948, 1960, 1972, 1984, 1999, 2011, 2023, 2035, 2050, 2062, 2074, 2086, 2125, 
    //2137, 2149, 2161, 2176, 2188, 2200, 2212, 2227, 2239, 2251, 2263, 2278, 2290, 2302, 2314, 2353, 2365, 2377, 2389, 2404, 2416, 2428, 2440, 2455, 2467, 2479, 2491, 2506, 2518, 2530, 2542, 2614, 2626, 
    //2638, 2650, 2665, 2677, 2689, 2701, 2716, 2728, 2740, 2752, 2791, 2803, 2815, 2827, 2842, 2854, 2866, 2878, 2893, 2905, 2917, 2929, 2944, 2956, 2968, 2980, 3019, 3031, 3043, 3055, 3070, 3082, 3094,
    // 3106, 3121, 3133, 3145, 3157, 3172, 3184, 3196, 3208, 3247, 3259, 3271, 3283, 3298, 3310, 3322, 3334, 3349, 3361, 3373, 3385, 3400, 3412, 3424, 3436, 
    
    

    bool goodBX = false;
    
    if (is2017) goodBX =  ( (treeVars.bxNb == 3) || (treeVars.bxNb == 46) || (treeVars.bxNb == 58) || (treeVars.bxNb == 70) || (treeVars.bxNb == 82) || (treeVars.bxNb == 121) || (treeVars.bxNb == 133) || (treeVars.bxNb == 145) || (treeVars.bxNb == 157) || (treeVars.bxNb == 172) || (treeVars.bxNb == 184) || (treeVars.bxNb == 196) || (treeVars.bxNb == 208) || (treeVars.bxNb == 223) || (treeVars.bxNb == 235) || (treeVars.bxNb == 247) || (treeVars.bxNb == 259) || (treeVars.bxNb == 274) || (treeVars.bxNb == 286) || (treeVars.bxNb == 298) || (treeVars.bxNb == 310) || (treeVars.bxNb == 349) || (treeVars.bxNb == 361) || (treeVars.bxNb == 373) || (treeVars.bxNb == 385) || (treeVars.bxNb == 400) || (treeVars.bxNb == 412) || (treeVars.bxNb == 424) || (treeVars.bxNb == 436) || (treeVars.bxNb == 451) || (treeVars.bxNb == 463) || (treeVars.bxNb == 475) || (treeVars.bxNb == 487) || (treeVars.bxNb == 502) || (treeVars.bxNb == 514) || (treeVars.bxNb == 526) || (treeVars.bxNb == 538) || (treeVars.bxNb == 577) || (treeVars.bxNb == 589) || (treeVars.bxNb == 601) || (treeVars.bxNb == 613) || (treeVars.bxNb == 628) || (treeVars.bxNb == 640) || (treeVars.bxNb == 652) || (treeVars.bxNb == 664) || (treeVars.bxNb == 679) || (treeVars.bxNb == 691) || (treeVars.bxNb == 703) || (treeVars.bxNb == 715) || (treeVars.bxNb == 730) || (treeVars.bxNb == 742) || (treeVars.bxNb == 754) || (treeVars.bxNb == 766) || (treeVars.bxNb == 826) || (treeVars.bxNb == 838) || (treeVars.bxNb == 850) || (treeVars.bxNb == 862) || (treeVars.bxNb == 877) || (treeVars.bxNb == 889) || (treeVars.bxNb == 901) || (treeVars.bxNb == 913) || (treeVars.bxNb == 928) || (treeVars.bxNb == 940) || (treeVars.bxNb == 952) || (treeVars.bxNb == 964) || (treeVars.bxNb == 1003) || (treeVars.bxNb == 1015) || (treeVars.bxNb == 1027) || (treeVars.bxNb == 1039) || (treeVars.bxNb == 1054) || (treeVars.bxNb == 1066) || (treeVars.bxNb == 1078) || (treeVars.bxNb == 1090) || (treeVars.bxNb == 1105) || (treeVars.bxNb == 1117) || (treeVars.bxNb == 1129) || (treeVars.bxNb == 1141) || (treeVars.bxNb == 1156) || (treeVars.bxNb == 1168) || (treeVars.bxNb == 1180) || (treeVars.bxNb == 1192) || (treeVars.bxNb == 1231) || (treeVars.bxNb == 1243) || (treeVars.bxNb == 1255) || (treeVars.bxNb == 1267) || (treeVars.bxNb == 1282) || (treeVars.bxNb == 1294) || (treeVars.bxNb == 1306) || (treeVars.bxNb == 1318) || (treeVars.bxNb == 1333) || (treeVars.bxNb == 1345) || (treeVars.bxNb == 1357) || (treeVars.bxNb == 1369) || (treeVars.bxNb == 1384) || (treeVars.bxNb == 1396) || (treeVars.bxNb == 1408) || (treeVars.bxNb == 1420) || (treeVars.bxNb == 1459) || (treeVars.bxNb == 1471) || (treeVars.bxNb == 1483) || (treeVars.bxNb == 1495) || (treeVars.bxNb == 1510) || (treeVars.bxNb == 1522) || (treeVars.bxNb == 1534) || (treeVars.bxNb == 1546) || (treeVars.bxNb == 1561) || (treeVars.bxNb == 1573) || (treeVars.bxNb == 1585) || (treeVars.bxNb == 1597) || (treeVars.bxNb == 1612) || (treeVars.bxNb == 1624) || (treeVars.bxNb == 1636) || (treeVars.bxNb == 1648) || (treeVars.bxNb == 1720) || (treeVars.bxNb == 1732) || (treeVars.bxNb == 1744) || (treeVars.bxNb == 1756) || (treeVars.bxNb == 1771) || (treeVars.bxNb == 1783) || (treeVars.bxNb == 1795) || (treeVars.bxNb == 1807) || (treeVars.bxNb == 1822) || (treeVars.bxNb == 1834) || (treeVars.bxNb == 1846) || (treeVars.bxNb == 1858) || (treeVars.bxNb == 1897) || (treeVars.bxNb == 1909) || (treeVars.bxNb == 1921) || (treeVars.bxNb == 1933) || (treeVars.bxNb == 1948) || (treeVars.bxNb == 1960) || (treeVars.bxNb == 1972) || (treeVars.bxNb == 1984) || (treeVars.bxNb == 1999) || (treeVars.bxNb == 2011) || (treeVars.bxNb == 2023) || (treeVars.bxNb == 2035) || (treeVars.bxNb == 2050) || (treeVars.bxNb == 2062) || (treeVars.bxNb == 2074) || (treeVars.bxNb == 2086) || (treeVars.bxNb == 2125) || (treeVars.bxNb == 2137) || (treeVars.bxNb == 2149) || (treeVars.bxNb == 2161) || (treeVars.bxNb == 2176) || (treeVars.bxNb == 2188) || (treeVars.bxNb == 2200) || (treeVars.bxNb == 2212) || (treeVars.bxNb == 2227) || (treeVars.bxNb == 2239) || (treeVars.bxNb == 2251) || (treeVars.bxNb == 2263) || (treeVars.bxNb == 2278) || (treeVars.bxNb == 2290) || (treeVars.bxNb == 2302) || (treeVars.bxNb == 2314) || (treeVars.bxNb == 2353) || (treeVars.bxNb == 2365) || (treeVars.bxNb == 2377) || (treeVars.bxNb == 2389) || (treeVars.bxNb == 2404) || (treeVars.bxNb == 2416) || (treeVars.bxNb == 2428) || (treeVars.bxNb == 2440) || (treeVars.bxNb == 2455) || (treeVars.bxNb == 2467) || (treeVars.bxNb == 2479) || (treeVars.bxNb == 2491) || (treeVars.bxNb == 2506) || (treeVars.bxNb == 2518) || (treeVars.bxNb == 2530) || (treeVars.bxNb == 2542) || (treeVars.bxNb == 2614) || (treeVars.bxNb == 2626) || (treeVars.bxNb == 2638) || (treeVars.bxNb == 2650) || (treeVars.bxNb == 2665) || (treeVars.bxNb == 2677) || (treeVars.bxNb == 2689) || (treeVars.bxNb == 2701) || (treeVars.bxNb == 2716) || (treeVars.bxNb == 2728) || (treeVars.bxNb == 2740) || (treeVars.bxNb == 2752) || (treeVars.bxNb == 2791) || (treeVars.bxNb == 2803) || (treeVars.bxNb == 2815) || (treeVars.bxNb == 2827) || (treeVars.bxNb == 2842) || (treeVars.bxNb == 2854) || (treeVars.bxNb == 2866) || (treeVars.bxNb == 2878) || (treeVars.bxNb == 2893) || (treeVars.bxNb == 2905) || (treeVars.bxNb == 2917) || (treeVars.bxNb == 2929) || (treeVars.bxNb == 2944) || (treeVars.bxNb == 2956) || (treeVars.bxNb == 2968) || (treeVars.bxNb == 2980) || (treeVars.bxNb == 3019) || (treeVars.bxNb == 3031) || (treeVars.bxNb == 3043) || (treeVars.bxNb == 3055) || (treeVars.bxNb == 3070) || (treeVars.bxNb == 3082) || (treeVars.bxNb == 3094) || (treeVars.bxNb == 3106) || (treeVars.bxNb == 3121) || (treeVars.bxNb == 3133) || (treeVars.bxNb == 3145) || (treeVars.bxNb == 3157) || (treeVars.bxNb == 3172) || (treeVars.bxNb == 3184) || (treeVars.bxNb == 3196) || (treeVars.bxNb == 3208) || (treeVars.bxNb == 3247) || (treeVars.bxNb == 3259) || (treeVars.bxNb == 3271) || (treeVars.bxNb == 3283) || (treeVars.bxNb == 3298) || (treeVars.bxNb == 3310) || (treeVars.bxNb == 3322) || (treeVars.bxNb == 3334) || (treeVars.bxNb == 3349) || (treeVars.bxNb == 3361) || (treeVars.bxNb == 3373) || (treeVars.bxNb == 3385) || (treeVars.bxNb == 3400) || (treeVars.bxNb == 3412) || (treeVars.bxNb == 3424) || (treeVars.bxNb == 3436) );  


    if (!is2017) goodBX =  ( (treeVars.bxNb == 4) || (treeVars.bxNb ==257) || (treeVars.bxNb ==336) ||(treeVars.bxNb ==415) ||(treeVars.bxNb ==811) ||(treeVars.bxNb ==1151) ||(treeVars.bxNb ==1230) ||(treeVars.bxNb ==1309) ||(treeVars.bxNb ==2045) ||(treeVars.bxNb ==2124) ||(treeVars.bxNb ==2302) ||(treeVars.bxNb ==2596) ||(treeVars.bxNb ==2939) ||(treeVars.bxNb ==3117) ||(treeVars.bxNb ==3097) ); 

    
    
    if (!goodBX) continue; 
    std::cout<<"  ------------------------------------------------------------------------------------ ------------------------------------------------------------------------------------ \n"<<std::endl;
    
    cout << "event number is " << entry << "  bxn is "<<treeVars.bxNb<<"\n" ;

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
      

      //if (tp > threshold) std::cout<<" tower number = "<<tower<<" with eta phi = "<<ieta <<"  "<<iphi <<"  "<<" and energy "<< tp << " emul max "<< maxOfTPEmul <<" at emul idx = "<<indexOfTPEmulMax<<std::endl;

      
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
      
      
    }//for (UInt_t tower = 0 ; tower < treeVars.nbOfTowers ; tower++) {
    
    /*(
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
    
    
    */
    
    
    

    // 1 refer to prefiring 
    // 2 refer to central bx
    // 3 refer to post firing 
    
    bool debug__ = false; 
    // commoin prin
    if (treeVars.nbOfL1IsoCands>0 || treeVars.nbOfL1NonisoCands > 0 || treeVars.nbOfL1preIsoCands >0 || treeVars.nbOfL1preNonisoCands > 0 || treeVars.nbOfL1postIsoCands > 0 || treeVars.nbOfL1postNonisoCands) {
      if(debug__) std::cout<<"number of nbOfL1IsoCands "<<treeVars.nbOfL1IsoCands
	       <<" nbOfL1NonisoCands "<<treeVars.nbOfL1NonisoCands
	       <<" nbOfL1preIsoCands "<<treeVars.nbOfL1preIsoCands
	       <<" nbOfL1preNonisoCands "<<treeVars.nbOfL1preNonisoCands
	       <<" nbOfL1postIsoCands "<<treeVars.nbOfL1postIsoCands
	       <<" nbOfL1postNonisoCands "<<treeVars.nbOfL1postNonisoCands
	       <<std::endl;
      
    }
    
    
    
    // prefiring isolated -2
    if (treeVars.nbOfL1preIsoCandsm2>0) {
      for (UInt_t iNL1preIsom2 = 0; iNL1preIsom2 < treeVars.nbOfL1preIsoCandsm2; iNL1preIsom2++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsom2
		 <<" eta = "<<treeVars.L1preIsoIetam2[iNL1preIsom2]
		 <<" phi = "<<treeVars.L1preIsoIphim2[iNL1preIsom2]
		 <<" energy = "<<treeVars.L1preIsoRankm2[iNL1preIsom2]
		 <<std::endl;
	if (treeVars.L1preIsoRankm2[iNL1preIsom2] > threshold){
	  ieta_vs_idx_L1Iso->Fill(1, treeVars.L1preIsoIetam2[iNL1preIsom2]);
	  ieta_vs_idx_L1Iso_binFix->Fill(1, (63./28.)*treeVars.L1preIsoIetam2[iNL1preIsom2] );
	}
      }
    }





    // prefiring isolated -1
    if (treeVars.nbOfL1preIsoCandsm1>0) {
      for (UInt_t iNL1preIsom1 = 0; iNL1preIsom1 < treeVars.nbOfL1preIsoCandsm1; iNL1preIsom1++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsom1
		 <<" eta = "<<treeVars.L1preIsoIetam1[iNL1preIsom1]
		 <<" phi = "<<treeVars.L1preIsoIphim1[iNL1preIsom1]
		 <<" energy = "<<treeVars.L1preIsoRankm1[iNL1preIsom1]
		 <<std::endl;
	if (treeVars.L1preIsoRankm1[iNL1preIsom1] > threshold) {
	  ieta_vs_idx_L1Iso->Fill(2, treeVars.L1preIsoIetam1[iNL1preIsom1]);
	  ieta_vs_idx_L1Iso_binFix->Fill(2, (63./28.)*treeVars.L1preIsoIetam1[iNL1preIsom1] );
	}
      }
    }



    // in time zero
    if (treeVars.nbOfL1preIsoCandszero>0) {
      for (UInt_t iNL1preIsozero = 0; iNL1preIsozero < treeVars.nbOfL1preIsoCandszero; iNL1preIsozero++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsozero
		 <<" eta = "<<treeVars.L1preIsoIetazero[iNL1preIsozero]
		 <<" phi = "<<treeVars.L1preIsoIphizero[iNL1preIsozero]
		 <<" energy = "<<treeVars.L1preIsoRankzero[iNL1preIsozero]
		 <<std::endl;
	if (treeVars.L1preIsoRankzero[iNL1preIsozero] > threshold){
	  ieta_vs_idx_L1Iso->Fill(3, treeVars.L1preIsoIetazero[iNL1preIsozero]);
	  ieta_vs_idx_L1Iso_binFix->Fill(3, (63./28.)*treeVars.L1preIsoIetazero[iNL1preIsozero] );
	}
      }
    }

    
    // postfiring isolated +1
    if (treeVars.nbOfL1preIsoCandsp1>0) {
      for (UInt_t iNL1preIsop1 = 0; iNL1preIsop1 < treeVars.nbOfL1preIsoCandsp1; iNL1preIsop1++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsop1
		 <<" eta = "<<treeVars.L1preIsoIetap1[iNL1preIsop1]
		 <<" phi = "<<treeVars.L1preIsoIphip1[iNL1preIsop1]
		 <<" energy = "<<treeVars.L1preIsoRankp1[iNL1preIsop1]
		 <<std::endl;
	if (treeVars.L1preIsoRankp1[iNL1preIsop1] > threshold) {
	  ieta_vs_idx_L1Iso->Fill(4, treeVars.L1preIsoIetap1[iNL1preIsop1]);
	  ieta_vs_idx_L1Iso_binFix->Fill(4, (63./28.)*treeVars.L1preIsoIetap1[iNL1preIsop1] );
	}
      }
    }


    // prefiring isolated +2
    if (treeVars.nbOfL1preIsoCandsp2>0) {
      for (UInt_t iNL1preIsop2 = 0; iNL1preIsop2 < treeVars.nbOfL1preIsoCandsp2; iNL1preIsop2++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsop2
		 <<" eta = "<<treeVars.L1preIsoIetap2[iNL1preIsop2]
		 <<" phi = "<<treeVars.L1preIsoIphip2[iNL1preIsop2]
		 <<" energy = "<<treeVars.L1preIsoRankp2[iNL1preIsop2]
		 <<std::endl;
	if (treeVars.L1preIsoRankp2[iNL1preIsop2] > threshold) {
	  ieta_vs_idx_L1Iso->Fill(5, treeVars.L1preIsoIetap2[iNL1preIsop2]);
	  ieta_vs_idx_L1Iso_binFix->Fill(5, (63./28.)*treeVars.L1preIsoIetap2[iNL1preIsop2] );
	}
      }
    }










    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // -----------------  NON Isolated -------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // prefiring isolated -2
    if (treeVars.nbOfL1preNonisoCandsm2>0) {
      for (UInt_t iNL1preIsom2 = 0; iNL1preIsom2 < treeVars.nbOfL1preNonisoCandsm2; iNL1preIsom2++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsom2
		 <<" eta = "<<treeVars.L1preNonisoIetam2[iNL1preIsom2]
		 <<" phi = "<<treeVars.L1preNonisoIphim2[iNL1preIsom2]
		 <<" energy = "<<treeVars.L1preNonisoRankm2[iNL1preIsom2]
		 <<std::endl;
	if (treeVars.L1preNonisoRankm2[iNL1preIsom2] > threshold){
	  ieta_vs_idx_L1NonIso->Fill(1, treeVars.L1preNonisoIetam2[iNL1preIsom2]);
	  ieta_vs_idx_L1NonIso_binFix->Fill(1, (63./28.)*treeVars.L1preNonisoIetam2[iNL1preIsom2]);
	}
      }
    }




    
    // prefiring isolated -1
    if (treeVars.nbOfL1preNonisoCandsm1>0) {
      for (UInt_t iNL1preIsom1 = 0; iNL1preIsom1 < treeVars.nbOfL1preNonisoCandsm1; iNL1preIsom1++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsom1
		 <<" eta = "<<treeVars.L1preNonisoIetam1[iNL1preIsom1]
		 <<" phi = "<<treeVars.L1preNonisoIphim1[iNL1preIsom1]
		 <<" energy = "<<treeVars.L1preNonisoRankm1[iNL1preIsom1]
		 <<std::endl;
	if (treeVars.L1preNonisoRankm1[iNL1preIsom1] > threshold) {
	  ieta_vs_idx_L1NonIso->Fill(2, treeVars.L1preNonisoIetam1[iNL1preIsom1]);
	  ieta_vs_idx_L1NonIso_binFix->Fill(2, (63./28.)*treeVars.L1preNonisoIetam1[iNL1preIsom1]);
	}
      }
    }



    // in time zero
    if (treeVars.nbOfL1preNonisoCandszero>0) {
      for (UInt_t iNL1preIsozero = 0; iNL1preIsozero < treeVars.nbOfL1preNonisoCandszero; iNL1preIsozero++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsozero
		 <<" eta = "<<treeVars.L1preNonisoIetazero[iNL1preIsozero]
		 <<" phi = "<<treeVars.L1preNonisoIphizero[iNL1preIsozero]
		 <<" energy = "<<treeVars.L1preNonisoRankzero[iNL1preIsozero]
		 <<std::endl;
	if (treeVars.L1preNonisoRankzero[iNL1preIsozero] > threshold){
	  ieta_vs_idx_L1NonIso->Fill(3, treeVars.L1preNonisoIetazero[iNL1preIsozero]);
	  ieta_vs_idx_L1NonIso_binFix->Fill(3, (63./28.)*treeVars.L1preNonisoIetazero[iNL1preIsozero]);
	}
      }
    }

    
    // postfiring isolated +1
    if (treeVars.nbOfL1preNonisoCandsp1>0) {
      for (UInt_t iNL1preIsop1 = 0; iNL1preIsop1 < treeVars.nbOfL1preNonisoCandsp1; iNL1preIsop1++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsop1
		 <<" eta = "<<treeVars.L1preNonisoIetap1[iNL1preIsop1]
		 <<" phi = "<<treeVars.L1preNonisoIphip1[iNL1preIsop1]
		 <<" energy = "<<treeVars.L1preNonisoRankp1[iNL1preIsop1]
		 <<std::endl;
	if (treeVars.L1preNonisoRankp1[iNL1preIsop1] > threshold){
	  ieta_vs_idx_L1NonIso->Fill(4, treeVars.L1preNonisoIetap1[iNL1preIsop1]);
	  ieta_vs_idx_L1NonIso_binFix->Fill(4, (63./28.)*treeVars.L1preNonisoIetap1[iNL1preIsop1]);
	}
      }
    }


    // prefiring isolated +2
    if (treeVars.nbOfL1preNonisoCandsp2>0) {
      for (UInt_t iNL1preIsop2 = 0; iNL1preIsop2 < treeVars.nbOfL1preNonisoCandsp2; iNL1preIsop2++) {
	if(debug__) std::cout<<" Isolated pre EGamma cadidates number "<<iNL1preIsop2
		 <<" eta = "<<treeVars.L1preNonisoIetap2[iNL1preIsop2]
		 <<" phi = "<<treeVars.L1preNonisoIphip2[iNL1preIsop2]
		 <<" energy = "<<treeVars.L1preNonisoRankp2[iNL1preIsop2]
		 <<std::endl;
	if (treeVars.L1preNonisoRankp2[iNL1preIsop2] > threshold) {
	  ieta_vs_idx_L1NonIso->Fill(5, treeVars.L1preNonisoIetap2[iNL1preIsop2]);
	  ieta_vs_idx_L1NonIso_binFix->Fill(5, (63./28.)*treeVars.L1preNonisoIetap2[iNL1preIsop2]);
	}
      }
    }


    /*



    // nominal isolated 
    if (treeVars.nbOfL1IsoCands>0) {
      for (int iNL1Iso = 0; iNL1Iso < treeVars.nbOfL1IsoCands; iNL1Iso++) {
	if(debug__) std::cout<<" Isolated EGamma cadidates number "<<iNL1Iso
		 <<" eta = "<<treeVars.L1IsoIeta[iNL1Iso]
		 <<" phi = "<<treeVars.L1IsoIphi[iNL1Iso]
		 <<" energy = "<<treeVars.L1IsoRank[iNL1Iso]
		 <<std::endl;
	if (treeVars.L1IsoRank[iNL1Iso] > threshold) ieta_vs_idx_L1Iso->Fill(2, treeVars.L1IsoIeta[iNL1Iso]);
      }
    }

    // postfiring isolated 
    if (treeVars.nbOfL1postIsoCands>0) {
      for (int iNL1postIso = 0; iNL1postIso < treeVars.nbOfL1postIsoCands; iNL1postIso++) {
	if(debug__) std::cout<<" Isolated post EGamma cadidates number "<<iNL1postIso
		 <<" eta = "<<treeVars.L1postIsoIeta[iNL1postIso]
		 <<" phi = "<<treeVars.L1postIsoIphi[iNL1postIso]
		 <<" energy = "<<treeVars.L1postIsoRank[iNL1postIso]
		 <<std::endl;
	if (treeVars.L1postIsoRank[iNL1postIso] > threshold) ieta_vs_idx_L1Iso->Fill(3, treeVars.L1postIsoIeta[iNL1postIso]);
      }
    }

    
    // prefiring non isolated 
    if (treeVars.nbOfL1preNonisoCands > 0) {
      for (int iNL1preNoniso = 0; iNL1preNoniso < treeVars.nbOfL1preNonisoCands; iNL1preNoniso++) {
	if(debug__) std::cout<<" Non Isolated EGamma pre cadidates number "<<iNL1preNoniso
		 <<" eta = "<<treeVars.L1preNonisoIeta[iNL1preNoniso]
		 <<" phi = "<<treeVars.L1preNonisoIphi[iNL1preNoniso]
		 <<" energy = "<<treeVars.L1preNonisoRank[iNL1preNoniso]
		 <<std::endl;
	if(treeVars.L1preNonisoRank[iNL1preNoniso]>threshold ) ieta_vs_idx_L1NonIso->Fill(1, treeVars.L1preNonisoIeta[iNL1preNoniso]);
      }
    }


    // nopminal non isolated 
    if (treeVars.nbOfL1NonisoCands > 0) {
      for (int iNL1Noniso = 0; iNL1Noniso < treeVars.nbOfL1NonisoCands; iNL1Noniso++) {
	if(debug__) std::cout<<" Non Isolated EGamma cadidates number "<<iNL1Noniso
		 <<" eta = "<<treeVars.L1NonisoIeta[iNL1Noniso]
		 <<" phi = "<<treeVars.L1NonisoIphi[iNL1Noniso]
		 <<" energy = "<<treeVars.L1NonisoRank[iNL1Noniso]
		 <<std::endl;
	if(treeVars.L1NonisoRank[iNL1Noniso]>threshold ) ieta_vs_idx_L1NonIso->Fill(2, treeVars.L1NonisoIeta[iNL1Noniso]);
      }
    }


    // postfiring non isolated 
    if (treeVars.nbOfL1postNonisoCands > 0) {
      for (int iNL1postNoniso = 0; iNL1postNoniso < treeVars.nbOfL1postNonisoCands; iNL1postNoniso++) {
	if(debug__) std::cout<<" Non Isolated EGamma post cadidates number "<<iNL1postNoniso
		 <<" eta = "<<treeVars.L1postNonisoIeta[iNL1postNoniso]
		 <<" phi = "<<treeVars.L1postNonisoIphi[iNL1postNoniso]
		 <<" energy = "<<treeVars.L1postNonisoRank[iNL1postNoniso]
		 <<std::endl;
	if(treeVars.L1postNonisoRank[iNL1postNoniso]>threshold ) ieta_vs_idx_L1NonIso->Fill(3, treeVars.L1postNonisoIeta[iNL1postNoniso]);
      }
    }
    */
    
    
    std::cout<< " \n \n ";
    
  }//for (int entry = 0 ; entry < treeentries ; ++entry) {
  
  TString outputfilename;
  TString cutval;
  cutval.Form("%d",threshold);
  if (is2017) outputfilename  = "PrefiringRateEE_2017data_FrontTrain_cutval_"+cutval+".root";
  if (!is2017) outputfilename = "PrefiringRateEE_2018data_FrontTrain_cutval_"+cutval+".root";
  
  TFile fout(outputfilename,"RECREATE");
  fout.cd();

  ieta_vs_iphi_TP->Write();
  ieta_vs_iphi_ETP->Write();
  
  ieta_vs_idx_L1Iso->Write();
  ieta_vs_idx_L1NonIso->Write();

  ieta_vs_idx_L1Iso_binFix->Write();
  ieta_vs_idx_L1NonIso_binFix->Write();
  
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

  
  for (int i=0; i<11; i++){
    TString postname;
    postname.Form("%d",i+18);
    std::cout<<" postname = "<<postname<<std::endl;
    idx_TP0_vec[i]->Write();
    idx_TP0_ETP_vec[i]->Write();
    
  }



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

