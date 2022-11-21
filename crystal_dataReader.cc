// crystal_dataReader.C - Macro to read the Giessen crystal test data into a ROOT tree format for follow up analysis.
//                        This is a heavily modified version of the PbWO4 crystal analysis macros used for ACCOS and 
//                        Giessen crystal tests, as written and modified in Genova (c. 2014) for use by E. Buchanan
//                        (crystal_ana.C) as a simple version for light transmission analysis, further modified for use
//                        in the second phase of Giessen crystal tests and in 2022 by a York summer student, M. Hawkins.
//
// Authors - R. De Vita, S. Fegan, E. Buchanan, M. Hawkins
//
// This version created on 18/11/22 in aftermath of summer project student work on the Giessen Crystal test data and codes
// Loads a csv data file and produces a root tree, with the intention of simplifying later analysis (see future macro) 
// csv file must contain the number of contained traces on the first line, and the name of the
// trace must have the form "PbWO-***-abc", where '***' is an  integer, and 'abc' is a trailing string indicating other information
// (these specifics of the csv file may vary)


#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <stdio.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <boost/tokenizer.hpp>


//#include "TCrystal.h"
#define Ncrystals 400
#define Nrad       25
#define Nann       17
#define Ntest     100
#define Nbroken    25
#define Ngie      400
#define Ngiewl   1151
#define maxScans 50

using namespace std;
using namespace boost;


void crystal_dataReader()
{


    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);

    gStyle->SetPadBorderMode(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadColor(10);

    gStyle->SetTitleFont(72, "X");
    gStyle->SetTitleFont(72, "Y");
    gStyle->SetTitleOffset(0.9, "X");
    gStyle->SetTitleOffset(1.2, "Y");
    gStyle->SetTitleSize(0.05, "X");
    gStyle->SetTitleSize(0.05, "Y");

    gStyle->SetLabelFont(72, "X");
    gStyle->SetLabelFont(72, "Y");
    gStyle->SetLabelFont(72, "Z");
    gStyle->SetPalette(1);
    gStyle->SetOptFit(111);
    gStyle->SetOptStat("nemri");
    //      gStyle->SetOptStat("");

    const Int_t colNum = 100;
    Int_t MyPalette[colNum];
    Double_t re[] = { 0.0, 0.5, 0.0, 1.0, 1.0, 0.5 };
    Double_t gr[] = { 0.0, 0.9, 1.0, 1.0, 0.0, 0.0 };
    Double_t bl[] = { 0.5, 1.0, 0.0, 0.0, 0.0, 0.0 };
    Double_t stop[] = { 0.0, 0.2, 0.5, 0.6, 0.8, 1.0 };
    Int_t FI = TColor::CreateGradientColorTable(6, stop, re, gr, bl, colNum);
    for (int i = 0; i < colNum; i++) MyPalette[i] = FI + i;
    gStyle->SetPalette(colNum, MyPalette);

    int    mystatus[Ncrystals][18] = { 0 };
    string myname[18] = { "AF","BF3","BF4","CF","AR","BR3","BR4","CR","L1","L2","L3","L4","LY","LT360","LT420","LT620","TTO","ALL" };
    float  myvalue[Ncrystals][18] = { 0 };
    float  mysigma[Ncrystals][18] = { 0 };
    int    mynumber[18] = { 0 };
    int    mytest[Ncrystals] = { 0 };
    int    mybroken[Ncrystals] = { 0 };


    // Create histos

    //exctracting results of Giessen test
    
    double gieID[Ngie] = { 0 };
    int gieSTAT[Ncrystals] = { 0 };
    double gieWL[Ngiewl];
    double giescale[Ngiewl];
    double gieratio[Ngiewl];
    double deltaK[Ngie];
    double gie_LTO[Ngie][Ngiewl] = { 0 };
    double gie_LTO_corr[Ngie][Ngiewl] = { 0 };
    double gie_LTO_ID[Ncrystals] = { 0 };
    double gie_LTO_stat[Ncrystals] = { 0 };
    double gie_LTO360_pre[Ncrystals] = { 0 };
    double gie_LTO360_post[Ncrystals] = { 0 };
    double gie_LTO360_acc[Ncrystals] = { 0 };
    double gie_LTO420_pre[Ncrystals] = { 0 };
    double gie_LTO420_post[Ncrystals] = { 0 };
    double gie_LTO420_acc[Ncrystals] = { 0 };
    double gie_LTO620_pre[Ncrystals] = { 0 };
    double gie_LTO620_post[Ncrystals] = { 0 };
    double gie_LTO620_acc[Ncrystals] = { 0 };
    double gie_DK[Ncrystals];
    double sic_DK[Ncrystals];
    double y[9][1151] = { 0 };
    double length = 0.2;
    FILE* gieFile;
    int  NgieFile = 0;


    cout << "\nReading Giessen data" << endl;
    //   for(int ifile=6; ifile<7; ifile++) {
    //     if(ifile==0) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX1_2_3_4_PROD.csv","r");
    //     if(ifile==1) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX5_6_PROD.csv","r");
    //     if(ifile==2) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX7_PROD.csv","r");
    //     if(ifile==3) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX7_149_irr_PROD_2.csv","r");
    //     if(ifile==4) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX8_PROD.csv","r");
    //     if(ifile==5) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/BOX9_10_PROD.csv","r");
    //     //gieFile = fopen ("/home/stuart/FT/Data/Crystal_337.csv","r");
    //     //if(ifile==6) gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/Box11_16-5__PROD.csv","r");
    //     if(ifile==6)

      //gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/Crystal_020.csv","r");
      //gieFile = fopen ("/home/stuart/FT/Crystals/Archive2/GiessenData/Crystal_023.csv","r");
      //gieFile = fopen ("/home/stuart/FTGiessenData/Crystal023.csv","r");
      //gieFile = fopen ("/home/stuart/FTGiessenData/Crystal020_023.csv","r");
      //gieFile = fopen ("/home/stuart/FT/GiessenIII/BOX_BTCP_PROD.csv","r");
    
    //gieFile = fopen("BOX7_PROD.csv", "r");
    //ifstream file("Box11-tuesday_PROD_truncated.csv", std::ifstream::in);
    ifstream file("/home/stuart/SideProjects/PWO/LT_Data/BOX1_2_3_4_PROD.csv", std::ifstream::in);
    
    //ifstream file3("BOX5_6_PROD.csv", std::ifstream::in);
    //ifstream file4("BOX8_PROD.csv", std::ifstream::in);
    //ifstream file5("BOX9_10_PROD.csv", std::ifstream::in);
    //file.open("BOX7_PROD.csv");
    //if (file.is_open()) {
        //cout << "File failed to open." << endl;
      //  return 1;
    //}
    int current_line = 0;
    string line;
    string word;
    istringstream iss;
    int nscan = 0;
    int nscan2 = 0;
    int crynum;
    int odd = 1;
    int counter = 0;
    int count = 0;
    int status;
    int crystal[161][161] = { 0 };
    std::pair <int, int> crystalInfo;
    vector < std::pair<int, int>> crystalData;
    

    

    while (!file.eof()) {
        
      //int counter = 0;
      current_line++;
      std::getline(file, line);
      //cout << line << endl;
      if (current_line == 1) {
	NgieFile = std::stoi(line);
	//cout << NgieFile << endl;
        
      }
      //if(file.is_open()) {
      //cout << "File failed to open." << endl;
      //  return 1;
      //}
      if (current_line == 2) {
	iss.str(line);
	while (iss.good()) {
	  iss >> word;
	  //cout << word <<endl;
	  nscan++;
	  char_separator<char> sep("_");
	  tokenizer<char_separator<char>> tokens(word, sep);
	  int Nword = 0;
	  for (const string& s : tokens) {
	    //cout << s << Nword << endl;
	    if (Nword == 1) {
	      crynum = std::stoi(s);
	      //cout << crynum << endl;
		cout << s << endl;
              
	    }
	    if (Nword == 2) {
	      //cout << s << endl;
	      if (s.compare("bef") == 0) {
		gieSTAT[nscan] = 0;
		cout << "bef" << endl;
		//gieSTAT[nscan] = status;
	      }
	      else if (s.compare("irr") == 0) {
		gieSTAT[nscan] = 1;
		cout << "irr" << endl;
		//gieSTAT[nscan] = status;
	      }
	    }
	    else {
	      //cout << "no thanks" << endl;
	    }
	    //cout << gieSTAT[nscan] << endl;

	    Nword++;
	  }
          
	  crystalInfo = make_pair(crynum, gieSTAT[nscan]);
	  crystalData.push_back(crystalInfo);
	  cout << crystalData.at(count).first << " " << crystalData.at(count).second << endl;
	  count++;
	}
      }
      //cout << crystalData.size() << endl;
      //cout << crystalData.at(1).first << endl;
      
      else {   
	odd = 1;
	//cout << odd << endl;
	char_separator<char> sep(" ");
	tokenizer<char_separator<char>> tokens(line, sep);
	
	for (const string& s : tokens) {
	  double gie_lto, gie_wl;
	  if (odd % 2 == 0) {
	    gie_lto = std::stod(s);
	    gie_LTO[(odd / 2) - 1][counter] = gie_lto;
	    // cout << odd << endl;
	  }
	  else {
	    gie_wl = std::stod(s);
	    gieWL[counter] = gie_wl;
	    //cout << odd << endl;
	  }
	  odd++;
	  
	}
	//cout << gie_LTO[(odd / 2) - 1][counter] << " " << gieWL[counter] << endl;
	
	counter++;
      }
    }


    for (int first = 0; first < NgieFile; first++){
      for (int j = 1; j < NgieFile; j++) {
	if (crystalData.at(first).first == crystalData.at(j).first) {
	  
	  // cout << crystalData.at(first).first << crystalData.at(first).second << crystalData.at(j).first << crystalData.at(j).second << endl;
	}
      }
      
    }

                
                
            
    
            
        
    
 
//     //int crystal[9] = { 145, 146, 136, 138, 139, 151, 162, 160, 156 };
//     //TGraph *deltas;
//     for (int iwl = 0; iwl < Ngiewl; iwl++) {//for every wavelength
        
        
//      y[0][iwl] = ((1 / length) * log(gie_LTO[0][iwl] / gie_LTO[9][iwl])); // compute delta k
//      y[1][iwl] = ((1 / length) * log(gie_LTO[1][iwl] / gie_LTO[10][iwl])); // compute delta k
//      y[2][iwl] = ((1 / length) * log(gie_LTO[2][iwl] / gie_LTO[11][iwl])); // compute delta k   
//      y[3][iwl] = ((1 / length) * log(gie_LTO[3][iwl] / gie_LTO[16][iwl])); // compute delta k
//      y[4][iwl] = ((1 / length) * log(gie_LTO[4][iwl] / gie_LTO[17][iwl])); // compute delta k    
//      y[5][iwl] = ((1 / length) * log(gie_LTO[5][iwl] / gie_LTO[14][iwl])); // compute delta k
//      y[6][iwl] = ((1 / length) * log(gie_LTO[6][iwl] / gie_LTO[12][iwl])); // compute delta k
//      y[7][iwl] = ((1 / length) * log(gie_LTO[7][iwl] / gie_LTO[13][iwl])); // compute delta k
//      y[8][iwl] = ((1 / length) * log(gie_LTO[8][iwl] / gie_LTO[15][iwl])); // compute delta k
    
    
//     }
    
    
//     TGraph* deltas[9];
//     for (int i = 0; i < 9; i++) {
//         deltas[i] = new TGraph(Ngiewl, gieWL, y[i]);
//         deltas[i]->GetXaxis()->SetTitle("Wavelength (nm)");
//         deltas[i]->GetYaxis()->SetTitle("DeltaK");
//         deltas[i]->SetTitle(Form("Crystal_%d", crystal[i]));
//         deltas[i]->GetYaxis()->SetRangeUser(-0.2, 1.5);


            
            


//     }
//     TCanvas* c90 = new TCanvas("c90", "Crystal 145 dk", 750, 1000);
    
//     deltas[0]->Draw("AP*");
//     c90->Print("crystal145dk.root");
//     TCanvas* c11a = new TCanvas("c11a", "Crystal 146 dk", 750, 1000);
//     deltas[1]->Draw("AP*");
//     c11a->Print("crystal146dk.pdf");
//     TCanvas* c11b = new TCanvas("c11b", "Crystal 136 dk", 750, 1000);
//     deltas[2]->Draw("AP*");
//     c11b->Print("crystal136dk.pdf");
//     TCanvas* c11c = new TCanvas("c11c", "Crystal 138 dk", 750, 1000);
//     deltas[3]->Draw("AP*");
//     c11c->Print("crystal138dk.pdf");
//     TCanvas* c11d = new TCanvas("c11d", "Crystal 139 dk", 750, 1000);
//     deltas[4]->Draw("AP*");
//     c11d->Print("crystal139dk.pdf");
//     TCanvas* c11e = new TCanvas("c11e", "Crystal 151 dk", 750, 1000);
//     deltas[5]->Draw("AP*");
//     c11e->Print("crystal151dk.pdf");
//     TCanvas* c11f = new TCanvas("c11", "Crystal 162 dk", 750, 1000);
//     deltas[6]->Draw("AP*");
//     c11f->Print("crystal162dk.pdf");
//     TCanvas* c11g = new TCanvas("c11g", "Crystal 160 dk", 750, 1000);
//     deltas[7]->Draw("AP*");
//     c11g->Print("crystal160dk.pdf");
//     TCanvas* c11h = new TCanvas("c11h", "Crystal 156 dk", 750, 1000);
//     deltas[8]->Draw("AP*");
//     c11h->Print("crystal156dk.pdf");
    
//     //deltakiwl[l] = new TGraph(Ngiewl, gieWL, delak[iwl]);

//     // compute delta k  



  
//   cout << "deltaK for 145 at 360 = " << ((1/length) * log (gie_LTO[0][1080]/gie_LTO[9][1080])) << endl; // compute delta k
//   cout << "deltaK for 146 at 360 = " << ((1 / length) * log(gie_LTO[1][1080] / gie_LTO[10][1080])) << endl; // compute delta k
//   cout << "deltaK for 136 at 360 = " << ((1 / length) * log(gie_LTO[2][1080] / gie_LTO[11][1080])) << endl; // compute delta k
//   cout << "deltaK for 138 at 360 = " << ((1 / length) * log(gie_LTO[3][1080] / gie_LTO[16][1080])) << endl; // compute delta k
//   cout << "deltaK for 139 at 360 = " << ((1 / length) * log(gie_LTO[4][1080] / gie_LTO[17][1080])) << endl; // compute delta k
//   cout << "deltaK for 151 at 360 = " << ((1 / length) * log(gie_LTO[5][1080] / gie_LTO[14][1080])) << endl; // compute delta k
//   cout << "deltaK for 162 at 360 = " << ((1 / length) * log(gie_LTO[6][1080] / gie_LTO[12][1080])) << endl; // compute delta k
//   cout << "deltaK for 160 at 360 = " << ((1 / length) * log(gie_LTO[7][1080] / gie_LTO[13][1080])) << endl; // compute delta k
//   cout << "deltaK for 156 at 360 = " << ((1 / length) * log(gie_LTO[8][1080] / gie_LTO[15][1080])) << endl; // compute delta k
  

  
//   cout << "deltaK for 145 at 420 = " << ((1/length) * log (gie_LTO[0][960]/gie_LTO[9][960])) << endl; // compute delta k
//   cout << "deltaK for 146 at 420 = " << ((1 / length) * log(gie_LTO[1][960] / gie_LTO[10][960])) << endl; // compute delta k
//   cout << "deltaK for 136 at 420= " << ((1 / length) * log(gie_LTO[2][960] / gie_LTO[11][960])) << endl; // compute delta k
//   cout << "deltaK for 138 at 420= " << ((1 / length) * log(gie_LTO[3][960] / gie_LTO[16][960])) << endl; // compute delta k
//   cout << "deltaK for 139 at 420= " << ((1 / length) * log(gie_LTO[4][960] / gie_LTO[17][960])) << endl; // compute delta k
//   cout << "deltaK for 151 at 420= " << ((1 / length) * log(gie_LTO[5][960] / gie_LTO[14][960])) << endl; // compute delta k
//   cout << "deltaK for 162 at 420= " << ((1 / length) * log(gie_LTO[6][960] / gie_LTO[12][960])) << endl; // compute delta k
//   cout << "deltaK for 160 at 420= " << ((1 / length) * log(gie_LTO[7][960] / gie_LTO[13][960])) << endl; // compute delta k
//   cout << "deltaK for 156 at 420= " << ((1 / length) * log(gie_LTO[8][960] / gie_LTO[15][960])) << endl; // compute delta k

//   cout << "deltaK for 145 at 620 = " << ((1 / length) * log(gie_LTO[0][560] / gie_LTO[9][560])) << endl; // compute delta k
//   cout << "deltaK for 146 at 620 = " << ((1 / length) * log(gie_LTO[1][560] / gie_LTO[10][560])) << endl; // compute delta k
//   cout << "deltaK for 136 at 620 = " << ((1 / length) * log(gie_LTO[2][560] / gie_LTO[11][560])) << endl; // compute delta k
//   cout << "deltaK for 138 at 620 = " << ((1 / length) * log(gie_LTO[3][560] / gie_LTO[16][560])) << endl; // compute delta k
//   cout << "deltaK for 139 at 620 = " << ((1 / length) * log(gie_LTO[4][560] / gie_LTO[17][560])) << endl; // compute delta k
//   cout << "deltaK for 151 at 620 = " << ((1 / length) * log(gie_LTO[5][560] / gie_LTO[14][560])) << endl; // compute delta k
//   cout << "deltaK for 162 at 620 = " << ((1 / length) * log(gie_LTO[6][560] / gie_LTO[12][560])) << endl; // compute delta k
//   cout << "deltaK for 160 at 620 = " << ((1 / length) * log(gie_LTO[7][560] / gie_LTO[13][560])) << endl; // compute delta k
//   cout << "deltaK for 156 at 620 = " << ((1 / length) * log(gie_LTO[8][560] / gie_LTO[15][560])) << endl; // compute delta k

//   cout << "light transmission ratio at 420 for crystal 145 irr " << gie_LTO[9][960] / gie_LTO[0][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 146 irr " << gie_LTO[10][960] / gie_LTO[1][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 136 irr " << gie_LTO[11][960] / gie_LTO[2][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 138 irr " << gie_LTO[16][960] / gie_LTO[3][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 139 irr " << gie_LTO[17][960] / gie_LTO[4][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 151 irr " << gie_LTO[14][960] / gie_LTO[5][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 162 irr " << gie_LTO[12][960] / gie_LTO[6][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 160 irr " << gie_LTO[13][960] / gie_LTO[7][960] << endl;
//   cout << "light transmission ratio at 420 for crystal 156 irr " << gie_LTO[15][960] / gie_LTO[8][960] << endl;
  
//   //wavelength[3] = {}
  
//   cout << "light transmission at 420 for crystal 156 irr " << gie_LTO[8][960] << endl;
  
  
//   //create canvasses
//   //TCanvas *canvasArray[Ncrystals];
//   TGraph *gr_LTO_bef[Ncrystals];
//   TGraph *gr_LTO_corr[Ncrystals];
//   TGraph *gr_LTO_irr[Ncrystals];
//   //TGraph *gr_LTO_ratio[Ncrystals];
//   //TGraph *gr_LTO_scale[Ncrystals];
//   //TGraph *gr_LTO_baseline[Ncrystals];
  
  
//   //cout << "Giessen Stat " << gieSTAT[cryNum][0] << endl;
  
//   for (int i = 0; i < Ncrystals; i++) {
 
//       gr_LTO_bef[i] = new TGraph(Ngiewl, gieWL, gie_LTO[i]);
      
//   }
//       for (int j = 0; j < 9; j++) {

          
//           //gr_LTO_bef[j]->SetTitle(Form("Crystal_%d", crystal[j]));


//           //gr_LTO_bef[i]->SetTitle(Form("Crystal_%i",i));
//           //gr_LTO_bef[i]->SetMarkerColor(i);
//           gr_LTO_bef[j]->GetXaxis()->SetTitle("Wavelength (nm)");
//           gr_LTO_bef[j]->GetYaxis()->SetTitle("LT (%)");


//       }

     
         
//   //Histograms
  
//   //   //     //c11->Divide(1,2); //divide the canvas into areas
//   //c11->cd(1);
  
//   gr_LTO_bef[0]->SetMarkerColor(1); //black
//   gr_LTO_bef[1]->SetMarkerColor(1); //red
//   gr_LTO_bef[2]->SetMarkerColor(1); //green
//   gr_LTO_bef[3]->SetMarkerColor(1); //blue
//   gr_LTO_bef[4]->SetMarkerColor(1); //yellow
//   gr_LTO_bef[5]->SetMarkerColor(1);
//   gr_LTO_bef[6]->SetMarkerColor(1); //black
//   gr_LTO_bef[7]->SetMarkerColor(1); //red
//   gr_LTO_bef[8]->SetMarkerColor(1); //green
//   gr_LTO_bef[9]->SetMarkerColor(2); //blue
//   gr_LTO_bef[10]->SetMarkerColor(2); //yellow
//   gr_LTO_bef[11]->SetMarkerColor(2);
//   gr_LTO_bef[12]->SetMarkerColor(2); //red
//   gr_LTO_bef[13]->SetMarkerColor(2); //green
//   gr_LTO_bef[14]->SetMarkerColor(2); //blue
//   gr_LTO_bef[15]->SetMarkerColor(2); //yellow
//   gr_LTO_bef[16]->SetMarkerColor(2);
//   gr_LTO_bef[17]->SetMarkerColor(2); //red
  

//   //gr_LTO_bef[24]->SetMarkerColor(1); //black
//   //gr_LTO_bef[25]->SetMarkerColor(2); //red
//   //gr_LTO_bef[26]->SetMarkerColor(3); //green

//   //gr_LTO_bef[27]->SetMarkerColor(5); //yellow
//   //gr_LTO_bef[28]->SetMarkerColor(6); //pink


  


//   TCanvas* c12 = new TCanvas("c12", "Crystal 145", 750, 1000);  //create canvas

//   gr_LTO_bef[0]->Draw("AP*");
//   gr_LTO_bef[9]->Draw("PSAME*");
//   c12->Print("crystal145.pdf");
//   TCanvas* c12a = new TCanvas("c12a", "Crystal 146", 750, 1000);  //create canvas

//   gr_LTO_bef[1]->Draw("AP*");
//   gr_LTO_bef[10]->Draw("PSAME*");
//   c12a->Print("crystal146.pdf");
//   TCanvas* c12b = new TCanvas("c12b", "Crystal 136", 750, 1000);  //create canvas

//   gr_LTO_bef[2]->Draw("AP*");
//   gr_LTO_bef[11]->Draw("PSAME*");
//   c12b->Print("crystal136.pdf");
//   TCanvas* c12c = new TCanvas("c12c", "Crystal 138", 750, 1000);  //create canvas

//   gr_LTO_bef[3]->Draw("AP*");
//   gr_LTO_bef[16]->Draw("PSAME*");
//   c12c->Print("crystal138.pdf");
//   TCanvas* c12d = new TCanvas("c12d", "Crystal 139", 750, 1000);  //create canvas

//   gr_LTO_bef[4]->Draw("AP*");
//   gr_LTO_bef[17]->Draw("PSAME*");
//   c12d->Print("crystal139.pdf");
//   TCanvas* c12e = new TCanvas("c12e", "Crystal 151", 750, 1000);  //create canvas

//   gr_LTO_bef[5]->Draw("AP*");
//   gr_LTO_bef[14]->Draw("PSAME*");
//   c12e->Print("crystal151.pdf");
//   TCanvas* c12f = new TCanvas("c12f", "Crystal 162", 750, 1000);  //create canvas

//   gr_LTO_bef[6]->Draw("AP*");
//   gr_LTO_bef[12]->Draw("PSAME*");
//   c12f->Print("crystal162.pdf");
//   TCanvas* c12g = new TCanvas("c12g", "Crystal 160", 750, 1000);  //create canvas

//   gr_LTO_bef[7]->Draw("AP*");
//   gr_LTO_bef[13]->Draw("PSAME*");
//   c12g->Print("crystal162.pdf");
//   TCanvas* c12h = new TCanvas("c12h", "Crystal 156", 750, 1000);  //create canvas

//   gr_LTO_bef[8]->Draw("AP*");
//   gr_LTO_bef[15]->Draw("PSAME*");
//   c12h->Print("crystal156.pdf");


  

  }
 
