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


void crystal_ana()
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
    // Crystal Size
    TH1F* hi_fwidth[4][2];
    for (int i = 0; i < 4; i++) {
        hi_fwidth[i][0] = new TH1F(Form("hi_fwidth_%i_0", i + 1), "", 200, 14.8, 15.2);
        hi_fwidth[i][0]->GetXaxis()->SetTitle("Width (mm)");
        hi_fwidth[i][0]->GetYaxis()->SetTitle("Nr. Crystals");
        hi_fwidth[i][1] = new TH1F(Form("hi_fwidth_%i_1", i + 1), "", 200, 0., 0.05);
        hi_fwidth[i][1]->GetXaxis()->SetTitle("(Max-Min)/2(mm)");
        hi_fwidth[i][1]->GetYaxis()->SetTitle("Nr. Crystals");
    }
    TH1F* hi_rwidth[4][2];
    for (int i = 0; i < 4; i++) {
        hi_rwidth[i][0] = new TH1F(Form("hi_rwidth_%i_0", i + 1), "", 200, 14.8, 15.2);
        hi_rwidth[i][0]->GetXaxis()->SetTitle("Width (mm)");
        hi_rwidth[i][0]->GetYaxis()->SetTitle("Nr. Crystals");
        hi_rwidth[i][1] = new TH1F(Form("hi_rwidth_%i_1", i + 1), "", 200, 0., 0.05);
        hi_rwidth[i][1]->GetXaxis()->SetTitle("(Max-Min)/2(mm)");
        hi_rwidth[i][1]->GetYaxis()->SetTitle("Nr. Crystals");
    }
    TH1F* hi_length[4][2];
    for (int i = 0; i < 4; i++) {
        hi_length[i][0] = new TH1F(Form("hi_length_%i_0", i + 1), "", 200, 199.8, 200.2);
        hi_length[i][0]->GetXaxis()->SetTitle("Width (mm)");
        hi_length[i][0]->GetYaxis()->SetTitle("Nr. Crystals");
        hi_length[i][1] = new TH1F(Form("hi_length_%i_1", i + 1), "", 200, 0., 0.05);
        hi_length[i][1]->GetXaxis()->SetTitle("(Max-Min)/2(mm)");
        hi_length[i][1]->GetYaxis()->SetTitle("Nr. Crystals");
    }
    TH1F* hi_ly_accos[2];
    hi_ly_accos[0] = new TH1F("hi_ly_accos_0", "", 200, 10., 20.);
    hi_ly_accos[0]->GetXaxis()->SetTitle("LY (p.e.)");
    hi_ly_accos[0]->GetYaxis()->SetTitle("Nr. Crystals");
    hi_ly_accos[1] = new TH1F("hi_ly_accos_1", "", 200, 0., 10.);
    hi_ly_accos[1]->GetXaxis()->SetTitle("Error (%)");
    hi_ly_accos[1]->GetYaxis()->SetTitle("Nr. Crystals");
    TH2F* hi_ly_accos_siccas_2d = new TH2F("hi_ly_accos_siccas_2d", "", 100, 10., 20., 100, 10., 20.);;
    hi_ly_accos_siccas_2d->GetXaxis()->SetTitle("SICCAS");
    hi_ly_accos_siccas_2d->GetYaxis()->SetTitle("ACCOS");
    TH1F* hi_ly_accos_siccas_1d = new TH1F("hi_ly_accos_siccas_1d", "", 200, -40., 40.);;
    hi_ly_accos_siccas_1d->GetXaxis()->SetTitle("ACCOS-SICCAS (%)");
    hi_ly_accos_siccas_1d->GetYaxis()->SetTitle("Nr. Crystals");
    TH2F* hi_ly_lt_accos = new TH2F("hi_ly_lt_accos", "", 100, 10., 50., 100, 10., 20.);;
    hi_ly_lt_accos->GetXaxis()->SetTitle("LT(360nm) (%)");
    hi_ly_lt_accos->GetYaxis()->SetTitle("LY (p.e.)");
    TH2F* hi_ly_lt_siccas = new TH2F("hi_ly_lt_siccas", "", 100, 10., 50., 100, 10., 20.);;
    hi_ly_lt_siccas->GetXaxis()->SetTitle("LT(360nm) (%)");
    hi_ly_lt_siccas->GetYaxis()->SetTitle("LY (p.e.)");

    TH1F* hi_trans_accos[4];
    TH2F* hi_trans_accos_siccas[4];
    for (int i = 0; i < 4; i++) {
        hi_trans_accos[i] = new TH1F(Form("hi_trans_accos_%i", i + 1), "", 200, 0., 100.);
        hi_trans_accos[i]->GetXaxis()->SetTitle("ACCOS LT (%)");
        hi_trans_accos[i]->GetYaxis()->SetTitle("Nr. Crystals");
        hi_trans_accos_siccas[i] = new TH2F(Form("hi_trans_accos_siccas_%i", i + 1), "", 200, 10., 90., 200, 10., 90.);
        hi_trans_accos_siccas[i]->GetXaxis()->SetTitle("SICCAS LT (%)");
        hi_trans_accos_siccas[i]->GetYaxis()->SetTitle("ACCOS LT (%)");
    }

    TH2F* hi_meantrTTO = new TH2F("hi_meantrTTO", "", 50, 300., 800., 200, 0., 100.);
    hi_meantrTTO->GetXaxis()->SetTitle("#lambda (nm)");
    hi_meantrTTO->GetYaxis()->SetTitle("Transverse Transmission (%)");
    TH2F* hi_errtrTTO = new TH2F("hi_errtrTTO", "", 50, 300., 800., 200, 0., 5.);
    hi_errtrTTO->GetXaxis()->SetTitle("#lambda (nm)");
    hi_errtrTTO->GetYaxis()->SetTitle("#sigma(Transverse Transmission) (%)");
    TH2F* hi_wl50trTTO = new TH2F("hi_wl50trTTO", "", 50, 0., 20., 200, 300., 400.);
    hi_wl50trTTO->GetXaxis()->SetTitle("x (cm)");
    hi_wl50trTTO->GetYaxis()->SetTitle("#lambda(50%) (nm)");
    TH1F* hi_errwl50trTTO = new TH1F("hi_errwl50trTTO", "", 200, 0., 5.);
    hi_errwl50trTTO->GetXaxis()->SetTitle("#delta(#lambda) (nm)");
    hi_errwl50trTTO->GetYaxis()->SetTitle("Nr. Crystals");

    TH1F* hi_dk_siccas = new TH1F("hi_dk_siccas", "", 120, 0., 1.2);
    hi_dk_siccas->GetXaxis()->SetTitle("dK(m^{-1})_{SICCAS}");
    hi_dk_siccas->GetYaxis()->SetTitle("Nr. Crystals");
    TH2F* hi_dk_lt_siccas = new TH2F("hi_dk_lt_siccas", "", 200, 20., 40., 120, 0., 1.2);
    hi_dk_lt_siccas->GetXaxis()->SetTitle("LT(360nm)_{SICCAS} (%)");
    hi_dk_lt_siccas->GetYaxis()->SetTitle("dK(m^{-1})_{SICCAS}");




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
    int  NgieFile;


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
	cout << line << endl;
        if (current_line == 1) {
            NgieFile = std::stoi(line);
            cout << NgieFile << endl;
            
        }
        //if(file.is_open()) {
            //cout << "File failed to open." << endl;
            //  return 1;
        //}
        if (current_line == 2) {
            iss.str(line);
            while (iss.good()) {
                iss >> word;
                nscan++;
                char_separator<char> sep("_");
                tokenizer<char_separator<char>> tokens(word, sep);
                int Nword = 0;
                for (const string& s : tokens) {
                    //cout << s << Nword << endl;
                    if (Nword == 1) {
                        crynum = std::stoi(s);
                        cout << crynum << endl;
                    
                    //}
                    if (Nword == 2) {
                        if (string(s) == "bef") {
                            gieSTAT[nscan] = 0;
                            //gieSTAT[nscan] = status;
                        }
                        else if (string(s) == "irr") {
                            gieSTAT[nscan] = 1;
                            //gieSTAT[nscan] = status;
                        }
                        else {
                            cout << "no thanks" << endl;
                        }
                        //cout << gieSTAT[nscan] << endl;
                        
                    }
                    Nword++;
                    
                }
                crystalInfo = make_pair(crynum, gieSTAT[nscan]);
                crystalData.push_back(crystalInfo);
                cout << crystalData.at(count).first << " " << crystalData.at(count).second << endl;
                count++;
		}
    
	    }
	}
        else 
        {   
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

                cout << crystalData.at(first).first << crystalData.at(first).second << crystalData.at(j).first << crystalData.at(j).second << endl;
            }
        }

    }
                
                
            
    
            
        
    
 
    //int crystal[9] = { 145, 146, 136, 138, 139, 151, 162, 160, 156 };
    //TGraph *deltas;
    for (int iwl = 0; iwl < Ngiewl; iwl++) {//for every wavelength
        
        
     y[0][iwl] = ((1 / length) * log(gie_LTO[0][iwl] / gie_LTO[9][iwl])); // compute delta k
     y[1][iwl] = ((1 / length) * log(gie_LTO[1][iwl] / gie_LTO[10][iwl])); // compute delta k
     y[2][iwl] = ((1 / length) * log(gie_LTO[2][iwl] / gie_LTO[11][iwl])); // compute delta k   
     y[3][iwl] = ((1 / length) * log(gie_LTO[3][iwl] / gie_LTO[16][iwl])); // compute delta k
     y[4][iwl] = ((1 / length) * log(gie_LTO[4][iwl] / gie_LTO[17][iwl])); // compute delta k    
     y[5][iwl] = ((1 / length) * log(gie_LTO[5][iwl] / gie_LTO[14][iwl])); // compute delta k
     y[6][iwl] = ((1 / length) * log(gie_LTO[6][iwl] / gie_LTO[12][iwl])); // compute delta k
     y[7][iwl] = ((1 / length) * log(gie_LTO[7][iwl] / gie_LTO[13][iwl])); // compute delta k
     y[8][iwl] = ((1 / length) * log(gie_LTO[8][iwl] / gie_LTO[15][iwl])); // compute delta k
    
    
    }
    
    
    TGraph* deltas[9];
    for (int i = 0; i < 9; i++) {
        deltas[i] = new TGraph(Ngiewl, gieWL, y[i]);
        deltas[i]->GetXaxis()->SetTitle("Wavelength (nm)");
        deltas[i]->GetYaxis()->SetTitle("DeltaK");
        deltas[i]->SetTitle(Form("Crystal_%d", crystal[i]));
        deltas[i]->GetYaxis()->SetRangeUser(-0.2, 1.5);


            
            


    }
    TCanvas* c90 = new TCanvas("c90", "Crystal 145 dk", 750, 1000);
    
    deltas[0]->Draw("AP*");
    c90->Print("crystal145dk.root");
    TCanvas* c11a = new TCanvas("c11a", "Crystal 146 dk", 750, 1000);
    deltas[1]->Draw("AP*");
    c11a->Print("crystal146dk.pdf");
    TCanvas* c11b = new TCanvas("c11b", "Crystal 136 dk", 750, 1000);
    deltas[2]->Draw("AP*");
    c11b->Print("crystal136dk.pdf");
    TCanvas* c11c = new TCanvas("c11c", "Crystal 138 dk", 750, 1000);
    deltas[3]->Draw("AP*");
    c11c->Print("crystal138dk.pdf");
    TCanvas* c11d = new TCanvas("c11d", "Crystal 139 dk", 750, 1000);
    deltas[4]->Draw("AP*");
    c11d->Print("crystal139dk.pdf");
    TCanvas* c11e = new TCanvas("c11e", "Crystal 151 dk", 750, 1000);
    deltas[5]->Draw("AP*");
    c11e->Print("crystal151dk.pdf");
    TCanvas* c11f = new TCanvas("c11", "Crystal 162 dk", 750, 1000);
    deltas[6]->Draw("AP*");
    c11f->Print("crystal162dk.pdf");
    TCanvas* c11g = new TCanvas("c11g", "Crystal 160 dk", 750, 1000);
    deltas[7]->Draw("AP*");
    c11g->Print("crystal160dk.pdf");
    TCanvas* c11h = new TCanvas("c11h", "Crystal 156 dk", 750, 1000);
    deltas[8]->Draw("AP*");
    c11h->Print("crystal156dk.pdf");
    
    //deltakiwl[l] = new TGraph(Ngiewl, gieWL, delak[iwl]);

    // compute delta k  



  
  cout << "deltaK for 145 at 360 = " << ((1/length) * log (gie_LTO[0][1080]/gie_LTO[9][1080])) << endl; // compute delta k
  cout << "deltaK for 146 at 360 = " << ((1 / length) * log(gie_LTO[1][1080] / gie_LTO[10][1080])) << endl; // compute delta k
  cout << "deltaK for 136 at 360 = " << ((1 / length) * log(gie_LTO[2][1080] / gie_LTO[11][1080])) << endl; // compute delta k
  cout << "deltaK for 138 at 360 = " << ((1 / length) * log(gie_LTO[3][1080] / gie_LTO[16][1080])) << endl; // compute delta k
  cout << "deltaK for 139 at 360 = " << ((1 / length) * log(gie_LTO[4][1080] / gie_LTO[17][1080])) << endl; // compute delta k
  cout << "deltaK for 151 at 360 = " << ((1 / length) * log(gie_LTO[5][1080] / gie_LTO[14][1080])) << endl; // compute delta k
  cout << "deltaK for 162 at 360 = " << ((1 / length) * log(gie_LTO[6][1080] / gie_LTO[12][1080])) << endl; // compute delta k
  cout << "deltaK for 160 at 360 = " << ((1 / length) * log(gie_LTO[7][1080] / gie_LTO[13][1080])) << endl; // compute delta k
  cout << "deltaK for 156 at 360 = " << ((1 / length) * log(gie_LTO[8][1080] / gie_LTO[15][1080])) << endl; // compute delta k
  

  
  cout << "deltaK for 145 at 420 = " << ((1/length) * log (gie_LTO[0][960]/gie_LTO[9][960])) << endl; // compute delta k
  cout << "deltaK for 146 at 420 = " << ((1 / length) * log(gie_LTO[1][960] / gie_LTO[10][960])) << endl; // compute delta k
  cout << "deltaK for 136 at 420= " << ((1 / length) * log(gie_LTO[2][960] / gie_LTO[11][960])) << endl; // compute delta k
  cout << "deltaK for 138 at 420= " << ((1 / length) * log(gie_LTO[3][960] / gie_LTO[16][960])) << endl; // compute delta k
  cout << "deltaK for 139 at 420= " << ((1 / length) * log(gie_LTO[4][960] / gie_LTO[17][960])) << endl; // compute delta k
  cout << "deltaK for 151 at 420= " << ((1 / length) * log(gie_LTO[5][960] / gie_LTO[14][960])) << endl; // compute delta k
  cout << "deltaK for 162 at 420= " << ((1 / length) * log(gie_LTO[6][960] / gie_LTO[12][960])) << endl; // compute delta k
  cout << "deltaK for 160 at 420= " << ((1 / length) * log(gie_LTO[7][960] / gie_LTO[13][960])) << endl; // compute delta k
  cout << "deltaK for 156 at 420= " << ((1 / length) * log(gie_LTO[8][960] / gie_LTO[15][960])) << endl; // compute delta k

  cout << "deltaK for 145 at 620 = " << ((1 / length) * log(gie_LTO[0][560] / gie_LTO[9][560])) << endl; // compute delta k
  cout << "deltaK for 146 at 620 = " << ((1 / length) * log(gie_LTO[1][560] / gie_LTO[10][560])) << endl; // compute delta k
  cout << "deltaK for 136 at 620 = " << ((1 / length) * log(gie_LTO[2][560] / gie_LTO[11][560])) << endl; // compute delta k
  cout << "deltaK for 138 at 620 = " << ((1 / length) * log(gie_LTO[3][560] / gie_LTO[16][560])) << endl; // compute delta k
  cout << "deltaK for 139 at 620 = " << ((1 / length) * log(gie_LTO[4][560] / gie_LTO[17][560])) << endl; // compute delta k
  cout << "deltaK for 151 at 620 = " << ((1 / length) * log(gie_LTO[5][560] / gie_LTO[14][560])) << endl; // compute delta k
  cout << "deltaK for 162 at 620 = " << ((1 / length) * log(gie_LTO[6][560] / gie_LTO[12][560])) << endl; // compute delta k
  cout << "deltaK for 160 at 620 = " << ((1 / length) * log(gie_LTO[7][560] / gie_LTO[13][560])) << endl; // compute delta k
  cout << "deltaK for 156 at 620 = " << ((1 / length) * log(gie_LTO[8][560] / gie_LTO[15][560])) << endl; // compute delta k

  cout << "light transmission ratio at 420 for crystal 145 irr " << gie_LTO[9][960] / gie_LTO[0][960] << endl;
  cout << "light transmission ratio at 420 for crystal 146 irr " << gie_LTO[10][960] / gie_LTO[1][960] << endl;
  cout << "light transmission ratio at 420 for crystal 136 irr " << gie_LTO[11][960] / gie_LTO[2][960] << endl;
  cout << "light transmission ratio at 420 for crystal 138 irr " << gie_LTO[16][960] / gie_LTO[3][960] << endl;
  cout << "light transmission ratio at 420 for crystal 139 irr " << gie_LTO[17][960] / gie_LTO[4][960] << endl;
  cout << "light transmission ratio at 420 for crystal 151 irr " << gie_LTO[14][960] / gie_LTO[5][960] << endl;
  cout << "light transmission ratio at 420 for crystal 162 irr " << gie_LTO[12][960] / gie_LTO[6][960] << endl;
  cout << "light transmission ratio at 420 for crystal 160 irr " << gie_LTO[13][960] / gie_LTO[7][960] << endl;
  cout << "light transmission ratio at 420 for crystal 156 irr " << gie_LTO[15][960] / gie_LTO[8][960] << endl;
  
  //wavelength[3] = {}
  
  cout << "light transmission at 420 for crystal 156 irr " << gie_LTO[8][960] << endl;
  
  
  //create canvasses
  //TCanvas *canvasArray[Ncrystals];
  TGraph *gr_LTO_bef[Ncrystals];
  TGraph *gr_LTO_corr[Ncrystals];
  TGraph *gr_LTO_irr[Ncrystals];
  //TGraph *gr_LTO_ratio[Ncrystals];
  //TGraph *gr_LTO_scale[Ncrystals];
  //TGraph *gr_LTO_baseline[Ncrystals];
  
  
  //cout << "Giessen Stat " << gieSTAT[cryNum][0] << endl;
  
  for (int i = 0; i < Ncrystals; i++) {
 
      gr_LTO_bef[i] = new TGraph(Ngiewl, gieWL, gie_LTO[i]);
      
  }
      for (int j = 0; j < 9; j++) {

          
          //gr_LTO_bef[j]->SetTitle(Form("Crystal_%d", crystal[j]));


          //gr_LTO_bef[i]->SetTitle(Form("Crystal_%i",i));
          //gr_LTO_bef[i]->SetMarkerColor(i);
          gr_LTO_bef[j]->GetXaxis()->SetTitle("Wavelength (nm)");
          gr_LTO_bef[j]->GetYaxis()->SetTitle("LT (%)");


      }

     
         
  //Histograms
  
  //   //     //c11->Divide(1,2); //divide the canvas into areas
  //c11->cd(1);
  
  gr_LTO_bef[0]->SetMarkerColor(1); //black
  gr_LTO_bef[1]->SetMarkerColor(1); //red
  gr_LTO_bef[2]->SetMarkerColor(1); //green
  gr_LTO_bef[3]->SetMarkerColor(1); //blue
  gr_LTO_bef[4]->SetMarkerColor(1); //yellow
  gr_LTO_bef[5]->SetMarkerColor(1);
  gr_LTO_bef[6]->SetMarkerColor(1); //black
  gr_LTO_bef[7]->SetMarkerColor(1); //red
  gr_LTO_bef[8]->SetMarkerColor(1); //green
  gr_LTO_bef[9]->SetMarkerColor(2); //blue
  gr_LTO_bef[10]->SetMarkerColor(2); //yellow
  gr_LTO_bef[11]->SetMarkerColor(2);
  gr_LTO_bef[12]->SetMarkerColor(2); //red
  gr_LTO_bef[13]->SetMarkerColor(2); //green
  gr_LTO_bef[14]->SetMarkerColor(2); //blue
  gr_LTO_bef[15]->SetMarkerColor(2); //yellow
  gr_LTO_bef[16]->SetMarkerColor(2);
  gr_LTO_bef[17]->SetMarkerColor(2); //red
  

  //gr_LTO_bef[24]->SetMarkerColor(1); //black
  //gr_LTO_bef[25]->SetMarkerColor(2); //red
  //gr_LTO_bef[26]->SetMarkerColor(3); //green

  //gr_LTO_bef[27]->SetMarkerColor(5); //yellow
  //gr_LTO_bef[28]->SetMarkerColor(6); //pink


  


  TCanvas* c12 = new TCanvas("c12", "Crystal 145", 750, 1000);  //create canvas

  gr_LTO_bef[0]->Draw("AP*");
  gr_LTO_bef[9]->Draw("PSAME*");
  c12->Print("crystal145.pdf");
  TCanvas* c12a = new TCanvas("c12a", "Crystal 146", 750, 1000);  //create canvas

  gr_LTO_bef[1]->Draw("AP*");
  gr_LTO_bef[10]->Draw("PSAME*");
  c12a->Print("crystal146.pdf");
  TCanvas* c12b = new TCanvas("c12b", "Crystal 136", 750, 1000);  //create canvas

  gr_LTO_bef[2]->Draw("AP*");
  gr_LTO_bef[11]->Draw("PSAME*");
  c12b->Print("crystal136.pdf");
  TCanvas* c12c = new TCanvas("c12c", "Crystal 138", 750, 1000);  //create canvas

  gr_LTO_bef[3]->Draw("AP*");
  gr_LTO_bef[16]->Draw("PSAME*");
  c12c->Print("crystal138.pdf");
  TCanvas* c12d = new TCanvas("c12d", "Crystal 139", 750, 1000);  //create canvas

  gr_LTO_bef[4]->Draw("AP*");
  gr_LTO_bef[17]->Draw("PSAME*");
  c12d->Print("crystal139.pdf");
  TCanvas* c12e = new TCanvas("c12e", "Crystal 151", 750, 1000);  //create canvas

  gr_LTO_bef[5]->Draw("AP*");
  gr_LTO_bef[14]->Draw("PSAME*");
  c12e->Print("crystal151.pdf");
  TCanvas* c12f = new TCanvas("c12f", "Crystal 162", 750, 1000);  //create canvas

  gr_LTO_bef[6]->Draw("AP*");
  gr_LTO_bef[12]->Draw("PSAME*");
  c12f->Print("crystal162.pdf");
  TCanvas* c12g = new TCanvas("c12g", "Crystal 160", 750, 1000);  //create canvas

  gr_LTO_bef[7]->Draw("AP*");
  gr_LTO_bef[13]->Draw("PSAME*");
  c12g->Print("crystal162.pdf");
  TCanvas* c12h = new TCanvas("c12h", "Crystal 156", 750, 1000);  //create canvas

  gr_LTO_bef[8]->Draw("AP*");
  gr_LTO_bef[15]->Draw("PSAME*");
  c12h->Print("crystal156.pdf");
  //new crystals
  //
  //
  //
  
  /*
  current_line = 0;
  
  nscan = 0;
  ifstream file2("BOX1_2_3_4_PROD.csv", std::ifstream::in);
  odd = 1;
  counter = 0;
  int counter2 = 0;
  while (!file2.eof()) {

      //int counter = 0;
      current_line++;
      getline(file2, line);
      if (current_line == 1) {
          NgieFile = std::stoi(line);
          cout << NgieFile << endl;

      }
      //if(file.is_open()) {
          //cout << "File failed to open." << endl;
          //  return 1;
      //}
      if (current_line == 2) {
          iss.str(line);
          while (iss.good()) {
              iss >> word;
              nscan++;
              char_separator<char> sep("_");
              tokenizer<char_separator<char>> tokens(word, sep);
              int Nword = 0;
              for (const string& s : tokens) {
                  //cout << s << Nword << endl;
                  if (Nword == 1) {
                      crynum = std::stoi(s);
                      cout << crynum << " " << counter2 << endl;
                      counter2++;
                  }
                  if (Nword == 2) {
                      cout << s << endl;
                  }
                  Nword++;

              }

          }
      }
      else
      {
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


  
  int crystal[39] = { 49, 37, 40, 41, 61, 56, 48, 53, 52 , 47, 215, 208, 231, 240, 206, 237, 233, 214, 236, 232, 244, 242, 259, 246, 265, 247, 252, 251, 254, 269,304, 322, 282, 276, 274, 268, 288, 298, 330 };
  //TGraph *deltas;
  for (int iwl = 0; iwl < Ngiewl; iwl++) {//for every wavelength


      y[0][iwl] = ((1 / length) * log(gie_LTO[0][iwl] / gie_LTO[11][iwl])); // compute delta k
      y[1][iwl] = ((1 / length) * log(gie_LTO[1][iwl] / gie_LTO[12][iwl])); // compute delta k
      y[2][iwl] = ((1 / length) * log(gie_LTO[2][iwl] / gie_LTO[10][iwl])); // compute delta k   
      y[3][iwl] = ((1 / length) * log(gie_LTO[3][iwl] / gie_LTO[15][iwl])); // compute delta k
      y[4][iwl] = ((1 / length) * log(gie_LTO[4][iwl] / gie_LTO[14][iwl])); // compute delta k    
      y[5][iwl] = ((1 / length) * log(gie_LTO[5][iwl] / gie_LTO[13][iwl])); // compute delta k
      y[6][iwl] = ((1 / length) * log(gie_LTO[6][iwl] / gie_LTO[16][iwl])); // compute delta k
      y[7][iwl] = ((1 / length) * log(gie_LTO[7][iwl] / gie_LTO[17][iwl])); // compute delta k
      y[8][iwl] = ((1 / length) * log(gie_LTO[8][iwl] / gie_LTO[19][iwl])); // compute delta k
      y[9][iwl] = ((1 / length) * log(gie_LTO[9][iwl] / gie_LTO[18][iwl])); // compute delta k
      y[10][iwl] = ((1 / length) * log(gie_LTO[20][iwl] / gie_LTO[38][iwl])); // compute delta k
      y[11][iwl] = ((1 / length) * log(gie_LTO[21][iwl] / gie_LTO[36][iwl])); // compute delta k   
      y[12][iwl] = ((1 / length) * log(gie_LTO[22][iwl] / gie_LTO[34][iwl])); // compute delta k
      y[13][iwl] = ((1 / length) * log(gie_LTO[23][iwl] / gie_LTO[43][iwl])); // compute delta k    
      y[14][iwl] = ((1 / length) * log(gie_LTO[24][iwl] / gie_LTO[63][iwl])); // compute delta k
      y[15][iwl] = ((1 / length) * log(gie_LTO[25][iwl] / gie_LTO[39][iwl])); // compute delta k
      y[16][iwl] = ((1 / length) * log(gie_LTO[26][iwl] / gie_LTO[42][iwl])); // compute delta k
      y[17][iwl] = ((1 / length) * log(gie_LTO[27][iwl] / gie_LTO[47][iwl])); // compute delta k
      y[18][iwl] = ((1 / length) * log(gie_LTO[28][iwl] / gie_LTO[64][iwl])); // compute delta k
      y[19][iwl] = ((1 / length) * log(gie_LTO[29][iwl] / gie_LTO[49][iwl])); // compute delta k
      y[20][iwl] = ((1 / length) * log(gie_LTO[30][iwl] / gie_LTO[52][iwl])); // compute delta k   
      y[21][iwl] = ((1 / length) * log(gie_LTO[31][iwl] / gie_LTO[66][iwl])); // compute delta k
      y[22][iwl] = ((1 / length) * log(gie_LTO[32][iwl] / gie_LTO[56][iwl])); // compute delta k    
      y[23][iwl] = ((1 / length) * log(gie_LTO[33][iwl] / gie_LTO[57][iwl])); // compute delta k
      y[24][iwl] = ((1 / length) * log(gie_LTO[35][iwl] / gie_LTO[61][iwl])); // compute delta k
      y[25][iwl] = ((1 / length) * log(gie_LTO[37][iwl] / gie_LTO[62][iwl])); // compute delta k
      y[26][iwl] = ((1 / length) * log(gie_LTO[40][iwl] / gie_LTO[67][iwl])); // compute delta k
      y[27][iwl] = ((1 / length) * log(gie_LTO[41][iwl] / gie_LTO[65][iwl])); // compute delta k   
      y[28][iwl] = ((1 / length) * log(gie_LTO[44][iwl] / gie_LTO[68][iwl])); // compute delta k
      y[29][iwl] = ((1 / length) * log(gie_LTO[45][iwl] / gie_LTO[69][iwl])); // compute delta k    
      y[30][iwl] = ((1 / length) * log(gie_LTO[46][iwl] / gie_LTO[70][iwl])); // compute delta k
      y[31][iwl] = ((1 / length) * log(gie_LTO[48][iwl] / gie_LTO[71][iwl])); // compute delta k
      y[32][iwl] = ((1 / length) * log(gie_LTO[50][iwl] / gie_LTO[72][iwl])); // compute delta k
      y[33][iwl] = ((1 / length) * log(gie_LTO[51][iwl] / gie_LTO[73][iwl])); // compute delta k
      y[34][iwl] = ((1 / length) * log(gie_LTO[53][iwl] / gie_LTO[74][iwl])); // compute delta k
      y[35][iwl] = ((1 / length) * log(gie_LTO[54][iwl] / gie_LTO[75][iwl])); // compute delta k    
      y[36][iwl] = ((1 / length) * log(gie_LTO[55][iwl] / gie_LTO[76][iwl])); // compute delta k
      y[37][iwl] = ((1 / length) * log(gie_LTO[58][iwl] / gie_LTO[78][iwl])); // compute delta k
      y[38][iwl] = ((1 / length) * log(gie_LTO[59][iwl] / gie_LTO[77][iwl])); // compute delta k
      y[39][iwl] = ((1 / length) * log(gie_LTO[60][iwl] / gie_LTO[79][iwl])); // compute delta k

  }

  fclose(gieFile);
  TGraph* deltas[39];
  for (int i = 0; i < 39; i++) {
      deltas[i] = new TGraph(Ngiewl, gieWL, y[i]);
      deltas[i]->GetXaxis()->SetTitle("Wavelength (nm)");
      deltas[i]->GetYaxis()->SetTitle("DeltaK");
      deltas[i]->SetTitle(Form("Crystal_%d", crystal[i]));
      deltas[i]->GetYaxis()->SetRangeUser(-0.2, 1.5);






  }
  TCanvas* c90 = new TCanvas("c90", "Crystal 49 dk", 750, 1000);

  deltas[0]->Draw("AP*");
  c90->Print("crystal49dk.pdf");
  TCanvas* c11a = new TCanvas("c11a", "Crystal 37 dk", 750, 1000);
  deltas[1]->Draw("AP*");
  c11a->Print("crystal37dk.pdf");
  TCanvas* c11b = new TCanvas("c11b", "Crystal 40 dk", 750, 1000);
  deltas[2]->Draw("AP*");
  c11b->Print("crystal40dk.pdf");
  TCanvas* c11c = new TCanvas("c11c", "Crystal 41 dk", 750, 1000);
  deltas[3]->Draw("AP*");
  c11c->Print("crystal41dk.pdf");
  TCanvas* c11d = new TCanvas("c11d", "Crystal 61 dk", 750, 1000);
  deltas[4]->Draw("AP*");
  c11d->Print("crystal61dk.pdf");
  TCanvas* c11e = new TCanvas("c11e", "Crystal 56 dk", 750, 1000);
  deltas[5]->Draw("AP*");
  c11e->Print("crystal56dk.pdf");
  TCanvas* c11f = new TCanvas("c11", "Crystal 48 dk", 750, 1000);
  deltas[6]->Draw("AP*");
  c11f->Print("crystal48dk.pdf");
  TCanvas* c11g = new TCanvas("c11g", "Crystal 53 dk", 750, 1000);
  deltas[7]->Draw("AP*");
  c11g->Print("crystal53dk.pdf");
  TCanvas* c11h = new TCanvas("c11h", "Crystal 52 dk", 750, 1000);
  deltas[8]->Draw("AP*");
  c11h->Print("crystal52dk.pdf");
  TCanvas* c50 = new TCanvas("c50", "Crystal 47 dk", 750, 1000);

  deltas[9]->Draw("AP*");
  c50->Print("crystal47dk.pdf");
  TCanvas* c51 = new TCanvas("c51", "Crystal 215 dk", 750, 1000);
  deltas[10]->Draw("AP*");
  c51->Print("crystal1215dk.pdf");
  TCanvas* c52 = new TCanvas("c52", "Crystal 208 dk", 750, 1000);
  deltas[11]->Draw("AP*");
  c52->Print("crystal208dk.pdf");
  TCanvas* c53 = new TCanvas("c53", "Crystal 231 dk", 750, 1000);
  deltas[12]->Draw("AP*");
  c53->Print("crystal231dk.pdf");
  TCanvas* c54 = new TCanvas("c54", "Crystal 240 dk", 750, 1000);
  deltas[13]->Draw("AP*");
  c54->Print("crystal240dk.pdf");
  TCanvas* c55 = new TCanvas("c55", "Crystal 206 dk", 750, 1000);
  deltas[14]->Draw("AP*");
  c55->Print("crystal237dk.pdf");
  TCanvas* c56 = new TCanvas("c56", "Crystal 237 dk", 750, 1000);
  deltas[15]->Draw("AP*");
  c56->Print("crystal237dk.pdf");
  TCanvas* c57 = new TCanvas("c57", "Crystal 233 dk", 750, 1000);
  deltas[16]->Draw("AP*");
  c57->Print("crystal233dk.pdf");
  TCanvas* c58 = new TCanvas("c58", "Crystal 214 dk", 750, 1000);
  deltas[17]->Draw("AP*");
  c58->Print("crystal214dk.pdf");
  TCanvas* c59 = new TCanvas("c59", "Crystal 236 dk", 750, 1000);

  deltas[18]->Draw("AP*");
  c59->Print("crystal236dk.pdf");
  TCanvas* c60 = new TCanvas("c60", "Crystal 232 dk", 750, 1000);
  deltas[19]->Draw("AP*");
  c60->Print("crystal232dk.pdf");
  TCanvas* c61 = new TCanvas("c61", "Crystal 244 dk", 750, 1000);
  deltas[20]->Draw("AP*");
  c61->Print("crystal244dk.pdf");
  TCanvas* c62 = new TCanvas("c62", "Crystal 242 dk", 750, 1000);
  deltas[21]->Draw("AP*");
  c62->Print("crystal242dk.pdf");
  TCanvas* c63 = new TCanvas("c63", "Crystal 259 dk", 750, 1000);
  deltas[22]->Draw("AP*");
  c63->Print("crystal259dk.pdf");
  TCanvas* c64 = new TCanvas("c64", "Crystal 246 dk", 750, 1000);
  deltas[23]->Draw("AP*");
  c64->Print("crystal246dk.pdf");
  TCanvas* c65 = new TCanvas("c65", "Crystal 265 dk", 750, 1000);
  deltas[24]->Draw("AP*");
  c65->Print("crystal265dk.pdf");
  TCanvas* c66 = new TCanvas("c66", "Crystal 247 dk", 750, 1000);
  deltas[25]->Draw("AP*");
  c66->Print("crystal247dk.pdf");
  TCanvas* c67 = new TCanvas("c67", "Crystal 252 dk", 750, 1000);
  deltas[26]->Draw("AP*");
  c67->Print("crystal252dk.pdf");
  TCanvas* c68 = new TCanvas("c68", "Crystal 251 dk", 750, 1000);

  deltas[27]->Draw("AP*");
  c68->Print("crystal251dk.pdf");
  TCanvas* c69 = new TCanvas("c69", "Crystal 254 dk", 750, 1000);
  deltas[28]->Draw("AP*");
  c69->Print("crystal254dk.pdf");
  TCanvas* c70 = new TCanvas("c70", "Crystal 269 dk", 750, 1000);
  deltas[29]->Draw("AP*");
  c70->Print("crystal269dk.pdf");
  TCanvas* c71 = new TCanvas("c71", "Crystal 304 dk", 750, 1000);
  deltas[30]->Draw("AP*");
  c71->Print("crystal304dk.pdf");
  TCanvas* c72 = new TCanvas("c72", "Crystal 322 dk", 750, 1000);
  deltas[31]->Draw("AP*");
  c72->Print("crystal322dk.pdf");
  TCanvas* c73 = new TCanvas("c73", "Crystal 282 dk", 750, 1000);
  deltas[32]->Draw("AP*");
  c73->Print("crystal282dk.pdf");
  TCanvas* c74 = new TCanvas("c74", "Crystal 276 dk", 750, 1000);
  deltas[33]->Draw("AP*");
  c74->Print("crystal276dk.pdf");
  TCanvas* c75 = new TCanvas("c75", "Crystal 274 dk", 750, 1000);
  deltas[34]->Draw("AP*");
  c75->Print("crystal274dk.pdf");
  TCanvas* c76 = new TCanvas("76", "Crystal 268 dk", 750, 1000);
  deltas[35]->Draw("AP*");
  c77->Print("crystal268dk.pdf");
  TCanvas* c77 = new TCanvas("c77", "Crystal 288 dk", 750, 1000);
  deltas[36]->Draw("AP*");
  c77->Print("crystal288dk.pdf");
  TCanvas* c78 = new TCanvas("c78", "Crystal 298 dk", 750, 1000);
  deltas[37]->Draw("AP*");
  c78->Print("crystal1298dk.pdf");
  TCanvas* c79 = new TCanvas("c79", "Crystal 330 dk", 750, 1000);
  deltas[38]->Draw("AP*");
  c79->Print("crystal330dk.pdf");*/
  //new crystals
  //
  //
/*
  current_line = 0;

  nscan = 0;
  ifstream file3("BOX5_6_PROD.csv", std::ifstream::in);
  odd = 1;
  counter = 0;
  int counter2 = 0;
  while (!file3.eof()) {

      //int counter = 0;
      current_line++;
      getline(file3, line);
      if (current_line == 1) {
          NgieFile = std::stoi(line);
          cout << NgieFile << endl;

      }
      //if(file.is_open()) {
          //cout << "File failed to open." << endl;
          //  return 1;
      //}
      if (current_line == 2) {
          iss.str(line);
          while (iss.good()) {
              iss >> word;
              nscan++;
              char_separator<char> sep("_");
              tokenizer<char_separator<char>> tokens(word, sep);
              int Nword = 0;
              for (const string& s : tokens) {
                  //cout << s << Nword << endl;
                  if (Nword == 1) {
                      crynum = std::stoi(s);
                      cout << crynum << " " << counter2 << endl;
                      counter2++;
                  }
                  if (Nword == 2) {
                      cout << s << endl;
                  }
                  Nword++;

              }

          }
      }
      else
      {
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


/*
  int crystal[20] = { 92, 126, 105, 107, 117, 115, 90, 96, 111, 364, 360, 102, 333, 350, 353, 351, 369, 355, 334, 357 };
  //TGraph *deltas;
  for (int iwl = 0; iwl < Ngiewl; iwl++) {//for every wavelength


      y[0][iwl] = ((1 / length) * log(gie_LTO[0][iwl] / gie_LTO[17][iwl])); // compute delta k
      y[1][iwl] = ((1 / length) * log(gie_LTO[1][iwl] / gie_LTO[18][iwl])); // compute delta k
      y[2][iwl] = ((1 / length) * log(gie_LTO[2][iwl] / gie_LTO[23][iwl])); // compute delta k   
      y[3][iwl] = ((1 / length) * log(gie_LTO[3][iwl] / gie_LTO[38][iwl])); // compute delta k
      y[4][iwl] = ((1 / length) * log(gie_LTO[4][iwl] / gie_LTO[39][iwl])); // compute delta k    
      y[5][iwl] = ((1 / length) * log(gie_LTO[5][iwl] / gie_LTO[40][iwl])); // compute delta k
      y[6][iwl] = ((1 / length) * log(gie_LTO[6][iwl] / gie_LTO[24][iwl])); // compute delta k
      y[7][iwl] = ((1 / length) * log(gie_LTO[7][iwl] / gie_LTO[22][iwl])); // compute delta k
      y[8][iwl] = ((1 / length) * log(gie_LTO[8][iwl] / gie_LTO[35][iwl])); // compute delta k
      y[9][iwl] = ((1 / length) * log(gie_LTO[9][iwl] / gie_LTO[36][iwl])); // compute delta k
      y[10][iwl] = ((1 / length) * log(gie_LTO[10][iwl] / gie_LTO[25][iwl])); // compute delta k
      y[11][iwl] = ((1 / length) * log(gie_LTO[11][iwl] / gie_LTO[26][iwl])); // compute delta k   
      y[12][iwl] = ((1 / length) * log(gie_LTO[12][iwl] / gie_LTO[27][iwl])); // compute delta k
      y[13][iwl] = ((1 / length) * log(gie_LTO[13][iwl] / gie_LTO[37][iwl])); // compute delta k    
      y[14][iwl] = ((1 / length) * log(gie_LTO[14][iwl] / gie_LTO[28][iwl])); // compute delta k
      y[15][iwl] = ((1 / length) * log(gie_LTO[15][iwl] / gie_LTO[30][iwl])); // compute delta k
      y[16][iwl] = ((1 / length) * log(gie_LTO[16][iwl] / gie_LTO[31][iwl])); // compute delta k
      y[17][iwl] = ((1 / length) * log(gie_LTO[19][iwl] / gie_LTO[32][iwl])); // compute delta k
      y[18][iwl] = ((1 / length) * log(gie_LTO[20][iwl] / gie_LTO[33][iwl])); // compute delta k
      y[19][iwl] = ((1 / length) * log(gie_LTO[21][iwl] / gie_LTO[34][iwl])); // compute delta k
      

  }

  fclose(gieFile);
  TGraph* deltas[20];
  for (int i = 0; i < 20; i++) {
      deltas[i] = new TGraph(Ngiewl, gieWL, y[i]);
      deltas[i]->GetXaxis()->SetTitle("Wavelength (nm)");
      deltas[i]->GetYaxis()->SetTitle("DeltaK");
      deltas[i]->SetTitle(Form("Crystal_%d", crystal[i]));
      deltas[i]->GetYaxis()->SetRangeUser(-0.2, 1.5);






  }
  TCanvas* d92 = new TCanvas("d90", "Crystal 92 dk", 750, 1000);

  deltas[0]->Draw("AP*");
  d90->Print("crystal92dk.pdf");
  TCanvas* d11a = new TCanvas("d11a", "Crystal 126 dk", 750, 1000);
  deltas[1]->Draw("AP*");
  d11a->Print("crystal126dk.pdf");
  TCanvas* d11b = new TCanvas("d11b", "Crystal 105 dk", 750, 1000);
  deltas[2]->Draw("AP*");
  d11b->Print("crystal40dk.pdf");
  TCanvas* d11c = new TCanvas("d11c", "Crystal 107 dk", 750, 1000);
  deltas[3]->Draw("AP*");
  d11c->Print("crystal107dk.pdf");
  TCanvas* d11d = new TCanvas("d11d", "Crystal 117 dk", 750, 1000);
  deltas[4]->Draw("AP*");
  d11d->Print("crystal117dk.pdf");
  TCanvas* d11e = new TCanvas("d11e", "Crystal 115 dk", 750, 1000);
  deltas[5]->Draw("AP*");
  d11e->Print("crystal115dk.pdf");
  TCanvas* d11f = new TCanvas("d11f", "Crystal 90 dk", 750, 1000);
  deltas[6]->Draw("AP*");
  d11f->Print("crystal90dk.pdf");
  TCanvas* d11g = new TCanvas("d11g", "Crystal 96 dk", 750, 1000);
  deltas[7]->Draw("AP*");
  d11g->Print("crystal96dk.pdf");
  TCanvas* d11h = new TCanvas("d11h", "Crystal 111 dk", 750, 1000);
  deltas[8]->Draw("AP*");
  d11h->Print("crystal111dk.pdf");
  TCanvas* d50 = new TCanvas("d50", "Crystal 364 dk", 750, 1000);

  deltas[9]->Draw("AP*");
  d50->Print("crystal364dk.pdf");
  TCanvas* d51 = new TCanvas("d51", "Crystal 360 dk", 750, 1000);
  deltas[10]->Draw("AP*");
  d51->Print("crystal360dk.pdf");
  TCanvas* d52 = new TCanvas("d52", "Crystal 102 dk", 750, 1000);
  deltas[11]->Draw("AP*");
  d52->Print("crystal102dk.pdf");
  TCanvas* d53 = new TCanvas("d53", "Crystal 333 dk", 750, 1000);
  deltas[12]->Draw("AP*");
  d53->Print("crystal333dk.pdf");
  TCanvas* d54 = new TCanvas("d54", "Crystal 350 dk", 750, 1000);
  deltas[13]->Draw("AP*");
  d54->Print("crystal350dk.pdf");
  TCanvas* d55 = new TCanvas("d55", "Crystal 353 dk", 750, 1000);
  deltas[14]->Draw("AP*");
  d55->Print("crystal353dk.pdf");
  TCanvas* d56 = new TCanvas("d56", "Crystal 351 dk", 750, 1000);
  deltas[15]->Draw("AP*");
  d56->Print("crystal351dk.pdf");
  TCanvas* d57 = new TCanvas("d57", "Crystal 369 dk", 750, 1000);
  deltas[16]->Draw("AP*");
  d57->Print("crystal369dk.pdf");
  TCanvas* d58 = new TCanvas("d58", "Crystal 355 dk", 750, 1000);
  deltas[17]->Draw("AP*");
  d58->Print("crystal355dk.pdf");
  TCanvas* d59 = new TCanvas("d59", "Crystal 334 dk", 750, 1000);

  deltas[18]->Draw("AP*");
  d59->Print("crystal334dk.pdf");
  TCanvas* d60 = new TCanvas("d60", "Crystal 357 dk", 750, 1000);
  deltas[19]->Draw("AP*");
  d60->Print("crystal357dk.pdf");
  
  
  */
//new crystals
//
//
//
/*
current_line = 0;

nscan = 0;
ifstream file4("BOX8_PROD.csv", std::ifstream::in);
odd = 1;
counter = 0;
int counter2 = 0;
while (!file4.eof()) {

    //int counter = 0;
    current_line++;
    getline(file4, line);
    if (current_line == 1) {
        NgieFile = std::stoi(line);
        cout << NgieFile << endl;

    }
    //if(file.is_open()) {
        //cout << "File failed to open." << endl;
        //  return 1;
    //}
    if (current_line == 2) {
        iss.str(line);
        while (iss.good()) {
            iss >> word;
            nscan++;
            char_separator<char> sep("_");
            tokenizer<char_separator<char>> tokens(word, sep);
            int Nword = 0;
            for (const string& s : tokens) {
                //cout << s << Nword << endl;
                if (Nword == 1) {
                    crynum = std::stoi(s);
                    cout << crynum << " " << counter2 << endl;
                    counter2++;
                }
                if (Nword == 2) {
                    cout << s << endl;
                }
                Nword++;

            }

        }
    }
    else
    {
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



int crystal[15] = { 139,22, 16, 12, 201, 229, 178, 200, 230, 11,31,24,34,1,21 };
//TGraph *deltas;
for (int iwl = 0; iwl < Ngiewl; iwl++) {//for every wavelength


    y[0][iwl] = ((1 / length) * log(gie_LTO[0][iwl] / gie_LTO[1][iwl])); // compute delta k
    y[1][iwl] = ((1 / length) * log(gie_LTO[3][iwl] / gie_LTO[18][iwl])); // compute delta k
    y[2][iwl] = ((1 / length) * log(gie_LTO[4][iwl] / gie_LTO[19][iwl])); // compute delta k   
    y[3][iwl] = ((1 / length) * log(gie_LTO[5][iwl] / gie_LTO[29][iwl])); // compute delta k
    y[4][iwl] = ((1 / length) * log(gie_LTO[6][iwl] / gie_LTO[30][iwl])); // compute delta k    
    y[5][iwl] = ((1 / length) * log(gie_LTO[7][iwl] / gie_LTO[31][iwl])); // compute delta k
    y[6][iwl] = ((1 / length) * log(gie_LTO[8][iwl] / gie_LTO[20][iwl])); // compute delta k
    y[7][iwl] = ((1 / length) * log(gie_LTO[9][iwl] / gie_LTO[21][iwl])); // compute delta k
    y[8][iwl] = ((1 / length) * log(gie_LTO[10][iwl] / gie_LTO[22][iwl])); // compute delta k
    y[9][iwl] = ((1 / length) * log(gie_LTO[11][iwl] / gie_LTO[26][iwl])); // compute delta k
    y[10][iwl] = ((1 / length) * log(gie_LTO[12][iwl] / gie_LTO[27][iwl])); // compute delta k
    y[11][iwl] = ((1 / length) * log(gie_LTO[13][iwl] / gie_LTO[23][iwl])); // compute delta k   
    y[12][iwl] = ((1 / length) * log(gie_LTO[14][iwl] / gie_LTO[24][iwl])); // compute delta k
    y[13][iwl] = ((1 / length) * log(gie_LTO[15][iwl] / gie_LTO[25][iwl])); // compute delta k    
    y[14][iwl] = ((1 / length) * log(gie_LTO[16][iwl] / gie_LTO[28][iwl])); // compute delta k


}

fclose(gieFile);
TGraph* deltas[15];
for (int i = 0; i < 15; i++) {
    deltas[i] = new TGraph(Ngiewl, gieWL, y[i]);
    deltas[i]->GetXaxis()->SetTitle("Wavelength (nm)");
    deltas[i]->GetYaxis()->SetTitle("DeltaK");



    deltas[i]->SetTitle(Form("Crystal_%d", crystal[i]));






}
TCanvas* c90 = new TCanvas("c90", "Crystal 139 dk", 750, 1000);
deltas[0]->Draw("AP*");
c90->Print("crystal139dk.pdf");
TCanvas* c11b = new TCanvas("c11b", "Crystal 22 dk", 750, 1000);
deltas[1]->Draw("AP*");
c11b->Print("crystal22dk.pdf");
TCanvas* c11c = new TCanvas("c11c", "Crystal 16 dk", 750, 1000);
deltas[2]->Draw("AP*");
c11c->Print("crystal16dk.pdf");
TCanvas* c11d = new TCanvas("c11d", "Crystal 12 dk", 750, 1000);
deltas[3]->Draw("AP*");
c11d->Print("crystal12dk.pdf");
TCanvas* c11e = new TCanvas("c11e", "Crystal 201 dk", 750, 1000);
deltas[4]->Draw("AP*");
c11e->Print("crystal201dk.pdf");
TCanvas* c11f = new TCanvas("c11", "Crystal 229 dk", 750, 1000);
deltas[5]->Draw("AP*");
c11f->Print("crystal229dk.pdf");
TCanvas* c11g = new TCanvas("c11g", "Crystal 178 dk", 750, 1000);
deltas[6]->Draw("AP*");
c11g->Print("crystal178dk.pdf");
TCanvas* c11h = new TCanvas("c11h", "Crystal 200 dk", 750, 1000);
deltas[7]->Draw("AP*");
c11h->Print("crystal139dk.pdf");
TCanvas* c11i = new TCanvas("c11i", "Crystal 230 dk", 750, 1000);
deltas[8]->Draw("AP*");
c11i->Print("crystal230dk.pdf");
TCanvas* c11j = new TCanvas("c11j", "Crystal 11 dk", 750, 1000);
deltas[9]->Draw("AP*");
c11j->Print("crystal11dk.pdf");
TCanvas* c11k = new TCanvas("c11k", "Crystal 31 dk", 750, 1000);
deltas[10]->Draw("AP*");
c11k->Print("crystal31dk.pdf");
TCanvas* c11l = new TCanvas("c11l", "Crystal 24 dk", 750, 1000);
deltas[11]->Draw("AP*");
c11l->Print("crystal24dk.pdf");
TCanvas* c11m = new TCanvas("c11m", "Crystal 34 dk", 750, 1000);
deltas[12]->Draw("AP*");
c11m->Print("crystal34dk.pdf");
TCanvas* c11n = new TCanvas("c11n", "Crystal 1 dk", 750, 1000);
deltas[13]->Draw("AP*");
c11n->Print("crystal1dk.pdf");
TCanvas* c11o = new TCanvas("c11o", "Crystal 36 dk", 750, 1000);
deltas[14]->Draw("AP*");
c11o->Print("crystal36dk.pdf");
*/
  //gr_LTO_bef[11]->Draw("PSAME*");
  //gr_LTO_bef[12]->Draw("PSAME*");
  //gr_LTO_bef[28]->Draw("PSAME*");

  
  //TGraph *scaled = new TGraph(iwl,gieWL,giescale);
  //TGraph *scaled = new TGraph(Ngiewl,gieWL,giescale);
  //scaled->SetMarkerColor(2);
  //scaled->Draw("PSAME*");
  

/*
  double gieSum[47];
  double gieAve[47];
  double d_LT_IR[4];
  double d_LT_420[4];
  double ratio_LT_IR[4];
  double ratio_LT_420[4];

  double d_LT_IR_023[4];
  double d_LT_420_023[4];
  double ratio_LT_IR_023[4];
  double ratio_LT_420_023[4];

  double a[2] = {0.5,0.9};
  double b[2] = {0.92,1.0};


  double c[8] = {0.995412, 0.972973, 0.984441, 0.964205, 0.98124, 0.936745, 0.977209, 0.953963};
  double d[8] = {0.865107, 0.765244, 0.850682, 0.748689, 0.82313, 0.526084, 0.783732, 0.710778};


  for(int i=0; i<47; i++) {//for every scan
    for(int iwl=80; iwl<200; iwl++) {//for every wavelength

    gieSum[i] += gie_LTO[i][iwl];

    //gie_LTO_corr[20][i] = 0.96 * gie_LTO[((int)gieSTAT[20])][i];

    }
    gieAve[i] = gieSum[i] / 120;

    //cout << gieSum[39] / 120 << endl;

  }

  for(int ii=0; ii<4; ii++){
    ratio_LT_IR[ii] =  gieAve[19+ii]/gieAve[18];
    d_LT_IR[ii] =  gieAve[18] - gieAve[19+ii];
    ratio_LT_420[ii] =  gie_LTO[19+ii][960]/gie_LTO[18][960];
    d_LT_420[ii] =  gie_LTO[18][960] - gie_LTO[19+ii][960];

    //if(ii<3){
      ratio_LT_IR_023[ii] =  gieAve[25+ii]/gieAve[24];
      d_LT_IR_023[ii] =  gieAve[24] - gieAve[25+ii];
      ratio_LT_420_023[ii] =  gie_LTO[25+ii][960]/gie_LTO[24][960];
      d_LT_420_023[ii] =  gie_LTO[24][960] - gie_LTO[25+ii][960];
      // }
  }


 

  cout << "light transmission ratio at 900ish for crystal 039 irr " << gieAve[45]/gieAve[44] << endl;

  cout << "expected light transmission ratio at 900ish for crystal 039 irr " << 0.84 + (0.1666 * (gie_LTO[45][960]/gie_LTO[44][960])) << endl;




  cout << "E_1 = " << ratio_LT_IR[0] << "  Diff_IR = " << d_LT_IR[0] << "  D_1 = " << ratio_LT_420[0] << "  Diff_420 = " << d_LT_420[0] << endl;
  cout << "E_1 = " << ratio_LT_IR_023[0] << "  Diff_IR = " << d_LT_IR_023[0] << "  D_1 = " << ratio_LT_420_023[0] << "  Diff_420 = " << d_LT_420_023[0] << endl;

  cout << "E_2 = " << ratio_LT_IR[1] << "  Diff_IR = " << d_LT_IR[1] << "  D_2 = " << ratio_LT_420[1] << "  Diff_420 = " << d_LT_420[1] << endl;
  cout << "E_2 = " << ratio_LT_IR_023[1] << "  Diff_IR = " << d_LT_IR_023[1] << "  D_2 = " << ratio_LT_420_023[1] << "  Diff_420 = " << d_LT_420_023[1] << endl;

  cout << "E_3 = " << ratio_LT_IR[2] << "  Diff_IR = " << d_LT_IR[2] << "  D_3 = " << ratio_LT_420[2] << "  Diff_420 = " << d_LT_420[2] << endl;
  cout << "E_3 = " << ratio_LT_IR_023[2] << "  Diff_IR = " << d_LT_IR_023[2] << "  D_3 = " << ratio_LT_420_023[2] << "  Diff_420 = " << d_LT_420_023[2] << endl;

  cout << "E_4 = " << ratio_LT_IR[3] << "  Diff_IR = " << d_LT_IR[3] << "  D_4 = " << ratio_LT_420[3] << "  Diff_420 = " << d_LT_420[3] << endl;
  cout << "E_4 = " << ratio_LT_IR_023[3] << "  Diff_IR = " << d_LT_IR_023[3] << "  D_4 = " << ratio_LT_420_023[3] << "  Diff_420 = " << d_LT_420_023[3] << endl;
  


  //TCanvas *c13=new TCanvas("c13","c13");  //create canvas
  //TGraph *diffCurve = new TGraph(4,d_LT_420,d_LT_IR);
  //diffCurve->SetTitle("Light Transmission Difference");
  //diffCurve->GetXaxis()->SetTitle("Difference from unirradiated at 420nm");
  //diffCurve->GetYaxis()->SetTitle("Difference from unirradiated at IR region (ave)");

  //diffCurve->Draw("AP*");


  //TCanvas *c13a=new TCanvas("c13a","c13a");  //create canvas
 // TGraph *diffCurve2 = new TGraph(4,d_LT_420_023,d_LT_IR_023);
  //diffCurve2->SetTitle("Light Transmission Difference");
  //diffCurve2->GetXaxis()->SetTitle("Difference from unirradiated at 420nm");
  //diffCurve2->GetYaxis()->SetTitle("Difference from unirradiated at IR region (ave)");
  //diffCurve2->SetMarkerColor(2); //red
  ///diffCurve2->Draw("PSAME*");

  //TCanvas *c14=new TCanvas("c14","c14");  //create canvas
  //TGraph *ratioPoints = new TGraph(2,a,b);
  //ratioPoints->GetXaxis()->SetTitle("Ratio at 420nm");
  //ratioPoints->GetYaxis()->SetTitle("Ratio at IR region (ave)");
  //ratioPoints->Draw("AP");

  //TGraph *ratioCurve = new TGraph(4,ratio_LT_420,ratio_LT_IR);
  //ratioCurve->SetTitle("Light Transmission Ratio");
  //ratioCurve->GetXaxis()->SetTitle("Ratio at 420nm");
  //ratioCurve->GetYaxis()->SetTitle("Ratio at IR region (ave)");

  //ratioCurve->Draw("PSAME*");

  //TCanvas *c14a=new TCanvas("c14a","c14a");  //create canvas
  //TGraph *ratioCurve2 = new TGraph(4,ratio_LT_420_023,ratio_LT_IR_023);
  //ratioCurve2->SetTitle("Light Transmission Ratio");
  //ratioCurve2->GetXaxis()->SetTitle("Ratio at 420nm");
  //ratioCurve2->GetYaxis()->SetTitle("Ratio at IR region (ave)");
  //ratioCurve2->SetMarkerColor(2); //red
  //ratioCurve2->Draw("PSAME*");

  //TCanvas *c15=new TCanvas("c15","c15");  //create canvas
  //TGraph *ratioP = new TGraph(8,d,c);
  //ratioP->GetXaxis()->SetTitle("Ratio at 420nm");
 // ratioP->GetYaxis()->SetTitle("Ratio at IR region (ave)");
 // ratioP->Draw("AP*");


//   /*
//   // saving histograms to file
//   TIter next(gDirectory->GetList());
//   TObject* obj;
//   TH1 *my_hi[1000];
//   int ihi=0;
//   while(obj= (TObject*)next()){
//   if(obj->InheritsFrom(TH1::Class())){
//   my_hi[ihi] = (TH1*)obj;
//   ihi++;
//   }
//   }
//   TFile myfile("crystals.root", "recreate");
//   for(int i=0; i<ihi; i++) my_hi[i]->Write();
//   myfile.Close();
//   */

  }
 
