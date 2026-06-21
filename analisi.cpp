#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>



using namespace std;
using namespace RooFit;


double findMin(TH1D* h) { //funzione che trova il bin minimo dell'istogramma in un intervallo di energie
  
  double Emin=3e4;
  double Emax=8e4; //intervallo in cui cercare, cambiare qui!!
  
  int nb = h->GetNbinsX();
  int binMin = -1;
  double minVal = 1e6;

  for (int b = 1; b <= nb; b++) {
      double center = h->GetBinCenter(b);
      if (center < Emin || center > Emax) continue;

      double c = h->GetBinContent(b);
      if (c < minVal) {
          minVal = c;
          binMin = b;
       }
   }

   return h->GetBinCenter(binMin);
}

double NoiseCut95(TF1* g) // si chiama così ormai, ma in realtà il target è impostato a  99%
{
    // 1. Integrale totale della Gaussiana nel range dell'istogramma (a 400 keV la gaussiana è approx nulla)
    double Itot = g->Integral(0, 400000);

    // 2. Target: 99% dell'integrale totale
    double target = 0.99 * Itot;

    // 3. Definiamo una funzione cumulativa (associa alla x l'integrale della funzione passata tra 0 e x)
    TF1* F = new TF1("F",
                     [g](double* x, double*) {
                         return g->Integral(0, x[0]);
                     },
                     0, 400000, 0);

    // 4. Trova E_cut tale che F(E_cut) = target
    double Ecut = F->GetX(target, 0, 400000);

    delete F;
    return Ecut;
}




void analisi(TString file= "filtered.root",  TString output="cut_output.root"){
  TFile *f = TFile::Open(file);

   // Recupera i due TTree
    TTree *tx = (TTree*) f->Get("filteredX");
    TTree *ty = (TTree*) f->Get("filteredY");

     
    vector<double>* xx=nullptr;
    vector<double>* zx=nullptr;
    vector<double>* Ex=nullptr;
   

    vector<double>* yy=nullptr;
    vector<double>* zy=nullptr;
    vector<double>* Ey=nullptr;
   

    tx->SetBranchAddress("xx", &xx);
    tx->SetBranchAddress("zx", &zx);
    tx->SetBranchAddress("Ex", &Ex);

    ty->SetBranchAddress("yy", &yy);
    ty->SetBranchAddress("zy", &zy);
    ty->SetBranchAddress("Ey", &Ey);

    
    int nX=tx->GetEntries();
    int nY=ty->GetEntries();

    gStyle->SetOptStat(0);

    
    
    //=========== Istogrammi===========
    TH1D *hEx = new TH1D("hEx", "Energia depositata, X;E_{dep} (eV);Conteggi", 200, 0, 400000);
    TH1D *hEy = new TH1D("hEy", "Energia depositata, Y;E_{dep} (eV);Conteggi", 200, 0, 400000);

    // Riempimento istogramma X
    for (int ev = 0; ev < nX; ev++) {
      tx->GetEntry(ev);
      for (double e : *Ex)
        hEx->Fill(e);
    }

    // Riempimento istogramma Y
    for (int ev = 0; ev < nY; ev++) {
        ty->GetEntry(ev);
        for (double e : *Ey)
	  hEy->Fill(e);
    }



    // Istogrammi per ogni piano:
    const int nplanes=5;
    int plane[nplanes]={0,250,500,750,1000};

    TH1D* hEx_plane[nplanes];
    TH1D* hEy_plane[nplanes];

    for(int p=0; p<5; p++){
      TString namex = Form("hEx_z%d", plane[p]);
      TString titlex = Form("Energia X, piano Z=%d;E_{dep} (eV);Conteggi", plane[p]);
      hEx_plane[p] = new TH1D(namex, titlex, 600, 0, 400000);

       TString namey = Form("hEy_z%d", plane[p]);
      TString titley = Form("Energia Y, piano Z=%d;E_{dep} (eV);Conteggi", plane[p]);
      hEy_plane[p] = new TH1D(namey, titley, 600, 0, 400000);
    }

    for (int ev = 0; ev < nX; ev++) {
    tx->GetEntry(ev);
    ty->GetEntry(ev);
    for (int i = 0; i < Ex->size(); i++) {  //ciclo per x
        double e = Ex->at(i);
        int z = (int) zx->at(i);

        for (int p = 0; p < nplanes; p++) {
            if (z == plane[p]) {
                hEx_plane[p]->Fill(e);
                break;
            }
        }
    }
    
    for (int i = 0; i < Ey->size(); i++) {  //ciclo per y
        double e = Ey->at(i);
        int z = (int) zy->at(i);

        for (int p = 0; p < nplanes; p++) {
            if (z == plane[p]) {
                hEy_plane[p]->Fill(e);
                break;
            }
        }
    }
}

    cout << "\n=== Entries per piano X ===" << endl;

    for (int p = 0; p < nplanes; p++) {
      cout << "Piano Z[" << plane[p] << "]: "
         << hEx_plane[p]->GetEntries()
         << " entries" << endl;

    }

    cout << "\n=== Entries per piano Y ===" << endl;

    for (int p = 0; p < nplanes; p++) {
      cout << "Piano Z[" << plane[p] << "]: "
         << hEy_plane[p]->GetEntries()
         << " entries" << endl;
    }


    //=============Fit noise=========================
    // facciamo il fit della gaussiana nel range (0,25000), per il quale "si vede a occhio" che l'andamento è gaussiano

    //==========Fit noise piano O X======================

    TCanvas* c_0x = new TCanvas("c_noise0x", "Fit rumore", 900, 600);
    c_0x->cd();
    
    TF1* g0x = new TF1("g0x", "gaus", 0, 25000);
    g0x->SetParameters(300, 0, 30000);



    // Fit binnato
    hEx_plane[0]->Fit(g0x, "R");   // "R" = rispetta il range

    double Emin0x = hEx_plane[0]->GetXaxis()->GetXmin();
    double Emax0x = hEx_plane[0]->GetXaxis()->GetXmax();

    TF1* gfull0x = new TF1("gfull0x", "gaus", Emin0x, Emax0x);
    gfull0x->SetParameters(g0x->GetParameter(0),
                         g0x->GetParameter(1),
                         g0x->GetParameter(2));
    gfull0x->SetLineColor(kRed);
    gfull0x->SetLineWidth(2);

    // Plot
    hEx_plane[0]->Draw();
    gfull0x->SetLineColor(kRed);
    gfull0x->Draw("same");


    //==========Fit noise piano 250 X======================
    TCanvas* c_250x = new TCanvas("c_noise250x", "Fit rumore", 900, 600);
    c_250x->cd();
    
    TF1* g250x = new TF1("g250x", "gaus", 0, 25000);
    g250x->SetParameters(300, 0, 30000);

    //g->SetParLimits(1, -2000, 2000);
    // g->SetParLimits(2, sigmaMin, sigmaMax);

    // Fit binnato
    hEx_plane[1]->Fit(g250x, "R");   // "R" = rispetta il range

    double Emin250x = hEx_plane[1]->GetXaxis()->GetXmin();
    double Emax250x = hEx_plane[1]->GetXaxis()->GetXmax();

    TF1* gfull250x = new TF1("gfull250x", "gaus", Emin250x, Emax250x);
    gfull250x->SetParameters(g250x->GetParameter(0),
                         g250x->GetParameter(1),
                         g250x->GetParameter(2));
    gfull250x->SetLineColor(kRed);
    gfull250x->SetLineWidth(2);

    // Plot
    hEx_plane[1]->Draw();
    gfull250x->SetLineColor(kRed);
    gfull250x->Draw("same");



     //==========Fit noise piano 1000 X======================
    TCanvas* c_1000x = new TCanvas("c_noise1000x", "Fit rumore", 900, 600);
    c_1000x->cd();
    
    TF1* g1000x = new TF1("g1000x", "gaus", 0, 25000);
    g1000x->SetParameters(700, 0, 10000);

    //g->SetParLimits(1, -2000, 2000);
    // g->SetParLimits(2, sigmaMin, sigmaMax);

    // Fit binnato
    hEx_plane[4]->Fit(g1000x, "R");   // "R" = rispetta il range

    double Emin1000x = hEx_plane[4]->GetXaxis()->GetXmin();
    double Emax1000x = hEx_plane[4]->GetXaxis()->GetXmax();

    TF1* gfull1000x = new TF1("gfull1000x", "gaus", Emin1000x, Emax1000x);
    gfull1000x->SetParameters(g1000x->GetParameter(0),
                         g1000x->GetParameter(1),
                         g1000x->GetParameter(2));
    gfull1000x->SetLineColor(kRed);
    gfull1000x->SetLineWidth(2);

    // Plot
    hEx_plane[4]->Draw();
    gfull1000x->SetLineColor(kRed);
    gfull1000x->Draw("same");




     //==========Fit noise piano O Y======================

    TCanvas* c_0y = new TCanvas("c_noise0y", "Fit rumore", 900, 600);
    c_0y->cd();
    
    TF1* g0y = new TF1("g0y", "gaus", 0, 25000);
    g0y->SetParameters(300, 0, 30000);

    //g->SetParLimits(1, -2000, 2000);
    // g->SetParLimits(2, sigmaMin, sigmaMax);

    // Fit binnato
    hEy_plane[0]->Fit(g0y, "R");   // "R" = rispetta il range

    double Emin0y = hEy_plane[0]->GetXaxis()->GetXmin();
    double Emax0y = hEy_plane[0]->GetXaxis()->GetXmax();

    TF1* gfull0y = new TF1("gfull0y", "gaus", Emin0y, Emax0y);
    gfull0y->SetParameters(g0y->GetParameter(0),
                         g0y->GetParameter(1),
                         g0y->GetParameter(2));
    gfull0y->SetLineColor(kRed);
    gfull0y->SetLineWidth(2);

    // Plot
    hEy_plane[0]->Draw();
    gfull0y->SetLineColor(kRed);
    gfull0y->Draw("same");


    //==========Fit noise piano 250 Y======================
    TCanvas* c_250y = new TCanvas("c_noise250y", "Fit rumore", 900, 600);
    c_250y->cd();
    
    TF1* g250y = new TF1("g250y", "gaus", 0, 25000);
    g250y->SetParameters(300, 0, 30000);

    //g->SetParLimits(1, -2000, 2000);
    // g->SetParLimits(2, sigmaMin, sigmaMax);

    // Fit binnato
    hEy_plane[1]->Fit(g250y, "R");   // "R" = rispetta il range

    double Emin250y = hEy_plane[1]->GetXaxis()->GetXmin();
    double Emax250y = hEy_plane[1]->GetXaxis()->GetXmax();

    TF1* gfull250y = new TF1("gfull250y", "gaus", Emin250y, Emax250y);
    gfull250y->SetParameters(g250y->GetParameter(0),
                         g250y->GetParameter(1),
                         g250y->GetParameter(2));
    gfull250y->SetLineColor(kRed);
    gfull250y->SetLineWidth(2);

    // Plot
    hEy_plane[1]->Draw();
    gfull250y->SetLineColor(kRed);
    gfull250y->Draw("same");



     //==========Fit noise piano 1000 X======================
    TCanvas* c_1000y = new TCanvas("c_noise1000y", "Fit rumore", 900, 600);
    c_1000y->cd();
    
    TF1* g1000y = new TF1("g1000y", "gaus", 0, 25000);
    g1000y->SetParameters(700, 0, 10000);

    //g->SetParLimits(1, -2000, 2000);
    // g->SetParLimits(2, sigmaMin, sigmaMax);

    // Fit binnato
    hEy_plane[4]->Fit(g1000y, "R");   // "R" = rispetta il range

    double Emin1000y = hEy_plane[4]->GetXaxis()->GetXmin();
    double Emax1000y = hEy_plane[4]->GetXaxis()->GetXmax();

    TF1* gfull1000y = new TF1("gfull1000y", "gaus", Emin1000y, Emax1000y);
    gfull1000y->SetParameters(g1000y->GetParameter(0),
                         g1000y->GetParameter(1),
                         g1000y->GetParameter(2));
    gfull1000y->SetLineColor(kRed);
    gfull1000y->SetLineWidth(2);

    // Plot
    hEy_plane[4]->Draw();
    gfull1000y->SetLineColor(kRed);
    gfull1000y->Draw("same");

    



    //Calcolo soglie per cut: eventualmente, con nuovi dati, CAMBIARE QUI.
    
    double cut0x=NoiseCut95(gfull0x);
    double cut250x=NoiseCut95(gfull250x);
    double cut1000x=NoiseCut95(gfull1000x);

    double cut0y=NoiseCut95(gfull0y);
    double cut250y=NoiseCut95(gfull250y);
    double cut1000y=NoiseCut95(gfull1000y);

    //Print energie di soglia
    cout<<" Piano 0: E_cutX = "<<cut0x<<" E_cutY = "<<cut0y<<endl;
    cout<<" Piano 250: E_cutX = "<<cut250x<<" E_cutY = "<<cut250y<<endl;
    cout<<" Piano 1000: E_cutX = "<<cut1000x<<" E_cutY = "<<cut1000y<<endl<<endl;




    
    

    TCanvas *c1 = new TCanvas("c1", "Edep X", 1200, 600);
    c1->Divide(2,1);
    c1->cd(1);
    gPad->SetLogy();
    hEx->SetLineColor(kBlue);
    hEx->SetFillColorAlpha(kBlue-9, 0.35);
    hEx->Draw();

    c1->cd(2);
    gPad->SetLogy();
    hEy->SetLineColor(kRed);
    hEy->SetFillColorAlpha(kRed-9, 0.35);
    hEy->Draw();

    
    TCanvas* c3 = new TCanvas("c3", "Edep per piano X", 1200, 800);
    c3->Divide(3,2);


    c3->cd(1);
    gPad->SetLogy();
    hEx_plane[0]->SetLineColor(kBlue);
    hEx_plane[0]->Draw();
    gfull0x->SetLineColor(kRed);
    gfull0x->Draw("same");

    c3->cd(2);
    gPad->SetLogy();
    hEx_plane[1]->SetLineColor(kBlue);
    hEx_plane[1]->Draw();
    gfull250x->SetLineColor(kRed);
    gfull250x->Draw("same");

    c3->cd(3);
    gPad->SetLogy();
    hEx_plane[2]->SetLineColor(kBlue);
    hEx_plane[2]->Draw();

    c3->cd(4);
    gPad->SetLogy();
    hEx_plane[3]->SetLineColor(kBlue);
    hEx_plane[3]->Draw();


    c3->cd(5);
    gPad->SetLogy();
    hEx_plane[4]->SetLineColor(kBlue);
    hEx_plane[4]->Draw();
    gfull1000x->SetLineColor(kRed);
    gfull1000x->Draw("same");

    


    

    TCanvas* c4 = new TCanvas("c4", "Edep per piano Y", 1200, 800);
    c4->Divide(3,2);


    c4->cd(1);
    gPad->SetLogy();
    hEy_plane[0]->SetLineColor(kRed);
    hEy_plane[0]->Draw();
    gfull0y->SetLineColor(kBlue);
    gfull0y->Draw("same");

    c4->cd(2);
    gPad->SetLogy();
    hEy_plane[1]->SetLineColor(kRed);
    hEy_plane[1]->Draw();
    gfull250y->SetLineColor(kBlue);
    gfull250y->Draw("same");

    c4->cd(3);
    gPad->SetLogy();
    hEy_plane[2]->SetLineColor(kRed);
    hEy_plane[2]->Draw();

    c4->cd(4);
    gPad->SetLogy();
    hEy_plane[3]->SetLineColor(kRed);
    hEy_plane[3]->Draw();

    c4->cd(5);
    gPad->SetLogy();
    hEy_plane[4]->SetLineColor(kRed);
    hEy_plane[4]->Draw();
    gfull1000y->SetLineColor(kBlue);
    gfull1000y->Draw("same");

    cout << "\n==================== TABELLA FIT X ====================\n";
    cout << left
     << setw(8)  << "Piano"
     << setw(16) << "Chi2/NDF"
     << setw(12) << "p-value"
     << setw(12) << "E_cut (eV)"
     << "\n--------------------------------------------------------\n";

    cout << setw(8)  << "0"
     << setw(16) << Form("%.2f/%d", g0x->GetChisquare(), g0x->GetNDF())
     << setw(12) << Form("%.3g", g0x->GetProb())
     << setw(12) << Form("%.1f", cut0x)
     << endl;

    cout << setw(8)  << "250"
     << setw(16) << Form("%.2f/%d", g250x->GetChisquare(), g250x->GetNDF())
     << setw(12) << Form("%.3g", g250x->GetProb())
     << setw(12) << Form("%.1f", cut250x)
     << endl;

    cout << setw(8)  << "1000"
     << setw(16) << Form("%.2f/%d", g1000x->GetChisquare(), g1000x->GetNDF())
     << setw(12) << Form("%.3g", g1000x->GetProb())
     << setw(12) << Form("%.1f", cut1000x)
     << endl;


    cout << "\n==================== TABELLA FIT Y ====================\n";
    cout << left
     << setw(8)  << "Piano"
     << setw(16) << "Chi2/NDF"
     << setw(12) << "p-value"
     << setw(12) << "E_cut (eV)"
     << "\n--------------------------------------------------------\n";

    cout << setw(8)  << "0"
     << setw(16) << Form("%.2f/%d", g0y->GetChisquare(), g0y->GetNDF())
     << setw(12) << Form("%.3g", g0y->GetProb())
     << setw(12) << Form("%.1f", cut0y)
     << endl;

    cout << setw(8)  << "250"
     << setw(16) << Form("%.2f/%d", g250y->GetChisquare(), g250y->GetNDF())
     << setw(12) << Form("%.3g", g250y->GetProb())
     << setw(12) << Form("%.1f", cut250y)
     << endl;

    cout << setw(8)  << "1000"
     << setw(16) << Form("%.2f/%d", g1000y->GetChisquare(), g1000y->GetNDF())
     << setw(12) << Form("%.3g", g1000y->GetProb())
     << setw(12) << Form("%.1f", cut1000y)
     << endl;




    

    //Istogrammi distribuzione del numero di hit di rumore (per la stima del valore di aspettazione)
    TH1I *h_noise_hits_X = new TH1I("h_noise_hits_X", "Hit sotto soglia in X;N_{hit};Conteggi", 10, 0, 10);
    TH1I *h_noise_hits_Y = new TH1I("h_noise_hits_Y", "Hit sotto soglia in Y;N_{hit};Conteggi", 10, 0, 10);

    for(int ev=0 ; ev < nX ; ev++){ //ciclo per riempire hist x
      tx->GetEntry(ev);
      int nhits = Ex->size();
      int n_hit_noise = 0;
      for(int ii=0 ; ii < nhits ; ii++){
	if((int)zx->at(ii) == 0 && Ex->at(ii)<cut0x){
	  n_hit_noise++;
	  continue;
	}
        
	if((int)zx->at(ii) == 250 && Ex->at(ii)<cut250x){
	  n_hit_noise++;
	  continue;
	}
        
	if((int)zx->at(ii) == 1000 && Ex->at(ii)<cut1000x){
	  n_hit_noise++;
	  continue;
	}
      }
      h_noise_hits_X->Fill(n_hit_noise);
    }

    for(int ev=0 ; ev < nY ; ev++){ //ciclo per riempire hist y
      ty->GetEntry(ev);
      int nhits = Ey->size();
      int n_hit_noise = 0;
      for(int ii=0 ; ii < nhits ; ii++){
	if((int)zy->at(ii) == 0 && Ey->at(ii)<cut0y){
	  n_hit_noise++;
	  continue;
	}
        
	if((int)zy->at(ii) == 250 && Ey->at(ii)<cut250y){
	  n_hit_noise++;
	  continue;
	}
        
	if((int)zy->at(ii) == 1000 && Ey->at(ii)<cut1000y){
	  n_hit_noise++;
	  continue;
	}
      }
      h_noise_hits_Y->Fill(n_hit_noise);
    }


    cout << "X noise entries = " << h_noise_hits_X->GetEntries() << endl;
    cout << "Y noise entries = " << h_noise_hits_Y->GetEntries() << endl;

    h_noise_hits_X->Print("all");
    
    // ===== Canvas unico con due istogrammi =====
    TCanvas* c_noise_hits = new TCanvas("c_noise_hits",
                                    "Distribuzione hit di rumore",
                                    1200, 600);
    c_noise_hits->Divide(2,1);

    // --- Plot X ---
    c_noise_hits->cd(1);
    h_noise_hits_X->SetLineColor(kBlue+1);
    h_noise_hits_X->SetLineWidth(2);
    h_noise_hits_X->SetFillColorAlpha(kBlue-8, 0.35);
    //h_noise_hits_X->SetStats(1);
    h_noise_hits_X->Draw();

    // --- Plot Y ---
    c_noise_hits->cd(2);
    h_noise_hits_Y->SetLineColor(kRed+1);
    h_noise_hits_Y->SetLineWidth(2);
    h_noise_hits_Y->SetFillColorAlpha(kRed-8, 0.35);
    //h_noise_hits_Y->SetStats(1);
    h_noise_hits_Y->Draw();

    // ===== Calcolo media e sigma =====
    double meanX  = h_noise_hits_X->GetMean();
    double sigmaX = h_noise_hits_X->GetStdDev();

    double meanY  = h_noise_hits_Y->GetMean();
    double sigmaY = h_noise_hits_Y->GetStdDev();

    // ===== Stampa risultati =====
    cout << "\n=== Statistiche hit di rumore ===" << endl;

    cout << "X-plane:" << endl;
    cout << "  Media = " << meanX  << endl;
    cout << "  Sigma = " << sigmaX << endl;

    cout << "\nY-plane:" << endl;
    cout << "  Media = " << meanY  << endl;
    cout << "  Sigma = " << sigmaY << endl;











    //=========== Nuovo tree eventi con taglio a basse energie, solo quelli con hit in almeno 4 piani ========

    TFile* fout = new TFile("cut_output.root", "RECREATE");

    TTree* tkx = new TTree("keptX", "Eventi X che passano il taglio");
    TTree* tky = new TTree("keptY", "Eventi Y che passano il taglio");

    vector<double> xxk, zxk;
    vector<double> yyk, zyk;
    vector<double> Exk, Eyk;

    tkx->Branch("xx", &xxk);
    tkx->Branch("zx", &zxk);
    tkx->Branch("Ex", &Exk);

    tky->Branch("yy", &yyk);
    tky->Branch("zy", &zyk);
    tky->Branch("Ey", &Eyk);

    

    for (int ev = 0; ev < nX; ev++) {
      tx->GetEntry(ev);
      ty->GetEntry(ev);

      xxk.clear();
      zxk.clear();
      yyk.clear();
      zyk.clear();
      Exk.clear();
      Eyk.clear();

    
    for (int i = 0; i < Ex->size(); i++) {  //ciclo per x
        double e = Ex->at(i);
        double z = zx->at(i);
	double x = xx->at(i);

      
	if ((int)z == plane[0]) {
	   if(e >= cut0x){
	     xxk.push_back(x);
	     zxk.push_back(z);
	     Exk.push_back(e);
	   }
	   continue;
         }
	if((int)z == plane[1]) {
	  if(e >= cut250x){
	    xxk.push_back(x);
	    zxk.push_back(z);
	    Exk.push_back(e);
	  }
	  continue;
        }
	if((int)z == plane[4]) {
	  if(e >= cut1000x){
	    xxk.push_back(x);
	    zxk.push_back(z);
	    Exk.push_back(e);
	  }
	  continue;
        }
	else{
	  xxk.push_back(x);
	  zxk.push_back(z);
	  Exk.push_back(e);
	}
    }

    for (int i = 0; i < Ey->size(); i++) {  //ciclo per y
        double e = Ey->at(i);
        double z = zy->at(i);
	double y = yy->at(i);

      
	if ((int)z == plane[0]) {
	   if(e >= cut0y){
	     yyk.push_back(y);
	     zyk.push_back(z);
	     Eyk.push_back(e);
	   }
	   continue;
         }
	if((int)z == plane[1]) {
	  if(e >= cut250y){
	    yyk.push_back(y);
	    zyk.push_back(z);
	    Eyk.push_back(e);
	  }
	  continue;
        }
	if((int)z == plane[4]) {
	  if(e >= cut1000y){
	    yyk.push_back(y);
	    zyk.push_back(z);
	    Eyk.push_back(e);
	  }
	  continue;
        }
	else{
	  yyk.push_back(y);
	  zyk.push_back(z);
	  Eyk.push_back(e);
	}
    }

    
    set<double> set_zx(zxk.begin(), zxk.end());
    set<double> set_zy(zyk.begin(), zyk.end());
    
    if(set_zx.size()>=4) //hit in ameno 4 piani distinti
      tkx->Fill();

    if(set_zy.size()>=4)
      tky->Fill();  

}




    
    //==== Print di roba ===


    cout<<"Numero eventi X="<<tkx->GetEntries()<<endl;
    cout<<"Numero eventi Y="<<tky->GetEntries()<<endl<<endl;
    
     for(int ev=0; ev<5; ev++){
      tkx->GetEntry(ev);
      cout << "\n=== Evento X " << ev << " ===" << endl;
      cout << "Ex.size() = " << Exk.size() << endl;
      for (size_t i = 0; i < Exk.size(); i++) {
         cout << "  Hit " << i
              << "  xx = " << xxk.at(i)
              << "  zx = " << zxk.at(i)
              << "  Ex = " << Exk.at(i)
	      << endl;
        }
    }

    
     //====== Controllo esistenza doppie hit nello stesso piano ===========
     int nxk=tkx->GetEntries();
     int nyk=tky->GetEntries();
     for(int ev=0; ev<nxk; ev++){
       tkx->GetEntry(ev);
       set<double> set_zx(zxk.begin(), zxk.end());

       if(set_zx.size()<zxk.size()){ 
	 cout<<"doppia hit per X trovta all'evento"<<ev<<endl;
	 break;
       }
     }
     for(int ev=0; ev<nyk; ev++){
       tky->GetEntry(ev);
       set<double> set_zy(zyk.begin(), zyk.end());

       if(set_zy.size()<zyk.size()){ 
	 cout<<"doppia hit per Y trovta all'evento"<<ev<<endl;
	 break;
       }
     }

     




     

    
    fout->Write();
    fout->Close();
    cout << "Dati scritti in " << output << endl;
    
    
    //f->Close();
}

  
