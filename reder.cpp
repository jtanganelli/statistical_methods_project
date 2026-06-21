#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <set>

using namespace std;

void reader(TString file, TString output="filtered.root") {

    // Apri il file ROOT
    TFile *f = TFile::Open(file);

    //Recupera Tree
    TTree *tree = (TTree*) f->Get("tree");
    
    double X[3][10];     // X[0]=x, X[1]=y, X[2]=z
    double Edep[2][10];  // Edep[0]=Ex, Edep[1]=Ey

    tree->SetBranchAddress("X", &X);
    tree->SetBranchAddress("Edep", &Edep);


    int nEvents = tree->GetEntries();
    cout << "Numero di eventi nel file: " << nEvents << endl;

    //stampa dati dei primi eventi
    for(int ev = 0 ; ev< 10 ; ev++ ){
      tree->GetEntry(ev);
      cout<<"===EVENTO "<<ev<<"==="<<endl;
      for(int ii = 0 ; ii < 10 ; ii++){
	cout<<"X="<<X[0][ii]
	    <<"     Y="<<X[1][ii]
	    <<"     Z="<<X[2][ii]<<endl;
	cout<<"     Ex="<<Edep[0][ii]
	    <<"     Ey="<<Edep[1][ii]<<endl<<endl;
      }
      cout<<endl;
    }
   
    // Togli i dati con valore di default
    TFile *fout = new TFile("filtered.root", "RECREATE");
    
    vector<double> xx, zx;
    vector<double> yy, zy;
    vector<double> Ex, Ey;

    TTree *tx = new TTree("filteredX", "Filtered hits X");
    TTree *ty = new TTree("filteredY", "Filtered hits Y");
    tx->Branch("xx",  &xx);
    tx->Branch("zx",  &zx);
    ty->Branch("yy",  &yy);
    ty->Branch("zy",  &zy);
    tx->Branch("Ex", &Ex);
    ty->Branch("Ey", &Ey);

    
    for(int ev=0; ev<nEvents; ev++){
      tree->GetEntry(ev);
      xx.clear();
      zx.clear();
      yy.clear();
      zy.clear();
      Ex.clear();
      Ey.clear();

      for(int h=0;h<10;h++){
	if(X[2][h]==-999)
	  continue;
	if(X[0][h]!=-999 && Edep[0][h]!=0){
	  xx.push_back(X[0][h]);
	  zx.push_back(X[2][h]);
	  Ex.push_back(Edep[0][h]);
	}
	if(X[1][h]!=-999 && Edep[1][h]!=0){
	  yy.push_back(X[1][h]);
	  zy.push_back(X[2][h]);
	  Ey.push_back(Edep[1][h]);
	}
      }
      tx->Fill();
      ty->Fill();

    }

    
    
      // print dei primi 5 eventi X
    for(int ev=0; ev<5; ev++){
      tx->GetEntry(ev);
      cout << "\n=== Evento X " << ev << " ===" << endl;
      cout << "Ex.size() = " << Ex.size() << endl;
      for (size_t i = 0; i < Ex.size(); i++) {
         cout << "  Hit " << i
              << "  xx = " << xx.at(i)
              << "  zx = " << zx.at(i)
              << "  Ex = " << Ex.at(i)
	      << endl;
        }
    }

    
    fout->Write();
    fout->Close();
    f->Close();
    cout << "Dati scritti in " << output << endl;
}
