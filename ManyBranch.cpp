#include<iostream>
#include "TTree.h"
#include "TMath.h"

int run(){

gROOT->cd();
TTree * t = new TTree("t", "a tree");
TRandom r;
Float_t A, B, C;
t->Branch("B1", &A, "A/F");
t->Branch("B2", &B, "B/F");
t->Branch("B3", &C, "C/F");

A = r.Rndm() * 10; 
B = r.Rndm() * 10; 
C = r.Rndm() * 10; 

t->Fill();
t->Scan()
  
return 0;

}
