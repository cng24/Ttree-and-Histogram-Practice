#include<iostream>
#include "TTree.h"
#include "TMath.h"

int run(){

gROOT->cd();
TTree * t = new TTree("t", "a tree");
TRandom r;
Float_t A, B, C;
t->Branch("Rand", &A, "A/F");
t->Branch("Ran", &B, "B/F");
t->Branch("Ra", &C, "C/F");

A = r.Rndm() * 10; 
B = r.Rndm() * 10; 
C = r.Rndm() * 10; 

t->Fill();
t->Scan()
  
return 0;

}
