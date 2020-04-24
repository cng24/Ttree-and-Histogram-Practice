# Ttree-and-Histogram-Practice

#making a Ttree in root with multiple Branches

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

#making a Ttree in root with multiple leaves

gROOT->cd();
TTree * t = new TTree("t", "a tree");
TRandom r;
Float_t A, B, C;
t->Branch("Rand", &A, "A/F",
                  &B, "B/F",
                  &C, "C/F");

A = r.Rndm() * 10; 
B = r.Rndm() * 10; 
C = r.Rndm() * 10; 

t->Fill();
t->Scan()



#making a general historgram

#make sure that -Y is put after ssh when signing into root for 
#the draw function on the histogram code

new TBrowser()

TH1D *h=new TH1D("h","h",100,0,100);
h-> Draw()



#making a historgram in a Ttree




#saving a file

TFile f("exampleTree", "recreate");
t->Write();
f->Close();

#check to see if it saved
f.ls() 


