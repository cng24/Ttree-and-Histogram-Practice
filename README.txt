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
Float_t Randos[3] = {0};

t->Branch("A", Randos, "TrialA[3]/F");
Randos[0] = r.Rndm() * 10;
Randos[1] = r.Rndm() * 10;
Randos[2] = r.Rndm() * 10;
t->Fill();
t->Scan()



#making a general historgram
#make sure that -Y is put after ssh when signing into root for 
#the draw function on the histogram code

new TBrowser() 

TH1D *h=new TH1D("h","h",100,0,100);
h-> Draw()




#making a historgram in a Ttree

gROOT->cd();
TTree *t = new TTree("t", "a tree");
TRandom r;
Float_t Randos[3] = {0};

t->Branch("A", Randos, "TrialA[3]/F");
Randos[0] = r.Rndm() * 10;
Randos[1] = r.Rndm() * 10;
Randos[2] = r.Rndm() * 10;
t->Fill();
t->Scan()

TH1D h1("h1","Histo from a Ttree",1000,0,10);
#name, title, number of bins, 

#HOW DO I FILL IT WITH THE ENTIRES FROM THE TTREE...https://root.cern/root/html530/TTree.html (potentially helpful?)
h1.Fill(A);
h1.Draw()

TH1D h1("h1","Histo from a Ttree",10,0,10);
h1.FillRandom("gaus",10000);
h1.Draw()




#saving a file

TFile f("exampleTree", "recreate");
t->Write();
f->Close();

#check to see if it saved
f.ls() 


