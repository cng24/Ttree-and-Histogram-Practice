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

#TH1D h1("h1","Histo from a Ttree",10,0,10);
#h1.FillRandom("gaus",10000);
#h1.Draw()
