
void LoopAll::myStatPhotonAnalysisRed(Util *ut, int jentry){

   if (PADEBUG)
      cout << "myStatPhotonAnalysRed START" << endl;

   counters[0]++;
	
   b_pho_n->GetEntry(jentry);
   b_pho_p4->GetEntry(jentry);
   b_pho_hoe->GetEntry(jentry);
   b_pho_sieie->GetEntry(jentry);

   Float_t best_mass = 0;
   for (int i=0; i<pho_n-1; i++) {
     TLorentzVector *pg1= (TLorentzVector *) pho_p4->At(i);
 
     for (int j=i+1; j<pho_n; j++) {
       TLorentzVector *pg2= (TLorentzVector *) pho_p4->At(j);
       TLorentzVector higgs = (*pg1) + (*pg2);
       Float_t mass = higgs.M();
       if (mass > best_mass){
        if ( pho_hoe[i]          < 0.05   && pho_hoe[j]         < 0.05
         && pho_sieie[i]         < 0.03    && pho_sieie[j]       < 0.03
         && pg1->Pt()            > 30      && pg2->Pt()          > 30
         && fabs(pg1->Eta())     < 2.4     && fabs(pg2->Eta())   < 2.4
         )
         best_mass = mass;
        }
       }
   }
   if (best_mass != 0){
     rooContainer->SetRealVar("mass",best_mass);
   }
 
   if(PADEBUG)
     cout<<"myStatPhotonAnalysisRed END"<<endl;

}
