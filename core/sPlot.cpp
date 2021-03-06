#include "sPlot.h"

#include <memory>
#include "TDirectory.h"

namespace HS{
  namespace FIT{


  
    Bool_t sPlot::Run(){
      cout<<"HS::FIT::sPlot::Do prelimanry fits "<<gDirectory->GetName()<<endl;

      if(FitManager::Run()==kFALSE) return kFALSE;
      
      //Note sPlot is much (10X) faster with tree store
      //Normal fit is 2X faster with vector...
      RooAbsData::setDefaultStorageType(RooAbsData::Tree);
      auto* dataset =dynamic_cast<RooDataSet*>( fCurrDataSet->emptyClone());
      dataset->append(*fCurrDataSet.get());
      RooAbsData::setDefaultStorageType(RooAbsData::Vector);
      
      auto *model=fCurrSetup->Model();
      
      ////////////////////////////////////////////////////////
       //sPlot
       cout<<"HS::FIT::sPlot::Run create sWeights "<<endl;
       fCurrSetup->Parameters().setAttribAll("Constant");
       
       fSPlot.reset(new RooStats::SPlot{"splot_fit", "sPlot Fit",
	     *dataset,model,fCurrSetup->Yields()});
       
       fCurrSetup->Parameters().setAttribAll("Constant",kFALSE);

       CreateWeights();
       
       // PlotDataModel();

       delete dataset;

       return kTRUE;
   }
    
    void sPlot::CreateWeights(){
      //Check that the fit was succesfull
      Double_t TotalYield=0;
      auto yields=fCurrSetup->Yields();
      for(Int_t iy=0;iy<yields.getSize();iy++)
    	TotalYield+=(dynamic_cast<RooRealVar*>(&yields[iy]))->getVal();
      
      if(TotalYield>0){ //got some weights
    	fWeights=std::make_shared<Weights>("HSsWeights");//initialise weights
    	fWeights->SetIDName(fCurrSetup->GetIDBranchName());
    	fWeights->SetTitle(fCurrSetup->GetName());
    	fWeights->SetFile(fCurrSetup->GetOutDir()+TString("Weights")+fCurrSetup->GetName()+fCurrSetup->GetTitle()+".root");
    	ExportWeights();
 	//fWeights->PrintWeight();
	fWeights->SortWeights();

      }
    
      else Warning("sPlot::sPlot()"," total weights 0, fit did not converge. Make sure the non-sweight fit to fix parameters was succesful. No weights will be assigned for these events");
      
    }

    void sPlot::ExportWeights(){
      cout<<"HSsPlot:ExportWeights "<<endl;
      const TString idname=fCurrSetup->GetIDBranchName();
      // cout<<"sPlot::ExportWeights()  "<<idname<<endl;
      const RooArgSet* vars=fCurrDataSet->get(0);
      Bool_t gotID=kFALSE;
      if(vars->find(idname))
	gotID=kTRUE;
      
      auto yields=fCurrSetup->Yields();
      
      Int_t NSpecies=yields.getSize();
      TVectorD eventW(NSpecies); //initialise weights vector
      for(Int_t iw=0;iw<NSpecies;iw++)//set name for each species, 
	fWeights->SetSpecies(TString(yields.at(iw)->GetName()).Remove(0,4));
      
      for(Long64_t ev=0;ev<fCurrDataSet->numEntries();ev++){//loop over events
	for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
	  eventW[iw]=fSPlot->GetSWeight(ev,yields.at(iw)->GetName());//get weight for this species
	if(gotID){//use ID from initial tree
	  const RooArgSet* vars=fCurrDataSet->get(ev);
	  fWeights->FillWeights((Long64_t)vars->getRealValue(idname),eventW);
	} //ID not defined just use entry number in dataset
	else fWeights->FillWeights(ev,eventW);
      }
    }
    
    weights_uptr sPlot::MergeWeights(){
      //in addition combine the weights into 1 and load them
      weights_uptr wts(new Weights("HSsWeights"));
      //Note the output file cannot contain the word Weights (because of Merge), hence Tweights!
      wts->Merge(SetUp().GetOutDir()+"/Weights",
		 SetUp().GetOutDir()+"/"+SetUp().GetName()+"Tweights.root",
		 "HSsWeights");
      //wts->Save();

      //reset to save and reopen
      wts.reset(new Weights{});
      wts->LoadSaved(SetUp().GetOutDir()+TString("Tweights.root"),"HSsWeights");
      
      return std::move(wts);
    }

    void sPlot::WeightedTree(){
      if(!fWeights.get()){
	if(Bins().GetSize()>0)
	  fWeights = MergeWeights();
	else{
	  fWeights.reset(new Weights());
	  fWeights->LoadSaved(SetUp().GetOutDir()+TString("Weights")+SetUp().GetName()+".root","HSsWeights");
	}
      }
      
      //Open tree file
      auto ftree=FiledTree::Read(Data().ParentTreeName(),
      				 Data().ParentName());
      //create a copy in a new file to append the weights to
      //Keep it in file as large trees can use too much memory
      fWeightedFiledTree=(FiledTree::CloneFull(ftree->Tree(),SetUp().GetOutDir()+"DataWeightedTree.root"));
  
      //delete original
      ftree.reset();
      //Add weights to tree
      fWeights->AddToTree(fWeightedFiledTree->Tree().get());
      fWeightedFiledTree->Tree()->SetBranchStatus("*",true);
    }

    void sPlot::DrawWeighted(const TString& var,const TString& wname,TString cut, const TString& opt){
      if(!fWeightedFiledTree.get())
	WeightedTree();

      if(!fWeightedFiledTree.get())
	return;

 
      if(cut==TString()) cut="1";
     
      fWeightedFiledTree->Tree()->Draw(var,wname+"*("+cut+")",opt);

    }
  
  }//namespace FIT 
}//namespace HS
