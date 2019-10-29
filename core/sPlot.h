////////////////////////////////////////////////////////////////
///
///Class:               sPlot
///Description:
///           

#pragma once


#include "FitManager.h"
#include "Setup.h"
#include "Weights.h"
#include <RooStats/SPlot.h>
 
namespace HS{
  namespace FIT{

    using splot_uptr = std::unique_ptr<RooStats::SPlot>;
    using weights_ptr =std::shared_ptr<Weights>;
    using weights_uptr =std::unique_ptr<Weights>;
    using tree_uptr =std::unique_ptr<TTree>;
    
    class sPlot  : public FitManager{
      
    public:
      sPlot()=default;
      sPlot(const sPlot&)=default;
      sPlot(sPlot&&)=default;
      ~sPlot() override =default;
      sPlot& operator=(const sPlot& other) = default;
      sPlot& operator=(sPlot&& other) = default;
      
      void Run() override;
 
      void Reset() override{
	FitManager::Reset();
        fSPlot.reset();
	fWeights.reset();
       }

      void CreateWeights();
      void ExportWeights();
      weights_uptr MergeWeights();
      void DrawWeighted(const TString& var,const TString& wname,TString cut="1",const TString& opt="");

      TTree* GetWeightedTree(){return fWeightedFiledTree->Tree().get();};
      void DeleteWeightedTree(){fWeightedFiledTree.reset();};
      
    protected:
      
      void WeightedTree();
      
    private:
      splot_uptr fSPlot; //!sPlot object
      weights_ptr fWeights;//!
      tree_uptr fWeightedTree;//!
      filed_uptr fWeightedFiledTree;//!

      ClassDefOverride(HS::FIT::sPlot,1);
    };
    
  }//namespace FIT
}//namespace HS

