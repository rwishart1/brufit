#include "PlotResults.h"
#include "RooHSEventsPDF.h"
#include "RooMcmc.h"
#include "FitManager.h"
#include <RooPlot.h>
#include <RooMsgService.h>
#include <TCanvas.h> 
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TLine.h>

namespace HS{
  namespace FIT{
    using namespace RooFit;

    CornerPlot::CornerPlot(Setup *setup, RooMcmc *mcmc)
    {
      fCanvases->SetName(TString("Corner Plot"));
      auto vars=setup->FitVars();

      auto tree = mcmc->GetTree();

      auto& pars = setup->ParsAndYields();
      std::cout<<"The Parameters are: "<<std::endl;
      Int_t Npars = pars.size();
      std::cout<<"The number of parameters is: "<<Npars<<std::endl;
      
      vector<Double_t> param(Npars);
      vector<Double_t> chain_mean(Npars);
      Int_t param_iter=0;

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters once to set values from tree
	  tree->SetBranchAddress(ipar->GetName(), &param[param_iter++]);
	}//close loop over setting params
      
      for(auto var : vars)
	{
      auto canName = "Corner Plot";
      auto canvas = new TCanvas(canName, canName);
      
      fCanvases->Add(canvas);
      canvas->Divide(Npars, Npars);
      //canvas->GetFrame()->SetBorderSize(0);
      
      gStyle->SetOptStat(0);
      gStyle->SetTitleY(1.0);
      gStyle->SetTitleTextColor(2);
      gStyle->SetTitleSize(0.2, "t");
      gStyle->SetLabelSize(0.125, "xy");
      gStyle->SetNdivisions(4, "xy");
      //gStyle->SetCanvasBorderMode(0);
      //gStyle->SetCanvasBorderSize(0);
      //gStyle->SetFrameBorderMode(0);
      //  gStyle->SetFrameBorderSize(0);
      gStyle->SetPadTopMargin(0);
      gStyle->SetPadBottomMargin(0);

      tree->UseCurrentStyle();
      Int_t counter=0;//Counter for new line on canvas      

      for (RooAbsArg* ipar : pars)
	{//Loop over parameters twice to draw corner plot
	 //fix the first parameter 
         //loop over a second parameter to draw the corresponding hists
	  Int_t int_counter =1; //Counter for across canvas

	  for (RooAbsArg* ipar2 : pars)
	    {//second loop 
	      
	      if(ipar==ipar2)
		{
		  auto can = canvas->cd(Npars*counter+int_counter); 
		  can->SetBorderSize(0);
		  can->SetTopMargin(0.0);
		  can->SetBottomMargin(0.1);
		  can->SetLeftMargin(0.19);
		  can->SetRightMargin(0.01);
		  TH1F* hist = new TH1F("hist", "hist", 100,-100000,100000);
		  TString DrawParInd1 = ipar->GetName();
		  TString DrawParInd = DrawParInd1 + ">>hist";
		  tree->Draw(DrawParInd);
		  Double_t mean = hist->GetMean();
		  Double_t sigma = hist->GetRMS();
		  TH1F *hist2 = new TH1F("hist2", "hist2", 100, mean-3*sigma, mean+3*sigma);
		  hist2->SetTitle(DrawParInd1);
		  TString DrawParFin = DrawParInd1 + ">>hist2";
		  tree->Draw(DrawParFin);

		  //Start Line  
		  /*
		  TLine *line = new TLine(hist->GetMean(),0,hist->GetMean(), hist->GetMaximum());
		  line->SetLineColor(kRed);
		  line->Draw(); 
		  */
		  //End Line   

		  //Start Mode Line
		  Int_t binmax = hist2->GetMaximumBin();
		  Double_t mode = hist2->GetBinCenter(binmax);
		  TLine *line = new TLine(mode, 0, mode, hist->GetMaximum());
		  line->SetLineColor(kRed);
		  line->Draw();
		  //End Mode Line
		  counter++;
		  
       		  break;
		}

	      else if(ipar!=ipar2)
		{
		  auto can = canvas->cd(Npars*counter+int_counter);
		  can->SetBorderSize(0);
		  can->SetTopMargin(0);
		  can->SetBottomMargin(0.1);
		  can->SetLeftMargin(0.19);
		  can->SetRightMargin(0.01);
		  /*
		  TH2F *hist = new TH2F("hist", "hist", 100, -100000, 100000, 100, -100000, 100000);
		  TString DrawPar = ipar->GetName();
		  TString DrawPar2 = ipar2->GetName();
		  TString title = DrawPar + ":" + DrawPar2;
		  TString Draw2D = DrawPar + ":" + DrawPar2 + ">>hist";
		  tree->Draw(Draw2D, "", "col");
		  Double_t meanX = hist->GetMean(1);
		  Double_t meanY = hist->GetMean(2);
		  Double_t rmsX = hist->GetRMS(1);
		  Double_t rmsY = hist->GetRMS(2);
		  TH2F *hist2 = new TH2F("hist2", "hist2", 100, meanX-3*rmsX, meanX+3*rmsX, 100, meanY-3*rmsY, meanY+3*rmsY);
		  hist2->SetTitle(title);
		  TString Draw2DFin = DrawPar + ":" + DrawPar2 + ">>hist2";
		  tree->Draw(Draw2DFin,"", "col");
		  */
		  //Start Line
		  /*
		  TLine* lineV = new TLine(meanX,meanY-3*rmsY, meanX, meanY+3*rmsY);
		  TLine* lineH = new TLine(meanX-3*rmsX, meanY, meanX+3*rmsX, meanY);
		  lineV->SetLineColor(kRed);
		  lineH->SetLineColor(kRed);
		  lineV->Draw();
		  lineH->Draw();
		  */
		  //End Line

		  //Start Mode Line One
		  /*
		  Int_t binmax = hist2->GetMaximumBin();
		  Int_t binmaxX, binmaxY, binmaxZ;
		  hist2->GetBinXYZ(binmax, binmaxX, binmaxY, binmaxZ);
		  Double_t modeX =  hist2->GetXaxis()->GetBinCenter(binmaxX);
		  Double_t modeY = hist2->GetYaxis()->GetBinCenter(binmaxY);
		  TLine *lineV = new TLine(modeX, meanY-3*rmsY, modeX, meanY+3*rmsY);
		  TLine *lineH = new TLine(meanX-3*rmsX, modeY, meanX+3*rmsX, modeY);
		  lineV->SetLineColor(kRed);
		  lineH->SetLineColor(kRed);
		  lineV->Draw();
		  lineH->Draw();
		  */ 
		  //End Mode Line One

		  //Start Mode Line Two
		  TH1F* histY = new TH1F("histY", "histY", 100,-100000,100000);
		  TString DrawParInd1Y = ipar->GetName();
		  TString DrawParIndY = DrawParInd1Y + ">>histY";
		  tree->Draw(DrawParIndY);
		  Double_t meanY = histY->GetMean();
		  Double_t sigmaY = histY->GetRMS();
		  TH1F *hist2Y = new TH1F("hist2Y", "hist2Y", 100, meanY-3*sigmaY, meanY+3*sigmaY);
		  TString DrawParFinY = DrawParInd1Y + ">>hist2Y";
		  tree->Draw(DrawParFinY);
		  Int_t binmaxY = hist2Y->GetMaximumBin();
		  Double_t modeY = hist2Y->GetBinCenter(binmaxY);

		  TH1F* histX = new TH1F("histX", "histX", 100,-100000,100000);
		  TString DrawParInd1X = ipar2->GetName();
		  TString DrawParIndX = DrawParInd1X + ">>histX";
		  tree->Draw(DrawParIndX);
		  Double_t meanX = histX->GetMean();
		  Double_t sigmaX = histX->GetRMS();
		  TH1F *hist2X = new TH1F("hist2X", "hist2X", 100, meanX-3*sigmaX, meanX+3*sigmaX);
		  TString DrawParFinX = DrawParInd1X + ">>hist2X";
		  tree->Draw(DrawParFinX);
		  Int_t binmaxX = hist2X->GetMaximumBin();
		  Double_t modeX = hist2X->GetBinCenter(binmaxX);

		  TH2F *hist2 = new TH2F("hist2", "hist2", 100, meanX-3*sigmaX, meanX+3*sigmaX, 100, meanY-3*sigmaY, meanY+3*sigmaY);
		  TString DrawPar = ipar->GetName();
		  TString DrawPar2 = ipar2->GetName();
		  TString title = DrawPar + ":" + DrawPar2;
		  TString Draw2D = DrawPar + ":" + DrawPar2 + ">>hist2";
		  hist2->SetTitle(title);
		  tree->Draw(Draw2D, "", "col");
		 
	      
		 
		  TLine *lineV = new TLine(modeX, meanY-3*sigmaY, modeX, meanY+3*sigmaY);
		  TLine *lineH = new TLine(meanX-3*sigmaX, modeY, meanX+3*sigmaX, modeY);
		  lineV->SetLineColor(kRed);
		  lineH->SetLineColor(kRed);
		  tree->Draw(Draw2D,"", "col");
		  lineV->Draw();
		  lineH->Draw(); 
		  
		  //End Mode Line Two
		  int_counter++;

		}
	      
	    }//End Second loop ipar2
	}//End First loop ipar
      canvas->Modified();
      canvas->Update();
      canvas->Draw();
      
      
	}//Loop over vars

    }//CornerPlot
  }//FIT
}//HS
