#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TBranch.h>
#include <TFile.h>
#include <TH3.h>
#include <TF1.h>
#include <TF3.h>
#include <TString.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TPad.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooCategory.h>
#include <RooArgList.h>
#include <RooBinning.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <TSystem.h>

#include <HWWLVJRooPdfs.h>

double cwwws[]        = {-12.,-6.,-2.,0.,2.,6.,12.};
double cwwws1[]        = {-12.,-6.,0.,6.,12.};
double cwwws2[]        = {-2.,0.,2.};
double ccws[]        = {-20.,-10.,-3.5,0.,3.5,10.,20.};
double ccws1[]        = {-20.,-10.,0.,10.,20.};
double ccws2[]        = {-3.5,0.,3.5};
double cbs[]        = {-60.,-30.,-10.,0.,10.,30.,60.};
double cbs1[]        = {-60.,-30.,0.,30.,60.};
double cbs2[]        = {-10.,0.,10.};
double vals[150][3];

double normSM;
RooWorkspace w("w","w");
RooWorkspace w2("w2","w2");


void interference(TString ch,TString sample)
{
    //vals gives the atgc values, SM for [0][i]
    vals[0][0]    = 0.;
    vals[0][1]    = 0.;
    vals[0][2]    = 0.;
    //set other atgc values of vals
    int count = 1; 
    for(int i = 0; i<5; i++)
        for(int j = 0; j<5; j++)
            for(int k = 0; k<5; k++)
            {
                if(cwwws1[i]==12 and ccws1[j]==20 and cbs1[k]==60)
                    continue;
                if(cwwws1[i]!=0 or ccws1[j]!=0 or cbs1[k]!=0)
                {
                    vals[count][0] = cwwws1[i];
                    vals[count][1] = ccws1[j];
                    vals[count][2] = cbs1[k];
                        count++;
                }
            }
    for(int i = 0; i<3; i++)
        for(int j = 0; j<3; j++)
            for(int k = 0; k<3; k++)
                if(cwwws2[i]!=0 or ccws2[j]!=0 or cbs2[k]!=0)
                {
                    vals[count][0] = cwwws2[i];
                    vals[count][1] = ccws2[j];
                    vals[count][2] = cbs2[k];
                    count++;
                }
    vector<TString> histonames;
    for(unsigned int i = 0; i<150; i++)
    {
        TString cwww_tmp    = vals[i][0] >= 0 ? ("cwww"+to_string(int(vals[i][0]))).c_str() : ("cwww_"+to_string(int(abs(vals[i][0])))).c_str();
        TString ccw_tmp        = vals[i][1] >= 0 ? ("ccw"+to_string(int(vals[i][1]))).c_str() : ("ccw_"+to_string(int(abs(vals[i][1])))).c_str();
        TString cb_tmp        = vals[i][2] >= 0 ? ("cb"+to_string(int(vals[i][2]))).c_str() : ("cb_"+to_string(int(abs(vals[i][2])))).c_str();
        if(vals[i][1] - int(vals[i][1]) != 0)
            //if it looks stupid but it works, it ain't stupid ;)
            ccw_tmp        = vals[i][1] >= 0 ? ("ccw"+to_string(int(vals[i][1]))+"_"+to_string(abs(int(10*(vals[i][1]-int(vals[i][1])))))).c_str() : ("ccw_"+to_string(int(abs(vals[i][1])))+"_"+to_string(abs(int(10*(vals[i][1]-int(vals[i][1])))))).c_str();
        TString name_tmp    = cwww_tmp + ccw_tmp + cb_tmp;
        histonames.push_back(name_tmp);
    }
    

    int binlo    = 900;
    int binhi    = 3500;
    int nbins    = 25;

    RooBinning bins(binlo,binhi);
    bins.addUniform(nbins,binlo,binhi);

    TFile * fileIn     = TFile::Open(TString("Input/"+sample+"atgc_tree.root"));

    TTreeReader reader("tree",fileIn);

    TTreeReaderValue<std::vector<double>> weights(reader,"weight");
    TTreeReaderValue<double> MWW_tree(reader,"MWW");
    TTreeReaderValue<int> channel_tree(reader,"channel");

    RooRealVar cwww("cwww","cwww",0,-15,15);        cwww.setConstant(true);
    RooRealVar ccw("ccw","ccw",0,-25,25);           ccw.setConstant(true);
    RooRealVar cb("cb","cb",0,-70,70);              cb.setConstant(true);
    RooRealVar MWW("MWW","MWW",2500,binlo,binhi);
    RooRealVar weight("weight","weight",1,0,500);
    RooCategory cat_cwww("cat_cwww","cat_cwww");
    RooCategory cat_ccw("cat_ccw","cat_ccw");
    RooCategory cat_cb("cat_cb","cat_cb");
    w.import(weight);

    RooArgSet parset(cwww,ccw,cb,weight);

    RooDataSet data3D("data3D","data3D",parset,"weight");

    for(unsigned int i = 0; i<150; i++)
    {
        RooDataHist hist(histonames[i],histonames[i],RooArgSet(MWW));
        w.import(hist);
    }

    RooDataHist hist_diff_cwww("hist_diff_cwww","hist_diff_cwww",RooArgSet(MWW));
    RooDataHist hist_diff_ccw("hist_diff_ccw","hist_diff_ccw",RooArgSet(MWW));
    RooDataHist hist_diff_cb("hist_diff_cb","hist_diff_cb",RooArgSet(MWW));
    RooDataHist hist_diff_cwww_ccw("hist_diff_cwww_ccw","hist_diff_cwww_ccw",RooArgSet(MWW));
    RooDataHist hist_diff_cwww_cb("hist_diff_cwww_cb","hist_diff_cwww_cb",RooArgSet(MWW));
    RooDataHist hist_diff_ccw_cb("hist_diff_ccw_cb","hist_diff_ccw_cb",RooArgSet(MWW));

    int channel = (ch=="el") ? 1 : 2;
    reader.SetEntry(0);
    int tmp = 0;    
    while(reader.Next())
    {
        MWW.setVal(*MWW_tree);
        if(*channel_tree==channel and MWW.getVal()>binlo)
        {
            for (int i = 0; i<150; i++)
            {
                double weight_tmp = (*weights)[i];
                w.data(histonames[i])->add(RooArgSet(MWW),weight_tmp);
            }
            //the order of weights can be found in order_of_weights.txt
            double weight_cwww_tmp        = (*weights)[112]-(*weights)[13];
            hist_diff_cwww.add(RooArgSet(MWW),weight_cwww_tmp);
            double weight_ccw_tmp        = (*weights)[72]-(*weights)[53];
            hist_diff_ccw.add(RooArgSet(MWW),weight_ccw_tmp);
            double weight_cb_tmp        = (*weights)[64]-(*weights)[61];
            hist_diff_cb.add(RooArgSet(MWW),weight_cb_tmp);
            double weight_cwww_ccw_tmp    = ((*weights)[122]-(*weights)[23]) - weight_cwww_tmp;//([cwww12,ccw20]-[cwww-12,ccw20]) - ([cwww12]-[cwww-12])
            hist_diff_cwww_ccw.add(RooArgSet(MWW),weight_cwww_ccw_tmp);
            double weight_cwww_cb_tmp    = ((*weights)[114]-(*weights)[110]) - weight_cb_tmp;//([cwww12,cb60]-[cwww12,cb-60]) - ([cb60]-[cb-60])
            hist_diff_cwww_cb.add(RooArgSet(MWW),weight_cwww_cb_tmp);
            double weight_ccw_cb_tmp     = ((*weights)[74]-(*weights)[70]) - weight_cb_tmp;//([ccw20,cb60]-[ccw20,cb-60]) - ([cb60]-[cb-60])
            hist_diff_ccw_cb.add(RooArgSet(MWW),weight_ccw_cb_tmp);
        }
        tmp++;
        if(tmp%5000==0)
            std::cout<<tmp<<std::endl;
    }


    gStyle->SetOptStat(0);
    TH1F* hist_cwwwccw         = (TH1F*)hist_diff_cwww_ccw.createHistogram("hist_cwww_ccw",MWW,RooFit::Binning(bins));
    hist_cwwwccw->Sumw2();hist_cwwwccw->SetTitle("");hist_cwwwccw->GetYaxis()->SetTitle("Arbitrary units");hist_cwwwccw->GetXaxis()->SetTitle("M_{WW}");hist_cwwwccw->SetMarkerColor(kBlack);
    TH1F* hist_ccwcb        = (TH1F*)hist_diff_ccw_cb.createHistogram("hist_ccw_cb",MWW,RooFit::Binning(bins));
    hist_ccwcb->Sumw2();hist_ccwcb->SetTitle("");hist_cwwwccw->GetYaxis()->SetTitle("Arbitrary units");hist_ccwcb->GetXaxis()->SetTitle("M_{WW}");hist_ccwcb->SetMarkerColor(kBlack);
    TH1F* hist_cwwwcb        = (TH1F*)hist_diff_cwww_cb.createHistogram("hist_cwww_cb",MWW,RooFit::Binning(bins));
    hist_cwwwcb->Sumw2();hist_cwwwcb->SetTitle("");hist_cwwwccw->GetYaxis()->SetTitle("Arbitrary units");hist_cwwwcb->GetXaxis()->SetTitle("M_{WW}");hist_cwwwcb->SetMarkerColor(kBlack);

    TCanvas can_cwwwccw("cwww_ccw","cwww_ccw",1);
    TCanvas can_ccwcb("ccw_cb","ccw_cb",1);
    TCanvas can_cwwwcb("cwww_cb","cwww_cb",1);

    can_cwwwccw.cd();
    hist_cwwwccw->Fit("expo");
    hist_cwwwccw->Draw("E1");
    can_cwwwccw.Draw();can_cwwwccw.Update();
    float slopeval_cwwwccw      = hist_cwwwccw->GetFunction("expo")->GetParameter(1);

    can_ccwcb.cd();
    TString fit_string = sample=="WZ"? "negexpo" : "expo";
    TF1 negexpo("negexpo","-exp([0]+[1]*x)",900,3500);
    hist_ccwcb->Fit(fit_string);
    hist_ccwcb->Draw("E1");
    can_ccwcb.Draw();can_ccwcb.Update();
    float slopeval_ccwcb        = hist_ccwcb->GetFunction(fit_string)->GetParameter(1);

    can_cwwwcb.cd();
    hist_cwwwcb->Fit(fit_string);
    hist_cwwwcb->Draw("E1");
    can_cwwwcb.Draw();can_cwwwcb.Update();
    std::cout<<"<ENTER>"<<std::endl;
    int pause;
    pause=getchar();


    //make complete signal model for plots
    w.factory("Exponential:SM_Pdf(MWW,a1[-0.001,-0.01,0.])");
    w.factory("Exponential:Pdf_ccw_lin(MWW,a33[-0.001,-0.01,0.])");
    w.factory("Exponential:Pdf_cb_lin(MWW,a44[-0.001,-0.01,0.])");
    w.factory("Exponential:Int_cwww_ccw(MWW,a5[-0.0001,-0.01,0.01])");
    w.factory("Exponential:Int_cwww_cb(MWW,a6[-0.001,-0.01,0.01])");
    w.factory("Exponential:Int_ccw_cb(MWW,a7[-0.001,-0.01,0.01])");
    w.factory("Exponential:Pdf_cwww(MWW,a2[-0.001,-0.01,0.])");
    w.factory("Exponential:Pdf_ccw(MWW,a3[-0.001,-0.01,0.])");
    w.factory("Exponential:Pdf_cb(MWW,a4[-0.0001,-0.1,0.])");


    RooAbsPdf * SM_Pdf          = w.pdf("SM_Pdf");
    RooAbsPdf * Pdf_ccw_lin     = w.pdf("Pdf_ccw_lin");
    RooAbsPdf * Pdf_cb_lin      = w.pdf("Pdf_cb_lin");
    RooAbsPdf * Int_cwww_ccw    = w.pdf("Int_cwww_ccw");
    RooAbsPdf * Int_cwww_cb     = w.pdf("Int_cwww_cb");
    RooAbsPdf * Int_ccw_cb      = w.pdf("Int_ccw_cb");
    RooAbsPdf * Pdf_cwww        = w.pdf("Pdf_cwww");
    RooAbsPdf * Pdf_ccw         = w.pdf("Pdf_ccw");
    RooAbsPdf * Pdf_cb          = w.pdf("Pdf_cb");

    RooRealVar N_SM("N_SM","N_SM",w.data("cwww0ccw0cb0")->sumEntries());        //hist0
    RooRealVar N__12("N_cwww__12","N_cwww__12",w.data("cwww_12ccw0cb0")->sumEntries());    //hist13, hist128
    RooRealVar N_12("N_cwww_12","N_cwww_12",w.data("cwww12ccw0cb0")->sumEntries());    //hist112, hist145
    RooRealVar N__20("N_ccw__20","N_ccw__20",w.data("cwww0ccw_20cb0")->sumEntries());    //hist53, hist134
    RooRealVar N_20("N_ccw_20","N_ccw_20",w.data("cwww0ccw20cb0")->sumEntries());    //hist72, hist139
    RooRealVar N__60("N_cb__60","N_cb__60",w.data("cwww0ccw0cb_60")->sumEntries());    //hist61, hist136
    RooRealVar N_60("N_cb_60","N_cb_60",w.data("cwww0ccw0cb60")->sumEntries());    //hist64, hist137
    RooRealVar N_12_20("N_cwww_ccw_12_20","N_cwww_ccw_12_20",w.data("cwww12ccw20cb0")->sumEntries());
    RooRealVar N_12_60("N_cwww_cb_12_60","N_cwww_cb_12_60",w.data("cwww12ccw0cb60")->sumEntries());
    RooRealVar N_20_60("N_ccw_cb_20_60","N_ccw_cb_20_60",w.data("cwww0ccw20cb60")->sumEntries());
    RooRealVar N_4norm("N_4norm","N_4norm",w.data("cwww_12ccw_20cb_60")->sumEntries());


    RooRealVar N2_tmp("N2_tmp","N2_tmp",(N_12.getVal()+N__12.getVal())/2-N_SM.getVal());
    RooRealVar N4_tmp("N4_tmp","N4_tmp",(N_20.getVal()+N__20.getVal())/2-N_SM.getVal());
    RooRealVar N5_tmp("N5_tmp","N5_tmp",(N_20.getVal()-N__20.getVal())/2);
    RooRealVar N6_tmp("N6_tmp","N6_tmp",(N_60.getVal()+N__60.getVal())/2-N_SM.getVal());
    RooRealVar N7_tmp("N7_tmp","N7_tmp",(N_60.getVal()-N__60.getVal())/2);
    RooRealVar N8_tmp("N8_tmp","N8_tmp",(N_12_20.getVal()+N_SM.getVal())-(N_12.getVal()+N_20.getVal()));
    RooRealVar N9_tmp("N9_tmp","N9_tmp",(N_12_60.getVal()+N_SM.getVal())-(N_12.getVal()+N_60.getVal()));
    RooRealVar N10_tmp("N10_tmp","N10_tmp",(N_20_60.getVal()+N_SM.getVal())-(N_20.getVal()+N_60.getVal()));

    RooArgList N_norm_list(N_SM,N2_tmp,N4_tmp,N5_tmp,N6_tmp,N7_tmp);
    N_norm_list.add(RooArgList(cwww,ccw,cb));
    N_norm_list.add(RooArgList(N8_tmp,N9_tmp,N10_tmp));

    TString cwww_f          = "+@1*(@6/12)**2";
    TString ccw_f           = "+@2*(@7/20)**2";
    TString cb_f            = "+@4*(@8/60)**2";
    //TString cwww_lin_f    = "+((@1-@2)/2)*(@7/12)";
    TString ccw_lin_f       = "+@3*(@7/20)";
    TString cb_lin_f        = "+@5*(@8/60)";
    TString cwww_ccw_f      = "+@9*(@6/12)*(@7/20)";
    TString cwww_cb_f       = "+@10*(@6/12)*(@8/60)";
    TString ccw_cb_f        = "+@11*(@7/20)*(@8/60)";
    RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+ccw_lin_f+cb_f+cb_lin_f+cwww_ccw_f+cwww_cb_f+ccw_cb_f,N_norm_list);

    RooFormulaVar N1("N1","@1/@0",RooArgList(Norm,N_SM));
    RooFormulaVar N2("N2","(@2*(@1/12)**2)/@0",RooArgList(Norm,cwww,N2_tmp));
    //RooFormulaVar N3("N3","(((@1-@2)/2)*(@3/12))/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
    RooFormulaVar N4("N4","(@2*(@1/20)**2)/@0",RooArgList(Norm,ccw,N4_tmp));
    RooFormulaVar N5("N5","(@2*(@1/20))/@0",RooArgList(Norm,ccw,N5_tmp));
    RooFormulaVar N6("N6","(@2*(@1/60)**2)/@0",RooArgList(Norm,cb,N6_tmp));
    RooFormulaVar N7("N7","(@2*(@1/60))/@0",RooArgList(Norm,cb,N7_tmp));
    RooFormulaVar N8("N8","(@3*(@1/12)*(@2/20))/@0",RooArgList(Norm,cwww,ccw,N8_tmp));
    RooFormulaVar N9("N9","(@3*(@1/12)*(@2/60))/@0",RooArgList(Norm,cwww,cb,N9_tmp));
    RooFormulaVar N10("N10","(@3*(@1/20)*(@2/60))/@0",RooArgList(Norm,ccw,cb,N10_tmp));



    //RooArgList N_list(N1,N2,N4,N6);
    RooArgList N_list(N1,N2,N4,N5,N6,N7);
    N_list.add(RooArgList(N8,N9,N10));
    //RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_ccw,*Pdf_cb);
    RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_ccw,*Pdf_ccw_lin,*Pdf_cb,*Pdf_cb_lin);
    Pdf_list.add(RooArgList(*Int_cwww_ccw,*Int_cwww_cb,*Int_ccw_cb));
    Pdf_list.Print();
    N_list.Print();
    RooAddPdf model1("model","model",Pdf_list,N_list);


    RooRealVar * a1         = w.var("a1");
    RooRealVar * a2         = w.var("a2"); a2->setConstant(true);
    //RooRealVar * a22    = w.var("a22"); a22->setConstant(true);
    RooRealVar * a3         = w.var("a3"); a3->setConstant(true);
    RooRealVar * a33        = w.var("a33");    a33->setConstant(true);
    RooRealVar * a4         = w.var("a4");     a4->setConstant(true);
    RooRealVar * a44        = w.var("a44");    a44->setConstant(true);
    RooRealVar * a5         = w.var("a5");    a5->setConstant(true);
    RooRealVar * a6         = w.var("a6");    a6->setConstant(true);
    RooRealVar * a7         = w.var("a7");    a7->setConstant(true);


//SM-fit
    cwww.setVal(0);    ccw.setVal(0); cb.setVal(0);     
    model1.fitTo(*w.data("cwww0ccw0cb0"));//hist0
    a1->setConstant(true);
//SM-interference-fits
    double N_SM_tmp_val = N_SM.getVal();//SM
    double N2_tmp_val = N2_tmp.getVal();//cwww quad
    double N4_tmp_val = N4_tmp.getVal();//ccw quad
    double N6_tmp_val = N6_tmp.getVal();//cb quad
    double N5_tmp_val = N5_tmp.getVal();//ccw lin
    double N7_tmp_val = N7_tmp.getVal();//cb lin
    N_SM.setVal(0);N2_tmp.setVal(0);N4_tmp.setVal(0);N6_tmp.setVal(0);    
    
    cwww.setVal(0);ccw.setVal(20); cb.setVal(0);
    a33->setConstant(false);
    model1.fitTo(hist_diff_ccw);
    a33->setConstant(true);
    
    cwww.setVal(0);ccw.setVal(0); cb.setVal(60);
    a44->setConstant(false);
    model1.fitTo(hist_diff_cb);
    a44->setConstant(true);


    //int cwww-ccw
    a5->setVal(slopeval_cwwwccw);
    a5->setConstant(true);
    a7->setVal(slopeval_ccwcb);
    a7->setConstant(true);
    //------------------------

    N_SM.setVal(N_SM_tmp_val);
    N2_tmp.setVal(N2_tmp_val);
    N4_tmp.setVal(N4_tmp_val);
    N5_tmp.setVal(N5_tmp_val);
    N6_tmp.setVal(N6_tmp_val);
    N7_tmp.setVal(N7_tmp_val);
//cwww-fit
    cwww.setVal(12); ccw.setVal(0); cb.setVal(0); 
    a2->setConstant(false);
    model1.fitTo(*w.data("cwww12ccw0cb0"));//hist13, hist128
    a2->setConstant(true);
//ccw-fit
    cwww.setVal(0);    ccw.setVal(20); cb.setVal(0); 
    a3->setConstant(false);
    model1.fitTo(*w.data("cwww0ccw20cb0"));//hist53, hist134
    a3->setConstant(true);
//cb-fit
    cwww.setVal(0); ccw.setVal(0); cb.setVal(60);
    a4->setConstant(false);
    model1.fitTo(*w.data("cwww0ccw0cb60"));//hist61, hist 136
    a4->setConstant(true);


    w2.import(model1);

    TFile * fileOut = new TFile("Output/genlevel_"+sample+"_"+ch+".root","RECREATE");
    w2.Write();
    fileOut->Close();





    ////////////////////
    /////SOME PLOTS/////
    ////////////////////
/*uncomment for plots
    //scaling for ccw and cd (for plots)
    double sum_SM    = w.data("cwww0ccw0cb0")->sumEntries();
    TH1F hist4scale_ccw("hist4scale_ccw","hist4scale_ccw",3,-30,30);
    hist4scale_ccw.SetBinContent(1,w.data("cwww0ccw_20cb0")->sumEntries()/sum_SM);
    hist4scale_ccw.SetBinContent(2,w.data("cwww0ccw0cb0")->sumEntries()/sum_SM);
    hist4scale_ccw.SetBinContent(3,w.data("cwww0ccw20cb0")->sumEntries()/sum_SM);
    hist4scale_ccw.Fit("pol2");
    TF1 * fitfunc_ccw    = hist4scale_ccw.GetFunction("pol2");
    RooRealVar par0_ccw("par0_ccw","par0_ccw",fitfunc_ccw->GetParameter(0));
    RooRealVar par1_ccw("par1_ccw","par1_ccw",fitfunc_ccw->GetParameter(1));
    RooRealVar par2_ccw("par2_ccw","par2_ccw",fitfunc_ccw->GetParameter(2));
    RooFormulaVar scaleshape_ccw("scaleshape_ccw","scaleshape_ccw","(@0+@1*@3+@2*@3**2)",RooArgList(par0_ccw,par1_ccw,par2_ccw,ccw));

    TH1F hist4scale_cb("hist4scale_cb","hist4scale_cb",3,-90,90);
    hist4scale_cb.SetBinContent(1,w.data("cwww0ccw0cb_60")->sumEntries()/sum_SM);
    hist4scale_cb.SetBinContent(2,w.data("cwww0ccw0cb0")->sumEntries()/sum_SM);
    hist4scale_cb.SetBinContent(3,w.data("cwww0ccw0cb60")->sumEntries()/sum_SM);
    hist4scale_cb.Fit("pol2");
    TF1 * fitfunc_cb    = hist4scale_cb.GetFunction("pol2");
    RooRealVar par0_cb("par0_cb","par0_cb",fitfunc_cb->GetParameter(0));
    RooRealVar par1_cb("par1_cb","par1_cb",fitfunc_cb->GetParameter(1));
    RooRealVar par2_cb("par2_cb","par2_cb",fitfunc_cb->GetParameter(2));
    RooFormulaVar scaleshape_cb("scaleshape_cb","scaleshape_cb","(@0+@1*@3+@2*@3**2)",RooArgList(par0_cb,par1_cb,par2_cb,cb));


    TCanvas cc1("smintccwpos","smintccwpos",1);
    cc1.cd();cc1.SetLogy();cc1.SetRightMargin(0.25);cc1.SetTickx();cc1.SetTicky();
    RooPlot * pp = MWW.frame(900,3500);
    w.data("cwww0ccw20cb0")->plotOn(pp,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins),RooFit::Name("ccw20"));
    cwww.setVal(0);ccw.setVal(20);cb.setVal(0);
    model1.plotOn(pp,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent),RooFit::Name("ccw20m"));
    w.data("cwww0ccw3_5cb0")->plotOn(pp,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins),RooFit::Name("ccw35"));
    cwww.setVal(0);ccw.setVal(3.5);cb.setVal(0);
    model1.plotOn(pp,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent),RooFit::Name("ccw35m"));
    w.data("cwww0ccw10cb0")->plotOn(pp,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins),RooFit::Name("ccw10"));
    cwww.setVal(0);ccw.setVal(10);cb.setVal(0);
    model1.plotOn(pp,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent),RooFit::Name("ccw10m"));
    w.data("cwww0ccw0cb0")->plotOn(pp,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));

    for(unsigned int i=0; i<10; i++)
    {
        ccw.setVal(2*i);
        if(i==0)
            model1.plotOn(pp,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
        else
            model1.plotOn(pp,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
        cout<<i<<" : "<<scaleshape_ccw.getVal()*sum_SM<<endl;
    }
    pp->GetYaxis()->SetRangeUser(0.002,5);
    pp->GetYaxis()->SetTitle("arb. units");
    pp->GetXaxis()->SetTitle("M_{WW} (GeV)");
    pp->SetTitle("");
    TLegend tt1(0.75,0.5,0.99,0.9);
    tt1.SetFillColor(0);
    tt1.SetFillStyle(0);
    tt1.SetBorderSize(0);
    tt1.SetHeader(ch);
    tt1.AddEntry(pp->getObject(0),"c_{W}/#Lambda^{2}=20 TeV^{-2}","le");
    tt1.AddEntry(pp->getObject(4),"c_{W}/#Lambda^{2}=10 TeV^{-2}","le");
    tt1.AddEntry(pp->getObject(2),"c_{W}/#Lambda^{2}=3.5 TeV^{-2}","le");
    tt1.AddEntry(pp->getObject(6),"c_{W}/#Lambda^{2}=0 TeV^{-2}","le");
    tt1.AddEntry(pp->getObject(8),"c_{W}/#Lambda^{2}=2,...,18 TeV^{-2}","l");
    pp->Draw();
    tt1.Draw();
    cc1.Draw();cc1.Update();

    TCanvas ccc2("smintccwneg","smintccwneg",1);
    ccc2.cd();ccc2.SetLogy();ccc2.SetRightMargin(0.25);ccc2.SetTickx();ccc2.SetTicky();
    RooPlot * pp2 = MWW.frame(900,3500);
    w.data("cwww0ccw_20cb0")->plotOn(pp2,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(-20);cb.setVal(0);
    model1.plotOn(pp2,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw_10cb0")->plotOn(pp2,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(-10);cb.setVal(0);
    model1.plotOn(pp2,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw_3_5cb0")->plotOn(pp2,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(-3.5);cb.setVal(0);
    model1.plotOn(pp2,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb0")->plotOn(pp2,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));

    ccw.setVal(-20);
    for(int i=0; i<10; i++)
    {
        ccw.setVal(2*i*(-1));
        if(i==0)
            model1.plotOn(pp2,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
        else
            model1.plotOn(pp2,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
        cout<<i<<" : "<<ccw.getVal()<<", "<<scaleshape_ccw.getVal()*sum_SM<<endl;
    }
    pp2->GetYaxis()->SetRangeUser(0.002,5);
    pp2->GetYaxis()->SetTitle("arb. units");
    pp2->GetXaxis()->SetTitle("M_{WW} (GeV)");
    pp2->SetTitle("");
    TLegend tt2(0.75,0.5,0.99,0.9);
    tt2.SetFillColor(0);
    tt2.SetFillStyle(0);
    tt2.SetBorderSize(0);
    tt2.SetHeader(ch);
    tt2.AddEntry(pp2->getObject(0),"c_{W}/#Lambda^{2}=-20 TeV^{-2}","le");
    tt2.AddEntry(pp2->getObject(2),"c_{W}/#Lambda^{2}=-10 TeV^{-2}","le");
    tt2.AddEntry(pp2->getObject(4),"c_{W}/#Lambda^{2}=-3.5 TeV^{-2}","le");
    tt2.AddEntry(pp2->getObject(6),"c_{W}/#Lambda^{2}=0 TeV^{-2}","le");
    tt2.AddEntry(pp2->getObject(8),"c_{W}/#Lambda^{2}=-2,...,-18 TeV^{-2}","l");
    pp2->Draw();
    tt2.Draw();
    ccc2.Draw();ccc2.Update();


    TCanvas ccc3("smintcbneg","smintcbneg",1);
    ccc3.cd();ccc3.SetLogy();ccc3.SetRightMargin(0.25);ccc3.SetTickx();ccc3.SetTicky();
    RooPlot * pp3 = MWW.frame(900,3500);
    w.data("cwww0ccw0cb_60")->plotOn(pp3,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(-60);
    model1.plotOn(pp3,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb_30")->plotOn(pp3,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(-30);
    model1.plotOn(pp3,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb_10")->plotOn(pp3,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(-10);
    model1.plotOn(pp3,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb0")->plotOn(pp3,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));

    cb.setVal(-20);
    for(int i=0; i<10; i++)
    {
        cb.setVal(6*i*(-1));
        if(i==0)
            model1.plotOn(pp3,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
        else
            model1.plotOn(pp3,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
        cout<<i<<" : "<<cb.getVal()<<", "<<scaleshape_cb.getVal()*sum_SM<<endl;
    }
    pp3->GetYaxis()->SetRangeUser(0.002,5);
    pp3->GetYaxis()->SetTitle("arb. units");
    pp3->GetXaxis()->SetTitle("M_{WW} (GeV)");
    pp3->SetTitle("");
    TLegend tt3(0.75,0.5,0.99,0.9);
    tt3.SetFillColor(0);
    tt3.SetFillStyle(0);
    tt3.SetBorderSize(0);
    tt3.SetHeader(ch);
    tt3.AddEntry(pp3->getObject(0),"c_{B}/#Lambda^{2}=-60 TeV^{-2}","le");
    tt3.AddEntry(pp3->getObject(2),"c_{B}/#Lambda^{2}=-30 TeV^{-2}","le");
    tt3.AddEntry(pp3->getObject(4),"c_{B}/#Lambda^{2}=-10 TeV^{-2}","le");
    tt3.AddEntry(pp3->getObject(6),"c_{B}/#Lambda^{2}=0 TeV^{-2}","le");
    tt3.AddEntry(pp3->getObject(8),"c_{B}/#Lambda^{2}=-6,...,-54 TeV^{-2}","l");
    pp3->Draw();
    tt3.Draw();
    ccc3.Draw();ccc3.Update();
    

    TCanvas ccc4("smintcbpos","smintcbpos",1);
    ccc4.cd();ccc4.SetLogy();ccc4.SetRightMargin(0.25);ccc4.SetTickx();ccc4.SetTicky();
    RooPlot * pp4 = MWW.frame(900,3500);
    w.data("cwww0ccw0cb60")->plotOn(pp4,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(60);
    model1.plotOn(pp4,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb30")->plotOn(pp4,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(30);
    model1.plotOn(pp4,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb10")->plotOn(pp4,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
    cwww.setVal(0);ccw.setVal(0);cb.setVal(10);
    model1.plotOn(pp4,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
    w.data("cwww0ccw0cb0")->plotOn(pp4,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));

    cb.setVal(20);
    for(int i=0; i<10; i++)
    {
        cb.setVal(6*i);
        if(i==0)
            model1.plotOn(pp4,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
        else
            model1.plotOn(pp4,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
        cout<<i<<" : "<<cb.getVal()<<", "<<scaleshape_cb.getVal()*sum_SM<<endl;
    }
    pp4->GetYaxis()->SetRangeUser(0.002,5);
    pp4->GetYaxis()->SetTitle("arb. units");
    pp4->GetXaxis()->SetTitle("M_{WW} (GeV)");
    pp4->SetTitle("");
    TLegend tt4(0.75,0.5,0.99,0.9);
    tt4.SetFillColor(0);
    tt4.SetFillStyle(0);
    tt4.SetBorderSize(0);
    tt4.SetHeader(ch);
    tt4.AddEntry(pp4->getObject(0),"c_{B}/#Lambda^{2}=60 TeV^{-2}","le");
    tt4.AddEntry(pp4->getObject(2),"c_{B}/#Lambda^{2}=30 TeV^{-2}","le");
    tt4.AddEntry(pp4->getObject(4),"c_{B}/#Lambda^{2}=10 TeV^{-2}","le");
    tt4.AddEntry(pp4->getObject(6),"c_{B}/#Lambda^{2}=0 TeV^{-2}","le");
    tt4.AddEntry(pp4->getObject(8),"c_{B}/#Lambda^{2}=6,...,54 TeV^{-2}","l");
    pp4->Draw();
    tt4.Draw();
    ccc4.Draw();ccc4.Update();
*/

    for(int i = 1; i<8; i++)
    {
        std::cout<<("a"+to_string(i)).c_str()<<": "<<w2.var(("a"+to_string(i)).c_str())->getVal()<<std::endl;
        if(i>2 and i<5)
            std::cout<<("a"+to_string(i)+to_string(i)).c_str()<<": "<<w2.var(("a"+to_string(i)+to_string(i)).c_str())->getVal()<<std::endl;
    }


    cwww.setVal(12); ccw.setVal(20); cb.setVal(60);
    std::cout<<"all params max: (negative signs compensate for negative fit functions)"<<std::endl;
    std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
    std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
    std::cout<<"N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
    std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;

    cwww.setVal(1.2); ccw.setVal(2); cb.setVal(6);
    std::cout<<"all params 1/10:"<<std::endl;
    std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
    std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
    std::cout<<"N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
    std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;

    exit(0);

}



