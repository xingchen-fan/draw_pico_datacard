#!/usr/bin/env python  
#script to generate efficiency plots and tables from datacards
import ROOT
import os

STR_CHI1PM = '#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[-0.2]{#scale[0.99]{^{#pm}}}}#kern[-1.3]{#scale[0.99]{_{1}}}'
STR_CHII = '#lower[-0.12]{#tilde{#chi}}#kern[+0.1]{#lower[0.2]{#scale[0.99]{^{0,#kern[+0.2]{#lower[-0.15]{#pm}}}}}}#kern[-4]{#scale[0.99]{_{i}}}'
STR_CHIJ = '#lower[-0.12]{#tilde{#chi}}#kern[+0.1]{#lower[0.2]{#scale[0.99]{^{0,#kern[+0.2]{#lower[0.15]{#mp}}}}}}#kern[-4]{#scale[0.99]{_{j}}}'
STR_CHI10 = '#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{1}}}'
STR_CHI20 = '#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{2}}}'
STR_CHI30 = '#lower[-0.12]{#tilde{#chi}}#kern[+0.2]{#lower[0.2]{#scale[0.99]{^{0}}}}#kern[-1.3]{#scale[0.99]{_{3}}}'
STR_XSOFT = 'X#lower[-0.2]{#scale[0.85]{_{soft}}}'
STR_MASS_ = 'm#kern[0.1]{#lower[-0.12]{_{'

def get_cn_higgsino_xsec(hig_mass):
  """function that returns CN Higgsino cross section times B(H->bb) in pb 
  given a Higgsino mass"""
  xsec = 0
  if (hig_mass == 125):
    xsec = 0.5824*0.5824*7.6022
  elif (hig_mass == 127):
    xsec = 0.5824*0.5824*7.6022
  elif (hig_mass == 150):
    xsec = 0.5824*0.5824*3.83231
  elif (hig_mass == 175):
    xsec = 0.5824*0.5824*2.26794
  elif (hig_mass == 200):
    xsec = 0.5824*0.5824*1.33562
  elif (hig_mass == 225):
    xsec = 0.5824*0.5824*0.860597
  elif (hig_mass == 250):
    xsec = 0.5824*0.5824*0.577314
  elif (hig_mass == 275):
    xsec = 0.5824*0.5824*0.400107
  elif (hig_mass == 300):
    xsec = 0.5824*0.5824*0.284855
  elif (hig_mass == 325):
    xsec = 0.5824*0.5824*0.20736
  elif (hig_mass == 350):
    xsec = 0.5824*0.5824*0.153841
  elif (hig_mass == 375):
    xsec = 0.5824*0.5824*0.116006
  elif (hig_mass == 400):
    xsec = 0.5824*0.5824*0.0887325
  elif (hig_mass == 425):
    xsec = 0.5824*0.5824*0.0686963
  elif (hig_mass == 450):
    xsec = 0.5824*0.5824*0.0537702
  elif (hig_mass == 475):
    xsec = 0.5824*0.5824*0.0424699
  elif (hig_mass == 500):
    xsec = 0.5824*0.5824*0.0338387
  elif (hig_mass == 525):
    xsec = 0.5824*0.5824*0.0271867
  elif (hig_mass == 550):
    xsec = 0.5824*0.5824*0.0219868
  elif (hig_mass == 575):
    xsec = 0.5824*0.5824*0.0179062
  elif (hig_mass == 600):
    xsec = 0.5824*0.5824*0.0146677
  elif (hig_mass == 625):
    xsec = 0.5824*0.5824*0.012062
  elif (hig_mass == 650):
    xsec = 0.5824*0.5824*0.00996406
  elif (hig_mass == 675):
    xsec = 0.5824*0.5824*0.00828246
  elif (hig_mass == 700):
    xsec = 0.5824*0.5824*0.00689981
  elif (hig_mass == 725):
    xsec = 0.5824*0.5824*0.00578355
  elif (hig_mass == 750):
    xsec = 0.5824*0.5824*0.0048731
  elif (hig_mass == 775):
    xsec = 0.5824*0.5824*0.00409781
  elif (hig_mass == 800):
    xsec = 0.5824*0.5824*0.00346143
  elif (hig_mass == 825):
    xsec = 0.5824*0.5824*0.0029337
  elif (hig_mass == 850):
    xsec = 0.5824*0.5824*0.0024923
  elif (hig_mass == 875):
    xsec = 0.5824*0.5824*0.00213679
  elif (hig_mass == 900):
    xsec = 0.5824*0.5824*0.00180616
  elif (hig_mass == 925):
    xsec = 0.5824*0.5824*0.00155453
  elif (hig_mass == 950):
    xsec = 0.5824*0.5824*0.00132692
  elif (hig_mass == 975):
    xsec = 0.5824*0.5824*0.00112975
  elif (hig_mass == 1000):
    xsec = 0.5824*0.5824*0.000968853
  elif (hig_mass == 1025):
    xsec = 0.5824*0.5824*0.000840602
  elif (hig_mass == 1050):
    xsec = 0.5824*0.5824*0.000731306
  elif (hig_mass == 1075):
    xsec = 0.5824*0.5824*0.000627083
  elif (hig_mass == 1100):
    xsec = 0.5824*0.5824*0.000538005
  elif (hig_mass == 1125):
    xsec = 0.5824*0.5824*0.00046747
  elif (hig_mass == 1150):
    xsec = 0.5824*0.5824*0.000405108
  elif (hig_mass == 1175):
    xsec = 0.5824*0.5824*0.000348261
  elif (hig_mass == 1200):
    xsec = 0.5824*0.5824*0.000299347
  elif (hig_mass == 1225):
    xsec = 0.5824*0.5824*0.000265935
  elif (hig_mass == 1250):
    xsec = 0.5824*0.5824*0.000240471
  elif (hig_mass == 1275):
    xsec = 0.5824*0.5824*0.000190411
  elif (hig_mass == 1300):
    xsec = 0.5824*0.5824*0.000160765
  elif (hig_mass == 1325):
    xsec = 0.5824*0.5824*0.000136272
  elif (hig_mass == 1350):
    xsec = 0.5824*0.5824*0.000111174
  elif (hig_mass == 1375):
    xsec = 0.5824*0.5824*9.74728e-05
  elif (hig_mass == 1400):
    xsec = 0.5824*0.5824*7.80263e-05
  elif (hig_mass == 1425):
    xsec = 0.5824*0.5824*6.96843e-05
  elif (hig_mass == 1450):
    xsec = 0.5824*0.5824*6.96962e-05
  elif (hig_mass == 1475):
    xsec = 0.5824*0.5824*4.98006e-05
  return xsec

def get_n1n2_higgsino_xsec(hig_mass):
  """function that returns N1N2 Higgsino cross section times B(H->bb) in pb 
  given a Higgsino mass"""
  xsec = 0
  if (hig_mass ==  125):
    xsec = 0.5824*0.5824*1.44725
  elif (hig_mass ==  127):
    xsec = 0.5824*0.5824*1.44725
  elif (hig_mass ==  150):
    xsec = 0.5824*0.5824*0.71514
  elif (hig_mass ==  175):
    xsec = 0.5824*0.5824*0.419059
  elif (hig_mass ==  200):
    xsec = 0.5824*0.5824*0.244213
  elif (hig_mass ==  225):
    xsec = 0.5824*0.5824*0.156286
  elif (hig_mass ==  250):
    xsec = 0.5824*0.5824*0.104252
  elif (hig_mass ==  275):
    xsec = 0.5824*0.5824*0.0719125
  elif (hig_mass ==  300):
    xsec = 0.5824*0.5824*0.0509994
  elif (hig_mass ==  325):
    xsec = 0.5824*0.5824*0.0369715
  elif (hig_mass ==  350):
    xsec = 0.5824*0.5824*0.0273286
  elif (hig_mass ==  375):
    xsec = 0.5824*0.5824*0.0205429
  elif (hig_mass ==  400):
    xsec = 0.5824*0.5824*0.0156691
  elif (hig_mass ==  425):
    xsec = 0.5824*0.5824*0.0120965
  elif (hig_mass ==  450):
    xsec = 0.5824*0.5824*0.00944017
  elif (hig_mass ==  475):
    xsec = 0.5824*0.5824*0.00743587
  elif (hig_mass ==  500):
    xsec = 0.5824*0.5824*0.00590757
  elif (hig_mass ==  525):
    xsec = 0.5824*0.5824*0.00473235
  elif (hig_mass ==  550):
    xsec = 0.5824*0.5824*0.0038167
  elif (hig_mass ==  575):
    xsec = 0.5824*0.5824*0.00309847
  elif (hig_mass ==  600):
    xsec = 0.5824*0.5824*0.00253015
  elif (hig_mass ==  625):
    xsec = 0.5824*0.5824*0.00207755
  elif (hig_mass ==  650):
    xsec = 0.5824*0.5824*0.00171418
  elif (hig_mass ==  675):
    xsec = 0.5824*0.5824*0.0014199
  elif (hig_mass ==  700):
    xsec = 0.5824*0.5824*0.00118113
  elif (hig_mass ==  725):
    xsec = 0.5824*0.5824*0.00098639
  elif (hig_mass ==  750):
    xsec = 0.5824*0.5824*0.000826366
  elif (hig_mass ==  775):
    xsec = 0.5824*0.5824*0.000694985
  elif (hig_mass ==  800):
    xsec = 0.5824*0.5824*0.000586211
  elif (hig_mass ==  825):
    xsec = 0.5824*0.5824*0.000495914
  elif (hig_mass ==  850):
    xsec = 0.5824*0.5824*0.000420556
  elif (hig_mass ==  875):
    xsec = 0.5824*0.5824*0.000361029
  elif (hig_mass ==  900):
    xsec = 0.5824*0.5824*0.000305935
  elif (hig_mass ==  925):
    xsec = 0.5824*0.5824*0.000262621
  elif (hig_mass ==  950):
    xsec = 0.5824*0.5824*0.00022285
  elif (hig_mass ==  975):
    xsec = 0.5824*0.5824*0.0001909
  elif (hig_mass ==  1000):
    xsec = 0.5824*0.5824*0.00016428
  elif (hig_mass ==  1025):
    xsec = 0.5824*0.5824*0.00014139
  elif (hig_mass ==  1050):
    xsec = 0.5824*0.5824*0.000121865
  elif (hig_mass ==  1075):
    xsec = 0.5824*0.5824*0.000105913
  elif (hig_mass ==  1100):
    xsec = 0.5824*0.5824*9.12469e-05
  elif (hig_mass ==  1125):
    xsec = 0.5824*0.5824*7.93058e-05
  elif (hig_mass ==  1150):
    xsec = 0.5824*0.5824*6.84561e-05
  elif (hig_mass ==  1175):
    xsec = 0.5824*0.5824*5.93602e-05
  elif (hig_mass ==  1200):
    xsec = 0.5824*0.5824*5.16263e-05
  elif (hig_mass ==  1225):
    xsec = 0.5824*0.5824*4.4906e-05
  elif (hig_mass ==  1250):
    xsec = 0.5824*0.5824*3.91587e-05
  elif (hig_mass ==  1275):
    xsec = 0.5824*0.5824*3.43135e-05
  elif (hig_mass ==  1300):
    xsec = 0.5824*0.5824*2.99353e-05
  elif (hig_mass ==  1325):
    xsec = 0.5824*0.5824*2.62223e-05
  elif (hig_mass ==  1350):
    xsec = 0.5824*0.5824*2.28072e-05
  elif (hig_mass ==  1375):
    xsec = 0.5824*0.5824*2.00393e-05
  elif (hig_mass ==  1400):
    xsec = 0.5824*0.5824*1.75031e-05
  elif (hig_mass ==  1425):
    xsec = 0.5824*0.5824*1.53144e-05
  elif (hig_mass ==  1450):
    xsec = 0.5824*0.5824*1.34572e-05
  elif (hig_mass ==  1475):
    xsec = 0.5824*0.5824*1.17047e-05
  return xsec

def signal_events_in_a(datacard_filename):
  """Get number of signal events in ABCD-A regions for a given resolved datacard"""
  a_region_yield = 0.0
  datacard = open(datacard_filename,'r')
  datacard_lines = datacard.read().split('\n')
  rates_line = datacard_lines[11].split()
  for rate_idx in range(len(rates_line)):
    #rows 9+12n and 11+12n are A region
    if (rate_idx % 12 == 9) or (rate_idx % 12 == 11):
      a_region_yield += float(rates_line[rate_idx])
  datacard.close()
  return a_region_yield

def signal_events_in_a_resboo(datacard_filename):
  """Get number of signal events in ABCD-A regions for a given combined datacard"""
  a_region_yield = 0.0
  datacard = open(datacard_filename,'r')
  datacard_lines = datacard.read().split('\n')
  rates_line = []
  for line in datacard_lines:
    if (line[0:4]=='rate'):
      rates_line = line.split()
  if (len(rates_line)==0):
    print('no rate line found in datacard '+datacard_filename)
  for rate_idx in range(len(rates_line)):
    #rows 1, 3, 5, 7, 9, and 11 are boosted A regions
    if (rate_idx < 12) and (rate_idx % 2 == 1):
      a_region_yield += float(rates_line[rate_idx])
    #rows 35+12n and 37+12n are resolved A region
    if (rate_idx >= 27) and ((rate_idx % 12 == 1) or (rate_idx % 12 == 11)):
      a_region_yield += float(rates_line[rate_idx])
  datacard.close()
  return a_region_yield

def signal_events_in_a_boo(datacard_filename):
  """Get number of signal events in boosted only ABCD-A regions for a given combined datacard"""
  a_region_yield = 0.0
  datacard = open(datacard_filename,'r')
  datacard_lines = datacard.read().split('\n')
  rates_line = []
  for line in datacard_lines:
    if (line[0:4]=='rate'):
      rates_line = line.split()
  if (len(rates_line)==0):
    print('no rate line found in datacard '+datacard_filename)
  for rate_idx in range(len(rates_line)):
    #rows 1, 3, 5, 7, 9, and 11 are boosted A regions
    if (rate_idx < 12) and (rate_idx % 2 == 1):
      a_region_yield += float(rates_line[rate_idx])
    #rows 35+12n and 37+12n are resolved A region
    #if (rate_idx >= 27) and ((rate_idx % 12 == 1) or (rate_idx % 12 == 11)):
    #  a_region_yield += float(rates_line[rate_idx])
  datacard.close()
  return a_region_yield

def higgsino_datacard_filename(higgsino_mass, lsp_mass):
  """Returns resolved datacard name given Higgsino and LSP mass"""
  higgsino_mass_str = str(higgsino_mass)
  if (higgsino_mass == 125):
    higgsino_mass_str = '127'
  return 'datacard-TChiHH_mChi-'+higgsino_mass_str+'_mLSP-'+str(lsp_mass)+'_Tune_2016,2017,2018_priority1_resolved.txt'

def higgsino_datacard_filename_resboo(higgsino_mass, lsp_mass, dimension='2D'):
  """Returns combined datacard name given Higgsino mass"""
  higgsino_mass_str = str(higgsino_mass)
  if (higgsino_mass == 125):
    higgsino_mass_str = '127'
  lsp_mass_str = str(lsp_mass)
  if (lsp_mass == 0):
    lsp_mass_str = '1'
  return dimension+'TChiHH'+higgsino_mass_str+'_LSP'+lsp_mass_str+'_Data_Combo.txt'

def decorate_graph(graph):
  """Adds CMS decorations to TGraph and saves it"""
  can = ROOT.TCanvas('can','can',1000,1000)
  can.cd()
  can.SetLogy()
  ROOT.gStyle.SetOptStat(0)
  can.SetMargin(0.12, 0.12, 0.12, 0.12);
  can.SetFillStyle(4000);
  hdummy = ROOT.TH1D('','',49,75,1300)
  #hdummy.GetYaxis().SetLabelSize(0.04)
  hdummy.GetYaxis().SetRangeUser(5.0e-6,8.0)
  #hdummy.GetYaxis().SetTitle('Resolved+Boosted Signal Efficiency')
  hdummy.GetYaxis().SetTitle('Resolved Signal Efficiency')
  #hdummy.GetYaxis().SetTitle('Boosted-Overlap Signal Efficiency')
  hdummy.GetYaxis().SetTitleSize(0.04)
  hdummy.GetYaxis().SetTitleOffset(1.4)
  #hdummy.GetYaxis().SetTitleOffset(0.565)
  hdummy.GetXaxis().SetLimits(75.0,1300.0)
  #hdummy.GetXaxis().SetLabelSize(0.04)
  hdummy.GetXaxis().SetTitle('m('+STR_CHI10+') [GeV]')
  hdummy.GetXaxis().SetTitleSize(0.04)
  hdummy.GetXaxis().SetTitleOffset(1.2)
  hdummy.Draw()
  graph.SetMarkerStyle(8)
  graph.SetMarkerSize(1.0)
  graph.Draw('P')
  lumilabel = ROOT.TLatex()
  lumilabel.SetTextSize(0.038)
  lumilabel.SetNDC(True)
  lumilabel.SetTextAlign(31)
  lumilabel.DrawLatex(1.0-0.12, 1.0-0.11,"#font[42]{137 fb^{-1} (13 TeV)}")
  labels = ROOT.TLatex()
  labels.SetNDC(False)
  labels.SetTextSize(0.04)
  labels.SetTextAlign(11)
  labels.DrawLatex(125, 2.5, '#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}')
  labels.SetTextSize(0.038)
  labels.DrawLatex(125, 0.7, '#font[42]{pp #rightarrow '+STR_CHII+'#kern[0.7]{'+STR_CHIJ+'}  #rightarrow '+STR_CHI10+'#kern[0.3]{'+STR_CHI10+'} + '+STR_XSOFT+'#rightarrow HH#tilde{G}#tilde{G} + '+STR_XSOFT+'}')
  labels.DrawLatex(125, 0.22, '#font[42]{'+STR_MASS_+STR_CHI20+'}}} = '+STR_MASS_+STR_CHI1PM+'}}} = '+STR_MASS_+STR_CHI10+'}}}#kern[0.5]{,} '+STR_MASS_+'#tilde{G}}}} = 1 GeV}')
  can.Print('plots/sigeff__onedim.pdf')

def decorate_graph2d(graph):
  """Adds CMS decorations to TGraph and saves it"""
  can = ROOT.TCanvas('can2','can2',1000,1000)
  can.cd()
  can.SetLogz()
  ROOT.gStyle.SetOptStat(0)
  can.SetMargin(0.13, 0.17, 0.15, 0.15);
  #can.SetFillStyle(4000);
  hdummy = ROOT.TH2D('hdummy2','',26,150,800,43,0,675)
  hdummy.GetZaxis().SetTitle('Resolved Signal Efficiency')
  #NOTE!! drawing multiple COLZ plots only works if the z-axis range is the same
  hdummy.GetZaxis().SetRangeUser(5.0e-6,1.0)
  hdummy.GetYaxis().SetTitle('m('+STR_CHI10+') [GeV]')
  hdummy.GetYaxis().SetTitleSize(0.04)
  hdummy.GetYaxis().SetTitleOffset(1.4)
  hdummy.GetXaxis().SetTitle('m('+STR_CHI20+') [GeV]')
  hdummy.GetXaxis().SetTitleSize(0.04)
  hdummy.GetXaxis().SetTitleOffset(1.2)
  hdummy.Draw('COLZ')
  graph.SetTitle('')
  graph.SetNpx(500)
  graph.SetNpy(500)
  graph.GetXaxis().SetLimits(150.0,800.0)
  graph.Draw('SAME COLZ')
  graph.GetHistogram().GetYaxis().SetTitle('m('+STR_CHI10+') [GeV]')
  graph.GetHistogram().GetYaxis().SetTitleSize(0.04)
  graph.GetHistogram().GetYaxis().SetTitleOffset(1.4)
  graph.GetHistogram().GetXaxis().SetTitle('m('+STR_CHI20+') [GeV]')
  graph.GetHistogram().GetXaxis().SetTitleSize(0.04)
  graph.GetHistogram().GetXaxis().SetTitleOffset(1.2)
  graph.GetHistogram().GetZaxis().SetTitle('Resolved Signal Efficiency')
  graph.GetHistogram().GetZaxis().SetTitleOffset(1.65)
  graph.GetHistogram().GetZaxis().SetRangeUser(5.0e-6,1.0)
  graph.GetHistogram().SetContour(500)
  graph.Draw('SAME COLZ')
  ROOT.gPad.RedrawAxis()
  lumilabel = ROOT.TLatex()
  lumilabel.SetTextSize(0.038)
  lumilabel.SetNDC(True)
  lumilabel.SetTextAlign(31)
  lumilabel.DrawLatex(1.0-0.17, 1.0-0.14,"#font[42]{137 fb^{-1} (13 TeV)}")
  labels = ROOT.TLatex()
  labels.SetNDC(False)
  labels.SetTextSize(0.04)
  labels.SetTextAlign(11)
  labels.DrawLatex(180, 625, '#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}')
  labels.SetTextSize(0.038)
  labels.DrawLatex(180, 575, '#font[42]{pp #rightarrow #kern[0.3]{'+STR_CHI20+'}#kern[0.3]{'+STR_CHI30+'}  #rightarrow HH#kern[0.3]{'+STR_CHI10+'}#kern[0.3]{'+STR_CHI10+'}}')
  labels.DrawLatex(180, 525, '#font[42]{'+STR_MASS_+STR_CHI20+'}}} = '+STR_MASS_+STR_CHI30+'}}}}')
  can.Print('plots/sigeff__twodim.pdf')

def make_efficiency_table(datacard_dir, higgsino_mass, lsp_mass):
  """Make latex table of yields from datacard"""
  datacard_filename = datacard_dir+'/'+higgsino_datacard_filename(higgsino_mass,lsp_mass)
  datacard = open(datacard_filename,'r')
  datacard_lines = datacard.read().split('\n')
  rates_line = datacard_lines[11].split()
  datacard.close()
  total_events = get_cn_higgsino_xsec(higgsino_mass)*137.35167*1000.0
  table_file = open('tables/signal_eff_'+str(higgsino_mass)+'.tex','w')
  table_file.write('\\documentclass[10pt,oneside]{report}\n')
  table_file.write('\\usepackage{graphicx,xspace,amssymb,amsmath,colordvi,colortbl,verbatim,multicol}\n')
  table_file.write('\\usepackage{multirow, rotating}\n\n')
  table_file.write('\\usepackage[active,tightpage]{preview}\n\n')
  table_file.write('\\usepackage{siunitx}\n')
  table_file.write('\\sisetup{round-mode = figures, round-precision=2}\n\n')
  table_file.write('\\renewcommand{\\arraystretch}{1.1}\n\n')
  table_file.write('\\begin{document}\n')
  table_file.write('\\begin{preview}\n')
  table_file.write('\\begin{tabular}{ cc c cccccc }\n')
  table_file.write(' & & \\multicolumn{7}{c}{\\large{Signal Efficiency[\%]}} \\\\\n')
  table_file.write('\\drmax & \\ptmiss [GeV] & Plane total& $N_{\\mathrm{b}}=2$ CSB & $N_{\\mathrm{b}}=2$ CSR & $N_{\\mathrm{b}}=3$ SB & $N_{\\mathrm{b}}=3$ SR & $N_{\\mathrm{b}}=4$ SB&  $N_{\\mathrm{b}}=4$ SR\\\\ \\hline\n')
  #write integrated yields
  table_file.write(' \\multicolumn{2}{c}{Integrated} ')
  integrated_row_yields = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  for rate_idx in range(len(rates_line)):
    if (rate_idx % 12 == 1):
      integrated_row_yields[1] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
    elif (rate_idx % 12 == 7):
      integrated_row_yields[2] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
    elif (rate_idx % 12 == 3):
      integrated_row_yields[3] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
    elif (rate_idx % 12 == 9):
      integrated_row_yields[4] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
    elif (rate_idx % 12 == 5):
      integrated_row_yields[5] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
    elif (rate_idx % 12 == 11):
      integrated_row_yields[6] += float(rates_line[rate_idx])
      integrated_row_yields[0] += float(rates_line[rate_idx])
  for abcd_cat_events in integrated_row_yields:
    table_file.write(' & '+str(round(abcd_cat_events/total_events*100.0,2)))
  table_file.write('\\\\\n')
  for drmax_idx in range(2):
    for met_idx in range(4):
      #write drmax column
      if (met_idx == 0):
        if (drmax_idx == 0):
          table_file.write('\\multirow{4}{*}{1.1--2.2}& ')
        else:
          table_file.write('\\multirow{4}{*}{$<$ 1.1}& ')
      else:
        table_file.write(' & ')
      #write met column
      if (met_idx == 0):
        table_file.write(' 150--200& ')
      elif (met_idx == 1):
        table_file.write(' 200--300& ')
      elif (met_idx == 2):
        table_file.write(' 300--400& ')
      else:
        table_file.write(' $>$ 400& ')
      bin_offset = 0
      if (drmax_idx == 0):
        bin_offset = 12
      bin_offset += 1
      bin_offset += 24*met_idx
      #write event numbers
      baseline_events = 0.0
      abcd_events = []
      for abcd_offset in [0,6,2,8,4,10]:
        baseline_events += float(rates_line[bin_offset+abcd_offset])
        abcd_events.append(float(rates_line[bin_offset+abcd_offset]))
      table_file.write(str(round(baseline_events/total_events*100.0,2)))
      for abcd_cat_events in abcd_events:
        table_file.write(' & '+str(round(abcd_cat_events/total_events*100.0,2)))
      table_file.write('\\\\\n')
  table_file.write('\\end{tabular}\n')
  table_file.write('\\end{preview}\n')
  table_file.write('\\end{document}')
  table_file.close()

if __name__ == "__main__":
  onedim_datacard_dir = 'datacards/tchihh_onedim_lumifix'
  #onedim_datacard_dir = 'datacards/combined_1DTChiHH'
  #onedim_datacard_dir = 'datacards/combined_2DTChiHH'
  twodim_datacard_dir = 'datacards/tchihh_twodim_lumifix'

  #make 1D plot
  #onedim_tchihh_efficiency = ROOT.TGraph(44)
  onedim_tchihh_efficiency = ROOT.TGraph(27)
  higgsino_mass = 150
  lsp_mass = 0
  point_idx = 0
  while higgsino_mass < 1200:
    #don't include non H->bb events in the denominator
    total_events = get_cn_higgsino_xsec(higgsino_mass)*137.35167*1000.0
    sr_events = signal_events_in_a(onedim_datacard_dir+'/'+higgsino_datacard_filename(higgsino_mass,lsp_mass))
    #sr_events = signal_events_in_a(datacard_dir+'/'+higgsino_datacard_filename(higgsino_mass,lsp_mass))
    onedim_tchihh_efficiency.SetPoint(point_idx, float(higgsino_mass), sr_events/total_events)
    #print('Adding point to graph, x='+str(higgsino_mass)+',y='+str(sr_events/total_events))
    higgsino_mass += 25
    point_idx += 1
  decorate_graph(onedim_tchihh_efficiency)

  #make 2D plot
  #twodim_tchihh_efficiency = ROOT.TH2D('twodim_tchihh_efficiency','',26,150,800,43,0,675)
  num_2d_points = len([name for name in os.listdir(twodim_datacard_dir) if os.path.isfile(twodim_datacard_dir+'/'+name)])
  twodim_tchihh_efficiency = ROOT.TGraph2D(num_2d_points)
  higgsino_mass = 150
  xbin = 1
  point_idx = 0
  while higgsino_mass < 1225:
    lsp_mass = 0
    ybin = 1
    #don't include non H->bb events in the denominator
    total_events = get_n1n2_higgsino_xsec(higgsino_mass)*137.35167*1000.0
    while lsp_mass < higgsino_mass:
      datacard_filename = twodim_datacard_dir+'/'+higgsino_datacard_filename(higgsino_mass,lsp_mass)
      if os.path.isfile(datacard_filename):
        sr_events = signal_events_in_a(datacard_filename)
        twodim_tchihh_efficiency.SetPoint(point_idx, float(higgsino_mass), float(lsp_mass), sr_events/total_events)
        #print('setting '+str(higgsino_mass)+', '+str(lsp_mass)+': '+str(sr_events/total_events))
        #twodim_tchihh_efficiency.SetBinContent(xbin,ybin,sr_events/total_events)
        point_idx += 1
      lsp_mass += 25
      ybin += 1
    higgsino_mass += 25
    xbin += 1
  decorate_graph2d(twodim_tchihh_efficiency)
  
  #make table
  make_efficiency_table(onedim_datacard_dir, 450, 0)
