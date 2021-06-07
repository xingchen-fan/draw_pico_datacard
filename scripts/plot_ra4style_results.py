#! /usr/bin/env python
import os
import math
from ROOT import *

#do_sig = True
#compile_table = True
#
#
# --- Official plot

plot_pulls = True

tag1 = '_nor4'
tag1_lbl = 'Pre-fit'
tag1_color = kPink+2

tag2 = '_r4' 
tag2_lbl = 'Post-fit'
tag2_color = kAzure+1

#tag_sig_c = '_r4_1900'
#
## --- Effect of systematics
## tag1 = '_nor4'
## tag1_lbl = 'R1--R3 fit'
#tag1_color = kBlack
#
## tag2 = '_nor4_nosys' 
## tag2_lbl = 'R1--R3 fit (no syst)'
#tag2_color = kGray+1
#
##         Bin names
##-------------------------------------
#bins = []
#for imj in ['lmj','hmj']:
#    for imet in ['lmet','mmet','hmet']:
#        for ir in ['r1','r2','r3','r4']:
#            if ir=='r1' or ir=='r3':
#                if imj=='lmj':
#                    bins.append('_'.join([ir,imet]))
#                else:
#                    continue
#            else:
#                for inb in ['lnb','mnb','hnb']:
#                    for inj in ['lnj','hnj']:
#                        bins.append('_'.join([ir,imet,inb, inj,imj]))
#
#nbins = len(bins)
## in the combine output they are ordered alphabetically
#bins_alpha = sorted(bins)

#         Reading the data
#-------------------------------------
#bkg_pre, ebkg_pre, pull_pre = [None]*nbins, [None]*nbins, [None]*nbins
#data, edata_up, edata_dn = [None]*nbins, [None]*nbins, [None]*nbins
#sig_nc, esig_nc = [None]*nbins, [None]*nbins
#
#file_pre = TFile("root/fitDiagnostics"+tag1+".root")
#bkg_in = file_pre.Get('shapes_fit_b/total_background')
#data_in = file_pre.Get('shapes_fit_b/total_data')
#sig_nc_in = file_pre.Get('shapes_prefit/total_signal')
#for ibin in range(0,nbins):
#    ibin_alpha = bins_alpha.index(bins[ibin])
#    x, y  = Double(0), Double(0)
#    data_in.GetPoint(ibin_alpha, x, y)
#    data[ibin] = y
#    edata_up[ibin] = data_in.GetErrorYhigh(ibin_alpha)
#    edata_dn[ibin] = data_in.GetErrorYlow(ibin_alpha)
#    bkg_pre[ibin] = bkg_in.GetBinContent(ibin_alpha+1)
#    ebkg_pre[ibin] = bkg_in.GetBinError(ibin_alpha+1)
#    pull_pre[ibin] = (data[ibin]-bkg_pre[ibin])/math.sqrt(bkg_pre[ibin]+ebkg_pre[ibin]*ebkg_pre[ibin])
#    sig_nc[ibin] = sig_nc_in.GetBinContent(ibin_alpha+1)
#    esig_nc[ibin] = sig_nc_in.GetBinError(ibin_alpha+1)
#file_pre.Close()
#
#bkg_post, ebkg_post, pull_post = [None]*nbins, [None]*nbins, [None]*nbins
#file_post = TFile("root/fitDiagnostics"+tag2+".root")
#bkg_in = file_post.Get('shapes_fit_b/total_background')
#for ibin in range(0,nbins):
#    ibin_alpha = bins_alpha.index(bins[ibin])
#    bkg_post[ibin] = bkg_in.GetBinContent(ibin_alpha+1)
#    ebkg_post[ibin] = bkg_in.GetBinError(ibin_alpha+1)
#    pull_post[ibin] = (data[ibin]-bkg_post[ibin])/math.sqrt(bkg_post[ibin]+ebkg_post[ibin]*ebkg_post[ibin])
#file_post.Close()
#
#sig_c, esig_c = [None]*nbins, [None]*nbins
#file_sig_c = TFile("root/fitDiagnostics"+tag_sig_c+".root")
#sig_c_in = file_sig_c.Get('shapes_prefit/total_signal')
#for ibin in range(0,nbins):
#    ibin_alpha = bins_alpha.index(bins[ibin])
#    sig_c[ibin] = sig_c_in.GetBinContent(ibin_alpha+1)
#    esig_c[ibin] = sig_c_in.GetBinError(ibin_alpha+1)
#file_sig_c.Close()
#
##         Making table
##-------------------------------------
#tab = []
#tab.append(open("tables/table_lmj"+tag1+"_vs"+tag2+".tex","w"))
#tab.append(open("tables/table_hmj"+tag1+"_vs"+tag2+".tex","w"))
#ncols = 4
#if (do_sig): ncols +=2
#
#tab_head = "\\begin{tabular}[tbp!]{ l cc"
#for i in range(ncols-3): tab_head += " r"
#tab_head += "} \n\\hline\\hline\n"
#tab_head += "${\\cal L}=137$ fb$^{-1}$ &"
#if do_sig:
#    tab_head += " T1tttt(2100,100) & T1tttt(1900,1250) &"
#tab_head += tag1_lbl+" & "+tag2_lbl+" & Obs. \\\\ \\hline\n"
## tab_head += tag1_lbl+" & Pull & "+tag2_lbl+" & Pull & Obs. \\\\ \\hline\n"
#for i in range(2): tab[i].write(tab_head)
#
#save_rows = { 'r1_lmet':'', 'r3_lmet':'', 'r1_mmet':'', 'r3_mmet':'', 'r1_hmet':'', 'r3_hmet':''}
#irow = 0
#for ibin in range(nbins):
#    tmp = bins[ibin].split("_")
#    ireg = tmp[0].replace("r","R")
#    imet = tmp[1].replace("lmet","$200<p_{\\rm T}^{\\text{miss}}\\leq350$ GeV")
#    imet = imet.replace("mmet","$350<p_{\\rm T}^{\\text{miss}}\\leq500$ GeV")
#    imet = imet.replace("hmet","$p_{\\rm T}^{\\text{miss}}> 500$ GeV")
#    inb,inj = '',''
#    if ireg=='R2' or ireg=='R4':
#        inb = tmp[2].replace("lnb","$N_{b} = 1$").replace("mnb","$N_{b} = 2$").replace("hnb","$N_{b} \\geq 3$")
#        if "hmet" in bins[ibin]:
#            inj = tmp[3].replace("lnj"," $6\\leq N_{jets} \\leq 7$").replace("hnj","$N_{jets} \\geq 8$")
#        else:
#            inj = tmp[3].replace("lnj"," $N_{jets} = 7$").replace("hnj","$N_{jets} \\geq 8$")
#    
#    itab = ibin/42
#    if (irow%14==0):
#        tab[itab].write("\\hline\n\\multicolumn{"+str(ncols)+"}{c}{"+imet+"}  \\\\ \\hline\n");
#
#    if ibin>41 and (ibin-42)%6==0: 
#        irow +=1
#        tmp_ = bins[ibin].split("_")
#        insert_bin = tmp_[0].replace("r2","r1").replace("r4","r3")+'_'+tmp_[1]
#        tab[itab].write(save_rows[insert_bin])
#        if '3' in insert_bin:
#            tab[itab].write("\\hline\n")
#
#
#    cols = []
#    if ireg=='R2' or ireg=='R4':
#        cols.append('{0:<30}'.format(ireg+": "+inb+", "+inj))
#    else:
#        cols.append('{0:<30}'.format(ireg))
#    if do_sig:
#        cols.append('{0:>10.1f}'.format(sig_nc[ibin], esig_nc[ibin]))
#        cols.append('{0:>10.1f}'.format(sig_c[ibin], esig_c[ibin]))
#    if (ireg=="R4"):
#        cols.append('{0:>20}'.format('${0:.1f} \\pm {1:.1f}$'.format(bkg_pre[ibin], ebkg_pre[ibin])))
#    else:
#        cols.append('{0:>20}'.format('${0:.0f} \\pm {1:.1f}$'.format(bkg_pre[ibin], ebkg_pre[ibin])))
#    # if ireg=="R4": 
#    #     cols.append('{0:>7.1f}'.format(pull_pre[ibin]))
#    # else: 
#    #     cols.append('')
#    cols.append('{0:>20}'.format('${0:.1f} \\pm {1:.1f}$'.format(bkg_post[ibin], ebkg_post[ibin])))
#    # if ireg=="R4": 
#    #     cols.append('{0:>7.1f}'.format(pull_post[ibin]))
#    # else: 
#    #     cols.append('')
#    cols.append('{0:>10.0f}'.format(data[ibin]))
#    tab[itab].write('&'.join(cols)+'\\\\\n')
#
#    if ibin<42 and ('1' in ireg or '3' in ireg): 
#        save_rows[bins[ibin]] = '&'.join(cols)+'\\\\\n'
#
#    if '3' in ireg: 
#        tab[itab].write("\\hline\n")
#
#    irow += 1
#
#for itab in tab:
#    itab.write("\\hline\\hline\n \\end{tabular}\n")
#    itab.close()
#
#
##         Making tables that can compile standalone
##----------------------------------------------------------
#fulltab = []
#
#for i in range(2):
#    tabname = tab[i].name.split("/")[-1]
#    fulltab.append(open("tables/full"+tabname,"w"))
#    with open("txt/header.tex") as head:
#        for line in head.readlines():
#            fulltab[i].write(line)
#    fulltab[i].write("\\begin{document}\n\\begin{preview}\n")
#    with open(tab[i].name) as body:
#        for line in body.readlines():
#            fulltab[i].write(line)
#    fulltab[i].write("\\end{preview}\n\\end{document}\n")
#    fulltab[i].close()
#    if (compile_table):
#        print "Converting "+fulltab[i].name+" ..."
#        os.system("pdflatex "+fulltab[i].name+" > /dev/null")
#
#if compile_table:
#    for i in range(2): 
#        print "open "+fulltab[i].name.replace(".tex",".pdf")

#currently manually enter data, later we will get it from combine
bkg_pre = [159.0, 91.4, 11.1, 2.8, 50.7, 27.7, 2.7, 2.9, 6.0, 1.97, 0.074, 0.88, 2.55, 0.96, 0.42, 0.69]
ebkg_pre = [12.0, 9.4, 3.4, 1.9, 8.6, 5.3, 1.3, 2.3, 1.6, 0.64, 0.083, 0.91, 0.84, 0.47, 0.42, 0.71]
bkg_post = [150.5, 92.1, 12.5, 2.8, 51.7, 32.6, 3.3, 1.34, 6.9, 2.18, 0.78, 0.52, 2.51, 1.33, 1.15, 0.79]
ebkg_post = [8.2, 6.6, 2.6, 1.2, 5.8, 4.1, 1.2, 0.88, 1.5, 0.63, 0.46, 0.48, 0.69, 0.55, 0.71, 0.58]
data = [141.0, 91.0, 14.0, 3.0, 51.0, 37.0, 4.0, 0.0, 9.0, 2.0, 4.0, 0.0, 1.0, 3.0, 1.0, 1.0]
bkg_pre += [33.8, 7.4, 1.68, 4.0, 0.89, 0.2]
ebkg_pre += [6.4, 1.9, 0.89, 1.4, 0.36, 0.17]
bkg_post += [37.1, 7.4, 1.53, 3.62, 0.75, 0.16]
ebkg_post += [4.4, 1.4, 0.64, 0.96, 0.25, 0.12]
data += [43.0, 7.0, 1.0, 4.0, 0.0, 0.0]
#edata_up = [math.sqrt(ibin) for ibin in data]
#edata_dn = [math.sqrt(ibin) for ibin in data]
edata_up = []
edata_dn = []
pull_pre = []
pull_post = []
#pulls as likelihood ratio significances
pull_pre = [-0.9, 0.0, 0.6, 0.1,   0.0, 1.1, 0.6, -1.9,   0.9, 0.0, 3.5, -1.1,   -1.0, 1.5, 0.6, 0.3,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#pull_post = []
for ibin in range(len(data)):
  #make Poisson errors
  alpha = 1.0 - 0.6827
  l = 0
  if (data[ibin] != 0):
    l = ROOT.Math.gamma_quantile(alpha/2,data[ibin],1.0)
  u =  ROOT.Math.gamma_quantile_c(alpha/2,data[ibin]+1,1.0)
  edata_up.append(u-data[ibin])
  edata_dn.append(data[ibin]-l)
  #make pulls
  w_pre = 1.0/ebkg_pre[ibin]**2
  w_pos = 1.0/ebkg_post[ibin]**2
  e_dat = edata_dn[ibin]
  if data[ibin] < bkg_pre[ibin]:
    e_dat = edata_up[ibin]
  w_dat = 1.0/e_dat**2
  e_com_pre = 1.0/math.sqrt(w_pre+w_dat)
  e_com_pos = 1.0/math.sqrt(w_pos+w_dat)
  l_com_pre = (bkg_pre[ibin]*w_pre+data[ibin]*w_dat)/(w_pre+w_dat)
  l_com_pos = (bkg_post[ibin]*w_pos+data[ibin]*w_dat)/(w_pos+w_dat)
  this_pull_pre = (data[ibin]-l_com_pre)/math.sqrt(e_dat**2-e_com_pre**2)
  this_pull_pos = (data[ibin]-l_com_pos)/math.sqrt(e_dat**2-e_com_pos**2)
  #pull_pre.append(this_pull_pre)
  #pull_post.append(this_pull_pos)
  #naive-pulls
  #pull_pre.append((data[ibin]-bkg_pre[ibin])/math.sqrt(bkg_pre[ibin]+ebkg_pre[ibin]*ebkg_pre[ibin]))
  #pull_post.append((data[ibin]-bkg_post[ibin])/math.sqrt(bkg_post[ibin]+ebkg_post[ibin]*ebkg_post[ibin]))
  #poisson significance
  #pull_pre.append(math.sqrt(2)*TMath.ErfInverse(-1.0+2.0*ROOT.Math.poisson_cdf(int(data[ibin]),bkg_pre[ibin])))
  pull_post.append(math.sqrt(2)*TMath.ErfInverse(-1.0+2.0*ROOT.Math.poisson_cdf(int(data[ibin]),bkg_post[ibin])))

print(pull_pre)
print(pull_post)

nbins = len(data)
nhbins = len(data)

#         Plotting the data
#-------------------------------------
miny, maxy = 0.051, 3000
htopdummy = TH1D("","",nhbins,0,nhbins+1)
grbkg_pre = TGraphErrors(nhbins)
grbkg_post = TGraphErrors(nhbins)
grdata = TGraphAsymmErrors(nhbins)
grpull_pre = TGraph(nhbins)
grpull_post = TGraph(nhbins)
i = 0
print(nbins)
for ibin in range(nbins):
    #if 'r4' not in bins[ibin]: continue
    i +=1
    grbkg_pre.SetPoint(i-1, i,bkg_pre[ibin])
    grbkg_pre.SetPointError(i-1, 0.5, ebkg_pre[ibin])    
    grpull_pre.SetPoint(i-1, i, pull_pre[ibin])

    grbkg_post.SetPoint(i-1, i, bkg_post[ibin])
    grbkg_post.SetPointError(i-1, 0.5, ebkg_post[ibin])
    grpull_post.SetPoint(i-1, i, pull_post[ibin])

    if ((not plot_pulls) and (data[ibin] == 0)):
      grdata.SetPoint(i-1, i, miny)
    else:
      grdata.SetPoint(i-1, i, data[ibin])
    grdata.SetPointEYhigh(i-1, edata_up[ibin])
    grdata.SetPointEYlow(i-1, edata_dn[ibin])
    #grdata.SetPointEXhigh(i, 0.)
    #grdata.SetPointEXlow(i, 0.)


# for ibin in range(nbins):
#     print bkg_pre[ibin], ebkg_pre[ibin]

can = TCanvas('c','c',1000,500)
can.cd()

an_region_size = 0.05
nb_title_size = 0.045
met_title_size = 0.04
met_range_size = 0.03
tpad_bottom = 0.3
tpad_bottom_margin = 0.
top_axis_label_size = 0.06
top_axis_title_size = 0.075
top_axis_title_offset = 0.5
if not plot_pulls:
  an_region_size = 0.041
  nb_title_size = 0.033
  met_title_size = 0.028
  met_range_size = 0.02
  tpad_bottom = 0.
  tpad_bottom_margin = 0.1
  top_axis_label_size = 0.04
  top_axis_title_size = 0.05
  top_axis_title_offset = 0.65

gStyle.SetOptStat(0)
top = TPad("top_pad", "top_pad", 0., tpad_bottom, 1., 1.)
can.SetMargin(0., 0., 0., 0.);
can.SetFillStyle(4000);

top.SetTopMargin(0.1)
top.SetBottomMargin(tpad_bottom_margin)
top.SetLeftMargin(0.1)
top.SetRightMargin(0.05)
top.SetFillStyle(4000);
top.Draw()

#     Top pad
# -----------------------------

top.cd()
top.SetLogy()
htopdummy.GetYaxis().SetLabelSize(top_axis_label_size)
htopdummy.GetYaxis().SetRangeUser(miny,maxy)
htopdummy.GetYaxis().SetTitle("Events / Bin")
htopdummy.GetYaxis().CenterTitle()
htopdummy.GetYaxis().SetTitleSize(top_axis_title_size)
htopdummy.GetYaxis().SetTitleOffset(top_axis_title_offset)
htopdummy.Draw()

if not plot_pulls:
  htopdummy.GetXaxis().SetLabelSize(top_axis_label_size)
  #htopdummy.GetXaxis().SetLabelOffset(0.02)
  htopdummy.GetXaxis().SetTitle("Bin number")
  htopdummy.GetXaxis().CenterTitle()
  htopdummy.GetXaxis().SetTitleSize(top_axis_title_size)
  htopdummy.GetXaxis().SetTitleOffset(top_axis_title_offset)

leg = TLegend(0.4, 0.9, 0.7, 0.98)
# leg.SetTextSize(30)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetNColumns(2)

grbkg_pre.SetLineColor(tag1_color)
grbkg_pre.SetLineWidth(1)
grbkg_pre.SetFillColor(0)
grbkg_pre.SetFillStyle(0)
grbkg_pre.Draw('2')
leg.AddEntry(grbkg_pre, tag1_lbl.replace("--","-"), "F")

grbkg_post.SetFillColor(tag2_color)
grbkg_post.SetFillStyle(3144)
grbkg_post.SetLineWidth(0)
grbkg_post.Draw('2')
leg.AddEntry(grbkg_post, tag2_lbl.replace("--","-"), "F")

grdata.SetMarkerStyle(20)
grdata.Draw('P')

leg.Draw()

cmslabel = TLatex()
cmslabel.SetTextSize(0.06)
cmslabel.SetNDC(kTRUE)
cmslabel.SetTextAlign(11)
cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS}")
# cmslabel.DrawLatex(top.GetLeftMargin()+0.005, 0.92,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}")
cmslabel.SetTextAlign(31)
cmslabel.DrawLatex(1-top.GetRightMargin()-0.005, 0.92,"#font[42]{137 fb^{-1} (13 TeV)}")

binlabel = TLatex()
binlabel.SetTextSize(an_region_size)
# binlabel.SetNDC(kTRUE)
binlabel.SetTextAlign(21)
binlabel.DrawLatex(4.5, 1600,"Resolved, High-#Delta R#lower[-0.1]{_{max}}")
binlabel.DrawLatex(12.5, 1600,"Resolved, Low-#Delta R#lower[-0.1]{_{max}}")
binlabel.DrawLatex(19.5, 1600,"Boosted")
binlabel.SetTextSize(nb_title_size)
ptmiss = "p#lower[-0.1]{_{T}}#kern[-0.25]{#scale[1.15]{#lower[0.2]{^{miss}}}}";

binlabel.DrawLatex(2.5, 800,"#font[52]{3b}")
binlabel.DrawLatex(6.5, 800,"#font[52]{4b}")
binlabel.DrawLatex(10.5, 800,"#font[52]{3b}")
binlabel.DrawLatex(14.5, 800,"#font[52]{4b}")
binlabel.DrawLatex(18, 800,"#font[52]{1H}")
binlabel.DrawLatex(21, 800,"#font[52]{2H}")

for i in range(4):
    binlabel.SetTextSize(met_title_size)
    binlabel.DrawLatex(2.5+i*4, 450,"#font[52]{"+ptmiss+" [GeV]}")
    binlabel.SetTextSize(met_range_size)
    binlabel.DrawLatex(1+i*4, 250,"#font[52]{150-200}")
    binlabel.DrawLatex(2+i*4, 250,"#font[52]{200-300}")
    binlabel.DrawLatex(3+i*4, 250,"#font[52]{300-400}")
    binlabel.DrawLatex(4+i*4, 250,"#font[52]{400+}")

for i in range(2):
    binlabel.SetTextSize(met_title_size)
    binlabel.DrawLatex(17+1+i*3, 450,"#font[52]{"+ptmiss+" [GeV]}")
    binlabel.SetTextSize(met_range_size)
    binlabel.DrawLatex(17+0+i*3, 250,"#font[52]{300-500}")
    binlabel.DrawLatex(17+1+i*3, 250,"#font[52]{500-700}")
    binlabel.DrawLatex(17+2+i*3, 250,"#font[52]{700+}")

#for i in range(2):
#    binlabel.DrawLatex(4+i*18, 350,"#font[52]{200 < "+ptmiss+"#leq 350 GeV}")
#    binlabel.DrawLatex(10+i*18, 350,"#font[52]{350 < "+ptmiss+"#leq 500 GeV}")
#    binlabel.DrawLatex(16+i*18, 350,"#font[52]{"+ptmiss+"#geq 500 GeV}")
#
#binlabel.SetTextSize(0.045)
#for i in range(6):
#    binlabel.DrawLatex(2+i*6, 150,"#font[52]{1b}")
#    binlabel.DrawLatex(4+i*6, 150,"#font[52]{2b}")
#    binlabel.DrawLatex(6+i*6, 150,"#font[52]{#geq3b}")

#binlabel.SetTextAlign(11)
#binlabel.SetTextSize(0.045)
#binlabel.DrawLatex(31.3, 50,"#font[52]{Low N#lower[-0.1]{_{jets}} (odd bin #)}")
#binlabel.DrawLatex(31.3, 25,"#font[52]{High N#lower[-0.1]{_{jets}} (even bin #)}")

a = TLine()
a.SetLineWidth(1)
a.SetLineStyle(3)
a.SetLineColor(kBlack)
a.DrawLine(8.5,miny,8.5,maxy)
a.DrawLine(16.5,miny,16.5,maxy)
a.DrawLine(4.5,miny,4.5,0.4*maxy)
a.DrawLine(12.5,miny,12.5,0.4*maxy)
a.DrawLine(19.5,miny,19.5,0.4*maxy)

#     Bottom pad
# -----------------------------
if plot_pulls:
  can.cd()
  bottom = TPad("bottom_pad", "bottom_pad", 0., 0., 1., 0.3)
  bottom.SetTopMargin(0.)
  bottom.SetBottomMargin(0.3)
  bottom.SetLeftMargin(0.1)
  bottom.SetRightMargin(0.05)
  bottom.SetFillStyle(4000);
  bottom.Draw()

  bottom.cd()
  hbotdummy = TH1D("","",nhbins,0,nhbins+1.0)
  hbotdummy.GetYaxis().SetRangeUser(-3.9,3.9)
  hbotdummy.GetYaxis().SetLabelSize(0.12)
  hbotdummy.GetYaxis().SetTitle("Sig.")
  hbotdummy.GetYaxis().CenterTitle()
  hbotdummy.GetYaxis().SetTitleSize(0.15)
  hbotdummy.GetYaxis().SetTitleOffset(0.2)
  
  hbotdummy.GetXaxis().SetLabelSize(0.12)
  hbotdummy.GetXaxis().SetLabelOffset(0.02)
  hbotdummy.GetXaxis().SetTitle("Bin number")
  hbotdummy.GetXaxis().CenterTitle()
  hbotdummy.GetXaxis().SetTitleSize(0.15)
  hbotdummy.GetXaxis().SetTitleOffset(0.9)
  
  hbotdummy.Draw()
  
  grpull_pre.SetLineColor(tag1_color)
  grpull_pre.SetLineWidth(1)
  grpull_pre.SetFillColor(0)
  grpull_pre.SetFillStyle(0)
  grpull_pre.Draw('B')
  
  grpull_post.SetFillColor(tag2_color)
  grpull_post.SetFillStyle(3144)
  #grpull_post.Draw('B')
  
  a = TLine()
  a.SetLineWidth(1)
  a.SetLineStyle(3)
  a.DrawLine(8.5, -3.9, 8.5, 3.9)
  a.DrawLine(16.5, -3.9, 16.5, 3.9)
  a.SetLineColor(kGray+1)
  a.DrawLine(4.5, -3.9, 4.5, 3.9)
  a.DrawLine(12.5, -3.9, 12.5, 3.9)
  a.DrawLine(19.5, -3.9, 19.5, 3.9)
  
  b = TLine()
  b.SetLineWidth(1)
  b.SetLineColor(kGray+1)
  b.SetLineStyle(1)
  b.DrawLine(0,0, nhbins+1,0)
  b.SetLineStyle(2)
  for i in [-1.,1.]:
      b.DrawLine(0,i, nhbins+1,i)
  b.SetLineStyle(3)
  for i in [-2.,2.]:
      b.DrawLine(0,i, nhbins+1,i)
  b.SetLineStyle(3)
  for i in [-3.,3.]:
      b.DrawLine(0,i, nhbins+1,i)

pname = 'plots/results_plot.pdf'
can.Print(pname)

print 'open', pname
