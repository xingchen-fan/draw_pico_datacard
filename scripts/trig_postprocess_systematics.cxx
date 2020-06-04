#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "TFile.h"
#include "TGraphAsymmErrors.h"

std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name);

int trigger_postprocessing() {
	TFile* f = TFile::Open("ntuples/trig_eff_old7.root");
	ofstream table_file;
	table_file.open("tables/trig_sys_table.tex");
	table_file << "\\begin{tabular}{lccccccc}\\hline \\hline \n";
        table_file << "E$_\\text{T}^\\text{miss}$ range& [150,160]& [160,180]& [180,200]& [200,225]& [225,250]& [250,300]& [300+]\\\\ \\hline \\hline \n";
	table_file << "Kinematic Variations& & & & & & & \\\\ \\hline \n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj3_ratio")),"$N_\\text{jets}>3$ (Nominal)") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj4_ratio")),"$N_\\text{jets}=4$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj5_ratio")),"$N_\\text{jets}=5$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb0_ratio")),"$N_\\text{b}=0$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb1_ratio")),"$N_\\text{b}=1$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb2_ratio")),"$N_\\text{b}\\geq2$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_amlow_ratio")),"$\\langle m\\rangle\\le 140$ ") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_amhigh_ratio")),"$\\langle m\\rangle > 140$ ") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxlow_ratio")),"$\\Delta R_\\text{max}\\le 2.2$") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxhigh_ratio")),"$\\Delta R_\\text{max}> 2.2$") << "\\hline \n";
	table_file << "Syst. uncertainty& ?& ?& ?& ?& ?& ?& ?\\\\ \\hline \\hline \n";
	table_file << "Ref Trigger ($p_{Tjet}>500$ GeV)& & & & & & & \\\\ \\hline \n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj3jet500_ratio")),"Ele27 (Nominal)") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_jetjet500_ratio")),"Jet500") << "\\hline \n";
	table_file << "Syst. uncertainty& ?& ?& ?& ?& ?& ?& ?\\\\ \\hline \\hline \n";
	table_file << "Data vs MC (MET120 only) & & & & & & & \\\\ \\hline \n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio")),"Data (nominal)") << "\n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio")),"MC") << "\\hline \n";
	table_file << "Syst. uncertainty& ?& ?& ?& ?& ?& ?& ?\\\\ \\hline \\hline \n";
	table_file << "Total syst. uncertainty& ?& ?& ?& ?& ?& ?& ?\\\\ \\hline \\hline \n";
	table_file << "\\end{tabular} \n";
        table_file.close();
	f->Close();
	return 0;
}

std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name) {
	std::string row_text = name;
	std::ostringstream str_stream;
	str_stream << std::fixed;
	double x, y;
        for (int i = 0; i<7; i++) {
		row_text = row_text + "& ";
		gr->GetPoint(i,x,y);
		y = y*100;
		str_stream << std::setprecision(1);
		str_stream << (y);
		row_text = row_text + str_stream.str();
		str_stream.str("");
		str_stream.clear();
		str_stream << std::setprecision(1);
		str_stream << (gr->GetErrorYlow(i)*100);
		row_text = row_text + "$_{-" + str_stream.str();
		str_stream.str("");
		str_stream.clear();
		str_stream << (gr->GetErrorYhigh(i)*100);
	       	row_text = row_text + "}^{+" + str_stream.str() + "}$";
		str_stream.str("");
		str_stream.clear();
	}
	row_text = row_text + "\\\\";
	return row_text;
}
