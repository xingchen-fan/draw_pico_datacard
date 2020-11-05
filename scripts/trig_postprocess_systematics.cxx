//script to generate systematic table
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include "TMath.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"

//helper functions
std::string make_systematic_table_row(TGraphAsymmErrors *gr, std::string name);
std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr);
std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table);
std::string make_systematic_table_row(std::vector<double> values, bool total=false);
std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values);

int trig_postprocess_systematics() {
	bool do_datamc = false;
	TFile* f = TFile::Open("ntuples/trig_eff_2018_full.root");
	int year = 2018;
	ofstream table_file;
	std::vector<std::vector<double>> variations_table;
	std::vector<double> variations_row;
	std::vector<std::vector<double>> reference_table;
	std::vector<double> reference_row;
	std::vector<std::vector<double>> datamc_table;
	std::vector<double> datamc_row;
	std::vector<std::vector<double>> systematics_rows;
	table_file.open("tables/trig_sys_table.tex");
	table_file << "\\begin{tabular}{lccccccc}\\hline \\hline \n";
        table_file << "E$_\\text{T}^\\text{miss}$ range& [150,160]& [160,180]& [180,200]& [200,225]& [225,250]& [250,300]& [300+]\\\\ \\hline \\hline \n";
	table_file << "\\multicolumn{8}{c}{Kinematic variations}\\\\ \\hline \n";
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj3_ratio")),"$N_\\text{jets}\\geq 4$ (Nominal)") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj3_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj4_ratio")),"$N_\\text{jets}=4$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj4_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj5_ratio")),"$N_\\text{jets}=5$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nj5_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb0_ratio")),"$N_\\text{b}=0$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb0_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb1_ratio")),"$N_\\text{b}=1$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb1_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb2_ratio")),"$N_\\text{b}\\geq2$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_nb2_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_amlow_ratio")),"$\\langle m\\rangle\\le$ 140 GeV, 300 $\\leq H_{T} <$ 400 GeV") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_amlow_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_amhigh_ratio")),"$\\langle m\\rangle >$ 140 GeV, 300 $\\leq H_{T} <$ 400 GeV") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_amhigh_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxlow_ratio")),"$\\Delta R_\\text{max}< 2.2$") << "\n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxlow_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxhigh_ratio")),"$\\Delta R_\\text{max}\\geq 2.2$") << "\\hline \n";
	variations_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_drmaxhigh_ratio"))));
	variations_row = get_systematic_table_variations(variations_table);
	systematics_rows.push_back(variations_row);
	table_file << make_systematic_table_row(variations_row) << "\n";
	table_file << "\\multicolumn{8}{c}{Reference trigger (for $p_\\text{Tjet}>500$ GeV)}\\\\ \\hline \n";
	if (year == 2016)
		table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_isoeljet500_ratio")),"Ele27 (Nominal)") << "\n";
	else
		table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_isoeljet500_ratio")),"Ele35 (Nominal)") << "\n";
	reference_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_isoeljet500_ratio"))));
	table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_jetjet500_ratio")),"Jet500") << "\\hline \n";
	reference_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_jetjet500_ratio"))));
	reference_row = get_systematic_table_variations(reference_table);
	systematics_rows.push_back(reference_row);
	table_file << make_systematic_table_row(reference_row) << "\n";
	if (do_datamc) {
		table_file << "\\multicolumn{8}{c}{Data compared to MC (MET120\\_HT60 trigger only)}\\\\ \\hline \n";
		table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio")),"Data (nominal)") << "\n";
		datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_met120_ratio"))));
		table_file << make_systematic_table_row(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio")),"MC") << "\\hline \n";
		datamc_table.push_back(get_systematic_table_values(static_cast<TGraphAsymmErrors*>(f->Get("hist_mc_ratio"))));
		datamc_row = get_systematic_table_variations(datamc_table);
		systematics_rows.push_back(datamc_row);
		table_file << make_systematic_table_row(datamc_row) << "\n";
	}
	table_file << make_systematic_table_row(add_variations_quadrature(systematics_rows), true) << "\n";
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

std::vector<double> get_systematic_table_values(TGraphAsymmErrors *gr) {
	std::vector<double> systematic_table_values;
	double x, y;
        for (int i = 0; i<7; i++) {
		gr->GetPoint(i,x,y);
		y = y*100;
		systematic_table_values.push_back(y);
	}
	return systematic_table_values;
}

std::vector<double> get_systematic_table_variations(std::vector<std::vector<double>> sys_table) {
	//0 is nominal
	std::vector<double> systematic_table_variations;
	for (unsigned int met_idx = 0; met_idx < sys_table[0].size(); met_idx++) {
		double max_difference = 0;
		for (unsigned int row = 1; row < sys_table.size(); row++) {
			if (TMath::Abs(sys_table[0][met_idx]-sys_table[row][met_idx]) > max_difference) {
				max_difference = TMath::Abs(sys_table[0][met_idx]-sys_table[row][met_idx]);
			}
		}
		systematic_table_variations.push_back(max_difference/sys_table[0][met_idx]);
	}
	return systematic_table_variations;
}

std::string make_systematic_table_row(std::vector<double> values, bool total) {
	std::string row_text = "Syst. uncertainty ";
	if (total) {
		row_text = "{\\bf Total Syst. uncertainty} ";
	}
	std::ostringstream str_stream;
	str_stream << std::fixed;
	double x, y;
        for (int i = 0; i<7; i++) {
		row_text = row_text + "& ";
		str_stream << std::setprecision(1);
		str_stream << (values[i]*100);
		row_text = row_text + str_stream.str();
		str_stream.str("");
		str_stream.clear();
	}
	row_text = row_text + "\\\\ \\hline \\hline";
	return row_text;
}

std::vector<double> add_variations_quadrature(std::vector<std::vector<double>> values) {
	std::vector<double> variations_quadrature;
	for (unsigned int met_idx = 0; met_idx < values[0].size(); met_idx++) {
		double quad_sum = 0;
		for (unsigned int row = 0; row < values.size(); row++) {
			quad_sum += values[row][met_idx]*values[row][met_idx];
		}
		variations_quadrature.push_back(TMath::Sqrt(quad_sum));
	}
	return variations_quadrature;
}
