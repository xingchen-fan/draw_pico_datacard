mkdir -p plots.trash
mv plots/*.pdf plots.trash
set -e // Exit when any command fails
mkdir -p  plots.an/qcd/2016
rm -f plots.an/qcd/2016/*.pdf
mkdir -p  plots.an/qcd/2017
rm -f plots.an/qcd/2017/*.pdf
mkdir -p  plots.an/qcd/2018
rm -f plots.an/qcd/2018/*.pdf
mkdir -p  plots.an/qcd/run2
rm -f plots.an/qcd/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2016
mv plots/*.pdf plots.an/qcd/2016
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2017
mv plots/*.pdf plots.an/qcd/2017
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2018
mv plots/*.pdf plots.an/qcd/2018
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2016,2017,2018
mv plots/*.pdf plots.an/qcd/run2

mkdir -p  plots.an/qcd.unblind/2016
rm -f plots.an/qcd.unblind/2016/*.pdf
mkdir -p  plots.an/qcd.unblind/2017
rm -f plots.an/qcd.unblind/2017/*.pdf
mkdir -p  plots.an/qcd.unblind/2018
rm -f plots.an/qcd.unblind/2018/*.pdf
mkdir -p  plots.an/qcd.unblind/run2
rm -f plots.an/qcd.unblind/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2016 --unblind
mv plots/*.pdf plots.an/qcd.unblind/2016
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2017 --unblind
mv plots/*.pdf plots.an/qcd.unblind/2017
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2018 --unblind
mv plots/*.pdf plots.an/qcd.unblind/2018
scons && ./run/higgsino/an_plot_syst_qcd.exe --year 2016,2017,2018 --unblind
mv plots/*.pdf plots.an/qcd.unblind/run2
