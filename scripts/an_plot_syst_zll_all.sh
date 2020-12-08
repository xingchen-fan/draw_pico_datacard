mkdir -p plots.trash
mv plots/*.pdf plots.trash
set -e // Exit when any command fails
mkdir -p plots.an/zll/2016
rm -f plots.an/zll/2016/*.pdf
mkdir -p plots.an/zll/2017
rm -f plots.an/zll/2017/*.pdf
mkdir -p plots.an/zll/2018
rm -f plots.an/zll/2018/*.pdf
mkdir -p plots.an/zll/run2
rm -f plots.an/zll/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2016
mv plots/*.pdf plots.an/zll/2016
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2017
mv plots/*.pdf plots.an/zll/2017
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2018
mv plots/*.pdf plots.an/zll/2018
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2016,2017,2018
mv plots/*.pdf plots.an/zll/run2

mkdir -p plots.an/zll.unblind/2016
rm -f plots.an/zll.unblind/2016/*.pdf
mkdir -p plots.an/zll.unblind/2017
rm -f plots.an/zll.unblind/2017/*.pdf
mkdir -p plots.an/zll.unblind/2018
rm -f plots.an/zll.unblind/2018/*.pdf
mkdir -p plots.an/zll.unblind/run2
rm -f plots.an/zll.unblind/run2/*.pdf
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2016 --unblind
mv plots/*.pdf plots.an/zll.unblind/2016
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2017 --unblind
mv plots/*.pdf plots.an/zll.unblind/2017
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2018 --unblind
mv plots/*.pdf plots.an/zll.unblind/2018
scons && ./run/higgsino/an_plot_syst_zll.exe --year 2016,2017,2018 --unblind
mv plots/*.pdf plots.an/zll.unblind/run2
