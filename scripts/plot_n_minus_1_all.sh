mkdir -p plots.trash
mv plots/*.pdf plots.trash
set -e // Exit when any command fails
mkdir -p plots.an/n_minus_1/2016
rm -f plots.an/n_minus_1/2016/*.pdf
mkdir -p plots.an/n_minus_1/2017
rm -f plots.an/n_minus_1/2017/*.pdf
mkdir -p plots.an/n_minus_1/2018
rm -f plots.an/n_minus_1/2018/*.pdf
mkdir -p plots.an/n_minus_1/run2
rm -f plots.an/n_minus_1/run2/*.pdf
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2016 
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2016 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2016 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2016 
mv plots/*.pdf plots.an/n_minus_1/2016
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2017 
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2017 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2017 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2017 
mv plots/*.pdf plots.an/n_minus_1/2017
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2018 
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2018 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2018 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2018 
mv plots/*.pdf plots.an/n_minus_1/2018
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year run2
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year run2
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year run2
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year run2
mv plots/*.pdf plots.an/n_minus_1/run2

mkdir -p plots.an/n_minus_1.unblind/2016
rm -f plots.an/n_minus_1.unblind/2016/*.pdf
mkdir -p plots.an/n_minus_1.unblind/2017
rm -f plots.an/n_minus_1.unblind/2017/*.pdf
mkdir -p plots.an/n_minus_1.unblind/2018
rm -f plots.an/n_minus_1.unblind/2018/*.pdf
mkdir -p plots.an/n_minus_1.unblind/run2
rm -f plots.an/n_minus_1.unblind/run2/*.pdf
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2016 --unblind
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2016 --unblind
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2016 --unblind
mv plots/*.pdf plots.an/n_minus_1.unblind/2016
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2017 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2017 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2017 --unblind 
mv plots/*.pdf plots.an/n_minus_1.unblind/2017
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2018 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2018 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2018 --unblind 
mv plots/*.pdf plots.an/n_minus_1.unblind/2018
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year run2 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year run2 --unblind 
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year run2 --unblind 
mv plots/*.pdf plots.an/n_minus_1.unblind/run2
