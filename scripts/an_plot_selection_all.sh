mkdir -p plots.trash
mv plots/*.pdf plots.trash
set -e // Exit when any command fails
mkdir -p plots.an/selection/2016
rm -f plots.an/selection/2016/*.pdf
mkdir -p plots.an/selection/2017
rm -f plots.an/selection/2017/*.pdf
mkdir -p plots.an/selection/2018
rm -f plots.an/selection/2018/*.pdf
mkdir -p plots.an/selection/run2
rm -f plots.an/selection/run2/*.pdf
scons && ./run/higgsino/an_plot_selection.exe --year 2016
mv plots/*.pdf plots.an/selection/2016
scons && ./run/higgsino/an_plot_selection.exe --year 2017
mv plots/*.pdf plots.an/selection/2017
scons && ./run/higgsino/an_plot_selection.exe --year 2018
mv plots/*.pdf plots.an/selection/2018
scons && ./run/higgsino/an_plot_selection.exe --year 2016,2017,2018
mv plots/*.pdf plots.an/selection/run2
