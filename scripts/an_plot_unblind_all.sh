mkdir -p plots.trash
mv plots/*.pdf plots.trash
mkdir -p tables.trash
mv tables/* tables.trash
mkdir -p plots.an/unblind
set -e // Exit when any command fails

scons&& run/higgsino/plot_kappas.exe --sample search --year run2 --scen mc --string_options search_ttbar_same_bins
scons&& run/higgsino/plot_kappas.exe --sample ttbar --unblind --year run2 --scen data --string_options search_ttbar_same_bins
scons&& run/higgsino/an_plot_syst_ttbar.exe --year run2 --unblind --string_options plot_planes,plot_btags
scons&& run/higgsino/table_sample_regions.exe --year run2 --no_signal --string_options split_mc_in_detail

mv plots/* plots.an/unblind
mv tables/table_search_regions.pdf plots.an/unblind
