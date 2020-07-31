mkdir -p tables.trash
mv tables/* tables.trash
set -e // Exit when any command fails
#mkdir -p plots.trash
#mv plots/*.pdf plots.trash
#mkdir -p plots.an/bkgest/2016
#rm -f plots.an/bkgest/2016/*.pdf
#mkdir -p plots.an/bkgest/2017
#rm -f plots.an/bkgest/2017/*.pdf
#mkdir -p plots.an/bkgest/2018
#rm -f plots.an/bkgest/2018/*.pdf
#mkdir -p plots.an/bkgest/run2
#rm -f plots.an/bkgest/run2/*.pdf
#scons && ./run/higgsino/an_plot_backgroundEstimation.exe --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {} 
#mv plots/*.pdf plots.an/bkgest/2016
#scons && ./run/higgsino/an_plot_backgroundEstimation.exe --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
#mv plots/*.pdf plots.an/bkgest/2017
#scons && ./run/higgsino/an_plot_backgroundEstimation.exe --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
#mv plots/*.pdf plots.an/bkgest/2018
#scons && ./run/higgsino/an_plot_backgroundEstimation.exe --year 2016,2017,2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
#mv plots/*.pdf plots.an/bkgest/run2

makePlot() {
  scons && $1 --year $2 --sample $3
  mv tables/*.tex tables.an/regions/$2
  cd tables.an/regions/$2
  pdflatex $4.tex
  cd -
  imgcat tables.an/regions/$2/$4.pdf
}


mkdir -p tables.an/regions/2016
rm -f tables.an/regions/2016/*
makePlot "./run/higgsino/table_sample_regions.exe" "2016" "search" "table_search_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2016" "ttbar" "table_ttbar_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2016" "zll" "table_zll_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2016" "qcd" "table_qcd_regions"

mkdir -p tables.an/regions/2017
rm -f tables.an/regions/2017/*
makePlot "./run/higgsino/table_sample_regions.exe" "2017" "search" "table_search_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2017" "ttbar" "table_ttbar_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2017" "zll" "table_zll_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2017" "qcd" "table_qcd_regions"

mkdir -p tables.an/regions/2018
rm -f tables.an/regions/2018/*
makePlot "./run/higgsino/table_sample_regions.exe" "2018" "search" "table_search_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2018" "ttbar" "table_ttbar_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2018" "zll" "table_zll_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "2018" "qcd" "table_qcd_regions"

mkdir -p tables.an/regions/run2
rm -f tables.an/regions/run2/*
makePlot "./run/higgsino/table_sample_regions.exe" "run2" "search" "table_search_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "run2" "ttbar" "table_ttbar_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "run2" "zll" "table_zll_regions"
makePlot "./run/higgsino/table_sample_regions.exe" "run2" "qcd" "table_qcd_regions"
