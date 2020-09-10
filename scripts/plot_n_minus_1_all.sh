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
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
mv plots/*.pdf plots.an/n_minus_1/2016
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
mv plots/*.pdf plots.an/n_minus_1/2017
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
mv plots/*.pdf plots.an/n_minus_1/2018
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year run2 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year run2 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year run2 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year run2 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
mv plots/*.pdf plots.an/n_minus_1/run2
