set -e // Exit when any command fails
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2016 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2017 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample search --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample ttbar  --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample zll    --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
scons && ./run/higgsino/plot_n_minus_1.exe --sample qcd    --year 2018 2>&1 | tee trash_log && cat trash_log | grep open | awk -F ' ' '{print $2}' | xargs -tI {} imgcat {}
