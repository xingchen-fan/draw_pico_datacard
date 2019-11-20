draw_pico
========

Repository inherited from: [https://github.com/ald77/ra4_draw](https://github.com/ald77/ra4_draw)

For examples of functionality see: [src/core/test.cxx](src/core/test.cxx)

## Setup

~~~~bash
git clone --recurse-submodules https://github.com/richstu/draw_pico
~~~~

## Higgsino useful commands

To plot overlap between the boosted and resolved analysis:

~~~~bash
./compile.py && ./run/higgsino/plot_regions.exe #boosted overlap in resolved regions
./compile.py && ./run/higgsino/plot_regions.exe --boo #resolved overlap in boosted regions
~~~~

Write all the datacards:

~~~~bash
./compile.py && ./run/higgsino/write_datacards.exe -p "127_1,150_1,175_1,200_1,225_1,250_1,275_1,300_1,325_1,350_1,375_1,400_1,425_1,450_1,475_1,500_1,525_1,550_1,575_1,600_1,625_1,650_1,675_1,700_1,725_1,750_1,775_1,800_1,825_1,850_1,875_1,900_1,925_1,950_1,975_1,1000_1,1125_1,1150_1,1175_1,1200_1,1225_1,1250_1,1275_1" -t resolved -o test/
~~~~

To do it for boosted, use `-t boosted`. Use option `--unblind` to include data. Then to get a limit interactively, e.g.:

~~~~bash
./compile.py && ./run/higgsino/scan_point.exe -f test/datacard-TChiHH_mChi-700_mLSP-1_Tune_2016_resolved.txt
~~~~

Run the full scan in the batch:

~~~~bash
./scripts/write_combine_cmds.py --card_dir test
~~~~

To combine outputs and make the limit plot

~~~~bash
cat test/scan_point*/limit*txt | sort >> resolved_test_limits.txt
./run/higgsino/plot_limit.exe -f resolved_test_limits.txt
~~~~


## Zgamma useful commands