draw_pico
========

Repository inherited from: [https://github.com/ald77/ra4_draw](https://github.com/ald77/ra4_draw)

For examples of functionality see: [src/core/test.cxx](src/core/test.cxx)

## Setup

To use the batch system:
~~~~bash
git clone --recurse-submodules https://github.com/richstu/draw_pico
source set_env.sh
~~~~

## Higgsino useful commands

To plot overlap between the boosted and resolved analysis:

~~~~bash
./run/higgsino/plot_regions.exe #boosted overlap in resolved regions
./run/higgsino/plot_regions.exe --boo #resolved overlap in boosted regions
~~~~

Write all the datacards:

~~~~bash
./run/higgsino/write_datacards.exe -t boosted -o boosted/
~~~~

For comparison of signal yields to other tables, remember to use `--recomet` such that the yield is not averaged with the one obtained using GenMET. To make a datacard for the boosted case, use `-t boosted`. Use option `--unblind` to include data. To run on a particular point add, e.g. `-p "700_1"`. Then to get a limit interactively:

~~~~bash
./run/higgsino/scan_point.exe -f test/datacard-TChiHH_mChi-700_mLSP-1_Tune_2016_resolved.txt
~~~~

Get limits for the full scan in the batch:

~~~~bash
./scripts/write_combine_cmds.py --card_dir test
~~~~

To combine outputs and make the limit plot

~~~~bash
cat test/scan_point*/limit*txt | sort >> resolved_test_limits.txt
./run/higgsino/plot_limit.exe -f resolved_test_limits.txt
~~~~


## Zgamma useful commands