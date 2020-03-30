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
./run/higgsino/write_datacards.exe -m N1N2 -o test
~~~~

For comparison of signal yields to other tables, remember to use `--recomet` such that the yield is not averaged with the one obtained using GenMET. To make a datacard for the boosted case, use `-t boosted`. Use option `--unblind` to include data. To run on a particular point add, e.g. `-p "700_1"`. Then to get a limit interactively:

~~~~bash
./run/higgsino/scan_point.exe -f test/datacard-TChiHH_mChi-700_mLSP-1_Tune_2016_resolved.txt
~~~~

Get limits for the full scan in the batch:

~~~~bash
./scripts/write_combine_cmds.py --card_dir test -m N1N2
~~~~

To combine outputs and make the limit plot

~~~~bash
cat test/scan_point*/limit*txt | sort >> resolved_test_limits.txt
./run/higgsino/plot_limit.exe -f resolved_test_limits.txt
~~~~

## Getting Higgsino cross-section (CN to N1N2) scale factors.

### Step 1. Download cross-section files

~~~~bash
scp -r lxplus:/afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV cross_section
~~~~

### Step 2. Fix bug in script

Set masses, xsecs, xsecUncs to 0 when initializing.

~~~~bash
cd cross_section
sed -i 's/std::vector<double>\* masses;/std::vector<double>\* masses=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecs;/std::vector<double>\* xsecs=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecUncs;/std::vector<double>\* xsecUncs=0;/' get_gaugino.C
~~~~

### Step 3. Run script for all mass points.

 model "CN" (mixing) or "N1N2" (no mixing).

~~~~bash
cd cross_section
../scripts/change_higgsino_cross_sections.py -i /net/cms29/cms29r0/pico/NanoAODv5/nano/2016/SMS-TChiHH_2D
~~~~

## Zgamma useful commands

