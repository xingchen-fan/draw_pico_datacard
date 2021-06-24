draw_pico
========

Repository inherited from: [https://github.com/ald77/ra4_draw](https://github.com/ald77/ra4_draw)

For examples of functionality see: [src/core/test.cxx](src/core/test.cxx)

## Setup

Many scripts require one to first source the `set_env.sh` script.

To use the batch system:
~~~~bash
git clone --recurse-submodules https://github.com/richstu/draw_pico
source set_env.sh
~~~~

## Higgsino useful commands

### To make datacards and get limits:

To plot overlap between the boosted and resolved analysis:

~~~~bash
./run/higgsino/plot_regions.exe #boosted overlap in resolved regions
./run/higgsino/plot_regions.exe --boo #resolved overlap in boosted regions
~~~~

Write all the datacards:

~~~~bash
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/tchihh_onedim/ -m CN -1 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/tchihh_twodim/ -m N1N2 -2 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_onedim/ -m T5HH -f -1 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion
~~~~

For comparison of signal yields to other tables, remember to use `--recomet` such that the yield is not averaged with the one obtained using GenMET. To make a datacard for the boosted case, use `-t boosted`. Use option `--unblind` to include data. To run on a particular point add, e.g. `-p "700_1"`. For computers with smaller amounts of memory, there may not be sufficient memory to generate all the 2D T5HH datacards at once. In this case, you can split the job in half with the `-s` flag, which sets a lower/upper bound on the gluino mass. 

~~~~bash
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion -s "-2100"
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion -s 2100
~~~~

Then to get a limit interactively:

~~~~bash
./run/higgsino/scan_point.exe -f datacards/tchihh_onedim/datacard-TChiHH_mChi-700_mLSP-0_Tune_2016,2017,2018_priority1_resolved.txt
~~~~

To get limits for the full scan in the batch, use the following script to generate batch commands.

~~~~bash
./scripts/write_combine_cmds.py --card_dir datacards/tchihh_onedim -m CN
./scripts/write_combine_cmds.py --card_dir datacards/tchihh_twodim -m N1N2
./scripts/write_combine_cmds.py --card_dir datacards/t5hh_twodim -m T5HH
./scripts/write_combine_cmds.py --card_dir datacards/t5hh_twodim -m T5HH
~~~~

These can be run locally with `./scripts/run_commands.py` or submitted to the UCSB batch system with `auto_submit_jobs.py`. To combine outputs and make the limit plot

~~~~bash
cat datacards/tchihh_onedim/scan_point*/limit*txt | sort >> tchihh_onedim_resolved_limits.txt
cat datacards/tchihh_twodim/scan_point*/limit*txt | sort >> tchihh_twodim_resolved_limits.txt
cat datacards/t5hh_onedim/scan_point*/limit*txt | sort >> t5hh_onedim_resolved_limits.txt
cat datacards/t5hh_twodim/scan_point*/limit*txt | sort >> t5hh_twodim_resolved_limits.txt
./run/higgsino/plot_limit.exe -f tchihh_onedim_resolved_limits.txt --drawData -t tchihh_onedim_resolved
./run/higgsino/limit_scan.exe -f tchihh_twodim_resolved_limits.txt -m N1N2 -t tchihh_twodim_resolved --unblind
./run/higgsino/plot_limit.exe -f t5hh_onedim_resolved_limits.txt --drawData -m T5HH -t t5hh_onedim_resolved
./run/higgsino/limit_scan.exe -f t5hh_twodim_resolved_limits.txt -m T5HH -t t5hh_twodim_resolved --unblind
~~~~

### To generate AN plots:

~~~~bash
./run/higgsino/an_plot_triggers.exe --unblind --year 2016 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year 2017 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year 2018 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year run2
./run/higgsino/an_plot_syst_ttbar.exe --unblind --year run2
./run/higgsino/an_plot_syst_zll.exe --unblind --year run2
./run/higgsino/an_plot_syst_qcd.exe --unblind --year run2
~~~~

### To generate paper plots:

~~~~bash
./run/higgsino/plot_search_unblind.exe -u -a -o plot_baseline,paper_style,plot_in_btags,plot_in_btags_with_met_split
./scripts/plot_ra4style_results.py
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
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -m CN
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -m N1N2
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -c
~~~~

## Zgamma useful commands

