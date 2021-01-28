#!/usr/bin/env python3
import subprocess
import argparse
import glob
import re

# The inputs needed to run this are at: /afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV
# N.B. Initiliaze the vectors connected to the branches to 0 to avoid seg faults in get_gaugino.C
# then run this inside the folder to avoid the need for further modifications

# ./higgsino_cross_sections.py -i xxx -m CN # Prints cross sections function
# ./higgsino_cross_sections.py -i xxx -c -m CN -tm N1N2 # Prints change CN to N1N2 NamedFunc

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('-i','--signal_path', required=True, default = "")
  parser.add_argument('-c','--change', action='store_true') # Print change or cross_section. Default is cross_section. 
  parser.add_argument('-m','--model', default = "CN") 
  parser.add_argument('-tm','--to_model', default = "N1N2") 
  args = parser.parse_args()

  # Find mass of higgsino
  signal_files = glob.glob(args.signal_path+'/*.root')
  chi_mass_list = set()
  for signal_file in signal_files:
    mChi = re.findall(r"mChi-\d+",signal_file)[0].replace('mChi-','')
    chi_mass_list.add(int(mChi))

  if args.change:
    print('const NamedFunc w_CNToN1N2("w_CNToN1N2", [](const Baby &b) -> NamedFunc::ScalarType{')
    print('  if(b.type() != 106000) return 1;')
  else:
    print('// Unit: xsec=pb, xsec_unc=%')
    if args.model == "CN":
      print('void higgsinoCrossSection(int hig_mass, double &xsec, double &xsec_unc) {')
    else:
      print('void higgsino2DCrossSection(int hig_mass, double &xsec, double &xsec_unc) {')

  # Print cross sections
  first_line = True
  for mass in sorted(chi_mass_list):
    if mass == 1500: continue
    result = subprocess.check_output('root -l -q \'get_gaugino.C("'+args.model+'","hino",'+str(mass)+')\'', shell=True)
    # print(result)
    result = result.decode()
    xsec = float(result.split(' is ')[-1].split(' [pb] ')[0])
    xsec_unc = float(result.split(' +/- ')[-1].split(' [rel')[0])
    if args.change:
      result = subprocess.check_output('root -l -q \'get_gaugino.C("'+args.to_model+'","hino",'+str(mass)+')\'', shell=True)
      # print(result)
      result = result.decode()
      to_xsec = float(result.split(' is ')[-1].split(' [pb] ')[0])
      to_xsec_unc = float(result.split(' +/- ')[-1].split(' [rel')[0])
      if first_line:  print('  if(b.mprod() =='+str(mass)+') return '+str(to_xsec)+'/'+str(xsec)+';')
      else: print('  else if(b.mprod() =='+str(mass)+') return '+str(to_xsec)+'/'+str(xsec)+';')
    else:
      if first_line: print('  if(hig_mass =={}) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, xsec, xsec_unc))
      else: print('  else if(hig_mass =={}) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, xsec, xsec_unc))
    first_line = False
  if args.change:
    print('  else return 0;')
    print('});')
  else:
    print('  else{ xsec = 0; xsec_unc = 0;}')
    print('}')

