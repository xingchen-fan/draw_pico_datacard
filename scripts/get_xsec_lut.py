#!/usr/bin/env python3

import subprocess

# The inputs needed to run this are at: /afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV
# N.B. Initiliaze the vectors connected to the branches to 0 to avoid seg faults in get_gaugino.C
# then run this inside the folder to avoid the need for further modifications

for mass in range(1500,1501,25):
  if mass==125: mass=127
  result = subprocess.check_output('root -l -q \'get_gaugino.C("CN","hino",'+str(mass)+')\'', shell=True)
  # print(result)
  result = result.decode()
  xsec = float(result.split(' is ')[-1].split(' [pb] ')[0])
  xsec_unc = float(result.split(' +/- ')[-1].split(' [rel')[0])
  print('else if(hig_mass =={}) {{ xsec = .5824*.5824*{}; xsec_unc = {}; return;}}'.format(mass, xsec, xsec_unc))


