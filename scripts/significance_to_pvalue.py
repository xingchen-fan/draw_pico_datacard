#!/usr/bin/env python
# http://www.physics.ucla.edu/~cousins/stats/combine_significance.cxx
import sys
import ROOT
import math
import argparse

parser = argparse.ArgumentParser(description='''Converts significance to a p-value. Shows the p-value when the significance is either one-sided or two-sided.''')
parser.add_argument('significance', help='Significance value to convet to a p-value')
args = parser.parse_args()

significance = float(args.significance)
pvalue = 0.5 - ROOT.Math.erf(significance/math.sqrt(2))/2
twotail_pvalue = 1.0 - ROOT.TMath.Erf(significance/math.sqrt(2))
print("significance: "+str(significance)+", (one-tail) pvalue: "+str(pvalue)+" (two-tail) pvalue: "+str(twotail_pvalue))
