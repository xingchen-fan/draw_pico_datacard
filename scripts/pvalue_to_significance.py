#!/usr/bin/env python
# http://www.physics.ucla.edu/~cousins/stats/combine_significance.cxx
import sys
import ROOT
import math
import argparse

parser = argparse.ArgumentParser(description='Converts a p-value to a one-tailed significance and a two-tailed significance')
parser.add_argument('pvalue', help='p-value to convert to a significance')
args = parser.parse_args()

pvalue = float(args.pvalue)
significance = ROOT.Math.normal_quantile_c(pvalue,1)
two_tail_significance = math.sqrt(2)*ROOT.TMath.ErfInverse(1 - pvalue);
print("pvalue: "+str(pvalue)+", (one-tail) significance: "+str(significance)+" (two-tail) significance: "+str(two_tail_significance))
