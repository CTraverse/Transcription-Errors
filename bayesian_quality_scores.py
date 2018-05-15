# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 15:49:06 2016

@author: chuck
"""
import math

#This script tests two methods of recalculating quality scores in sequences from the CirSeq method
##CirSeq generates linked repeats of the same mRNA fragment and each base in reach repeat has its own quality score
##Does a Bayesian approach give better results than the original method?

nucs = ['A', 'T', 'C', 'G', 'N']

Rep1 = ['A', '=']
Rep2 = ['T', '?']
Rep3 = ['T', 'I']

##Bayes' method
def probability(base, call, quality):
    if base == call:
        return (1 - 10**(float(ord(quality)-33)/-10))
    else:
        return (10**(float(ord(quality)-33)/-10))/3


post_dist = []
for base in nucs:

    consensus_numerator = probability(base,Rep1[0],Rep1[1]) * probability(base,Rep2[0],Rep2[1]) * probability(base,Rep3[0],Rep3[1])
    consensus_denominator = 0
    for nuc in nucs:
        consensus_denominator += probability(nuc,Rep1[0],Rep1[1]) * probability(nuc,Rep2[0],Rep2[1]) * probability(nuc,Rep3[0],Rep3[1])
    post_dist.append([base, float(consensus_numerator / consensus_denominator)])
    
max_prob = max(post_dist, key = lambda x:x[1])

print(post_dist)

print(max_prob, chr(int(round(-10*math.log10(1-max_prob[1])))+33), int(round(-10*math.log10(1-max_prob[1]))))

##########


###Acevedo method
if Rep1[0] == Rep2[0] == Rep3[0]:
    ReCalculatedQualityScore = int(round((float(ord(Rep1[1]) + ord(Rep2[1]) + ord(Rep3[1]) - 99)/3.0)))
    
elif Rep1[0] == Rep2[0] and Rep3[0] != Rep1[0]:
    RawProbability = (10**(float(ord(Rep1[1])-33)/-10.0))*(10**(float(ord(Rep2[1])-33)/-10.0))*(1-(1/3.0)*(10**(float(ord(Rep3[1])-33)/-10.0)))
    ReCalculatedQualityScore = int(round(-10*math.log10(RawProbability)/3.0))
    print(RawProbability)
    
elif Rep1[0] == Rep3[0] and Rep2[0] != Rep1[0]:
    RawProbability = (10**(float(ord(Rep1[1])-33)/-10.0))*(10**(float(ord(Rep3[1])-33)/-10.0))*(1-(1/3.0)*(10**(float(ord(Rep2[1])-33)/-10.0)))
    ReCalculatedQualityScore = int(round(-10*math.log10(RawProbability)/3.0))
    print(RawProbability)
    
elif Rep2[0] == Rep3[0] and Rep2[0] != Rep1[0]:
    RawProbability = (10**(float(ord(Rep2[1])-33)/-10.0))*(10**(float(ord(Rep3[1])-33)/-10.0))*(1-(1/3.0)*(10**(float(ord(Rep1[1])-33)/-10.0)))
    ReCalculatedQualityScore = int(round(-10*math.log10(RawProbability)/3.0))
    print(RawProbability)
else:
    "Won't work"

print("Acevedo", ReCalculatedQualityScore)