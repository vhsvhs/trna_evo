#
# This script assumes that you've already run runme.py
#

DATADIR = "DATA" 
SUMMARYDIR = DATADIR + "/summaries"
MSADIR = DATADIR + "/alignments"

import math, os, re, sys
import scipy as sc

def get_ac_from_name(name):
    #print name, name.__len__()
    return name[ (name.__len__() )-3: ]

def sem(values):
    return ( sc.var(values) / math.sqrt( float(values.__len__()) ) )

def find_acswitch_file(species):
    for f in os.listdir( SUMMARYDIR ):
        if f.__contains__("monophyly"):
            continue
        f = re.sub(".acswitches.txt", "", f)
        if species.__contains__(f):
            ret = SUMMARYDIR + "/" + f + ".acswitches.txt"
            #print "24:", species, ret
            return ret
    
       
def find_alignment_file(species):
    for f in os.listdir( MSADIR ):
        if f.__contains__("reduced"):
            continue
        f = re.sub(".struct.muscle.fasta", "", f)
        if species.__contains__(f):
            ret = MSADIR + "/" + f + ".struct.muscle.fasta"
            #print "37:", species, ret
            return ret


def get_canonical_species(species):
    for f in os.listdir( MSADIR ):
        if f.__contains__("reduced"):
            continue
        f = re.sub(".struct.muscle.fasta", "", f)
        if species.__contains__(f):
            return f

codon_usage_db = sys.argv[1]

maxusage = 0
switch_usages = [] # an array of codon usages for switched trnas
nonswitch_usages = [] # an array of codon usages for trnas that did not switch

fout = open(DATADIR + "/summary.codon.usage.txt", "w")
fout.write("Species\tMean codon usage for nonswitched tRNAs\tMean codon usage for switched tRNAs\n")

fin = open(codon_usage_db, "r")
acs = [] # an array of all anti-codon triplets
for ac in acs:
    ac_count[ac] = 0        
lines = fin.readline().split("\r")    
for l in lines:
    if l.startswith("Species"):
        l = l.strip()
        tokens = l.split("\t")
        acs = tokens[1: tokens.__len__()-1 ]
        print acs
    else:
        # scan codon usage for this species. . .
        species = l.split("\t")[0]
        species = re.sub(" ", ".", species)
        
        tokens = l.split("\t")[1:]
        ac_usage = {}
        for i in range(0, acs.__len__()):
            ac_usage[ acs[i] ] = float(tokens[i])
            if float(tokens[i]) > maxusage:
                maxusage = float(tokens[i])
        
        this_switch_usages = [] # codon usage data for this species
        this_nonswitch_usages = []
        
        # scan AC switches for this species. . .
        switched_trnas = []

        acswitchfile = find_acswitch_file( species )        
        if acswitchfile == None:
            print "I can't find AC switch summary for", species
            continue
        if False == os.path.exists( acswitchfile ):
            print "I can't find", acswitchfile
            continue
                    
        switchfin = open(acswitchfile, "r")
        for l in switchfin.readlines():
            if l.startswith("No detected"):
                continue
            if l.__len__() < 2:
                continue
            tokens = l.split()
            #if tokens[0].__contains__("Non"):
            thisac = get_ac_from_name( tokens[3] )
            switched_trnas.append(tokens[3])
            switch_usages.append( ac_usage[thisac] )
            this_switch_usages.append( ac_usage[thisac] )
        switchfin.close()
        
        # scan all ACs for this species. . .
        if switched_trnas.__len__() > 0:
            apath = find_alignment_file( species )
            allfin = open(apath, "r")
            for l in allfin:
                if l.startswith(">"):
                    thislabel = re.sub(">", "", l.strip() )
                    if thislabel not in switched_trnas and thislabel.__contains__("trna"):
                        thisac = get_ac_from_name( thislabel )
                        if thisac in ac_usage:
                            nonswitch_usages.append( ac_usage[thisac] )
                            this_nonswitch_usages.append( ac_usage[thisac] )
            line = get_canonical_species(species) + "\t%.3f"%sc.mean(this_nonswitch_usages) + "\t%.3f"%sem(this_nonswitch_usages) + "\t%.3f"%sc.mean(this_switch_usages) + "\t%.3f"%sem(this_switch_usages)
            fout.write(line + "\n")
            print line
fin.close()
fout.close()
print sc.mean(nonswitch_usages), sc.mean(switch_usages)