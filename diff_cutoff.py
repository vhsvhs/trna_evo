import re, os

DDIFF_CUTOFFS = [0.0, 0.02, 0.04, 0.1]

species_countns = {}
species_countsyn = {}
species_nunique = {}

def parse_species():
    fin = open("DATA/summary.species.txt", "r")
    for l in fin.xreadlines():
        tokens = l.split()
        if tokens.__len__() == 8:
            species = tokens[0]
            species_nunique[species] = int(tokens[2])
            species_countns[species] = 0
            species_countsyn[species] = 0
    fin.close()

def parse_summary(path, cutoff):
    fin = open(path, "r")
    for l in fin.xreadlines():
        if l.split().__len__() == 6:
            tokens = l.split()
            species = tokens[0]
            type = tokens[1]
            d_diff = float(tokens[4])
            if d_diff <= cutoff:
                if type == "SY":
                    species_countsyn[species] += 1
                else:
                    species_countns[species] += 1
    fin.close()

parse_species()
for cutoff in DDIFF_CUTOFFS:
    parse_summary("DATA/all.acswitches.txt", cutoff)
    fout = open("DATA/all.acswitches." + cutoff.__str__() + ".txt", "w")
    fout.write("Species\tN_unique\tN_nonsyn.\tN_syn.\n")
    for species in species_countns:
        fout.write(species + "\t" + species_nunique[species].__str__() + "\t" + species_countns[species].__str__() + "\t" + species_countsyn[species].__str__() + "\n")
    fout.close()