import urllib

from argparser import *
ap = ArgParser(sys.argv)

from plot_tools import *

species_list_path = ap.getArg("--switch_list")

BASE_URL = "http://www.kazusa.or.jp"

KNOWN_CODONS = []

def query_species(species_name):
    tokens = species_name.split()
    name = re.sub(" ", "+", species_name)
    url = BASE_URL + "/codon/cgi-bin/spsearch.cgi?species=" + name + "&c=s"
    #print url
    f = urllib.urlopen(url)
    content = f.read()
    lines = content.split("\n")

    option_url = {} # key = species name of the option, value = URL for data
    option_ngenes = {} # key = species name of the option, value = the number of genes used to calculate this species' codon usage
    
    for ii in range(0, lines.__len__() ):
        l = lines[ii]
        if l.startswith("<I>") and False == l.startswith("<I><A"): # this is a line with a query result...
            this_name = l.split(">")[1].split("<")[0]
            if this_name.startswith("mitochondrion") or this_name.startswith("chloroplast"):
                continue
            this_ngenes = int( l.split("]:")[1].split("<")[0] )
            this_url = BASE_URL + lines[ii-1].split("\"")[1]
            option_url[ this_name ] = this_url
            option_ngenes[ this_name ] = this_ngenes
        
    if option_url.__len__() == 0:
        print "No match for", species_name
        return False
    
    #
    # What is the best option?
    #
    max_ngenes = 0
    best_option = None
    for option in option_ngenes:
        if option_ngenes[option] > max_ngenes:
            max_ngenes = option_ngenes[option]
            best_option = option
    
    if best_option.startswith("mitochondrion") or best_option.startswith("chloroplast"):
        print "No match for", species_name
        return False
            
    print "Match:", species_name, "-->", best_option
    
    return (option_ngenes[best_option], get_codon_usage_from_url( option_url[best_option] ) )
    
    #for option in option_url:
    #    print option, option_url[option], option_ngenes[option]
            
def get_codon_usage_from_url(url):
    """Returns the codon usage table at a specific URL."""    
    codon_usage = {}
    
    f = urllib.urlopen(url)
    content = f.read()
    lines = content.split("\n")
    foundit = False
    for l in lines:
        if l.startswith("</PRE>"):
            foundit = False
            
        if foundit:
            tokens = l.split(")")
            for t in tokens:
                st = t.split()
                if st.__len__() > 2:
                    codon = re.sub(" ", "", st[0])
                    if codon not in KNOWN_CODONS:
                        KNOWN_CODONS.append( codon )
                    usage = re.sub(" ", "", st[1])
                    usage = re.sub("\(", "", usage)
                    codon_usage[codon] = float( usage )
                elif st.__len__() > 1: # the weird case where st[0] and st[1] are fused.
                    if st[0].__contains__("("):
                        #print st
                        partial = re.sub(" ", "", st[0])
                        codon = partial[0:3]
                        if codon not in KNOWN_CODONS:
                            KNOWN_CODONS.append( codon )
                        #print "81:", codon
                        usage = re.sub("\(", "", partial[3:])
                        #print "83:", usage
                        codon_usage[codon] = float( usage )
                    elif st[1].__contains__("("):
                        #print st
                        codon = re.sub(" ", "", st[0])
                        if codon not in KNOWN_CODONS:
                            KNOWN_CODONS.append( codon )
                        #print "81:", codon
                        usage = st[1].split("(")[0]
                        #print "83:", usage
                        codon_usage[codon] = float( usage )                        
    
        if l.startswith("<PRE>"):
            foundit = True
    
    return codon_usage

def get_species_list():
    fin = open(species_list_path, "r")
    lines = fin.readlines()
    lines = lines[1:]
    species_list = []
    for l in lines:
        tokens = l.split()
        species_name  = re.sub("\.", " ", tokens[0])
        species_name = re.sub("-", ". ", species_name )
        if species_name not in species_list:
            species_list.append( species_name )
    return species_list


def write_codon_usage_table(data, ngenes):
    """data[species][codon] = usage"""
    species = data.keys()
    species.sort()
    
    KNOWN_CODONS.sort()

    fout = open("codon_usages.txt", "w")
    fout.write("SPECIES\tNgenes\t")
    fout.write( "\t".join(KNOWN_CODONS) )
    fout.write("\tTotal")
    fout.write("\n")
    
    for sp in species:
        # first sum the usage proportions.
        sum_usage = 0.0
        for c in KNOWN_CODONS:
            if c in data[sp]:
                sum_usage += data[sp][c]
        
        # then write the usages.
        fout.write(sp + "\t")
        fout.write(ngenes[sp].__str__() + "\t")
        for c in KNOWN_CODONS:
            if c in data[sp]:
                fout.write( "%.4f"%( data[sp][c]/1000.0) )
            else:
                fout.write("n/a")
            fout.write("\t")
        fout.write((sum_usage/1000.0).__str__())
        fout.write("\n")
    fout.close()
    
    all_ngenes = []
    for sp in species:
        all_ngenes.append( ngenes[sp] )
    print "\n. I found", species.__len__(), "species."
    print "\n. Mean Ngenes", mean(all_ngenes), "sd=", sd(all_ngenes)

#################################
#
# main

species_list = get_species_list()
no_match_species = []

species_codon_usage = {}
species_ngenes = {}
for sp in species_list:
    x = query_species(sp)
    if x == False:
        no_match_species.append( sp )
    else:
        ngenes = x[0]
        usages = x[1]
        species_codon_usage[sp] = usages
        species_ngenes[sp] = ngenes
        
write_codon_usage_table(species_codon_usage, species_ngenes)
fout = open("codon_skipped_species.txt", "w")
for sp in no_match_species:
    fout.write(sp + "\n")
fout.close()
