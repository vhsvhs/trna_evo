"""
This script translates the database spreadsheet into a list of switch.

--in PATH, where PATH is a list of switches, probably all.acswitches.txt

--restrict X, where X can be BAC or EUK or nothing
"""


from argparser import *
ap = ArgParser(sys.argv)

from codons import *

def get_anticodon_from_trna_name(name):
    """Returns the anticodon from the trna name"""
    return name[ name.__len__()-3: ]

def get_aa_from_trna_name(name):
    return name[ name.__len__()-6:name.__len__()-3 ]

def get_switch_list(restrict = None, datatype="codon"):
    fin = open( ap.getArg("--in") )
    switches = [] # array of tuples (from, to, ns/sy)
    lines = fin.readlines()[1:]
    fin.close()
    for l in lines:
        if l.__len__() > 5:
            tokens = l.split()
            type = tokens[2]
            kingdom = tokens[1] 
            if restrict != None:
                if kingdom != restrict:
                    continue           
            if datatype == "codon":
                from_trna = get_anticodon_from_trna_name( tokens[3] )
                to_trna = get_anticodon_from_trna_name( tokens[4] )
                switches.append(  (from_trna, to_trna, type) )
            elif datatype == "aa":
                from_trna = get_aa_from_trna_name( tokens[3] )
                to_trna = get_aa_from_trna_name( tokens[4] )                
                switches.append(  (from_trna, to_trna, type) )
            
    count_ns = 0
    count_sy = 0
    for s in switches:
        if s[2] == "NS":
            count_ns += 1
        elif s[2] == "SY":
            count_sy += 1
    print "\n. For", restrict
    print ". I found", switches.__len__(), "anticodon switches."
    print ".", count_ns, "non-synonymous switches."
    print ".", count_sy, "synonymous switches."
    return switches

def write_switch_list( switches, keyword=None ):
    """Writes a simple spreadsheet listing only the switches and their type (NS or SY).
    This spreadsheet excludes other information, such as species name and phylogenetic distances."""
    if keyword == None:
        keyword = ""
    else:
        keyword = "." + keyword
    fout = open( ap.getArg("--in") + keyword +  ".short.txt", "w")
    for s in switches:
        fout.write( "\t".join(s) + "\n")
    fout.close()

def write_aa_switch_list( switches, keyword=None ):
    """Writes a simple spreadsheet listing only the switches and their type (NS or SY).
    This spreadsheet excludes other information, such as species name and phylogenetic distances."""
    if keyword == None:
        keyword = ""
    else:
        keyword = "." + keyword
    fout = open( ap.getArg("--in") + keyword +  ".aa.short.txt", "w")
    for s in switches:
        fout.write( "\t".join(s) + "\n")
    fout.close()
    
def translate_switch_list( switches ):
    """Translates the switch list to amino acids.
    Returns a matrix of a.a. switches for nonsynonymous switches."""
    ts = []
    for s in switches:
        if s[2] == "NS":
            ts.append(   (CODONS_AA[ s[0] ], CODONS_AA[ s[1] ], s[2])   )
    return ts
    
def count_switch_freqs(switches):
    """switches = array of tuples (from_trna, to_trna, type)"""
    m = {} # key = from tRNA, value = hash, key = to tRNA, value = count
    count_ns = 0
    for s in switches:
        count_ns += 1
        if s[0] not in m:
            m[ s[0] ] = {}
        if s[1] not in m[ s[0] ]:
            m[ s[0] ][ s[1] ] = 0
        m[ s[0] ][ s[1] ] += 1
            
    # Normalize
    #for a in m:
    #    for b in m[a]:
    #        m[a][b] = float(m[a][b]) / count_ns
    return m

def print_codon_freq_matrix( m, keyword=None ):
    KNOWN_CODONS = CODONS_AA.keys()
    KNOWN_CODONS.sort()
    
    if keyword == None:
        keyword = ""
    else:
        keyword = "." + keyword
    fout = open(ap.getArg("--in") + keyword + ".fmat.txt", "w")
    fout.write(".\t")
    fout.write( "\t".join(KNOWN_CODONS) )
    fout.write("\n")
    for c in KNOWN_CODONS:
        fout.write(c + "\t")
        for d in KNOWN_CODONS:
            if c in m:
                if d in m[c]:
                    fout.write( m[c][d].__str__() + "\t")
                else:
                    fout.write("0\t")
            else:
                fout.write("0\t")
        fout.write("\n")
    fout.close()


def print_aa_freq_matrix(m, keyword=None):
    """ m is a translated freq matrix."""
    AAS = []
    for c in CODONS_AA:
        if CODONS_AA[c] not in AAS:
            AAS.append(CODONS_AA[c])
    
    AAS.sort()

    if keyword == None:
        keyword = ""
    else:
        keyword = "." + keyword
    fout = open(ap.getArg("--in") + keyword + ".aa.fmat.txt", "w")
    fout.write(".\t")
    fout.write( "\t".join(AAS) )
    fout.write("\n")
    
    for c in AAS:
        fout.write(c + "\t")
        for d in AAS:
            if c in m:
                if d in m[c]:
                    fout.write( m[c][d].__str__() + "\t")
                else:
                    fout.write("0\t")
            else:
                fout.write("0\t")
        fout.write("\n")
    fout.close() 
        

#################################
#
# main
#
kingdoms = ["EUK", "BAC", None]
for k in kingdoms:
    switches = get_switch_list(restrict = k, datatype="codon")
    write_switch_list(switches, keyword=k)
    m = count_switch_freqs(switches)
    print_codon_freq_matrix(m, keyword=k)
    
    aa_switches = get_switch_list(restrict = k, datatype="aa")
    write_aa_switch_list(aa_switches, keyword=k)
    am = count_switch_freqs(aa_switches)
    print_aa_freq_matrix(am, keyword=k)


        