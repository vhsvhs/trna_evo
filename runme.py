#########################################################
#
# tRNA anti-codon switching analysis
#
# To what extent is anti-codon switching prevalent
# across the tree of life?
#
# April 2013
# Victor Hanson-Smith
# victor.hanson-smith@ucsf.edu
#
############################################################
#
# To use this script, follow these instructions:
#
# 1. Install the following software:
#   a. tRNAscan-SE, for structurally aligning tRNA sequences
#   b. Muscle, for aligning sequences
#   c. RAxML, for inferring maximum likelihood phylogenies
#   d. DendroPy, a python package for phylogenetics
#   e. (optional) mpirun, part of the message-passing interface (MPI) suite.
#   f. (optional) mpi_dispatch, a lightweight program for distributing 
#      jobs across nodes in a distributed-memory multiprocessor.
#
# 2. Configure the global parameters in this script (see below)
#    to point to the locations of software on your machine.
#
# 3. Invoke this script, as follows:
#   
#    $> python runme.py PATH
#
#    . . . where PATH points to a fasta-formatted file containing 
#    the database of tRNA sequences.
#    It is important that the sequence names in this fasta file
#    adhere to the naming convention used by the Genomic
#    tRNA database (http://gtrnadb.ucsc.edu).
#
#    Specifically, use the following format:
#    >Genus_species.<trna name>-AlaTGC (start-stop) Ala (TGC)
#
#    An good example, from http://gtrnadb.ucsc.edu:
#    >Acaryochloris_marina_MBIC11017_chr.trna66-AlaTGC (1408891-1408816)  Ala (TGC) 76 bp  Sc: 99.24
#
# 4. Collect the output from this script.  This script may take
#    several hours to complete.  Output will be written to the
#    directory you specify below (in the global parameters).
#
# OUTPUT
#
# For each species in the database, this script produces 
# the following output files:
#
# << documentation needed here! >>
#
#
# ARCHITECTURE OF THIS SCRIPT
#
# After these comments, the script is divided into the following section:
# 1. Global parameters, configure for your system
# 2. Data structures and primitive functions
# 3. "Main" which unites all the primitives into an analysis pipeline.
#
##########################################################

#######################################################
#
# ---> Configure these global parameters for your system. . .
#
LOGDIR = "LOGS"
DATADIR = "DATA"         # All output data will be written to this folder.
TEMPDIR = "TEMP"
SEQDIR = DATADIR + "/raw_sequences"
RAXMLDIR = DATADIR + "/raxml_output"
TREEDIR = DATADIR + "/trees"
STRUCTDIR = DATADIR + "/tRNAscan_output"
MSADIR = DATADIR + "/alignments"
SUMMARYDIR = DATADIR + "/summaries"

ALLDIRS = [LOGDIR, DATADIR, TEMPDIR, SEQDIR, RAXMLDIR, TREEDIR, STRUCTDIR, MSADIR, SUMMARYDIR]

TRNASCAN = "trnascan-SE" # Points to the executable: http://selab.janelia.org/tRNAscan-SE/
PSUEDOGENE_CUTOFF = 40

MUSCLE = "muscle"    # Points to the executable for Muscle, a program for multiple sequence alignment
RAXML = "raxml -T 2" # Points to the exectuable for RAxML, a program for ML phylogenetic inference

USE_MPI = True  # if False, then the following two variables can be ignored. . .
MPIRUN = "mpirun -np 29 --machinefile hosts.txt"
MPIDISPATCH = "/common/bin/mpi_dispatch"

#######################################################
#
# You shouldn't need to modify anything below here. . . 
#
import os, sys, re
import cPickle as pickle
from dendropy import Tree, treecalc

#
# Global data structures that will be filled with data
# during the analysis. . .
#

#
species_trna_seq = {} # key = species, value = hashtable: key = tRNA name, value = tRNA sequence (cleaned, with anti-codon changed to NNN).  This hashtable contains only unique sequences; redundant sequences are discarded rather than added to this hashtable.
species_alltrnanames = {} # key = species, value = list of all tRNA names in the databse for this species.  Note that some tRNAs will not be in species_trna_seq because their sequences are redundant. 
species_trna_dups = {} # [species][trna name] = a list of other tRNA names that have identical sequences to this tRNA.
species_countreject = {} # key = species, value = count of rejected tRNAs from database.
species_nscount = {} # key = species, value = count of putative nonsynonymous anticodon-switching events.
species_scount = {} # key = species, value = count of putative synonymous anticodon-switching events.
species_ac_nscount = {} # key = species, value = hash; key = codon, value = count of anticodon shifts to this codon
species_ac_scount = {} # same as species_ac_nscount, but count nonsynonymous anticodon shifts only.
ac_count = {} # key = anticodon, value = count of tRNAs among all species.

##################################

def pickle_globals():
    p = [species_trna_seq, species_alltrnanames, species_trna_dups, species_countreject, species_nscount, species_scount, species_ac_nscount, species_ac_scount, ac_count]
    pickle.dump( p, open( TEMPDIR + "/save.p", "w" ) )
    
def unpickle_globals():
    p = pickle.load( open( TEMPDIR + "/save.p", "r" ) )
    species_trna_seq = p[0]
    species_alltrnanames = p[1]
    species_trna_dups = p[2]
    species_countreject = p[3]
    species_nscount = p[4]
    species_scount = p[5]
    species_ac_nscount = p[6]
    species_ac_scount = p[7]
    ac_count = p[8]

def remove_problem_chars(x):
    #print x
    """Removes characters within taxa names that are problematic."""
    chars = ["\)", "\(", "\[", "\]", "\:"]
    for c in chars:
        x = re.sub(c, "", x)
    return x

def line_to_species(line):
    """Input: fasta taxa line from Genomic tRNA database all-trnas.fa.  Output: species name"""
    line = line.strip()
    s = re.sub(">", "", line)
    tokens = s.split()[0].split(".")[0].split("_")
    species = ""
    for i in range(0, tokens.__len__()):
        t = tokens[i]
        if t.__contains__("scaffold") or t.__contains__("Scaffold"):
            break
        elif t.startswith("chr") or t.startswith("Chr"):
            break
        elif t.__contains__("Contig") or t.__contains__("contig"):
            break
        elif t.startswith("GL") or t.startswith("ABH") or t.startswith("ADFV"):
            break
        elif t.__contains__("plasmid"):
            break
        elif t.startswith("JH") or t.startswith("ACBE") or t.startswith("WH") or t.startswith("NA") or t.startswith("ACF") or t.startswith("Ultra") or t.startswith("ultra"):
            break
        if t.__contains__("random"):
            break
        if i >=2 and t == tokens[0]:
            break
        alldigits = True
        for c in t:
            if c.isalpha():
                alldigits = False
                break
        if alldigits:
            break
        species += t + "."
    species = species[0: species.__len__()-1 ]
    species = remove_problem_chars(species)
    return species

def line_to_name(line):
    """Input: fasta taxa line from Genomic tRNA database all-trnas.fa. Output: name of tRNA"""
    line = line.strip()
    s = re.sub(">", "", line)
    name = s.split()[0].split(".")[1]
    name = remove_problem_chars( name )
    return name

def line_to_anticodon(line):
    """Input: fasta taxa line. Output: the anti-codon preference of this tRNA."""
    tokens = line.split()
    ac = tokens[3]
    ac = re.sub("\(", "", ac)
    ac = re.sub("\)", "", ac)
    return ac

def get_ac_from_name(name):
    #print name, name.__len__()
    return name[ (name.__len__() )-3: ]

def get_aa_from_name(name):
    return name[ (name.__len__() )-6 : (name.__len__())-3]

def is_taxa_good(taxaline):
    """This method rejects tRNA sequences, for a variety of reasons. . ."""
    if taxaline.__contains__("chrM"): # chromosomal
        return  False
    if taxaline.__contains__("???"): # unknown anticodon
        return False
    if float(taxaline.split()[7]) < PSUEDOGENE_CUTOFF: # Sc score < 40, so it's probably a psuedogene.
        return False
    return True

def split_and_clean_database(path):
    # path points to a FASTA-formatted file of tRNA sequences, with sequence name formats
    # used by http://gtrnadb.ucsc.edu
    
    fin = open(path, "r")
    print "\n. OK, I found the tRNA database at", path
    print ". I'm reading the tRNA sequences.  This may take a while . . ."
    currtaxa = None
    currseq = ""
    lines = fin.readlines()
    fin.close()

    lout = open(LOGDIR + "/" + path + ".rejects", "w")

    i = -1
    while(i < lines.__len__()-1):
        i += 1
        line = lines[i]
        line = line.strip()
        #for line in fin.xreadlines():
        if line.startswith(">"): # then we found a sequence title
            species = line_to_species(line)
            if species not in species_countreject:
                species_countreject[species] = 0
            
            if False == is_taxa_good(line):
                lout.write(line + "\n")
                species_countreject[species] += 1
            if is_taxa_good(line):
               species = line_to_species(line)
               if species.__len__() < 1: # skip empty line species.
                   continue
               trna = line_to_name(line)
               thisac = line_to_anticodon(line)
               if thisac in ac_count:
                   ac_count[thisac] += 1
               else:
                   ac_count[thisac] = 1
               seq = ""
               j = i+1
               while (j < lines.__len__()):
                   if False == lines[j].startswith(">"):
                       seq += lines[j].strip()
                   else:
                       break
                   j += 1
               #print species, trna

               if species not in species_alltrnanames:
                   species_alltrnanames[species] = []
               species_alltrnanames[species].append( trna )

               if species not in species_trna_seq:
                   species_trna_seq[species] = {}

               if species not in species_trna_dups:
                   species_trna_dups[species] = {}

               # Is this tRNA sequence unique for this species?
               # ...examine all the sequences we currently have for this species:
               found_dup = False

               for n in species_trna_seq[species]:
                   if species_trna_seq[species][n] == seq and get_ac_from_name(n) == thisac:
                       found_dup = True
                       # record the duplicate:
                       #if n not in species_trna_dups[species]:
                       #    species_trna_dups[species][n] = []
                       #print name, "is a dup of", n
                       species_trna_dups[species][n].append( trna )
                       break
               # If yes, save this tRNA sequence
               if found_dup == False:
                   species_trna_seq[species][trna] = seq
                   #print name, "is not a dup."
                   species_trna_dups[species][trna] = []
                   #print species, name, currseq
    lout.close()

    for species in species_trna_seq:
        print "\n. OK, I found", species_alltrnanames[species].__len__(), "total,", species_trna_seq[species].__len__(), "unique tRNA sequences in", species

def write_fasta_for_species(species):
    """This method writes the contents of species_trna_seq for species
    to a FASTA-formatted file."""
    species_n_t = {} # a copied version of species_trna_seq, but with redundant tRNAs removed.
    fastapath = SEQDIR + "/" + species + ".fa" 
    fout = open(fastapath, "w")
    for name in species_trna_seq[species]:
        fout.write(">" + name + "\n")
        fout.write(species_trna_seq[species][name] + "\n")
    fout.close()
    return fastapath
    
def write_trnascan_commands(species_fasta):
    commands = []
    species_list = species_fasta.keys()
    species_list.sort()
    for species in species_list:
        infile = STRUCTDIR + "/" + species + ".struct.txt "
        os.system("rm " + infile)
        c = TRNASCAN + " -f " + infile + " " + species_fasta[species]
        commands.append(c)
    spath = "trnascan.commands.sh"
    fout = open(spath, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    return spath

def trnascan_to_fasta(species):
    """Converts the output from tRNA-Scan-SE to a FASTA file."""
    opath = STRUCTDIR + "/" + species + ".struct.txt"
    if False == os.path.exists(opath):
        return
    fin = open(opath, "r")
    lines = fin.readlines()
    fin.close()
    name_seq = {}
    currname = None
    currseq = None
    acstart = 0
    acstop = 0
    for l in lines:
        if l.__contains__("Length:"):
            currname = l.strip().split(".")[0]
        elif l.startswith("Typ"):
            acstart = int(l.split()[5].split("-")[0])
            acstop = acstart + 3
        elif l.startswith("Seq:"):
            rawseq = l.strip().split()[1]
            editseq = rawseq[0:(acstart-1)] + "NNN" + rawseq[(acstop-1):]
            #print currname, "replaced", rawseq[acstart-1:acstop-1], "with NNN"
            seq = ""
            for c in editseq:
                if c.isupper():
                    seq += c
            name_seq[currname] = seq
            currname = None
            acstart = 0
            acstop = 0
    fout = open(STRUCTDIR + "/" + species + ".struct.fasta", "w")
    for name in name_seq:
        fout.write(">" + name + "\n")
        fout.write(name_seq[name] + "\n")
    fout.close()


def run_muscle():
    commands = []
    for species in species_trna_seq:
        fapath = STRUCTDIR + "/" + species + ".struct.fasta"
        if False == os.path.exists(fapath):
            continue
        command = MUSCLE + " -in " + fapath
        command += " -out " + MSADIR + "/" + species + ".struct.muscle.fasta"
        commands.append( command )
    spath = "muscle_commands.sh"
    fout = open(spath, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    #exit()
    print "\n. OK, I'm aligning the tRNA sequences, using the commands in", spath
    if USE_MPI:
        os.system(MPIRUN + " " + MPIDISPATCH + " " + spath)
    else:
        os.system("source " + spath)

def run_raxml():
    commands = []
    for species in species_trna_seq:
        msapath = MSADIR + "/" + species + ".struct.muscle.fasta"
        if os.path.exists(msapath):
            #print msapath
            command = RAXML + " -m GTRCAT -c 4 -n " + species + " -s " + msapath + " -p 12345"
            commands.append(command)
        else:
            print "\n. I can't find", msapath, " -- I'm skipping it."
    spath = "raxml_commands.sh"
    fout = open(spath, "w")
    for c in commands:
        fout.write(c + "\n")
    fout.close()
    #exit()
    print "\n. OK, I'm running RAxML, using the commands in", spath
    if USE_MPI:
        os.system(MPIRUN + " " + MPIDISPATCH + " " + spath)
    else:
        os.system("source " + spath)


def count_trna_types(species):
    trna_count = {}
    for trna in species_alltrnanames[species]:
        ac = get_ac_from_name(trna)
        if ac not in trna_count:
            trna_count[ac] = 1
        else:
            trna_count[ac] += 1
    return trna_count

def pretty_print_trees():
    print "\n. OK, I'm reformatting the RAxML results for nice printing..."
    """Reformats the phylogeny, such that each taxon label looks like this:
    trna12-AlaTCT[6/7]
    . . . where 6 is the number of sequences collapsed into this sequence, and 7 is the number of total tRNAs in the databse."""
    species_list = species_trna_seq.keys()
    species_list.sort()
    for species in species_list:
        #print species_trna_dups[species]
        treepath = RAXMLDIR + "/RAxML_result." + species
        if False == os.path.exists( treepath ):
            continue
        newtreepath = TREEDIR + "/" + species + ".tree"
        t = Tree()
        t.read_from_path(treepath, "newick")
        print " -->", treepath
        trna_count = count_trna_types(species)
        #print trna_count
        newts = t.__str__()
        for taxon in t.taxon_set:
            #print "372:", taxon.label
            thisac = get_ac_from_name(taxon.label)
            count_this_type = trna_count[thisac]
            count_dups = species_trna_dups[species][taxon.label].__len__() + 1
            if count_dups <= 1:
                count_dups = ""
            else:
                count_dups = "(" + count_dups.__str__() + ")"

            newts = re.sub( taxon.label, (taxon.label + count_dups + "[" + count_this_type.__str__()+ "]"), newts)
        fout = open(newtreepath, "w")
        fout.write( newts + "\n" )
        fout.close()


def debug411(species):
    """This methos id depricated."""
    # for debugging:
    countns = 0
    counts = 0
    for ac in species_ac_nscount[species]:
        countns += species_ac_nscount[species][ac]
    for ac in species_ac_scount[species]:
        counts += species_ac_scount[species][ac]
    print "491: verify: ns =", countns, "s=", counts, "sum=", (countns + counts)
    # end debugging

def asses_monophyly(t):
    """t is a DendroPy Tree."""
    """This function returns a hashtable, where key = anticodon preference X,
    value = the number of tRNAs with a.c. other than X that must be invoked to make
    the X clade monophyletic. """
    
    t.is_rooted = False
    t.update_splits()
    
    # First, sort the leaf nodes by their anticodon preference.
    ac_labels = {} # key = a.c., value = list of Node objects
    for i, t1 in enumerate(t.taxon_set):
        thisac = get_ac_from_name( t1.label )
        if thisac not in ac_labels:
            ac_labels[ thisac ] = []
        ac_labels[ thisac ].append( t1.label )
    
    # Next, find the MRCA for each set of a.c. nodes
    for ac in ac_labels:
        mrca = t.mrca(taxon_labels=ac_labels[ac])
        print ac, mrca

    # to-do: count the number of nodes descendant from mrsa.

def find_anticodon_switches():
    print "\n. OK, I'm searching for switched anticodons. . ."
    species_list = species_trna_seq.keys()
    species_list.sort()
    #print "504:", species_list
    allpath = DATADIR + "/all.acswitches.txt"
    allout = open(allpath, "w")
    allout.write("Species\tswitch type\testimated from\tto\tto ac\tto aa\td\n")
    #
    # FOR EACH SPECIES. . . 
    #
    for species in species_list:
        print species
        rpath = SUMMARYDIR + "/" + species + ".acswitches.txt"
        treepath = RAXMLDIR + "/RAxML_result." + species
        if os.path.exists(treepath):
            species_nscount[species] = 0
            species_scount[species] = 0
            species_ac_nscount[species] = {}
            species_ac_scount[species] = {}
            fout = open(rpath, "w") # a summary of found ac switches will be written here.
            t = Tree()
            t.read_from_path(treepath, "newick")
            print "\n. Calculating all pairwise distances between sequences on tree:", treepath
            pdm = treecalc.PatristicDistanceMatrix(t) # matrix of pairwise distances between taxa
            
            asses_monophyly(t)
            
            #
            # FOR EACH tRNA SEQUENCE. . .
            # Goal: for each tRNA find the min. distance to another tRNA with the same
            # anticodon, then find the min. distance to another tRNA that is of a different
            # anticodon type.
            print "\."
            for i, t1 in enumerate(t.taxon_set):                
                min2same = None
                min2diff = None
                closest_diff = None # taxon label of closest same-anti-codon tRNA sequence to sequence t1.
                myac = get_ac_from_name(t1.label)
                myaa = get_aa_from_name(t1.label)
                if myaa == "Met":
                    continue
                for t2 in t.taxon_set:
                    if t1 == t2:
                        continue
                    thisac = get_ac_from_name(t2.label)
                    thisaa = get_aa_from_name(t2.label)
                    d = pdm(t1, t2)
                    if myac == thisac:
                        if min2same == None:
                            min2same = d
                        elif min2same > d:
                            min2same = d
                    elif myac != thisac:
                        if min2diff == None:
                            min2diff = d
                            closest_diff = t2.label
                        elif min2diff > d:
                            min2diff = d
                            closest_diff = t2.label
                if min2same == None:
                    min2same = 0.0 # in the event of singletons
                if min2same > min2diff and min2diff != None: # . . . then we've identified an anticodon shift:
                    thataa = get_aa_from_name(closest_diff)
                    if thataa  == myaa: # synonymous shift
                        species_scount[species] += 1
                        fout.write("Synonymous" + " " + closest_diff + " -> " + t1.label + " d:" + min2diff.__str__()  + "\n")
                        allout.write(species + "\tSY\t" + closest_diff + "\t" + t1.label + "\t" + myac + "\t" + myac + "\t%.4f"%min2diff + "\n")                                     
                        print "  . Syn." + " " + closest_diff + " -> " + t1.label + " d:" + min2diff.__str__()
                        if myac not in species_ac_scount[species]:
                            species_ac_scount[species][myac] = 1
                        else:
                            species_ac_scount[species][myac] += 1
                    elif thataa != myaa and thataa != "Met": # nonsynonymous shift
                        species_nscount[species] += 1
                        fout.write("Nonsynonymous" + " " + closest_diff + " -> " + t1.label + " d:" + min2diff.__str__()  + "\n")
                        allout.write(species + "\tNS\t" + closest_diff + "\t" + t1.label + "\t" + myac + "\t" + myac + "\t%.4f"%min2diff + "\n")   
                        print "  . Nonsyn." + " " + closest_diff + " -> " + t1.label + " d:" + min2diff.__str__()
                        if myac not in species_ac_nscount[species]:
                            species_ac_nscount[species][myac] = 1
                        else:
                            species_ac_nscount[species][myac] += 1
            
            if species_nscount[species] == 0 and species_scount[species] == 0:
                fout.write("No detected switched anitcodons for " + species + "\n") 
            fout.close()
            print ".", species, "has", species_nscount[species], "putative nonsynonymously switched anticodons."    
            print ".", species, "has", species_scount[species], "putative synonymously switched anticodons."        
        else:
            print ". I skipped species", species, "because I can't find the ML tree."
    allout.close()


def write_summaries():
    """Writes a tab-seprated text file with a summary of the numer of tRNAs in each species, and 
    the number of putative anticodon switching events."""
    fout = open(DATADIR + "/summary.species.txt", "w")
    header = "Species\tN tRNAs\tN unique\tN rejected\ttRNAs\tN nonsyn. switches\tProportion of total\tN syn. switches\tProportion of total\n"
    fout.write(header)

    species_list = species_trna_seq.keys()
    species_list.sort()
    for species in species_list:
        line = species + "\t"
        # N tRNA in species:
        if species in species_trna_seq:
            line += species_alltrnanames[species].__len__().__str__() + "\t"
        else:
            line += "None\t"

        # N unique tRNAs in species:
        if species in species_trna_seq:
            line += species_trna_seq[species].__len__().__str__() + "\t"
        else:
            line += "None\t"

        # N rejected
        if species in species_countreject:
            line += species_countreject[species].__str__() + "\t"
        else:
            line += "None\t"
        
        # N nonsynonymously switched anticodons in species:
        if species in species_nscount:
            line += species_nscount[species].__str__() + "\t"
        else:
            line += "None\t"       

        # Proportion of total tRNAs, nonsynonymously switched anticodons in species:
        if species in species_nscount:
            line += "%.3f"%(float(species_nscount[species])/species_alltrnanames[species].__len__()) + "\t"
        else:
            line += "None\t"       

        # N synonymously switched anticodons in species:
        if species in species_scount:
            line += species_scount[species].__str__() + "\t"
        else:
            line += "None\t"        

        # Proportion of total tRNAs, synonymously switched anticodons in species:
        if species in species_scount:
            line += "%.3f"%(float(species_scount[species])/species_alltrnanames[species].__len__()) + "\t"
        else:
            line += "None\t"        
        fout.write(line + "\n")
    fout.close()

    """Writes a tab-seperated file with the number of anticodon changes to each codon."""
    all_ac = []
    ac_s = {}
    ac_ns = {}
    for species in species_ac_nscount:
        for ac in species_ac_nscount[species]:
            if ac not in all_ac:
                all_ac.append(ac)
            if ac not in ac_ns:
                ac_ns[ac] = 0
            ac_ns[ac] += species_ac_nscount[species][ac]
        for ac in species_ac_scount[species]:
            if ac not in all_ac:
                all_ac.append(ac)
            if ac not in ac_s:
                ac_s[ac] = 0
            ac_s[ac] += species_ac_scount[species][ac]

    fout = open(DATADIR + "/summary.codons.txt", "w")
    fout.write("Anticodon\t N syn.\t N nonsyn.\n")
    for ac in all_ac:
        sstr = "0"
        if ac in ac_s:
            sstr = ac_s[ac].__str__()
        nsstr = "0"
        if ac in ac_ns:
            nsstr = ac_ns[ac].__str__()
        #print ac, sstr, nsstr
        fout.write(ac + "\t" + sstr.__str__() + "\t" + nsstr.__str__() + "\n") 
    fout.close()

#############################################
#
# Main. . .
#

#
# ** 0. Setup the workspace. . .
#
for d in ALLDIRS:
    if os.path.exists(d) == False:
        os.system("mkdir " + d)

#
# ** 1. Read the database, write an individual FASTA file for each species.
#    This step is required, because it initializes global dictionaries that
#    are used in later steps.
#
split_and_clean_database(sys.argv[1]) # If downloaded from the tRNA database, then sys.argv[1] will be "all-trnas.fa"
pickle_globals()
#exit()
"""
species_fasta = {}
species_list = species_trna_seq.keys()
species_list.sort()
for species in species_list:
    species_fasta[species] = write_fasta_for_species(species)

#
# 2. Use tRNAscan-SE to identify introns and align the tRNA sequences. . .
#
spath = write_trnascan_commands(species_fasta)
print "\n. OK, I'm running tRNAscan.  This may take a while. . ."
if USE_MPI:
    os.system(MPIRUN + " " + MPIDISPATCH + " " + spath)
else:
    os.system("source " + spath)
print "\n. OK, I'm parsing the results from tRNAscan. . ."
for species in species_list:
    trnascan_to_fasta( species )

#
# 3. Re-align the tRNA sequences and build ML phylogenies. . .
#
print "\n. OK, I'm aligning the tRNA sequences. . ."
run_muscle()
run_raxml()
os.system("mv ./RAxML* ./" + RAXMLDIR + "/") # Move the RAxML results into the data folder
"""

# 3b. Reformat the RAxML phylogeny for printing. . .
#pretty_print_trees()
#exit()

#
# 4. *** Scan the ML phylogenies for switched anti-codons. . .
#
#unpickle_globals()
find_anticodon_switches()
write_summaries()
