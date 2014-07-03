MIN_BIN_LENGTH = 0.0
MAX_BIN_LENGTH = 1.0
BIN_SIZE =    0.05
TINY = 0.000001

PLOTDIR = "PLOTS"

import math, re, os

# set is an array of floats
def mean(set):
    if set.__len__() == 0:
        return None
    sum = 0.0
    for x in set:
        sum += x
    return sum / float( set.__len__() )

# standard deviation
def sd(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() ) )

# calculates variance
def var(set):
    avg = mean(set)
    if avg == None:
        return None
    sumofsquares = 0.0
    for x in set:
        sumofsquares += (x - avg)**2
    return math.sqrt( float( sumofsquares ) / float( set.__len__() - 1 ) ) 

def stderr(set):
    return (sd(set) / math.sqrt( set.__len__() ) )

def calculate_bins(vals):    
    max_bin = None
    bins = {} # key = floor, value = count
    
    for v in vals:
        if max_bin == None:
            max_bin = v
        elif max_bin < v:
            max_bin = v
    MAX_BIN_LENGTH = max_bin
    
    i = MIN_BIN_LENGTH
    while(i < MAX_BIN_LENGTH):
        bins["%.4f"%i] = 0.0
        i += BIN_SIZE
    
    for b in vals:
        if b < MIN_BIN_LENGTH or b > MAX_BIN_LENGTH:
            continue
        
        bin = MIN_BIN_LENGTH
        while (bin <= b):
            bin += BIN_SIZE
        print bin, BIN_SIZE
        bins["%.4f"%(bin - BIN_SIZE)] += 1
    
    normalized_bins = {}
    for b in bins.keys():
        normalized_bins[b] = float(bins[b]) / vals.__len__()
    return normalized_bins

#
# data[xgroup][series] = value
#
def barplot1(data, xlab, ylab, filekeyword):
    pointsets = data.keys()
    pointsets.sort()
    
    finalset = pointsets[ pointsets.__len__()-1 ]
    tablepath = PLOTDIR + "/barplot.table." + filekeyword + ".txt"
    fout = open(tablepath, "w")
    for p in pointsets:
        if p != finalset:
            fout.write(p.__str__() + "\t")
        else:
            fout.write(p.__str__() )
    fout.write("\n")
    
    maxy = 0.0
    for p in pointsets:
        if data[p] > maxy:
            maxy = data[p]
        if data[p] < TINY:
            data[p] = TINY
        
        if p != finalset:
            fout.write( data[p].__str__() + "\t")
        else:
            fout.write( data[p].__str__() )            
    fout.write("\n")
    fout.close()
    maxy = maxy * 1.2
      
    pdfpath = PLOTDIR + "/barplot." + filekeyword + ".pdf"
    cranstr = "pdf(\"" + pdfpath + "\", width=8, height=4);\n"    
    cranstr += "bars <- read.table(\"" + tablepath + "\", header=T, sep=\"\\t\")\n"
    
    print pointsets

    cranstr += "pointsets <- c("
    for p in pointsets:
        cranstr += (p).__str__() + ","
    cranstr = re.sub(",$", "", cranstr)
    cranstr += ");\n"
    
    cranstr += "barx = barplot(as.matrix(bars), beside=TRUE, col=c('black'), xlab='" + xlab + "', ylab='" + ylab + "', ylim=range(0.0," + maxy.__str__() + "), names.arg=pointsets);\n"
    cranpath = PLOTDIR + "/barplot." + filekeyword + ".cran"
    fout = open(cranpath, "w")
    fout.write( cranstr )
    fout.close()
    
    os.system("r --no-save < " + cranpath)
    
def histo_from_summary(path):
    diff_values = [] # min. distance to tRNA with different anticodon preference
    same_values = [] # min. distance to tRNA with same anticodon preference
    ratio_values = [] # ratio of d_diff / d_same
    fin = open(path, "r")
    for l in fin.xreadlines():
        if l.split().__len__() == 6:
            diff_values.append( float(l.split()[4]) )
            same_values.append( float(l.split()[5]) )
            ratio_values.append( float(l.split()[5]) - float(l.split()[4]) )
    fin.close()
    bins = calculate_bins( diff_values )
    barplot1(bins, "d_diff", "frequency", "all.acswitches.d_diff")
    bins = calculate_bins( same_values )
    barplot1(bins, "d_same", "frequency", "all.acswitches.d_same")
    bins = calculate_bins( ratio_values )
    barplot1(bins, "d_same - d_diff", "frequency", "all.acswitches.s-d")
    
#histo_from_summary("DATA/all.acswitches.txt")