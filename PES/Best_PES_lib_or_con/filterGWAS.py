import argparse
import sys
import csv
from copy import copy

def check_args(args=None):
    """Parse arguments"""
    parser = argparse.ArgumentParser(description="Filter a GWAS by a list of genomic regions.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gwas', help="GWAS summary statistics file containing at least the chromosome and position of the variants. Accepted chromosome codes: 1-22 for autosomal chromosomes, X or 23 for the X chromosome, Y or 24 for the Y chromosome and M, MT or 25 for mitochondrial DNA.", required=True)
    parser.add_argument('-s', '--sep', help="Delimiter character for GWAS summary statistics. Default = tab ('\t').", default="\t")
    parser.add_argument('-f', '--filter', help="Tab-separated filter list, containing three columns: chromosome, start position, and end position. Expects 1-based inclusive coordinates, relative to the +/sense strand (i.e. start < end). An optional fourth column including strand information ('+'/'-') can be included, otherwise +/sense strand is assumed.", required=True)
    parser.add_argument('-H', '--noheader', help="Indicate if summary statisitcs lack a header line.", action="store_true")
    parser.add_argument('-c', '--chr', help="Chromosome column number in GWAS file.", required=True)
    parser.add_argument('-b', '--bp', help="Position column number in GWAS file.", required=True)
    parser.add_argument('-u', '--upstream', help="Size of additional upstream region to include. Default = 0.", default="0")
    parser.add_argument('-d', '--downstream', help="Size of additional downstream region to include. Default = 0", default="0")
    parser.add_argument('-o', '--output', help="Output file name for filtered GWAS.", required=True)
    return(parser.parse_args(args))

def filterGwas(gwasFile, sep, noheader, filterFile, chrCol, bpCol, upstream, downstream):
    try:
        c = int(chrCol) - 1
        b = int(bpCol) - 1
        u = int(upstream)
        d = int(downstream)
    except:
        print("ERROR: --chr, --bp, --upstream, and --downstream expect integer arguments.")
        sys.exit()
    
    sep = sep.replace("\\t", "\t")

    filterDict = {}
    for i in range(1, 26):
        filterDict[i] = []
    gwas = []
    with open(filterFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            chrom = line[0]
            try:
                chrom = int(chrom)
            except ValueError:
                if chrom == 'X':
                    chrom = 23
                elif chrom == 'Y':
                    chrom = 24
                elif chrom == 'M' or chrom == 'MT':
                    chrom = 25
                else:
                    print("WARNING: invalid chromosome code encountered. Skipping:", "\t".join(line))
                    continue
            if len(line) == 4:
                if line[3] == '-':
                    tmp = u
                    u = d
                    d = u
                elif line[3] != '+':
                    print("WARNING: invalid strand code encountered. Skipping:", "\t".join(line))
                    continue
            start = int(line[1])
            end = int(line[2])
            if start < end:
                start = start - u
                end = end + d + 1
            else:
                print("WARNING: start position greater or equal to end position. Skipping:", "\t".join(line))
            filterDict[chrom].append(range(start, end))
    
    with open(gwasFile, 'r') as f:
        filteredGwasList = []
        reader = csv.reader(f, delimiter=sep)
        if not noheader:
            next(reader)
        for line in reader:
            chrom = copy(line[c])  # chr column
            try:
                chrom = int(chrom)
            except ValueError:
                if chrom == 'X':
                    chrom = 23
                elif chrom == 'Y':
                    chrom = 24
                elif chrom == 'M' or chrom == 'MT':
                    chrom = 25
                else:
                    print("WARNING: invalid chromosome code encountered. Skipping:", "\t".join(line))
                    continue
            bp = int(line[b])  # bp column
            if any(bp in i for i in filterDict[chrom]):
                filteredGwasList.append(line)
    return(filteredGwasList)


if __name__ == "__main__":
    arguments = check_args(sys.argv[1:])
    newGwas = filterGwas(arguments.gwas, arguments.sep, arguments.noheader, arguments.filter, arguments.chr, arguments.bp, arguments.upstream, arguments.downstream)
    newGwas = ['\t'.join([str(i) for i in j]) for j in newGwas]
    with open(arguments.output, 'w') as o:
      o.write('\n'.join(newGwas) + '\n')
