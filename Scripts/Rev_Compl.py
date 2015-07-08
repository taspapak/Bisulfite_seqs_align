#!/usr/bin/python
import argparse
import sys
import os
from itertools import groupby

#Author:Papakonstantinou Anastasios

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

#Command-line Parameters
def arg_parsing ():

    parser = argparse.ArgumentParser(description='Reverse complement of DNA sequence.')
    parser.add_argument('-i','--input', help='Input file name(Fasta format)',required=True)
    parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()

    return args

#Check if input is fasta
def check(fasta):

     #If the input file is in another directory
     if "/" in fasta:
         inp_name=fasta.rsplit('/', 1)[1]
      
     #If the input file is in the same directory
     else:
         inp_name=fasta

     
     #Input a fasta
     if inp_name[-5:] == "fasta":
        os.system('cat ' + fasta + '> ./Temp | mv Temp ' + inp_name)

     else:
        print "\nThere is no file with the specified name or this in not a correct file format!\n"
        sys.exit()
    
     return inp_name

#Parse fasta sequences
def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    headerlist=[]
    seqlist=[]

    fh = open(fasta_name)
    
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        headerlist.append(header)
        seqlist.append(seq)
    #Returns a list of all sequences
    return headerlist,seqlist

#Get the reverse complement
def reverse_complement(seq):
    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)

    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)

    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)

    return bases

#Apply conversion for all sequences
def reverse_all(whole):
    
    list1=[]
    
    for i in whole:
      conv=reverse_complement(i)
      list1.append(conv)

    return list1


#Split string into equal chunks to write in file with nice indentation
def chunks(l, n):

    chunklist=[]
    
    for i in xrange(0, len(l), n):
        chunklist.append(l[i:i+n])

    return chunklist

#Write converted sequences to file
def write2file(bisul,head,out_f):

    f = open(out_f,"w")

    for seq in range(len(bisul)):

        chunk=chunks(bisul[seq],127)
           
        f.write(">")
        f.write(head[seq]+"\n")
     
        for i in chunk:
           f.write(i+"\n")
       
    print "\nConversion Complete!\n"
    f.close()
   
#Call methods
def main():
    
    arguments=arg_parsing()

    inp_file=arguments.input
    out_file=arguments.output
   
    fasta=check(inp_file) 
    headerlist,seqlist=fasta_iter(fasta)
    revseq=reverse_all(seqlist)
    write2file(revseq,headerlist,out_file)

if __name__ == "__main__":
    main()


