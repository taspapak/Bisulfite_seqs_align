#!/usr/bin/python
from itertools import groupby
import commands
import random
import argparse
import os
import sys

#Author:Papakonstantinou Anastasios

#Command-line Parameters
def arg_parsing ():
    
 parser = argparse.ArgumentParser(description='Blast search from reads in Fastq,Fasta or from unmapped reads in map file:\n \
 a)Conversion of input To Fasta.\n \
 b)Selection of random user-defined sequences from the Fasta file.\n \
 c)Writing of the queries in a file and BLAST them  again the specified nt databases(Unconv_nt,FW_Conv_nt,RC_Conv_nt).\n \
 d)Returns a file with the best hits for each query from all the databases.\n \
 ')

 requiredNamed=parser.add_argument_group("required arguments")

 requiredNamed.add_argument('-i','--input', help='Input file name(gz format)',required=True)
 requiredNamed.add_argument('-o','--output',help='Output file name of Blast results', required=True)

 parser.add_argument('-n','--number',   help='Number of sequences to extract(Default=1000)',type=int,required=False)
 parser.add_argument('-w','--wordsize', help='Word_Size(Default=24)',type=int,required=False)
 parser.add_argument('-e','--evalue',   help='E-value threshold(Default=1e-06)',type=float,required=False)

 args = parser.parse_args()
 
 return args


#Check input and converts fastq to fasta directly
def Conver(fastq):

    #If the input file is compressed
       if  fastq[-2:] == "gz":

          #If the input file is in another directory
          if "/" in fastq:
              inp_name=fastq[:-3].rsplit('/', 1)[1]

          else:
              inp_name=fastq[:-3]

       # Handle different inputs

          #Input a fastq
          if inp_name[-5:] == "fastq":
              os.system('zcat ' + fastq + '| paste - - - - ' + '| sed "s/^@/>/g"' + '  | cut -f1-2 | tr "\t" "\n"         > ./' + inp_name)

          #Input a map
          elif inp_name[-3:] == "map":
              os.system('zcat ' + fastq + '| paste -' + '''| awk ' $5=="-" {print ">"$1"\\n"$2}' '''  + ' > ./' + inp_name)

          #Input a fasta
          elif inp_name[-5:] == "fasta":
              os.system('zcat ' + fastq + ' > ./' + inp_name)

          else:
             print "\nThere is no file with the specified name or this in not a correct file format!"
             sys.exit()

       #If the input file is not compressed
       elif fastq[-2:] != "gz":

          #If the input file is in another directory
          if "/" in fastq:
              inp_name=fastq.rsplit('/', 1)[1]

          else:
              inp_name=fastq

          #Input a fastq
          if inp_name[-5:] == "fastq":
             os.system('cat ' + fastq + '| paste - - - - ' + '| sed "s/^@/>/g"' + '  | cut -f1-2 | tr "\t" "\n"         > ./Temp | mv Temp ' + inp_name)

          #Input a map
          elif inp_name[-3:] == "map":
             os.system('cat ' + fastq + '| paste -' + '''| awk ' $5=="-" {print ">"$1"\\n"$2}' '''  + ' > ./Temp | mv Temp ' + inp_name)

          #Input a fasta
          elif inp_name[-5:] == "fasta":
              os.system('cat ' + fastq + '> ./Temp | mv Temp ' + inp_name)

          else:
             print "\nThere is no file with the specified name or this in not a correct file format!"
             sys.exit()


       return inp_name


#Gets fasta file from above and parses fasta sequences
def fasta_iter(fasta_name):
     
    fh = open(fasta_name)
    
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header,seq

#Random sequences to be picked up
def pick_random(iterable, samplesize,fastq):

    results = []

    for i, v in enumerate(iterable):
        r = random.randint(0, i)
        if r < samplesize:
            if i < samplesize:
                results.insert(r, v) # add first samplesize items in random order
            else:
                results[r] = v # at a decreasing rate, replace random items

    if len(results) < samplesize:
        print "\nYou specified larger number that the existing number of sequences that is " + str(len(results)) + "!\n"
        #os.system("rm "+fastq)
        sys.exit()

    #If the input if a fastq (check last char)
    if fastq[-1] == "q":
         rand_name = str(samplesize)+"Random_" + fastq[:-6] + ".fasta"

    #If the input if a map (check last char)
    elif fastq[-1] == "p":
         rand_name = str(samplesize)+"Random_" + fastq[:-4] + ".fasta"

    #If the input if a fasta (check last char)
    elif fastq[-1] == "a":
         rand_name = str(samplesize)+"Random_" + fastq    

    f=open(rand_name,"w")
    
    for i in results:
        f.write(">")

        f.write(str(i[0])+"\n")
        f.write(str(i[1])+"\n")
    
    f.close()

    os.system("rm " +fastq)

    return rand_name
    
#Blast commands 
def blast(my_query,eval,wsize):

 tf = open ("task.file","w")

  
 #Blast against nt databases
 tf.write(
"/apps/BLAST+/2.2.28/bin/blastn -task megablast -db /project/production/apapakon/nt/DB/Unconverted_nt_db/Unconverted_nt_db   -outfmt '6 qseqid stitle qcovhsp pident bitscore evalue' -query " + my_query + " -max_target_seqs 1 -word_size " + str(wsize) + " -reward 1 -penalty -2   -evalue  " + str(eval) + "  -num_threads 8 -out Unconv_Res -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n"

"/apps/BLAST+/2.2.28/bin/blastn -task megablast -db /project/production/apapakon/nt/DB/FW_Converted_nt_db/FW_Converted_nt_db -outfmt '6 qseqid stitle qcovhsp pident bitscore evalue' -query " + my_query+ "  -max_target_seqs 1 -word_size " + str(wsize) + " -reward 1 -penalty -2   -evalue  " + str(eval) + "  -strand plus -num_threads 8 -out FW_conv_Res -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n"

"/apps/BLAST+/2.2.28/bin/blastn -task megablast -db /project/production/apapakon/nt/DB/RC_Converted_nt_db/RC_Converted_nt_db -outfmt '6 qseqid stitle qcovhsp pident bitscore evalue'  -query " + my_query + " -max_target_seqs 1 -word_size " + str(wsize) +" -reward 1 -penalty -2  -evalue  " + str(eval) + "  -strand plus -num_threads 8 -out RC_conv_Res -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n")
 
 tf.close()
 
 #Blast against Homo Sapiens databases
 '''
 tf.write(
  
  "/apps/BLAST+/2.2.28/bin/blastn -task megablast -db /project/production/apapakon/Homo_Sapiens/DB/Unconverted_db/Unconv_HS_Genome       -outfmt '6 qseqid stitle qcovhsp pident sstrand  bitscore evalue'  -query " + my_query + " -max_target_seqs 1  -word_size " + str(wsize) + "  -reward 1 -penalty -2  -evalue " + str(eval) + " -strand plus  -num_threads 4 -out Unconv_Res  -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n"

 "/apps/BLAST+/2.2.28/bin/blastn  -task megablast -db /project/production/apapakon/Homo_Sapiens/DB/UNMETH_Conv_db/Unmeth_HS_Genome       -outfmt '6 qseqid stitle qcovhsp pident sstrand  bitscore evalue'  -query " + my_query+ "  -max_target_seqs 1  -word_size " + str(wsize) + "  -reward 1 -penalty -2  -evalue " + str(eval) + "  -num_threads 4 -out FW_conv_Res -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n"

 "/apps/BLAST+/2.2.28/bin/blastn -task megablast -db /project/production/apapakon/Homo_Sapiens/DB/METH_Conv_db/Meth_HS_Genome            -outfmt '6 qseqid stitle qcovhsp pident sstrand  bitscore evalue'  -query " + my_query + " -max_target_seqs 1  -word_size " + str(wsize) + "  -reward 1 -penalty -2  -evalue " + str(eval) + "  -num_threads 4 -out RC_conv_Res -best_hit_overhang 0.25  -best_hit_score_edge 0.05\n")
 tf.close()

 '''

 f = open("megablast.cmd","w")

 f.write("#!/bin/bash\n"
 "# @ job_name         = doblast\n"
 "# @ initialdir       = .\n"
 "# @ class            = lowprio\n"
 "# @ output           = /project/production/apapakon/Program/stout/blast_out1000\n"
 "# @ error            = /project/production/apapakon/Program/stout/blast_err1000 \n"
 "# @ total_tasks      = 3\n"
 "# @ wall_clock_limit = 04:00:00\n"
 "## number of openmp processes\n"
 "# @ cpus_per_task = 8\n"
 "# @ tasks_per_node = 1\n"

 "\n"

 "date1=$(date +'%s')\n"

 "/apps/GREASY/latest/bin/greasy task.file\n"

 "date2=$(date +'%s')\n"

 "diff=$(($date2-$date1))\n"

 "echo $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed.\n")

 f.close()

 string = commands.getoutput("mnsubmit megablast.cmd")
 print "\n"
 print (string+ "(Blast jobs)")

 return string.rsplit(" ",1)[1]

#Joins results from the 3 databases,sorts them by query name and evalue and finally keeps the best hits for each query
def manage_res(blast_res,blast_id,my_query):

 f = open("filtering.cmd","w")

 f.write( "#!/bin/bash\n"
          "# @ job_name         = filter\n"
          "# @ initialdir       = .\n"
          "# @ output           = /project/production/apapakon/Program/stout/filt_blast_out1000\n"
          "# @ error            = /project/production/apapakon/Program/stout/filt_blast_err1000 \n"
          "# @ total_tasks      = 1\n"
          "# @ wall_clock_limit = 01:00:00\n"
          "## number of openmp processes\n"
          "# @ cpus_per_task = 1\n"
          "# @ tasks_per_node = 1\n"
          "\n"

          "awk '!x[$1]++' Unconv_Res  > un_tophits\n"
          "awk '!x[$1]++' FW_conv_Res > fw_tophits\n"
          "awk '!x[$1]++' RC_conv_Res > rc_tophits\n"


          "cat un_tophits fw_tophits rc_tophits > all_tophits\n"
          "awk '{print $NF,$0}' all_tophits | sort -k1,1g | sort --stable -k2,2 | cut -f2- -d' ' > All_sorted\n"
          "awk '!x[$1]++' All_sorted | awk '{print $NF,$0}' | sort -k1,1g | cut -f2- -d' ' > "+blast_res +"\n"
          "sed -i '1i# Query:"+my_query+ "' "+blast_res+ "\n"
          "sed -i '2i# Databases:FW_Conv_nt,RC_Conv_nt,Unconv_nt\\n' " +blast_res+ "\n"
          "sed -i '4i# Fields:query id,subject title,%hsp coverage,%identity,bitscore,evalue' "+blast_res+ "\n"
          "awk '$1 == \"#\" {print $1}' "+ blast_res + " | awk '{if ($2==\"PREDICTED:\" || $2==\"TPA:\") print $2,$3,$4; else print $2,$3;}'  "+ blast_res + "  > Scripts/Total_Organisms\n"
          
          #For nt (FW VS RC)
          "awk '$1 == \"#\" {print $1}' "+ blast_res + " | grep \"FW_Conv\" "+ blast_res + " | awk '{if ($2==\"PREDICTED:\" || $2==\"TPA:\") print $2,$3,$4; else print $2,$3;}'  > Scripts/FW_Organisms\n"
          "awk '$1 == \"#\" {print $1}' "+ blast_res + " | grep \"RC_Conv\" "+ blast_res + " | awk '{if ($2==\"PREDICTED:\" || $2==\"TPA:\") print $2,$3,$4; else print $2,$3;}'  > Scripts/RC_Organisms\n"

          #For Human(Meth VS Umeth)
          #"awk '$1 == \"#\" {print $1}' "+ blast_res + " | grep \"_METH\"   "+ blast_res + " | awk '{if ($2==\"PREDICTED:\" || $2==\"TPA:\") print $2,$3,$4; else print $2,$3;}'  > Scripts/FW_Organisms\n"
          #"awk '$1 == \"#\" {print $1}' "+ blast_res + " | grep \"_UNMETH\" "+ blast_res + " | awk '{if ($2==\"PREDICTED:\" || $2==\"TPA:\") print $2,$3,$4; else print $2,$3;}'  > Scripts/RC_Organisms\n"
          "Scripts/./Align_Report.py \n"

          "rm un_tophits fw_tophits rc_tophits all_tophits\n"
          "rm Unconv_Res FW_conv_Res RC_conv_Res\n")
          #"rm task.file greasy*\n")


 f.close()


 string=commands.getoutput("mnsubmit -dep afterok:"+blast_id+" filtering.cmd"  )
 print (string+ "(Filtering of Results(Depends on blast)")

 os.system("rm filtering.cmd")
 os.system("rm megablast.cmd")

#Call methods
def main():
    
    arguments=arg_parsing()

    file_input=arguments.input
    blast_res=arguments.output
    num_seqs=arguments.number
    e_value= arguments.evalue
    word_size = arguments.wordsize

    if num_seqs is None:
      num_seqs=1000

    if e_value  is None:
      e_value=1e-06

    if word_size is None:
      word_size=24


    fasta=Conver(file_input)
    sequences=fasta_iter(fasta)
    query=pick_random(sequences,num_seqs,fasta)

    blast_job_id=blast(query,e_value,word_size)
    manage_res(blast_res,blast_job_id,query)


if __name__ == "__main__":
    main()

