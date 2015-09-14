#This is a version that is meant to work with Fastq files
import Bio
from Bio import SeqUtils
from Bio import SeqIO
import sys


#for records and adapter, should be sys.argv[1 and 2]
fastqfile = sys.argv[1] #This is the input fasta file
adapter = sys.argv[2] #This is the input adapter as a string
keepreads = sys.argv[3] #True or false, this will determine whether or not reads are kept. If true, it will keep reads that do not have the adapter in it. If false, it will get rid of those reads.
removeadapters = sys.argv[4] #True or false, if this is True, the adapters will be removed. If true, removes the adapters from the sequences. If false, it keeps them.
end_defn = sys.argv[5] #If 5, the primer is removed from the 5' end of the sequence. If 3, then it is removed from the 3' end of the sequence.
adapter_name = sys.argv[6] #This is the name of the adapter that you can put into the output text file.


# Here is the command for the test:  python cut_degen.py 'test.fastq' 'GAACWAYWYCT' 'True' 'True'  '5' 'test'

keepreads = str(keepreads)
removeadapters = str(removeadapters)
fastqfile=str(fastqfile)
end_defn = str(end_defn)

fh = open(fastqfile, mode='r+')
len_adapter = len(adapter)
count_adapter_found = 0
count_adapter_not_found = 0
total_seq_count = 0

parsed = SeqIO.parse(fh, format="fastq")

output_fh_name = "output.fastq"

if fastqfile=="test3prime.fastq":
    output_fh_name="output2.fastq"

output_fh = open(output_fh_name, mode='w+')

output_text_name = "output.txt"
if fastqfile=="test3prime.fastq":
    output_text_name="output2.txt"
output_text_fh = open(output_text_name, mode='w+')


for record in parsed:
    try:
        sequence = str(record.seq)
        search = SeqUtils.nt_search(sequence, adapter) #This will search the
        index = int(search[1]) #If it finds the adapter, is the starting index from which it was found.
        adapter_start = index
        adapter_end = index+len_adapter
        count_adapter_found +=1
        total_seq_count+=1
        if removeadapters == "True": #if the value is true, it removes the adapters from the sequences.
            if end_defn=="5":
                record = record[adapter_end:] #If a 5' adapter, you remove adapter from beginning
            elif end_defn=="3":
                record = record[:adapter_start] #If it is a 3' adapter, you remove the adapter at the end
        elif removeadapters == "False": #if the value is false, it does not remove the adapters from the sequences.
            record = record
        SeqIO.write(record, output_fh, format="fastq") #No matter what, write the reads.
    except IndexError:
        count_adapter_not_found+=1
        total_seq_count+=1
        record = record
        if keepreads=="True":
            SeqIO.write(record, output_fh, format="fastq")
        elif keepreads=="False":
            pass
        else:
            pass

output_fh.close()

percent_cut = 100*(float(count_adapter_found)/float(total_seq_count))


output_text_fh.write("The total number of sequences that were analyzed was %i.\n\n"%total_seq_count)
output_text_fh.write("Adapter was found and removed for %i sequences (%i%% of total).\n\n"%(count_adapter_found, percent_cut))

if keepreads =="True":
    output_text_fh.write("Sequences that did not contain the adapter were kept.\n\n")
elif keepreads=="False":
    output_text_fh.write("Sequences that did not contain the adapter were removed from the dataset.\n\n")
if removeadapters=="True":
    output_text_fh.write("The adapters were removed from the dataset.\n\n")
elif removeadapters=="False":
    output_text_fh.write("The adapters were not removed from the dataset.\n\n")
if end_defn=="5":
    output_text_fh.write("Adapters were removed from the 5\' end.\n\n")
elif end_defn=="3":
    output_text_fh.write("Adapters were removed from the 3\'end.\n\n")

output_text_fh.write("The name of the adapter that was removed was named %s, and had the sequence %s.\n\n"%(adapter_name,adapter))
output_text_fh.close()