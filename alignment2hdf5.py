import numpy as np
import h5py
import math
import os
import glob
import re

            
def _read_fasta(fp):
    """Function generator that returns headers, sequences, and index one by one
    given a fasta file (single or multiline)
    """
    #define globals
    header = None 
    seq = ""
    index = 0
    
    #iterate line by line in file
    for line in fp:
        #remove returns and new lines
        line = line.strip()
        
        #if line start with > it consider it a header,
        #but if header already exists, it means that previous
        #sequence and header is already in variable, in this case
        #just yield header and sequence and continue.
        #after that populate again header and clean seq
        if line.startswith(">"):
            if header: 
                yield (header, seq, index)
            header = line.split()[0][1:] #just maintain first element in multiinfo headers
            seq =  ""
            index += 1
        
        #if line does not start with > it means that it is just a
        #sequence, so concatenate this with previous seqs (for multiline fastas)
        else:
            seq += line
    
    #at the end of the file (last line), just check if the header is 
    #non empty and yield the last hearder and sequence
    if header: 
        yield (header, seq, index)
        


def _build_hdf5(phy, phynames, phymap, scaffold_names, scaffold_lengths, output):
    """Build hdf5 with the information provided by other functions"""
 
    with h5py.File(output, 'w') as h:
            h["phy"] = np.asarray(phy, dtype="u1")
            h["phymap"] = np.asarray(phymap)
            h["scaffold_names"] = np.asarray(scaffold_names, dtype="S")
            h["scaffold_lengths"] = np.asarray(scaffold_lengths)
            h["phymap"].attrs["reference"] = "converted-with-alignment2hdf5"
            h["phymap"].attrs["phynames"] = np.asarray(phynames, dtype="S")
            h["phymap"].attrs["columns"] = [b"chroms", b"phy0", b"phy1", b"pos0", b"pos1",]
        


def split_fasta_to_hdf5(fasta, number_loci=1, output=None):
    """
    Convert single and multiline fasta alignment into HFD5 file compatible with ipyrad tools and Superbpp. 
    It is capable of spliting sequences into multiple loci.
    
    This function assumes matrix is square (same length for all samples)
    
    Giving the following fasta and setting number_loci=2:
        >sample_1
        ACGGCACGTAAAGTA
        >sample_2
        ACTGCACGTAG
        GGTA
    
    It will create a HDF5 files with 2 loci:
    sample_1: ACGGCAC and GTAAAGTA
    sample_2: ACTGCAC and GTAGGGTA

    Parameters:
    fasta (string)
    Path to fasta file to be converted.

    number_loci (integer)
    Number of final loci that will be produced. Default: 1

    output (string)
    Name of the file and extension of the final produced file. By default it is used the same name of the fasta file.
    """
    # create empty arrays
    phymap = []
    scaffold_names = []
    scaffold_lengths = []
    phynames = []
    phy = []
    
    with open(fasta) as fp:
        
        previous_length = None
        
        for header, seq, index in _read_fasta(fp):
            
            #code base in their unicode form and append to phy. Missing data is encoded as Zero (0) (apparentely this is incorrect)
            # phy.append([0 if base in ["N","-","?"] else ord(base) for base in seq.strip().upper()])
            #code base in their unicode form, but maintaining N as valid char
            phy.append([0 if base in ["-","?"] else ord(base) for base in seq.strip().upper()])
            
            #append name of sequence to phynames
            phynames.append(header)
            
            #detect non-square matrix
            if previous_length and len(seq) != previous_length:
                raise Exception(f"Looks like this matrix is not square, >{header} has a different length")
            previous_length = len(seq)
            
            # if is the first sequence create phymap and scaffold dataset
            #this assumes that matrix is perfectly square
            if index == 1:

                # get length
                length = len(seq.strip())

                ## if number_loci is larger than length set error 
                if number_loci > length:
                    raise Exception(f"Impossible to get the number of loci requested ({number_loci}), the number is larger than bases in the alignment ({length})")

                length_loci = length / number_loci
                adjusted_length_loci = math.floor(length_loci)

                # split each locus one by one
                for idx_locus in range(number_loci):
                    start = idx_locus*adjusted_length_loci
                    end = start + adjusted_length_loci


                    # if length is not divisible, include the remainder bases in the last locus
                    if idx_locus == range(number_loci)[-1] and end < length:
                        end += length%number_loci


                    # fill phymap, scaffold_lengths, and scaffold_names 
                    phymap.append([idx_locus + 1, start, end, 0, end])
                    scaffold_names.append(f"loc-{idx_locus + 1}")
                    scaffold_lengths.append(end-start)

        #create hdf5 file
        #define default hdf5 path if is not provided
        if not output:
            path = os.path.dirname(fasta)
            base = os.path.basename(fasta)
            name = os.path.splitext(base)[0]
            output = os.path.join(path, name + ".hdf5")
        _build_hdf5(phy, phynames, phymap, scaffold_names, scaffold_lengths, output)
        
        


def multiple_fastas_to_hdf5(list_fastas, output="alignment.hdf5"):
    """
    Convert multiple gene alignments (one fasta per gene or locus) into one HFD5 file compatible with ipyrad tools and Superbpp. 

    Giving two fastas for each loci:
    gene1.fna
        >sample_1
        ACGGCAC
        >sample_2
        ACTGCAC
        >sample_3
        ACTGCAA
        
    gene2.fna
        >sample_1
        GTAAAGTA
        >sample_2
        GTAGGGTA
    
    It will create a HDF5 file with the follow information:
    sample_1: ACGGCACGTAAAGTA
    sample_2: ACTGCACGTAGGGTA
    sample_3: ACTGCAANNNNNNNN
    
    Parameters:
    list_fastas (string)
    Path where fastas are located, it can be something like: "./aligned/*.fna" or a list ["gene1.fa", "gene2.fa, gene3.fa"]
    
    output (string)
    Name of the file and extension of the final produced file. Default: "alignment.hdf5"

    """
    #load the first and create a dictionary with sample name as key and start appending genes
    #this is an alignment, so it is expected be square
  
    dataset = {}
    previous_length = 0
    length = 0

    # create empty arrays
    phymap = []
    scaffold_names = []
    scaffold_lengths = []
    phynames = []
    phy = []

    if type(list_fastas) == str:
        fastas = glob.glob(list_fastas)

    # do a initial iteration over fastas to complete list of samples
    all_samples = set()
    for fasta in fastas:
        with open(fasta) as fp:
            for line in fp:
                if line.startswith(">"):
                    all_samples.add(line.split()[0][1:]) #add sample to set 

    # iterate over all fastas to create a unique matrix 
    for index_fasta, fasta in enumerate(fastas):
        #get only name for fasta
        base = os.path.basename(fasta)
        name = os.path.splitext(base)[0]

        with open(fasta) as fp:
            #declare set to catch any sample with the locus in current fasta
            samples_with_locus = set()

            #read fasta 
            for header, seq, index in _read_fasta(fp):
                if header in dataset:
                    dataset[header] += seq
                else:
                    dataset[header] = seq

                #add sample with locus
                samples_with_locus.add(header)


                # if is the first sequence create phymap and scaffold dataset
                #this assumes that matrix is perfectly square
                #each locus can have variable length but all samples in each locus should be the same
                #it is an aligment after all.
                if index == 1:
                    # get length
                    length = len(seq.strip())

                    # calculate start and end point
                    start = previous_length
                    end = previous_length + length
                    previous_length = end

                    # fill phymap, scaffold_lengths, and scaffold_names 
                    phymap.append([index_fasta + 1, start, end, 0, end])
                    scaffold_names.append(name)
                    scaffold_lengths.append(length)


        #when fasta file reach to the end check the missing samples to populate with zeros. 
        missing_samples = all_samples - samples_with_locus
        for sample in missing_samples:
            if sample in dataset:
                dataset[sample] += "N"*length
            else:
                dataset[sample] = "N"*length


    #use dataset to populate phynames and phy
    for sample in dataset:
        phynames.append(sample)
        phy.append([0 if base in ["-","?"] else ord(base) for base in dataset[sample].strip().upper()])
    
    #create hdf5 file
    _build_hdf5(phy, phynames, phymap, scaffold_names, scaffold_lengths, output)




def _read_nexus(nexus):
    """Function generator that returns block and content of the block by lines
    given a nexus file (sequential or interleaved)
    """
    block = None
    keywords = ["charpartition", "taxpartition"]

    for line in nexus:
        if len(line) > 1:
            line = line.strip().lower()

            #Detect start comments 
            if line.startswith("["):
                block = "comment"

            #Detect if comment close, even in the same line or future lines
            #and remove block
            if line.endswith("]"):
                block = None
                continue

            #Skip all info if block is a comment
            if block == "comment":
                continue

            #compatibility with execute command, extract info from other file
            #mostly used to separate command nexus from matrix nexus
            if line.startswith("execute"):
                data_nexus = line.split()[1].replace(";","")
                with open(data_nexus, "r") as fp:
                    for ex_line in _read_nexus(fp):
                        yield ex_line

            ##
            ## Identify blocks
            ##

            #Detect the keyword matrix, enter in seqs block skip to next line
            if line.startswith("matrix"):
                block = "seqs"
                continue
 
            #Detect the keywords enter in proper block
            #and remove keywords, names and operands from actual lines
            for keyword in keywords:
                if line.startswith(keyword):
                    block = keyword
                    line = re.sub(f"{keyword}.*=", "", line)

            ## Yield info inside blocks (removing single semicolons)                
            if block and len(line) > 2:
                yield block, line

            #close block when semicolon is present
            if line.endswith(";"):
                block = None



def nexus_to_hdf5(nexus, output=None, extract_other_info=True):
    """
    Convert nexus alignment into HFD5 file compatible with ipyrad tools and Superbpp.
    
    This function assumes matrix is square (same length for all samples), and requires
    data and charpartition blocks. Other blocks will be ignored.  Data block could have
    a sequential or interleaved DNA matrix. However 'matchchar' parameter is not supported.
    This function can use 'execute' command in the nexus file to load data and charpartition 
    blocks from other nexus files.

    Example of file.nex:
        #NEXUS
        [This is an example of nexus file]

        Begin data;
            Dimensions ntax=6 nchar=48;
            Format datatype=nucleotide gap=- missing=?;
            Matrix
        a1    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
        a2    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
        b1    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
        b2    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
        c1    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
        c2    CTGATTTACATGTCAGATGTTTTTACTAGTTCCCAACAGTTTCTCATG
            ;
        End;

        [charpartition block is requiered]
        charpartition lociset =
        1: 1-10,
        2: 11-20,
        3: 21-30,
        4: 31-40,
        5: 41-48;
        end;


    Parameters:
    fasta (string)
    Path to fasta file to be converted.

    extract_other_info (boolean)
    Extract taxpartition block and convert into an IMAP file. Default: True

    output (string)
    Name of the file and extension of the final produced file. By default it is used the same name of the nexus file.
    """
    dataset = {}
    imap = []

    # create empty arrays
    phymap = []
    scaffold_names = []
    scaffold_lengths = []
    phynames = []
    phy = []
    number_loci = 0

    with open(nexus) as fp:
        for block, line in _read_nexus(fp):
            #fill dataset dict with sequences
            if block == "seqs":
                name_sample, seq = line.split()
                if name_sample in dataset:
                    dataset[name_sample] += seq
                else:
                    dataset[name_sample] = seq
                

            if block == "charpartition":
                # idx_locus + 1, start, end, 0, end]
                number_loci += 1
                locus, range = (i.strip() for i in line.split(":"))
                start, end = re.sub(",|;|\s", "", range).split("-") #remove any remaining separator and split the string
                start = int(start) - 1 
                end = int(end)
                phymap.append([number_loci, start, end, 0, end])
                scaffold_names.append(f"loc-{locus}")
                scaffold_lengths.append(end-start)
                

            if block == "taxpartition" and extract_other_info:
                species, individuals = (i.strip() for i in line.split(":"))
                individuals = re.sub(",|;", "", individuals)
                individuals = individuals.split()
                for individual in individuals:
                    imap.append(f"{species}\t{individual}")


    #use dataset to populate phynames and phy
    for sample in dataset:
        phynames.append(sample)
        phy.append([0 if base in ["-","?"] else ord(base) for base in dataset[sample].strip().upper()])
    
    
    #create hdf5 file
    #define default hdf5 path if is not provided
    if not output:
        path = os.path.dirname(nexus)
        base = os.path.basename(nexus)
        name = os.path.splitext(base)[0]
        output = os.path.join(path, name)
    _build_hdf5(phy, phynames, phymap, scaffold_names, scaffold_lengths, f"{output}.hdf5")

    # save imap if requested
    if extract_other_info:
        with open(f"{output}.popfile.txt", "w") as imap_file:
            imap_file.writelines(f"{i}\n" for i in imap)






if __name__ == "__main__":
    # #test with simple nexus
    # with open("nexus.nex") as fp:
    #     for line in _read_nexus(fp):
    #         print(line)

    # #test with nexus that is linked to a data matrix in another file
    # with open("nexuslink.nex") as fp:
    #     for line in _read_nexus(fp):
    #         print(line)

    # #test with interleaving
    # with open("nexusinterleaved.nex") as fp:
    #     for line in _read_nexus(fp):
    #         print(line)


    nexus_to_hdf5("./dev/nexus.nex")
    # nexus_to_hdf5("nexuslink.nex")
    # nexus_to_hdf5("nexusinterleaved.nex")