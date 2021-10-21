import numpy as np
import h5py
import math
import os
import glob
import re


# def convert(fasta, nloci=None, ns=None, hdf5=None, quiet=False):
#     """
#     DEPRECATED: Expanded compatibility with multi and single line fastas
#     Convert fasta alignments into HFD5 file compatible with ipyrad tools. 
#     It is capable of spliting sequences into multiple loci by size or by Ns chains.
    
#     fasta (path)
#     Fasta file with alignment to convert
    
#     nloci (integer)
#     Number of loci to split secuence.
#     Example: given a a sequence like: ACTGCACGTAGGGTA and nloci=2
#     fasta2hdf5 will consider two loci: ACTGCAC and GTAGGGTA
    
#     ns (integer)
#     Number of Ns that denotate the end of a locus and the begining of another.
#     Example: given a sequence like: ACTGCACGNNNNNNNNNNTAGGGTA and ns=10
#     fasta2hdf5 will consider two loci: ACTGCACG and TAGGGTA
    
#     hdf5 (path)
#     Path to save the resulting HDF5 file. Default=Same location and name as fasta
#     """
    
#     print("""-----------------------------------------------------------------------------------------------------------
#     DEPRECATION NOTICE: This function is deprecated, use instead: split_fasta_to_hdf5() or multiple_fastas_to_hdf5()
#     -----------------------------------------------------------------------------------------------------------
#     """)
    
#     if nloci and ns:
#         raise Exception("Only one mode is allowed not both. 1) arbitrarily split the sequence in N loci or 2) Ns as locus separator")
#     elif not nloci and not ns:
#         raise Exception("Define the method to delimitate loci from sequences with nloci OR ns")
        
#     #define default hdf5 path
#     if not hdf5:
#         path = os.path.dirname(fasta)
#         base = os.path.basename(fasta)
#         name = os.path.splitext(base)[0]
#         hdf5 = os.path.join(path, name + ".hdf5")
        
#     with open(fasta) as f:
        
#         phynames = []
#         phy = []

#         for idx, line in enumerate(f):
#             #if line is a header extract the first element before space as name
#             if line[0] == ">":
#                 phynames.append(line.split()[0][1:])

#             # else extract the sequence info
#             else:
#                 #Mode arbitrary n loci
#                 if nloci:
#                     # if is the first sequence create phymap and scaffold dataset
#                     if idx == 1:
                       
#                         # create empty arrays
#                         phymap = []
#                         scaffold_names = []
#                         scaffold_lengths = []

#                         # get length
#                         length = len(line.strip())

#                         ## if nloci is provided 
#                         if nloci > length:
#                             raise Exception(f"Impossible to get the number of loci requested ({nloci}), the number is larger than bases in the alignment ({length})")

#                         length_loci = length / nloci
#                         adjusted_length_loci = math.floor(length_loci)

#                         # split each locus one by one
#                         for idx_locus in range(nloci):
#                             start = idx_locus*adjusted_length_loci
#                             end = start + adjusted_length_loci


#                             # if length is not divisible, include the remainder bases in the last locus
#                             if idx_locus == range(nloci)[-1] and end < length:
#                                 end += length%nloci


#                             # fill phymap, scaffold_lengths, and scaffold_names 
#                             phymap.append([idx_locus + 1, start, end, 0, end])
#                             scaffold_names.append(f"loc-{idx_locus + 1}")
#                             scaffold_lengths.append(end-start)

#                     # prepare phy, for now add sequence by sequence to the file. 
#                     phy.append([0 if base in ["N","-","?"] else ord(base) for base in line.strip().upper()])
                    
                    
#                     # certainly this will fill the memory try somethng like:
#                     #   def append(self, values):
#                     #     with h5py.File(self.datapath, mode='a') as h5f:
#                     #         dset = h5f[self.dataset]
#                     #         dset.resize((self.i + 1, ) + shape)
#                     #         dset[self.i] = [values]
#                     #         self.i += 1
#                     #         h5f.flush()

#                 #Mode loci separated by NNNNN chains
#                 if ns:
#                     # if is the first sequence create phymap and scaffold dataset
#                     if idx == 1:

#                         # create empty arrays
#                         phymap = []
#                         scaffold_names = []
#                         scaffold_lengths = []

#                         #get location of loci                       
#                         for idx_locus, locus in enumerate(re.finditer("[^=]+", line.strip().upper().replace("N"*ns,"="))):
#                             start = locus.start() - idx_locus
#                             end = locus.end() - idx_locus

                            
#                             phymap.append([idx_locus + 1, start, end, 0, end])
#                             scaffold_names.append(f"loc-{idx_locus + 1}")
#                             scaffold_lengths.append(end-start)

#                     phy.append([0 if base in ["N","-","?"] else ord(base) for base in line.strip().upper().replace("N"*ns,"")])



#         with h5py.File(hdf5, 'w') as h:
#             h["phy"] = np.asarray(phy, dtype="u1")
#             h["phymap"] = np.asarray(phymap)
#             h["scaffold_names"] = np.asarray(scaffold_names, dtype="S")
#             h["scaffold_lengths"] = np.asarray(scaffold_lengths)
#             h["phymap"].attrs["reference"] = "imported-from-fasta"
#             h["phymap"].attrs["phynames"] = np.asarray(phynames, dtype="S")
#             h["phymap"].attrs["columns"] = [b"chroms", b"phy0", b"phy1", b"pos0", b"pos1",]
        
#         if not quiet: 
#             print(f"HDF5 file saved at: {hdf5}")




            
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
        line = line.rstrip()
        
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
        

def _read_nexus(fp):
    pass


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
    Parameters:
    list_fastas (string)
    Path where fastas are located, it can be something like: "./aligned/*.fna" or a list ["gene1.fa", "gene2.fa, gene3.fa"]
    
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




def nexus_to_hdf5(nexus, output=None, extract_other_info=True, quiet=False, force=False):
    """aaa"""
    #formats: nloci, 
  
    pass



if __name__ == "__main__":
    multiple_fastas_to_hdf5("*.FNA")
