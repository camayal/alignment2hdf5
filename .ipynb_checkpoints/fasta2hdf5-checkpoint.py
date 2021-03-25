import numpy as np
import h5py
import math
import os


def convert(fasta, nloci=1, ns=None, hdf5=None):
    """
    Convert fasta alignments into HFD5 file compatible with ipyrad tools. 
    It is capable of spliting sequences into multiple loci by size or by Ns chains.
    
    fasta (path)
    Fasta file with alignment to convert
    
    nloci (integer)
    Number of loci to split secuence. Default=1
    Example: given a a sequence like: ACTGCACGTAGGGTA and nloci=2
    fasta2hdf5 will consider two loci: ACTGCAC and GTAGGGTA
    
    ns (integer)
    Number of Ns that denotate the end of a locus and the begining of another.
    Example: given a sequence like: ACTGCACGNNNNNNNNNNTAGGGTA and ns=10
    fasta2hdf5 will consider two loci: ACTGCACG and TAGGGTA
    
    hdf5 (path)
    Path to save the resulting HDF5 file. Default=Same location and name as fasta
    """

#     if nloci and ns:
#         raise Exception("Only one mode is allowed not both. 1) arbitrarily split the sequence in N loci or 2) Ns as locus separator")
#     elif nloci == None and ns == None:
#         raise Exception("Define the method to delimitate loci from sequences with nloci OR ns")
        
    #define default hdf5 path
    if not hdf5:
        path = os.path.dirname(fasta)
        base = os.path.basename(fasta)
        name = os.path.splitext(base)[0]
        hdf5 = os.path.join(path, name + ".hdf5")
        
    with open(fasta) as f:
        
        phynames = []
        phy = []

        for idx, line in enumerate(f):
            #if line is a header extract the first element before space as name
            if line[0] == ">":
                phynames.append(line.split()[0][1:])

            # else extract the sequence info
            else:
                if nloci:
                    # if is the first sequence create phymap and scaffold dataset
                    if idx == 1:
                       
                        # create empty arrays
                        phymap = []
                        scaffold_names = []
                        scaffold_lengths = []

                        # get length
                        length = len(line.strip())

                        ## if nloci is provided 
                        if nloci > length:
                            raise Exception(f"Impossible to get the number of loci requested ({nloci}), the number is larger than bases in the alignment ({length})")

                        length_loci = length / nloci
                        adjusted_length_loci = math.floor(length_loci)

                        # split each locus one by one
                        for idx_locus in range(nloci):
                            start = idx_locus*adjusted_length_loci
                            end = start + adjusted_length_loci


                            # if length is not divisible, include the remainder bases in the last locus
                            if idx_locus == range(nloci)[-1] and end < length:
                                end += length%nloci


    #                         print(f"seq:{line.rstrip()[start:end]}, start: {start}, end: {end}")

                            # fill phymap, scaffold_lengths, and scaffold_names 
                #             print("proto-phymap: ",[idx_locus + 1, start, end, 0, end])
                            phymap.append([idx_locus + 1, start, end, 0, end])
                #             print(f"loc-{idx_locus + 1}")
                            scaffold_names.append(f"loc-{idx_locus + 1}")
                #             print(end-start)
                            scaffold_lengths.append(end-start)

                        #add them to the hdf5 file

                    # prepare phy, for now add sequence by sequence to the file. 
                    phy.append([0 if base in ["N","-","?"] else ord(base) for base in line.strip().upper()])
                    # certainly this will fill the memory try somethng like:
#                       def append(self, values):
#                         with h5py.File(self.datapath, mode='a') as h5f:
#                             dset = h5f[self.dataset]
#                             dset.resize((self.i + 1, ) + shape)
#                             dset[self.i] = [values]
#                             self.i += 1
#                             h5f.flush()


                if ns:
                    pass
                    #get sequence



        with h5py.File(hdf5, 'w') as h:
            h["phy"] = np.asarray(phy)
            h["phymap"] = np.asarray(phymap)
            h["scaffold_names"] = np.asarray(scaffold_names, dtype="S")
            h["scaffold_lengths"] = np.asarray(scaffold_lengths)
            h["phymap"].attrs["reference"] = "imported-from-fasta"
            h["phymap"].attrs["phynames"] = np.asarray(phynames, dtype="S")
            h["phymap"].attrs["columns"] = [b"chroms", b"phy0", b"phy1", b"pos0", b"pos1",]
        
        print(f"HDF5 file saved at: {hdf5}")



    