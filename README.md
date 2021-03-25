# fasta2hdf5

Simple converter for fasta alignments into HDF5 (Hierarchical Data Format) used in some downstream analyses in [ipyrad](https://github.com/dereneaton/ipyrad) and other software like [superBPP](https://github.com/eaton-lab/superbpp).

This converter can split the sequence in similar sized loci (using the parameter `nloci` and the number of loci wanted). 

### Usage as module:
```
import fasta2hdf5
fasta2hdf5.convert("test/simple.fasta", nloci=7)
```

### Usage as CLI script:
```
ToDo
```


Additionally, it can split loci using multiple Ns that separate them (using the parameter `ns` and the number of Ns that separate each loci).

### Usage as module:
```
import fasta2hdf5
fasta2hdf5.convert("test/nchains.fasta", ns=30)
```
### Usage as CLI script
```
ToDo
```

ToDo: Multiple loci in individual fasta files.
