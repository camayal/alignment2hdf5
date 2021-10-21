# alignment2hdf5

Simple converter for FASTA and Nexus alignments into HDF5 (Hierarchical Data Format) used in some downstream analyses in [ipyrad](https://github.com/dereneaton/ipyrad) and other software like [superBPP](https://github.com/eaton-lab/superbpp).

## This converter has several modes:

### 1. Convert multiple fasta files in a single folder into a hdf5 file.

Genes or loci must be saved in individual fasta files. They could not include all samples. For example:

`gene1.fna`
```text
>sample_1
ACGGCAC
>sample_2
ACTGCAC
>sample_3
ACTGCAA
``` 
`gene2.fna`
```text
>sample_1
GTAAAGTA
>sample_2
GTAGGGTA
```

#### Usage as module:
```python
import alignment2hdf5
alignment2hdf5.multiple_fastas_to_hdf5("./test/genes/*.FNA", output="./test/alignment.hdf5")
```

#### Usage as CLI script:
```
ToDo
```

### 2. Split a fasta file into multiple loci having the same length and convert it into a hdf5 file.
#### Usage as module:

Fasta files can be single-lined or multi-lined (interleaved), for example:

`simple.fa` 
```text
>sample_1 single line
ACGGCACGTAAAGTA
>sample_2 multiline
ACTGCACGTAG
GGTA
```

```python
import alignment2hdf5
alignment2hdf5.split_fasta_to_hdf5("./test/simple.fasta", number_loci=4, output="./test.simple.hdf5")
```

#### Usage as CLI script:
```
ToDo
```

### 3. Convert nexus file into a hdf5 file.
#### Usage as module:
```
import alignment2hdf5
ToDo
```

#### Usage as CLI script:
```
ToDo
```