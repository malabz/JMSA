[website](http://lab.malab.cn/~cjt/MSA/soft.html)

### Usage

```shell
usage: java -jar  [-m] mode [-i] path [-o] path

  necessary arguments: 
    -i  Input file path (nucleotide sequences in fasta format)
    -m  three align option (1.Tree 2.Center 3.Cluster)
         1. Tree     accurate but slowest
         2. Center   rough but fastest
         3. Cluster  accurate and faster

  optional arguments: 
    -o  Output file path (default: infile.mode)
```

