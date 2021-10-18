## Multiple Sequence Alignment for DNA/RNA


[website](http://lab.malab.cn/~cjt/MSA/soft.html)

### Usage

```shell
usage: java -jar jmsa-**.jar  [-m] mode [-i] path [-o] path

  necessary arguments: 
    -i  Input file path (nucleotide sequences in fasta format)
    -m  three align option (1.Tree 2.Center 3.Cluster)
         1. Tree     accurate but slowest
         2. Center   rough but fastest
         3. Cluster  accurate and faster

  optional arguments: 
    -o  Output file path (default: infile.mode)
```

### Change Log
---

- 18-10-2021, version 0.3-alpha
  
  add cluster Tree mode, and fix some problem in FM index

- 3-9-2021, version 0.2-alpha
  
  fix clustering problems

- 12-5-2021, version 0.1-alpha
  
  inital version
  

### Development Environment

- JDK 1.8

- Intellij IDEA (Maven)
