## About this software
This software solves the [Seriation](http://www.jstatsoft.org/v25/i03) problem finding a suitable linear order for a set of proteins. The result is a list of proteins ordered in one dimension such that functionally associated proteins are closer.

![Figure 1](figure/F1.png)

## Authors
The software was developed by [Felipe Kuentzer](http://lattes.cnpq.br/1979213773480902), in collaboration with 
Douglas G. Ávila, Alexandre Pereira, Gabriel Perrone, Samoel da Silva, [Alexandre Amory](http://lattes.cnpq.br/2609000874577720), and [Rita de Almeida](http://lattes.cnpq.br/4672766298301524).

The version provided here was modified by [Clovis Ferreira dos Reis](http://lattes.cnpq.br/5487049518249525) to improve the textual feedback and to avoid bugs like:
* Duplication of identifiers on the ordering output.
* Segmentation fault while reading an input file containing many nodes.

## Download and compilation
Compilation requires GCC. To compile this software invoke the following commands on the shell:
<pre>
> wget https://github.com/arthurvinx/seriation/archive/master.zip
> unzip master.zip
> cd seriation-master/
> gcc ordering1D.c -o ordering1D -lm
</pre>

## How to use
To execute the software invoke this command on the shell:

<pre>
> ./ordering1D f=[absolute path to association file]
</pre>

Parameters list:
<pre>
> ./ordering1D

An association file name is necessary! No default!

Parameters list:
        f=Association file
        i=Number of isothermal steps
        m=Number of Monte Carlo steps
        c=Cooling factor
        a=Alpha value
        p=Percentual energy for initial temperature
        s=Random seed
</pre>

Parameters default values:
<pre>
i=100
m=2000
c=0.5
a=1.0
p=0.0001
</pre>

## Input
The input is a text file describing an undirected protein-protein interaction (PPI) network. This repository
contains an [example file](data/362663.protein.links.900.v11.0.txt) from *Escherichia coli*. In this example, the nodes are labeled by ENSEMBL Peptide IDs.

Protein-protein interaction network data can be downloaded from [STRING](https://string-db.org/). You may choose to download the information with the subscores per channel and tune your filters. The input must be a file containing two columns, no header, with rows composed by the IDs of two proteins that interact with each other.

## Outputs
Two text files will be saved in the association file directory, one containing the prefix "energy_" detailing the ordering process, and one containing the prefix "ordering_" (this will be your ordered list). The lower the final energy, the better the ordered list. To improve the outputs, I suggest to increase the number of Monte Carlo steps to 20000.

This repository contains an [example of the output produced by this software](output/ordering_362663.protein.links.900.v11.0.txt)
for the *Escherichia coli* PPI network.

## License
The source code is distributed under the terms of the GNU General Public License v3 [GPL](http://www.gnu.org/copyleft/gpl.html).

## How to cite this software
If you are using this package on your research please cite:

* [Kuentzer, F. A. et al. (2014). Optimization and analysis of seriation algorithm for ordering protein networks.
IEEE International Conference on Bioinformatics and Bioengineering, 231-237.](https://doi.org/10.1109/BIBE.2014.43)

## Similar softwares
* [Seriation R Package](http://www.jstatsoft.org/v25/i03), available at [CRAN](http://cran.r-project.org/web/packages/seriation/index.html).