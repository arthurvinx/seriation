## About this software
This software solves the [Seriation](http://www.jstatsoft.org/v25/i03) problem finding a suitable linear order for a set of proteins. The result is a list of proteins ordered in one dimension such that functionally associated proteins are closer.

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
Type "./ordering1D" to show the options:

<pre>
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

To execute the software type:

<pre>
./ordering1D f=<absolute path>
</pre>

The input file is a textual file describing an undirected network of nodes. In our [example](data/362663.protein.links.900.v11.0.txt) the nodes are labeled by ENSEMBL Peptide IDs.

For a quicker test you can execute smaller dataset, like the Escherichia Coli.

> ./cfm-seriation f=Escherichia_coli.dat

The results are save in the input file directory.

## Inputs

The input file is a textual file describing an undirected network of nodes (in our examples the nodes are labeled by ENSEMBL Peptide IDs). Example:
<pre>
362663.ECP_0002	362663.ECP_4228
362663.ECP_0002	362663.ECP_0939
362663.ECP_4228 362663.ECP_0002
362663.ECP_0939 362663.ECP_0002
</pre>

This repository contains an example input file from Escherichia coli named [362663.protein.links.900.v11.0.txt](362663.protein.links.900.v11.0.txt).

## Outputs

The output is a text file with the order of the network nodes. Example:

<pre>
Protein	dim1
Z5822	0
Z5823	1
Z2911	2
Z2910	3
Z2909	4
Z4123	5
Z4124	6
Z3105	7
Z3106	8
...
</pre>

The following image represents the Homo sapiens network with a *random ordering*.

![initial](initial.png)

The next image represents the Homo sapiens network **'seriated'**.

![final](final.png)

## License
The source code is distributed under the terms of the GNU General Public License v3 [GPL](http://www.gnu.org/copyleft/gpl.html).

## How to cite this software
If you are using this package on your research please cite:

* [Kuentzer, F. A. et al. (2014). Optimization and analysis of seriation algorithm for ordering protein networks.
IEEE International Conference on Bioinformatics and Bioengineering, 231-237.](https://doi.org/10.1109/BIBE.2014.43)

## Similar softwares
* [Seriation R Package](http://www.jstatsoft.org/v25/i03), available at [CRAN](http://cran.r-project.org/web/packages/seriation/index.html).