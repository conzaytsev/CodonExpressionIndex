# CodonExpressionIndex
Python module for analysis of codon influence on protein expression.

## Requirements
`Python 3`, `scipy`

## Dataset
This model requires protein expression data for specific organism.
Dataset should be presented as a .csv file containing three columns: protein ID, number of protein copies per cell, gene sequence.
First line of the table is reserved for the header.

Example:

    Protein ID;Protein copies;Gene sequence
    P69776;38022.530250746298;ATGAAAGCTACTAAACTGGTACTGGGCGCGGTAATCCTGGGTTCTACTCTGCTGGCAGGTTGCTCCAGCAACGCTAAAATCGATCAGCTGTCTTCTGACGTTCAGACTCTGAACGCTAAAGTTGACCAGCTGAGCAACGACGTGAACGCAATGCGTTCCGACGTTCAGGCTGCTAAAGATGACGCAGCTCGTGCTAACCAGCGTCTGGACAACATGGCTACTAAATACCGCAAGTAA
