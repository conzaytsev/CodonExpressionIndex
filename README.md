# Codon Expression Index
Python module for analysis of codon influence on protein expression.

## Requirements
`Python >3.9`, `scipy`

## Citation
Zaytsev K, Bogatyreva N, Fedorov A. Link Between Individual Codon Frequencies and Protein Expression: Going Beyond Codon Adaptation Index. Int J Mol Sci. 2024 Oct 29;25(21):11622. https://doi.org/10.3390/ijms252111622

## Dataset
This model requires protein expression data for a specific organism.
Dataset should be presented as a .csv file containing three columns: protein ID, number of protein copies per cell, gene sequence.
First line of the table is reserved for the header.


Example:

    Protein ID;Protein copies;Gene sequence
    P69776;38022.530250746298;ATGAAAGCTACTAAACTGGTACTGGGCGCGGTAATCCTGGGTTCTACTCTGCTGGCAGGTTGCTCCAGCAACGCTAAAATCGATCAGCTGTCTTCTGACGTTCAGACTCTGAACGCTAAAGTTGACCAGCTGAGCAACGACGTGAACGCAATGCGTTCCGACGTTCAGGCTGCTAAAGATGACGCAGCTCGTGCTAACCAGCGTCTGGACAACATGGCTACTAAATACCGCAAGTAA
    P0AG51;4376.7908013480401;ATGGCAAAGACTATTAAAATTACTCAAACCCGCAGTGCAATCGGTCGTCTGCCGAAACACAAGGCAACGCTGCTTGGCCTGGGTCTGCGTCGTATTGGTCACACCGTAGAGCGCGAGGATACTCCTGCTATTCGCGGTATGATCAACGCGGTTTCCTTCATGGTTAAAGTTGAGGAGTAA
    P0A7Q6;3802.6818971891998;ATGAAAGTTCGTGCTTCCGTCAAGAAATTATGCCGTAACTGCAAAATCGTTAAGCGTGATGGTGTCATCCGTGTGATTTGCAGTGCCGAGCCGAAGCATAAACAGCGCCAAGGCTGA

The provided dataset for E. Coli ATCC 25922 is based on the experimental data by Jacek R. Wiśniewski and Dariusz Rakus (Quantitative analysis of the Escherichia coli proteome, https://doi.org/10.1016/j.dib.2014.08.004).

## Installation
From PyPI:

    pip install cei
    
For the latest development version:
    
    pip install git+https://github.com/conzaytsev/CodonExpressionIndex.git

## Quickstart
Importing the module:

    import cei

To use the Codon Expression Index model:

    model = cei.CodonExpressionIndex(path_to_dataset)
        
To use the Codon Productivity model:

    model = cei.CodonProductivity(path_to_dataset)
        
In order to use the default dataset for E. Coli ATCC 25922 leave the brackets empty.
        
    model = cei.CodonExpressionIndex()
    model = cei.CodonProductivity()

Codon scores for each of the models as a dict:

    model.scores
        
Predict number of protein copies per cell based on the nucleotide sequence:

    model.predict('ATG...')
