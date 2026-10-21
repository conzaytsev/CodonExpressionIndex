from math import log, sqrt, isnan, exp
from scipy.stats import kendalltau, linregress
from importlib.resources import files
from . import datasets


def kendallr(query, target):
    result = kendalltau(query, target)[0]
    if isnan(result):
        return 0
    else:
        return result


class CodonAdaptationIndex:
    def __init__(self, path_to_dataset='ATCC_25922.csv', fraction = 0.125):
        with open(path_to_dataset, 'r') as file:
            data = file.readlines()

        genes = []
        quantities = []
        for i in range(len(data)):
            data_string = data[i].replace('\n', '').split(';')
            if len(data_string[2]) % 3 == 0:
                genes.append(data_string[2])
                quantities.append(log(float(data_string[1]) + 1))

        codon_dictionary = [['TTT', 'TTC'],
                            ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                            ['ATT', 'ATC', 'ATA'],
                            ['ATG'],
                            ['GTT', 'GTC', 'GTA', 'GTG'],
                            ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                            ['CCT', 'CCC', 'CCA', 'CCG'],
                            ['ACT', 'ACC', 'ACA', 'ACG'],
                            ['GCT', 'GCC', 'GCA', 'GCG'],
                            ['TAT', 'TAC'],
                            ['CAT', 'CAC'],
                            ['CAA', 'CAG'],
                            ['AAT', 'AAC'],
                            ['AAA', 'AAG'],
                            ['GAT', 'GAC'],
                            ['GAA', 'GAG'],
                            ['TGT', 'TGC'],
                            ['TGG'],
                            ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                            ['GGT', 'GGC', 'GGA', 'GGG']]

        codons = []
        for line in codon_dictionary:
            codons.extend(line)

        reverse_dictionary = {}
        for i in range(len(codon_dictionary)):
            for codon in codon_dictionary[i]:
                reverse_dictionary[codon] = i

        codon_count = {codon: 0 for codon in codons}
        for gene in genes[: min(int(len(genes) * fraction), 70)]:
            for j in range(0, len(gene) - 3, 3):
                if gene[j: j + 3] in codons:
                    codon_count[gene[j: j + 3]] += 1

        for codon in codons:
            if codon_count[codon] == 0:
                codon_count[codon] = 1

        rscu = {codon: codon_count[codon] * len(codon_dictionary[reverse_dictionary[codon]]) / sum(codon_count[synonymous_codon] for synonymous_codon in codon_dictionary[reverse_dictionary[codon]]) for codon in codons}
        scores = {codon: log(rscu[codon] / max(rscu[synonymous_codon] for synonymous_codon in codon_dictionary[reverse_dictionary[codon]])) for codon in codons}
        self.scores = scores

        train_scores = []
        for gene in genes:
            gene_score = 0
            codon_counter = 0
            for i in range(0, len(gene), 3):
                if gene[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                    gene_score += scores[gene[i: i + 3]]
                    codon_counter += 1
            train_scores.append(exp(gene_score / codon_counter))

        slope, intercept, r, p, std_err = linregress(train_scores, quantities)
        self.slope = slope
        self.intercept = intercept
        self.predictions = [self.predict(gene) for gene in genes]


    def predict(self, sequence):
        gene_score = 0
        codon_counter = 0
        for i in range(0, len(sequence), 3):
            if sequence[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                gene_score += self.scores[sequence[i: i + 3]]
                codon_counter += 1

        return exp(gene_score / codon_counter)


class CodonExpressionIndex:
    def __init__(self, path_to_dataset='ATCC_25922.csv'):
        with open(path_to_dataset, 'r') as file:
            data = file.readlines()

        genes = []
        quantities = []
        for i in range(len(data)):
            if i > 0:
                data_string = data[i].replace('\n', '').split(';')
                if len(data_string[2]) % 3 == 0:
                    genes.append(data_string[2])
                    quantities.append(log(float(data_string[1]) + 1))

        codon_dictionary = [['TTT', 'TTC'],
                            ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                            ['ATT', 'ATC', 'ATA'],
                            ['ATG'],
                            ['GTT', 'GTC', 'GTA', 'GTG'],
                            ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                            ['CCT', 'CCC', 'CCA', 'CCG'],
                            ['ACT', 'ACC', 'ACA', 'ACG'],
                            ['GCT', 'GCC', 'GCA', 'GCG'],
                            ['TAT', 'TAC'],
                            ['CAT', 'CAC'],
                            ['CAA', 'CAG'],
                            ['AAT', 'AAC'],
                            ['AAA', 'AAG'],
                            ['GAT', 'GAC'],
                            ['GAA', 'GAG'],
                            ['TGT', 'TGC'],
                            ['TGG'],
                            ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                            ['GGT', 'GGC', 'GGA', 'GGG']]

        triplets = []
        for item in codon_dictionary:
            triplets.extend(item)

        reverse_dictionary = {}
        for i in range(len(codon_dictionary)):
            for triplet in codon_dictionary[i]:
                reverse_dictionary[triplet] = i

        distribution = {triplet: [0 for i in range(len(genes))] for triplet in triplets}
        for i in range(len(genes)):
            sequence = genes[i]
            for j in range(0, len(sequence) - 3, 3):
                if sequence[j: j + 3] in triplets:
                    distribution[sequence[j: j + 3]][i] += 1

        for i in range(len(genes)):
            gene_sum = 0
            for triplet in triplets:
                gene_sum += distribution[triplet][i]

            for triplet in triplets:
                distribution[triplet][i] /= gene_sum

        scores = {triplet: float(3 * kendallr(distribution[triplet], quantities) * sqrt(len(genes) * (len(genes) - 1)) / sqrt(2 * (2 * len(genes) + 5))) for triplet in triplets}
        self.scores = scores

        train_scores = []
        for gene in genes:
            gene_score = 0
            codon_counter = 0
            for i in range(0, len(gene), 3):
                if gene[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                    gene_score += scores[gene[i: i + 3]]
                    codon_counter += 1
            train_scores.append(gene_score / codon_counter)

        slope, intercept, r, p, std_err = linregress(train_scores, quantities)
        self.slope = slope
        self.intercept = intercept
        self.expression = quantities


    def predict(self, test_sequence):
        gene_score = 0
        codon_counter = 0
        for i in range(0, len(test_sequence), 3):
            if test_sequence[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                gene_score += self.scores[test_sequence[i: i + 3]]
                codon_counter += 1

        return self.slope * gene_score / codon_counter + self.intercept


class CodonProductivity:
    def __init__(self, path_to_dataset='ATCC_25922.csv'):
        with open(path_to_dataset, 'r') as file:
            data = file.readlines()

        genes = []
        quantities = []
        for i in range(len(data)):
            if i > 0:
                data_string = data[i].replace('\n', '').split(';')
                genes.append(data_string[2])
                quantities.append(float(data_string[1]))

        codon_dictionary = [['TTT', 'TTC'],
                            ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
                            ['ATT', 'ATC', 'ATA'],
                            ['ATG'],
                            ['GTT', 'GTC', 'GTA', 'GTG'],
                            ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
                            ['CCT', 'CCC', 'CCA', 'CCG'],
                            ['ACT', 'ACC', 'ACA', 'ACG'],
                            ['GCT', 'GCC', 'GCA', 'GCG'],
                            ['TAT', 'TAC'],
                            ['CAT', 'CAC'],
                            ['CAA', 'CAG'],
                            ['AAT', 'AAC'],
                            ['AAA', 'AAG'],
                            ['GAT', 'GAC'],
                            ['GAA', 'GAG'],
                            ['TGT', 'TGC'],
                            ['TGG'],
                            ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                            ['GGT', 'GGC', 'GGA', 'GGG']]

        triplets = []
        for item in codon_dictionary:
            triplets.extend(item)

        reverse_dictionary = {}
        for i in range(len(codon_dictionary)):
            for triplet in codon_dictionary[i]:
                reverse_dictionary[triplet] = i

        aa_copies = {triplet: 0 for triplet in triplets}
        triplet_copies = {triplet: 0 for triplet in triplets}
        for i in range(len(genes)):
            sequence = genes[i]
            expression = quantities[i]
            for j in range(0, len(sequence), 3):
                if sequence[j: j + 3] in triplets:
                    aa_copies[sequence[j: j + 3]] += expression
                    triplet_copies[sequence[j: j + 3]] += 1

        scores = {}
        for triplet in triplets:
            if triplet_copies[triplet] == 0:
                scores[triplet] = 0
            else:
                scores[triplet] = aa_copies[triplet] / triplet_copies[triplet]
        self.scores = scores

        train_scores = []
        for gene in genes:
            gene_score = 0
            codon_counter = 0
            for i in range(0, len(gene), 3):
                if gene[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                    gene_score += scores[gene[i: i + 3]]
                    codon_counter += 1
            train_scores.append(gene_score / codon_counter)

        for i in range(len(quantities)):
            quantities[i] = log(float(quantities[i]) + 1)

        slope, intercept, r, p, std_err = linregress(train_scores, quantities)
        self.slope = slope
        self.intercept = intercept
        self.expression = quantities


    def predict(self, test_sequence):
        gene_score = 0
        codon_counter = 0
        for i in range(0, len(test_sequence), 3):
            if test_sequence[i: i + 3] not in ['TAA', 'TAG', 'TGA']:
                gene_score += self.scores[test_sequence[i: i + 3]]
                codon_counter += 1

        return self.slope * gene_score / codon_counter + self.intercept
