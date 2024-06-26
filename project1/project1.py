import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from graphviz import Digraph
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

def getSequences(filename):

    sequences = []
    with open(filename, "r") as file:
        for i, line in enumerate(file):
            if i % 2 == 1: # every other line
                sequences.append(line.strip())
    return sequences

# zip pairs the elements at each location together for the different sequences
# there is a value of 1 at the location if they are not equal
# the return value is the amount of mismatches between the two sequences
def getHammingDistances(sequence1, sequence2):
    # I'm not sure how we should handle indels, so they're just treated as
    # another character. an indel in both sequences is a match not a mismatch
    return sum(x != y for x, y in zip(sequence1, sequence2))

# uses getHammingDistances to compare the alignments to each other and return
# a matrix of hamming distances from each sequence to each sequence
def getDistanceMatrix(sequences):

    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))

    for i in range(num_sequences):
        for j in range(i, num_sequences):
            distance = getHammingDistances(sequences[i], sequences[j])
            distance_matrix[j,i] = distance
            distance_matrix[i,j] = distance

    return distance_matrix.astype(np.int64)

# recursive function to graph the tree to print
def add_node(node):
    # if the node doesn't have any children
    if node.is_terminal():
        tree.node(str(node.name))
    else:
        tree.node(str(node.name))

        # clades[0] and [1] are the left and right children respectively
        # add nodes for each of the children of the current node
        add_node(node.clades[0])
        add_node(node.clades[1])

        # add an edge from primary node to left child
        tree.edge(node.name, node.clades[0].name, label=str(round(node.clades[0].branch_length, 3)))
        # add an edge from primary node to right child
        tree.edge(node.name, node.clades[1].name, label=str(round(node.clades[1].branch_length, 3)))


calculator = DistanceCalculator('identity')
sequences = getSequences("msa.pir")
seqs = [SeqRecord(Seq(sequences[0]), id="seq1"), SeqRecord(Seq(sequences[1]), id="seq2"),
        SeqRecord(Seq(sequences[2]), id="seq3"),SeqRecord(Seq(sequences[3]), id="seq4"),
        SeqRecord(Seq(sequences[4]), id="seq5"),SeqRecord(Seq(sequences[5]), id="seq6")]
# function to input alignments in so they're of the right form to be read in by the bioPhylo functions
alignment = MultipleSeqAlignment(seqs)

# normalizes distances to all be between 0 and 1,
# sequences are more similar for low values, less for higher values
matrix = calculator.get_distance(alignment)
print(matrix)

constructor = DistanceTreeConstructor()
upgma_tree = constructor.upgma(matrix)

tree = Digraph()

# start at the root of the tree and add each node to the graph to print
add_node(upgma_tree.root)
tree.render('tree', format='png')

# write the tree to a gv file
with open('tree1.gv', 'w') as file:
    file.write(tree.source)

# biophylo function to get the nodes from the tree
nodes = upgma_tree.get_nonterminals() + upgma_tree.get_terminals()
internalMatrix = [[0.0 for x in nodes] for y in nodes]

for i, node1 in enumerate(nodes):
    for j, node2 in enumerate(nodes):
        # calcualte the distance from each node to each other node in the tree
        internalMatrix[i][j] = upgma_tree.distance(node1, node2)

labels = [node.name for node in nodes]

# format and print out the table of the distance matrix with all the inner nodes included
print('\t' + '\t'.join(labels))
for i, data in enumerate(internalMatrix):
    print(labels[i] + '\t' + '\t'.join('{:.3f}'.format(x) for x in data))

# distance matrix from hamming distances, not used in tree building
#print(getDistanceMatrix(getSequences("msa.pir")))