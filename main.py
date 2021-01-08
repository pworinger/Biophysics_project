import numpy as np
import networkx as nx
import time
from model import Model
from network_visualisation import draw


m = Model(perturbed_gene=10)

# print(m.graph.nodes())
# print(m.graph.edges())
#
# print(nx.get_node_attributes(m.graph, "expression_level"))
# print(nx.get_edge_attributes(m.graph, "alpha"))
# print(nx.get_edge_attributes(m.graph, "beta"))
#
# print(m.graph.node["STE12"]["expression_level"])
# for n in nx.all_neighbors(m.graph, "STE12"):
#     print(n)
# print(m.graph.adj["STE12"]["SRY1"]["alpha"])

# print(type(m.alpha))
# print(m.alpha)

print(nx.number_of_nodes(m.graph))
print(nx.number_of_edges(m.graph))
# print(nx.get_edge_attributes(m.graph, "alpha"))
# print(nx.get_edge_attributes(m.graph, "beta"))
# print(nx.get_edge_attributes(m.graph, "weight"))
# print(nx.find_cycle(m.graph, source="SWI5", orientation="original"))
# print(nx.find_cycle(m.graph, source="SWI5", orientation="original"))
# print(nx.find_cycle(m.graph, source="SWI5", orientation="original"))
# print(nx.find_cycle(m.graph, source="SWI5", orientation="original"))
# print(nx.find_cycle(m.graph, source="SWI5", orientation="original"))

# for gene in ['CLN3', 'CLN1', 'CLN2', 'CDH1',  'CDC20', 'CLB5', 'CLB6', 'SIC1', 'CLB1', 'CLB2', 'MCM1', 'SWI5']:
#     print(gene)
#     for gene2 in nx.neighbors(m.graph, gene):
#         print("     ", gene2)
#         print(list(nx.all_simple_paths(m.graph, source=gene2, target=gene, cutoff=4)))

import pickle

cycles3 = []
for gene in nx.neighbors(m.graph, "SWI5"):
    cycles3.append(list(nx.all_simple_paths(m.graph, source=gene, target="SWI5", cutoff=3)))

with open('cycles3.pkl', 'wb') as f:
    pickle.dump(cycles3, f)

# del cycles3
#
# with open('cycles3.pkl', 'rb') as f:
#     myCycles3 = pickle.load(f)
#     print(len(myCycles3))
#     print([len(l) for l in myCycles3])
#     print(myCycles3)
# del myCycles3
#
# cycles4 = []
# for gene in (m.graph).neighbors("SWI5"):
print(len(list((m.graph).successors("SWI5"))))
print(len(list((m.graph).predecessors("SWI5"))))
print(nx.is_directed(m.graph))
#     cycles4.append(list(nx.all_simple_paths(m.graph, source=gene, target="SWI5", cutoff=4)))
#
# with open('cycles4.pkl', 'wb') as f:
#     pickle.dump(cycles4, f)
#
# del cycles4

# cyclesB = []
# cyclesE = []
# t = time.time()
# for gene in (m.graph).successors("SWI5"):
#     #print(gene)
#     cyclesB.append(list(nx.all_simple_paths(m.graph, source=gene, target="SWI5", cutoff=3)))
# print(time.time() - t)
#
# t = time.time()
# for gene in (m.graph).predecessors("SWI5"):
#     #print(gene)
#     cyclesE.append(list(nx.all_simple_paths(m.graph, source="SWI5", target=gene, cutoff=3)))
# print(time.time() - t)

print("55555555555555555555555555555555555555555555555555555555555555555555555")

cycles5 = []
for gene in (m.graph).predecessors("SWI5"):
    print(gene)
    cycles5.append(list(nx.all_simple_paths(m.graph, source="SWI5", target=gene, cutoff=5)))

with open('cycles5.pkl', 'wb') as f:
    pickle.dump(cycles5, f)

del cycles5

#print(list(nx.simple_cycles(m.graph)))

#draw(m.graph, "test")

#-----------------testing_evolve------------------------------------------

# temp = np.copy(m.y)


# m.evolve(10)

# print(np.concatenate((temp.reshape(len(temp),1), m.y.reshape(len(m.y),1)), axis=1)[0:100,:])

#----------------testing_cycle--------------------------------------------
# temp = np.copy(m.y)
#
# m.cycle(0.01, num_step = int(5e4), track=346)
#
# print(np.concatenate((temp.reshape(len(temp),1), m.y.reshape(len(m.y),1)), axis=1)[0:10,:])
#
# m.visualize(0.01)


def rank_cycles(loops):
    flat_list = [item for sublist in loops for item in sublist]
    scores = []
    for sublist in flat_list:
        # score measures the strength of the loop
        score = 1
        previous = 'SWI5'
        for node in sublist:
            score = score * m.graph.get_edge_data(previous, node, default = 0)['alpha']
            previous = node
        scores.append(score)
    scores.sort(reverse=True, key=abs)
    return scores