import numpy as np
import networkx as nx
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

draw(m.graph, "test")

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
