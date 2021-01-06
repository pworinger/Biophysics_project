import networkx as nx
from model import Model

m = Model()

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

print(type(m.alpha))
print(m.alpha)