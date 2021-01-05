import networkx as nx
from model import Model

m = Model()

print(m.data.nodes())
print(m.data.edges())

print(nx.get_node_attributes(m.data, "expression_level"))
print(nx.get_edge_attributes(m.data, "alpha"))
print(nx.get_edge_attributes(m.data, "beta"))

print(m.data.node["STE12"]["expression_level"])
for n in nx.all_neighbors(m.data, "STE12"):
    print(n)
print(m.data.adj["STE12"]["SRY1"]["alpha"])