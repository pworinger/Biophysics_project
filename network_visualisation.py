# import networkx as nx
# import pydot
# import math
#
# def draw(G, file_name):
#
#     width = {k: math.log(1E7 * v) for k, v in nx.get_edge_attributes(G, "weight").items()}
#     nx.set_edge_attributes(G, width, name="penwidth")
#
#     def get_color(b):
#         if b:
#             return "green"
#         else:
#             return "red"
#     edge_color = {k: get_color(v) for k, v in nx.get_edge_attributes(G, "sign").items()}
#     nx.set_edge_attributes(G, edge_color, name="color")
#
#     nx.set_node_attributes(G, "filled", name="style")
#     nx.set_node_attributes(G, "15", name="fontsize")
#     #nx.set_node_attributes(G, True, name="fixedsize")
#     nx.set_node_attributes(G, 1, name="width")
#     nx.set_node_attributes(G, 1, name="height")
#
#     nx.set_edge_attributes(G, 6, name="len")
#
#     nx.drawing.nx_pydot.write_dot(G, file_name)
#     graph_a = pydot.graph_from_dot_file(file_name)
#     graph_a[0].write_svg(+ file_name + ".svg", prog="neato")
#

# from model import Model
#
# nodes = ['SWI5', 'DAL5', 'ASH1', 'PCL5', 'REE1', 'PCL5', 'UGA1', 'HO']
# nodes2 = ['SWI5', 'PRP9', 'CTR1', 'REE1', 'FRE1', 'PCL5']
# m = Model(nodes_subgraph=nodes2)
#
# draw(m.graph, "network2")
