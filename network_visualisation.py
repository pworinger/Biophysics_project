import networkx as nx
import pydot
import math

def draw(G, path_fig):

    width = {k: math.log(1E7 * v) for k, v in nx.get_edge_attributes(G, "weight").items()}
    nx.set_edge_attributes(G, width, name="penwidth")

    def get_color(b):
        if b:
            return "green"
        else:
            return "red"
    edge_color = {k: get_color(v) for k, v in nx.get_edge_attributes(G, "sign").items()}
    nx.set_edge_attributes(G, edge_color, name="color")

    nx.set_node_attributes(G, "filled", name="style")
    nx.set_node_attributes(G, "10", name="fontsize")
    #nx.set_node_attributes(G, True, name="fixedsize")
    nx.set_node_attributes(G, 0.7, name="width")
    nx.set_node_attributes(G, 0.7, name="height")

    nx.set_edge_attributes(G, 3, name="len")

    nx.drawing.nx_pydot.write_dot(G, path_fig)
    graph_a = pydot.graph_from_dot_file(path_fig)
    graph_a[0].write_svg("images/" + path_fig + ".svg", prog="neato")
    print("done")