import pandas as pd
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import math
import pydot

class Model():
    def __init__(self, path_edge_list="idea_ode_coefficients.tsv", path_node_id="idea_wide_format_data.txt", nodes_subgraph=[]):
        full_graph = self.create_full_graph(path_edge_list, path_node_id)
        if nodes_subgraph !=[]:
            self.graph = full_graph.subgraph(nodes_subgraph)
        else:
            self.graph = full_graph
        self.alpha = nx.adjacency_matrix(self.graph, nodelist=None, weight="alpha")
        self.beta = nx.adjacency_matrix(self.graph, nodelist=None, weight="beta")
        self.y = np.ones(self.graph.number_of_nodes())

    # this function load the data from the tsv files whose paths are provided as arguments
    # and convert the dataframe containing edge list to a networkx object,
    # putting alpha and beta values as attributes of the edges
    def create_full_graph(self, path_edge_list, path_node_id):
        df = pd.read_csv(path_edge_list, sep="\t")
        df.columns = ["cause", "effect", "alpha", "beta"]
        df["weight"] = abs(df["alpha"]) + df["beta"]*df["beta"]
        df["sign"] = df["alpha"] + df["beta"] > 0
        G = nx.from_pandas_edgelist(df, source="cause", target="effect", edge_attr=["alpha", "beta", "weight", "sign"], create_using=nx.DiGraph())
        self.index_to_name = pd.read_csv(path_node_id, sep="\t", usecols=["GENE"])
        self.index_to_name.to_dict()["GENE"]
        self.name_to_index = pd.read_csv(path_node_id, sep="\t", usecols=["GENE"]).reset_index().set_index("GENE").to_dict()["index"]
        #G = nx.relabel_nodes(G, self.name_to_index, copy=True)
        nx.set_node_attributes(G, self.index_to_name, name="GENE")
        nx.set_node_attributes(G, self.name_to_index, name="index")
        return G

    def evolve(self, dt):
        # want to multiply all edges going to zero with y
        # those edges correspond to alpha[:,0]
        # therefor we can transpose alpha and multiply it with y
        sum_ = self.alpha @ (self.y - 1) + self.y * (self.beta @ self.y) - self.beta @ np.ones(len(self.y))
        # if we remove the log we can avoid dividing by y
        self.y = sum_*dt + self.y

    # this function apply an impulsive perturbation to the perturbed_gene
    # and let the system evolve, while monitoring the expression levels for
    # the genes listed in track argument
    def cycle(self, dt, perturbed_gene, track = [0], num_step = None):
        perturbed_gene_index = self.graph.node[perturbed_gene]['index']
        if num_step:
            self.y[perturbed_gene_index] = 1.01
            self.dt = dt
            self.time_series = np.empty((num_step, len(track)))

            self.evolve(dt)
            self.time_series[0, :] = self.y[track]

            self.y[perturbed_gene_index] = 1
            for i in range(1,num_step):
                if (i+1) % int(num_step/10) ==0:
                    print(str((i+1)/num_step*100) + '% are done')
                self.evolve(dt)
                self.time_series[i, :] = self.y[track]
        else:
            print("")
            raise NotImplementedError

    # this function is used to create a plot of the data collected when running cycle function
    def visualize(self, legends):
        plt.figure(figsize=(10, 7))
        t = np.linspace(0, len(self.time_series) * self.dt - self.dt, len(self.time_series))
        plt.xlabel('time [minutes]')
        plt.ylabel('relative expresssion []')
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', "gold", "black"]
        for i in range(self.time_series.shape[1]):
            plt.plot(t, (self.time_series[:, i]), color=colors[i], )
        plt.ylim((1 - 0.001, 1 + 0.001))
        plt.legend(legends)
        plt.xlim((0, 1850))
        plt.xticks(np.arange(0, 1900, 100))
        plt.savefig('expression_levels.png', dpi=300)
        plt.xlim((0, 1850))
        plt.ylim((0.999, 1.001))
        plt.show()

    # this function is used to create a plot of the data collected when running cycle function,
    # after applying log to the expression levels and using log scale on y-axis
    # to enable to see oscillations despite fast divergence
    def visualize_double_log(self, legends):
        t = np.linspace(0, len(self.time_series)*self.dt-self.dt, len(self.time_series))
        plt.figure(figsize=(10, 7))
        plt.xlabel('time [minutes]')
        plt.ylabel('log relative expresssion []')
        colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', "gold", "black"]
        for i in range(self.time_series.shape[1]):
            plt.plot(t, np.log(self.time_series[:,i]), color=colors[i])
        plt.legend(legends)
        plt.yscale('log')
        plt.xlim((0, 1850))
        plt.ylim((1e-20, 10 ))
        plt.xticks(np.arange(0, 1900, 100))
        plt.savefig('expression_levels_double_log.png', dpi=300)
        plt.show()

    # this function is used to create a .svg image representing the network
    # a file in dot language is first created, and then used to generate the image
    def draw(self, file_name):
        G = self.graph

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
        nx.set_node_attributes(G, "15", name="fontsize")
        # nx.set_node_attributes(G, True, name="fixedsize")
        nx.set_node_attributes(G, 1, name="width")
        nx.set_node_attributes(G, 1, name="height")

        nx.set_edge_attributes(G, 6, name="len")

        nx.drawing.nx_pydot.write_dot(G, file_name)
        graph_a = pydot.graph_from_dot_file(file_name)
        graph_a[0].write_svg(file_name + ".svg", prog="neato")
        print("file '" + file_name+ ".svg' was successfully created")