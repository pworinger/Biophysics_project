import pandas as pd
import networkx as nx
import numpy as np
import math

class Model():
    def __init__(self, path_edge_list="idea_ode_coefficients.tsv", path_node_id="idea_wide_format_data.txt"):
        full_graph = self.create_full_graph(path_edge_list, path_node_id)
        self.graph = self.create_subgraph(full_graph)
        self.alpha = nx.adjacency_matrix(self.graph, nodelist=None, weight="alpha")
        self.beta = nx.adjacency_matrix(self.graph, nodelist=None, weight="beta")
        self.y = self.initial_expression_level(self.graph.number_of_nodes())

    def create_full_graph(self, path_edge_list, path_node_id):
        df = pd.read_csv(path_edge_list, sep="\t")
        df.columns = ["cause", "effect", "alpha", "beta"]
        G = nx.from_pandas_edgelist(df, source="cause", target="effect", edge_attr=["alpha", "beta"], create_using=nx.DiGraph())
        self.index_to_name = pd.read_csv(path_node_id, sep="\t", usecols=["GENE"])
        # self.index_to_name = self.index_to_name[self.index_to_name["GENE"].isin(G.nodes())]
        self.index_to_name.to_dict()["GENE"]
        self.name_to_index = pd.read_csv(path_node_id, sep="\t", usecols=["GENE"]).reset_index().set_index("GENE").to_dict()["index"]
        G = nx.relabel_nodes(G, self.name_to_index, copy=True)
        nx.set_node_attributes(G, self.index_to_name, name="GENE")
        return G

    lst = ['CLN3', 'CLN1', 'CLN2', 'CDH1', 'SWI5', 'CDC20', 'CLB5', 'CLB6', 'SIC1', 'CLB1', 'CLB2', 'MCM1']

    def create_subgraph(self, full_graph, list_of_start_genes=lst, method="include_neighbours_at_max_distance", max_distance=1):
        index_lst = [self.name_to_index[gene] for gene in list_of_start_genes]
        if method == "include_neighbours_at_max_distance":
            acc = set()
            for gene in index_lst:
                if max_distance == 1:
                    acc = acc.union(set(nx.all_neighbors(full_graph, gene)))
                else:
                    acc = acc.union(set(nx.ego_graph(full_graph, gene, radius=max_distance, center=True, undirected=True, distance=None).nodes()))
            return full_graph.subgraph(list(acc))
        else:
            raise NotImplementedError

    def initial_expression_level(self, size):
        return np.random.rand(size)

    # def evolve(self, dt):
    #     temp = np.copy(self.y)
    #     for i in range(len(self.y)):
    #         sum_ = sum(self.alpha[:,i]*(temp-1) + self.beta[:,i]*(temp[i]*temp-1))
    #         self.y[i] = math.exp(dt*sum_/temp[i] + math.log(temp[i]))

    def evolve(self, dt):
        # want to multiply all edges going to zero with y
        # those edges correspond to alpha[:,0]
        # therefor we can transpose alpha and multiply it with y
        sum_ = self.alpha.T @ (self.y - 1) + self.y * (self.beta.T @ self.y) - self.beta.T @ np.ones(len(self.y))
        self.y = np.exp( dt*sum_/self.y) * self.y
