import pandas as pd
import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt

class Model():
    def __init__(self, path_edge_list="idea_ode_coefficients.tsv", path_node_id="idea_wide_format_data.txt", perturbed_gene = 0):
        full_graph = self.create_full_graph(path_edge_list, path_node_id)
        self.graph = self.create_subgraph(full_graph)
        self.alpha = nx.adjacency_matrix(self.graph, nodelist=None, weight="alpha")
        self.beta = nx.adjacency_matrix(self.graph, nodelist=None, weight="beta")
        self.y = np.ones(self.graph.number_of_nodes())
        self.y[perturbed_gene] = np.random.uniform(0.8, 1.2)
        # self.y = self.initial_expression_level(self.graph.number_of_nodes())

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

    # def initial_expression_level(self, size, seed = 0):
    #     np.random.seed(seed=seed)
    #     return np.random.uniform(1e-3, 1 - 1e-3, size)

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
        # self.y = np.exp(dt*sum_/self.y) * self.y
        # if we remove the log we can avoid dividing by y
        self.y = sum_*dt + self.y
        # in the following I looked at what happened it I don't transpose alpha and beta
        # sum_ = self.alpha @ (self.y - 1) + self.y * (self.beta @ self.y) - self.beta @ np.ones(len(self.y))
        # self.y = np.exp(dt*sum_/self.y) * self.y
        
        
    def cycle(self, dt, track = 0, num_step = None):
        if num_step:
            self.time_series = np.empty(num_step)
            for i in range(num_step):
                if (i+1) % int(num_step/10) ==0:
                    print(str((i+1)/num_step*100) + '% are done')
                self.evolve(dt)
                self.time_series[i] = self.y[track]
        else:
            raise NotImplementedError

    def visualize(self, dt):
        # if self.time_series in dir():
        t = np.linspace(0, len(self.time_series)*dt-dt, len(self.time_series))
        plt.plot(t, self.time_series)
        plt.show()
