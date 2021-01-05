import pandas as pd
import networkx as nx

class Model():
    def __init__(self, path_edge_list="idea_ode_coefficients.tsv"):
        self.full_graph = self.create_full_graph(path_edge_list)
        self.data = self.create_subgraph()
        self.initialize_expression_level()

    def create_full_graph(self, path_edge_list):
        df = pd.read_csv(path_edge_list, sep="\t")
        return nx.from_pandas_edgelist(df, source="cause", target="effect", edge_attr=["linear", "quad"])

    lst = ['CLN3', 'CLN1', 'CLN2', 'CDH1', 'SWI5', 'CDC20', 'CLB5', 'CLB6', 'SIC1', 'CLB1', 'CLB2', 'MCM1']

    def create_subgraph(self, list_of_start_genes=lst, method="include_neighbours_at_max_distance", max_distance=1):
        if method == "include_neighbours_at_max_distance":
            acc = set()
            for gene in list_of_start_genes:
                if max_distance == 1:
                    acc = acc.union(set(nx.all_neighbors(self.full_graph, gene)))
                else:
                    acc = acc.union(set(nx.ego_graph(self.full_graph, gene, radius=max_distance, center=True, undirected=True, distance=None).nodes()))
            return self.full_graph.subgraph(list(acc))
        else:
            raise NotImplementedError

    def initialize_expression_level(self):
        nx.set_node_attributes(self.data, 0, name="expression_level") # TODO use a dict-like argument to specify different value for each node

    def evolve(self):
        raise NotImplementedError