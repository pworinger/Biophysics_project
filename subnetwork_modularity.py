def grow_subnetwork(network, seeds):
    G = network.to_undirected(reciprocal=False, as_view=False)
    subnetwork_nodes = set()
    for seed in seeds:
        subnetwork_nodes.add(seed)
    new_node = True
    first_step = True
    while new_node:
        if first_step:
            first_step = False
            Wint = 0
            Wext = 0
            for seed in seeds:
                for neighbour, dict in G[seed].items():
                    if neighbour not in subnetwork_nodes:
                        Wext += dict["weight"]
                    else:
                        Wint += dict["weight"]
            Wint = Wint / 2
            w_tot = Wint + Wext
            M = Wint / (1 + w_tot * w_tot)
            neighbours = set(G[seed].keys())
        else:
            Wint = stored_wint
            Wext = stored_wext
            M = stored_M

        best_node = 0
        for candidate in neighbours:
            # compute modularity
            add_wint = 0
            add_wext = 0
            for neigh in G[candidate]:
                if neigh in subnetwork_nodes:
                    add_wint += G[candidate][neigh]["weight"]
                else:
                    add_wext += G[candidate][neigh]["weight"]
            new_w_int = Wint + add_wint
            new_w_ext = Wext - add_wint + add_wext
            new_w_tot = new_w_int + new_w_ext
            new_M = new_w_int / (1 + new_w_tot * new_w_tot)
            if new_M >= M:
                best_node = candidate
                stored_wint = new_w_int
                stored_wext = new_w_ext
                stored_M = new_M

        if best_node != 0:
            neighbours.remove(best_node)
            neighbours = neighbours.union(set(G[best_node].keys()) - subnetwork_nodes)
            subnetwork_nodes.add(best_node)
        else:
            new_node = False
        #print(subnetwork_nodes)
        print("current modularity: ", stored_M)
    return subnetwork_nodes, M

# ################################################################################
#
# from model import Model
#
# m = Model()
#
# genes_12 = ['CLN3', 'SWI5', 'CLN1', 'CLN2', 'CDH1',  'CDC20', 'CLB5', 'CLB6', 'SIC1', 'CLB1', 'CLB2', 'MCM1']
#
# subN, M = grow_subnetwork(m.graph, genes_12)
#
# print(len(list(subN)))
#
# import pickle
#
# with open('subN_with12.pkl', 'wb') as f:
#     pickle.dump(subN, f)