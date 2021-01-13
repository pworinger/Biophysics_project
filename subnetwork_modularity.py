def grow_subnetwork(network, seeds):
    # change network to undirected graph
    G = network.to_undirected(reciprocal=False, as_view=False)

    # initialize the subnetwork with the set of nodes provided as argument
    subnetwork_nodes = set()
    for seed in seeds:
        subnetwork_nodes.add(seed)
    new_node = True
    first_step = True

    # algorithm terminates when no new node can be added whitout reducing modularity
    while new_node:
        # at first step, Wint and Wext have to be computed from scratch
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
            Wint = Wint / 2 # each edge inside the network was counted twice
            Wtot = Wint + Wext
            M = Wint / (1 + Wtot * Wtot)
            neighbours = set(G[seed].keys())
        else:
            Wint = stored_Wint
            Wext = stored_Wext
            M = stored_M

        best_node = 0

        # for each neighbour of current subnetwork
        # compute modularity by updating Wint and Wext
        # and select the node that increase modularity the most
        for candidate in neighbours:
            add_wint = 0
            add_wext = 0
            for neigh in G[candidate]:
                if neigh in subnetwork_nodes:
                    add_wint += G[candidate][neigh]["weight"]
                else:
                    add_wext += G[candidate][neigh]["weight"]
            new_Wint = Wint + add_wint
            new_Wext = Wext - add_wint + add_wext
            new_Wtot = new_Wint + new_Wext
            new_M = new_Wint / (1 + new_Wtot * new_Wtot)
            if new_M >= M:
                stored_Wint = new_Wint
                stored_Wext = new_Wext
                stored_M = new_M
                best_node = candidate

        # if a valid node is found, it is added to the network
        if best_node != 0:
            neighbours.remove(best_node)
            neighbours = neighbours.union(set(G[best_node].keys()) - subnetwork_nodes)
            subnetwork_nodes.add(best_node)
        else:
            new_node = False
        #print(subnetwork_nodes)
        print("current modularity: ", stored_M)
    return subnetwork_nodes, M