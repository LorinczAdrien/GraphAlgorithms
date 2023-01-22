# GraphAlgorithms
 A selection of graph algorithms on directed, undirected graphs.

# Important
- Documentation, functions names are all written in Hungarian.

# Algorithms
- General graph algorithms (Implemented in parent class 'Graf')
    - BFS                                       - szelessegi_bejaras
    - DFS                                       - melysegi_bejaras
    - Topological sort                          - topologiai_rendezes
    - Moore's shortest path algorithm           - Moore_SP
    - Dijkstra's shortest path algorithm        - Dijsktra_SP
    - Bellman-Ford shortest path algorithm      - Bellman_Ford_SP
    - Circle detection                          - van_kor
    - Negative circle detection                 - van_negativ_kor
    - Is graph regular                          - regularis

- Undirected graph algorithms
    - Prim's MST algorithm                      - Prim_MST
    - Kruskal's MST algorithm                   - Kruskal_MST
    - Boruvka's MST algorithm                   - Boruvka_MST
    - Reverse-delete MST algorithm              - Forditott_torles_MST
    - Biconnected components                    - Tarjan_BiConnect
    - TSP 2 approximation                       - TSP_2_approximate
    - BFS shortest path between two nodes       - BFS_SP

- Directed graph algorithms
    - Transitive closure                        - tranzitiv_lezaras or tranzitiv_lezaras_matrix
    - +- strongly connected components          - plusz_minusz_algoritmus
    - Kosaraju's strongly connected components  - Kosaraju_algoritmusa
    - Tarjan's StrongConnect algorithm          - Tarjan_StrongConnect
    - Johnson's SP between all nodes            - Johnson
    - Floyd-Warshall SP between all nodes       - Floyd-Warshall
    - Edmonds-Karp maximum flow                 - Edmonds_Karp
    - Pump maximum flow algorithm               - Pumpalo_algoritmus
    - Critical path algorithm                   - Kritikus_ut_masodik_modell

- 'GridGraph' algorithms
    - BFS shortest path between two nodes       - BFS_SP
    - A* shortest path between two nodes        - A_stat_SP

- Directed tree algorithms
    - Center node                               - center_node
    - Lowest Common Ancestor                    - LCA

- Undirected tree algorithms
    - Center node                               - center_node

- Binary tree algorithms
    - Preorder, Inorder, Postorder traversals   - preorder_bejaras, inorder_bejaras, postorder_bejaras
    - Prufer code algorithm                     - Prufer_kodolas

# Other
- The '<<' operator is overloaded for every class -> It displays the relevant information stored in class.
- The '+', '+=' operators are overloaded on: Graf, IranyitatlanGraf, IranyitottGraf -> Add edges, a list of edges or even a whole graph to the current graph.
- The '-', '-=' operators are overloaded on: Graf, IranyitatlanGraf, IranyitottGraf -> Remove edges, a list of edges or even a whole graph from the current graph.