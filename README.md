# mat_comp
```
from ortools.linear_solver import pywraplp
import os
import networkx as nx 
import random
import matplotlib.pyplot as plt
from itertools import combinations

# Função para ler arquivo
def lerArquivo(arquivo): 
    with open(arquivo, 'r') as file:
        linhas = file.readlines()
    
    # 1ª linha: número de vértices, número de centrais, número de arestas
    nv, nc, na = map(int, linhas[0].split())
    
    centrais = []
    arestas = []
    
    # Linhas seguintes: centrais
    for i in range(1, nc + 1): #4, i = 0
        vertice, grau_minimo = map(int, linhas[i].split())
        centrais.append((vertice, grau_minimo))
    
    # Linhas seguintes: arestas
    for i in range(nc + 1, nc + na + 1):
        u, v, custo = map(int, linhas[i].split())
        arestas.append((u, v, custo))
        
    #criar grafo para esses parametros passados
    G = nx.Graph()  
    G.add_nodes_from(range(1, nv + 1))  # Add vértices
    for u, v, custo in arestas:
        G.add_edge(u, v, weight=custo)  # Add arestas e seus pesos

    #print(f"nv: {nv}")
    #print(f"nc: {nc}")
    #print(f"na: {na}")
    #print(f"centrais: {centrais}c")
    #print(f"arestas: {arestas}a")
    
    if not nx.is_connected(G):
        print("O grafo dado não é conexo.")
        
    return nv, nc, G, centrais, arestas

# Função para resolver o problema e retirar soluções com ciclos (1ª forma)
# Método Eliminação de ciclos
def eliminarCiclos(nv, nc, G, centrais, arestas):

#    print(f"N°vert: {nv}")
  #  print(f"N°centrais: {nc}")
  #  print(f"Centrais: {centrais}c")
   # print(f"Arestas: {arestas}a")
    
    # Inicializar o solver
    solver = pywraplp.Solver.CreateSolver('SCIP')
    solver.SetTimeLimit(1800 * 1000) #1800 seg

    T = nx.Graph()  # inicializar grafo vazio
    T.add_nodes_from(G.nodes())  # tem mesmas vértices
    if not nx.is_connected(G):
        print("graph given is not conexo")
        return None

    # criar um x binário para cada aresta (u, v)
    x = {}  # dicionario
    for (u, v, custo) in arestas:  # preencher com (u1, v1): x1
        x[(u, v)] = solver.IntVar(0, 1, f'x_{u}_{v}') #.format(u,v)?
        
    arestas_centrais = [(u, v, custo) for (u, v, custo) in arestas if u in [i for i, _ in centrais] and v in [i for i, _ in centrais]]
    
    # Função objetivo: minimizar o custo total das arestas selecionadas
    print("Definindo função obj...")
    solver.Minimize(solver.Sum(custo * x[(u, v)] for (u, v, custo) in arestas))

    # Restrição 1: Soma de todas as variáveis x(c,c) = c - 1 (ps: c = número de centrais)
    print("Definindo restrições para centrais...")
    solver.Add(solver.Sum(x[u,v] for (u, v, _) in arestas_centrais) == (nc - 1))
    
    # Restrição 2: evitar ciclo (subtours)
    #criar S, eliminar ciclos -> restrições limita n° de arestas desse S  (n° arestas de S <= |S| - 1)
    print("Adicionando restrições para evitar ciclos...")
    for size in range(2, nc+1): #size entre 2 e nc - 1 (todos subconj possíveis)
        print(f"Tamanho: {size}")
        for subconj in combinations(range(1,nc+1), size): #todas combinações com vértices de size 0 a nc
            #(0,1,2,3)
            # i = 0, j assume 1, 2, 3 -> (0,1) (0,2) (0,3)
            # i = 1, j assume 2, 3 -> (1,2) (1,3)
            # i = 2, j assume 3 -> (2,3)    
            #todas combinações possíveis: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
            #calcular só as arestas existentes no Grafo original
            #print(f"Subconjunto considerado: {subconj}")
            arestas_subconj = [(u, v) for (u, v, _) in arestas_centrais if u in subconj and v in subconj and u < v] 
            atual = solver.Sum(x[(u, v)] for (u, v) in arestas_subconj if (u, v) in x or (v, u) in x)
            solver.Add(atual <= size - 1)          
            #print(f"Arestas no subconjunto {subconj}: todos são...{arestas_subconj} Valor da soma das arestas no subconjunto {atual}")
            
    # Restrição 3: Cada vértice central deve obedecer seu grau mínimo
    print("Definindo restrições de grau mínimo...")
    for i, grau_min in centrais: 
        solver.Add(solver.Sum(x[(u, v)] for (u, v, _) in arestas if u == i or v == i) >= grau_min)
    
    # Restrição 4: todos verts t(ñ centrais) tem exatamente uma aresta entre central e ñ central
    print("Definindo restrições para vértices não centrais...")
    for t in range(1, nv+1):
        if t not in [i for i, _ in centrais]:
            solver.Add(solver.Sum(x[(u, v)] for (u, v, custo) in arestas if (u == t and v in [i for i, _ in centrais]) or (v == t and u in [i for i, _ in centrais])) == 1)
    
    # Resolver o modelo
    status = solver.Solve()

    # Preencher o caminho da solução com arestas selecionadas
    if status == pywraplp.Solver.OPTIMAL:
        caminho_solucao = []
        for (u, v, custo) in arestas:
            if x[(u, v)].solution_value() == 1:
                caminho_solucao.append((u, v, custo))
                T.add_edge(u, v, weight=custo)
    
        print("Caminho da solução (arestas da árvore geradora mínima):")
        for aresta in caminho_solucao:
            print(aresta)
        
        print('Solution:')
        print('Objective value =', solver.Objective().Value()) #custo total
        print('Problem solved in %f milliseconds' % solver.wall_time())
        print('Problem solved in %d iterations' % solver.iterations())
    else:
        print('The problem does not have an optimal solution.')

    return T

# Ler arquivos e resolver o problema
pasta = 'arqsFabio/Instancias'
arquivos = [f for f in os.listdir(pasta) if f.endswith('.txt')]

# Exibir opções de arquivos para o usuário
print("Opções:")
for i, x in enumerate(arquivos):
    print(f"{i + 1}: {x}")

# Escolher um arquivo específico
teste = int(input("Escolha o número do arquivo que deseja usar: ")) - 1
arquivo_selecionado = os.path.join(pasta, arquivos[teste])

# Nome do arquivo escolhido
print(f"\nArquivo escolhido: {arquivo_selecionado}\n")
nv, nc, G, centrais, arestas = lerArquivo(arquivo_selecionado)

# Resolver o problema usando a primeira forma
eliminarCiclos(nv, nc, G, centrais, arestas)

```
