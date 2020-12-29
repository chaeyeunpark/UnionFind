import matplotlib.pyplot as plt
import sys
import json

def plot_lattice(ax, L):
    for i in range(L):
        ax.hlines(i, -0.5, L-0.5, color='gray')
        ax.vlines(i, -0.5, L-0.5, color='gray')

    ax.set_xlim([-0.5, L-0.5])
    ax.set_ylim([-0.5, L-0.5])

def to_vertex_index(L, row, col):
    return (((row + L) % L)) * L + (col + L) % L;

def vertex_to_coord(L, vidx):
    return (vidx // L, vidx % L)

def to_edge(L, edge_index):
    row, col = divmod(edge_index, 2*L)

    if col < L: # horizontal edge
        u = to_vertex_index(L, row, col)
        v = to_vertex_index(L, row, col+1)

    else:
        u = to_vertex_index(L, row, col - L)
        v = to_vertex_index(L, row+1, col - L)

    return (min(u,v), max(u,v))

def lattice_to_plot(coord):
    return (coord[1], coord[0])

def edge_index_to_plot(L, edge_index):
    edge = to_edge(L, edge_index)
    if (edge[1] - edge[0]) == 1 or (edge[1] - edge[0]) == L:
        return [[lattice_to_plot(vertex_to_coord(L, edge[0])),
                lattice_to_plot(vertex_to_coord(L, edge[1]))]]
    elif edge[1] - edge[0] == L-1:
        x, y = lattice_to_plot(vertex_to_coord(L, edge[1]))
        return [[(x,y), (x+0.5,y)], [(-0.5,y), (0,y)]]
    else:
        x, y = lattice_to_plot(vertex_to_coord(L, edge[1]))
        return [[(x,y), (x,y+0.5)], [(x,-0.5), (x,0)]]

def plot_errors(ax, L, errors, error_type):
    from matplotlib import collections  as mc

    edges = []
    for n in range(2*L*L):
        if errors[n] % 2 != 0:
            edges += edge_index_to_plot(L, n)

    if error_type == 'Z':
        color = 'blue'
    else:
        color = 'red'

    lc = mc.LineCollection(edges, colors=color, linewidths=2)
    ax.add_collection(lc)

def plot_symdrome(ax, L, syndromes, error_type):
    if error_type == 'Z':
        xs = []
        ys = []
        for n in syndromes:
            coord = vertex_to_coord(L,n)
            x, y = lattice_to_plot(coord)
            xs.append(x)
            ys.append(y)
        ax.scatter(xs, ys, marker='x', s=200, linewidths=0.2)
    else:
        xs = []
        ys = []
        for n in syndromes:
            coord = vertex_to_coord(L,n)
            x, y = lattice_to_plot(coord)
            xs.append(x+0.5)
            ys.append(y+0.5)
        ax.scatter(xs, ys, marker='o', s=200, linewidths=0.2)


def symdrome_locations(L, syndromes):
    return [i for i in range(L*L) if syndromes[i] == 1]

if __name__ == '__main__':
    fig, ax = plt.subplots(figsize=(8,8))

    with open(sys.argv[1]) as f:
        error_json = json.load(f)

    error_type = sys.argv[2]
    assert(error_type == 'X' or error_type == 'Z')

    if error_type == 'Z':
        errors = error_json['z_errors']
        syndromes = error_json['syndromes_for_z']
    else:
        errors = error_json['x_errors']
        syndromes = error_json['syndromes_for_x']

    L = error_json['L']
    plot_lattice(ax, L)
    plot_errors(ax, L, errors, error_type)

    plot_symdrome(ax, L, symdrome_locations(L,syndromes), error_type)

    plt.show()
