import matplotlib
from Bio.SubsMat import MatrixInfo as scoringMatrices
from Bio.Alphabet import IUPAC
from Bio import pairwise2
import itertools
import functools
import collections
import random
import numpy as np
import networkx as nx
from sklearn.neighbors import KDTree
import matplotlib.pyplot as mpl
matplotlib.use('Agg')


# TODO Decide if you want it 2d or scalable ...
# TODO state nomenclature for functions, camelcase or underscores

class SymetricDict(collections.OrderedDict):
    def __init__(self, keys, *args, **kwargs):
        """
        A symetric dictionary, view usage examples for a a quick demonstration ...

        Parameters
        ----------
        keys : keys that one wants for the dimensions
        args : placeholder for extra arguments to be passed to other functions
        kwargs : placeholder for keyword atguments apssed to further functions

        Examples
        --------
        >>> foo = SymetricDict(['A', 'B'])
        >>> foo
        SymetricDict([('A', OrderedDict([('A', None), ('B', None)])), ('B', OrderedDict([('A', None), ('B', None)]))])
        >>> foo['A', 'B'] = 2 # this replaces the value in both ['A']['B'] and ['B']['A']
        >>> print(foo)
        SymetricDict([('A', OrderedDict([('A', None), ('B', 2)])), ('B', OrderedDict([('A', 2), ('B', None)]))])
        >>> foo['A', 'C'] = 3 # adds the keys when absent
        >>> print(foo)
        SymetricDict([('A', OrderedDict([('A', None), ('B', 2), ('C', 3)])), \
('B', OrderedDict([('A', 2), ('B', None)])), ('C', OrderedDict([('A', 3)]))])
        >>> print(foo['A', 'C']) # and well ...
        3
        >>> print(foo['C', 'A'])
        3
        >>> print(foo['C']['A']) # If you want to you can use the standard double square parens notation ...
        3
        """

        super(SymetricDict, self).__init__()
        elems = list(set(keys))
        elems.sort()
        for x in elems:
            for y in elems:
                self[x, y] = None

    def __getitem__(self, key):
        if (len(key) is 1) or (type(key) is str):
            val = collections.OrderedDict.__getitem__(self, key)
        else:
            val = functools.reduce(collections.OrderedDict.__getitem__, key, self)
        return val

    def __setitem__(self, key, val):
        if len(key) is 2:
            for keys in itertools.permutations(key, 2):
                for x in keys:
                    if x not in self:
                        collections.OrderedDict.__setitem__(self, x, collections.OrderedDict())
                collections.OrderedDict.__setitem__(self[keys[0]], keys[1], val)
        else:
            print("Dont knwo how to deal with that ...")
            print(key)

    def __delitem__(self, key):
        collections.OrderedDict.__delitem__(self, key)

    def __iter__(self):
        return collections.OrderedDict.__iter__(self)

    def __len__(self):
        return collections.OrderedDict.__len__(self)

    def check_lengths(self):
        second_dim_lens = list([set([len(x) for x in self.values()])])
        assert len(second_dim_lens) == 1, 'More than one length for the second dimension of the object'
        assert len(self) == second_dim_lens[0], 'Dim 1 and 2 have differing lengths'
        return True

    # TODO implement __instance_check__

    def closest_neighbors(self):
        all_top_values = [
            sorted(range(len(i)), key=lambda a: i[a])[-2:]  # TODO, get a way for them not to return themselves ...
            for i in [list(x.values()) for x in self.values()]
        ]
        return all_top_values


"""
TODO: make something to check that
- length of each dict is the same
- all dimensions have the same names and are same length
- fix bug when symdict is modified by doing
- Make this an extension of an ordered dict ...
>>> dict[x][y] = foo # does not modify dict[y][x]
"""


def build_2d_dict(combinations):
    """
    Gets a nested iterable object and returns a nested dictionary with the keys being the elements in the original
    nested sequence.

    NOTE: Has been deprecated, use a SymetricDict for this task :)

    Parameters
    ----------
    combinations : a nested iterable object

    Returns
    -------
    out : dict
        a nested dictionary, modifies the order of the first index to be alphabetical

    Examples
    --------
    >>> my_list = list(itertools.combinations(['a', 'b','c'], 2))
    >>> print(my_list)
    [('a', 'b'), ('a', 'c'), ('b', 'c')]
    >>> build_2d_dict(itertools.combinations(['a', 'b','c'], 2))
    {'a': {'b': None, 'c': None}, 'b': {'a': None, 'c': None}, 'c': {'a': None, 'b': None}}

    """
    combinations = list(combinations)
    unique_first = list(set([x[0] for x in combinations]))
    unique_first.sort()
    my_nested = {}
    for x in unique_first:
        if x not in my_nested:
            my_nested[x] = {}
        for _, y in [i for i in combinations if i[0] == x]:
            my_nested[x][y] = None
            if y not in my_nested:
                my_nested[y] = {}
            my_nested[y][x] = my_nested[x][y]
    return my_nested


def getPermutations(letters, n, sort=True, drop_repeats=True):
    """

    Parameters
    ----------
    letters : iterable object of whose we are going to get the permutations
    n : number of elements to permute
    sort : should we sort the input letters (and therefore the output ...)
    drop_repeats : should we drop the repeated elements

    Returns
    -------
    A list with the permutations posible

    Examples
    --------
    >>> getPermutations("ABA", 2, sort = True, drop_repeats = False) # Note that it outputs repeated elements
    ['AA', 'AA', 'AB', 'AA', 'AA', 'AB', 'BA', 'BA', 'BB']
    >>> getPermutations("ABA", 2, sort = True, drop_repeats = True)
    ['AA', 'AB', 'BA', 'BB']
    """
    if drop_repeats:
        letters = set(letters)
    if sort:
        letters = list(letters)
        letters.sort()
    my_permutations = ["".join(x) for x in itertools.product(letters, repeat=n)]
    return my_permutations


class BlosumDictGenerator(SymetricDict):
    # TODO change name ... this is not really a generator
    def __init__(self, sequences, scoringmatrx=scoringMatrices.blosum62):
        """

        Parameters
        ----------
        sequences : iterable object
            An iterable object whose elements will be used as keys
        scoringmatrx : dict
            a scoring matrix for the comparissosns that will be carried out, defaults to scoringMatrices.blosum62

        Examples
        --------

        >>> something = BlosumDictGenerator('PEPTIDE')
        >>> print(something['D', 'E'])
        None
        >>> x = something ('D', 'E')
        >>> print(x)
        2.0
        >>> print(something['D', 'E'])
        2.0

        """
        SymetricDict.__init__(self, sequences)
        self.scoringmatrx = scoringmatrx
        self.saved_calculations = 0
        self.calculations = 0

    def __call__(self, seq1, seq2, verbose=False):
        if self[seq1, seq2] is None:
            if verbose:
                print("Calculating from {} to {}".format(seq1, seq2))
            self[seq1, seq2] = pairwise2.align.localdx(seq1, seq2, self.scoringmatrx, score_only=True)
        else:
            if verbose:
                print("Skipping calculation from {} to {}".format(seq1, seq2))
            self.saved_calculations += 1
        return self[seq1, seq2]

    # TODO add a way to count the number of NAs missing
    # TODO add a way to force the calculation of all the matrix ...


def Build_blosum_distance_dict(sequences, verbose=True, scoring_matrx=scoringMatrices.blosum62):
    """
    Gets a list of strings interpretable as sequences with the given dictionary and returns a nested dictionary with
    the distances between all pairs (as calculated by the smith-waterman algorithm)

    Parameters
    ----------
    sequences : iterable
        a list of strings interpretable as sequences
    verbose : logical
        wether the function should print progress output
    scoring_matrx: dict
        a scoring matrix object as implemented by biopython (defaults to blosum62)

    Returns
    -------
    out : dict
        a nested dictionary of the elements whose keys are the provided sequences and the values are the score

    Examples
    --------
    >>> Build_blosum_distance_dict(['PEP', 'PTY', 'TIDE'], verbose = False)
    SymetricDict([('PEP', OrderedDict([('PEP', 19.0), ('PTY', 7.0), ('TIDE', 5.0)])), \
('PTY', OrderedDict([('PEP', 7.0), ('PTY', 19.0), ('TIDE', 5.0)])), \
('TIDE', OrderedDict([('PEP', 5.0), ('PTY', 5.0), ('TIDE', 20.0)]))])
    """

    my_dict = SymetricDict(sequences)

    for x in my_dict.keys():
        for y in my_dict[x].keys():
            if my_dict[x, y] is not None:
                continue
            if verbose:
                print("Calculating from {} to {}".format(x, y))
            my_dict[x, y] = pairwise2.align.localdx(x, y, scoring_matrx, score_only=True)

    return my_dict
# TODO add something to incorporate your alphabet
# TODO make a way to lazily load the scores on demand and cache them ...


@functools.singledispatch
def getNmers(arg, n, fill_missing=False, default_alphabet=IUPAC.IUPACProtein.letters):
    """

    Parameters
    ----------


    Returns
    -------
    object

    Examples
    --------
    >>> getNmers('ABAB', 2, fill_missing = False, default_alphabet = 'ABC')
    OrderedDict([('AB', 1), ('BA', 1)])
    >>> getNmers('ABAB', 2, fill_missing = True, default_alphabet = 'ABC')
    OrderedDict([('AA', 0), ('AB', 1), ('AC', 0), ('BA', 1), ('BB', 0), ('BC', 0), ('CA', 0), ('CB', 0), ('CC', 0)])
    >>> getNmers('ABAB', 2, fill_missing = True, default_alphabet = 'AB')
    OrderedDict([('AA', 0), ('AB', 1), ('BA', 1), ('BB', 0)])
    >>> getNmers(['BAB', 'ABA'], 2, fill_missing = False, default_alphabet = 'AB')
    OrderedDict([('BAB', OrderedDict([('AB', 1), ('BA', 1)])), ('ABA', OrderedDict([('AB', 1), ('BA', 1)]))])
    >>> # Note that it processes each independently
    """
    # TODO add a way to check your alphabet and remove underscores (borders of sequences)
    my_nmers = [arg[i:i + n] for i in range(len(arg) - n + 1)]
    my_nmers = {x: x.count(x) for x in my_nmers}

    if fill_missing:
        all_posibles = getPermutations(default_alphabet, n)
        for k in all_posibles:
            my_nmers.setdefault(k, 0)

    my_nmers = collections.OrderedDict(sorted(my_nmers.items()))
    return my_nmers


@getNmers.register(list)
def _(arg, *args, **kwargs):
    my_nmer_list = collections.OrderedDict((x, getNmers(x, *args, **kwargs)) for x in arg)
    return my_nmer_list


def BuildNmerDistanceDict(n, alphabet=IUPAC.IUPACProtein.letters, fistn=0):
    """

    Parameters
    ----------
    n : length of the nmers that will be generated
    fistn : optional argument to specify that only certain number of combinations should be used

    Returns
    -------
    A nested dictionary with the combinations and scores

    Examples
    --------
    >>> BuildNmerDistanceDict(1, "PET")
    SymetricDict([('E', OrderedDict([('E', 5.0), ('P', 0.0), ('T', 0.0)])), \
('P', OrderedDict([('E', 0.0), ('P', 7.0), ('T', 0.0)])), \
('T', OrderedDict([('E', 0.0), ('P', 0.0), ('T', 5.0)]))])


    TODO
    ----
    Add an argument to stablish the alphabet being used
    Use our symetrical dict object as an output
    """
    perms = getPermutations(alphabet, n)
    # print(str(len(perms))) this prints the number of distances that will be calculated
    if fistn is not 0:
        perms = perms[0:fistn]
    distances = Build_blosum_distance_dict(perms, verbose=False)
    return distances


# TODO check why BuildNmerDistanceDict(2, 'PEP') would be different than (2, 'PEP')


def nmer_knn(sequences, nmer_size, k=2, default_alphabet=IUPAC.IUPACProtein.letters):
    """

    Parameters
    ----------
    sequences : sequences to be used to in splitting to nmers and clustered
    nmer_size : size of the nmers to build
    k : Number of neighbors to return
    default_alphabet : The default alpabet to be used

    Returns
    -------
    A nested list with the indexes of the nearest neighbors for each sequence provided

    Examples
    --------
    >>> testing_windows = [ "ABCDE", "DEFG", "HIJK", "JKHI", "ABKI" ]
    >>> foo = nmer_knn(testing_windows, 2, 2, default_alphabet='ABCDEFGHIJK')
    >>> print(foo)
    [[0 1]
     [1 0]
     [2 3]
     [3 2]
     [4 0]]

    """
    # TODO check for the performance of multiple distance metrics in this space
    # TODO check for the performance of multiple nmer sizes

    nmer_list = getNmers(sequences, nmer_size, fill_missing=True, default_alphabet=default_alphabet)
    myarray = np.array([
        [int(i) for i in elem.values()]
        for elem in nmer_list.values()
    ])
    kdt = KDTree(myarray, leaf_size=20, metric='manhattan')
    nn = kdt.query(myarray, k=k, return_distance=False)

    return nn


def knn_to_graph(nn_indeces, labels, labelname='Sequence'):
    """

    Parameters
    ----------
    nn_indeces
    labels
    labelname

    Returns
    -------

    """
    g = nx.Graph()
    g.add_nodes_from([x for x, _ in enumerate(nn_indeces)])

    edges_list = []
    for x in nn_indeces:
        for y in x[1:]:
            edges_list.append((x[0], y))

    g.add_edges_from(edges_list)

    for x in g.nodes:
        g.nodes[x][labelname] = labels[x]

    return g


def default_graph_plotting(g, labelname='Sequence', colour='r', filename='myfig.png'):
    """

    Parameters
    ----------
    g : A graph object
    labelname : Node atributes to be used as labels
    colour : Colour specifications of the nodes
    filename : file name to which the plot will be saved
    """
    pos = nx.spring_layout(g, k=1/3.5, scale=0.5)
    mpl.figure(num=None, figsize=(10, 10), dpi=150, facecolor='w', edgecolor='k')
    nx.draw(g, pos, node_color=colour, labels=nx.get_node_attributes(g, labelname), font_size=12, with_labels=True)
    if filename is not None:
        mpl.savefig(filename)


def shuffle_sequence(seq):
    newseq = ''.join(random.sample(seq, len(seq)))
    return newseq


def invert_sequence(seq):
    newseq = seq[::-1]
    return newseq


@functools.singledispatch
def flatten(iterable):
    yield iterable


@flatten.register(list)
@flatten.register(set)
@flatten.register(tuple)
def _(iterable):
    for i in iterable:
        yield from flatten(i)


"""
# Dummy Example

testing_windows = [
    "ABCDE", "DEFG", "HIJK",
    "JKHI","ABKI",]

foo = [
    getNmers(x, 2, fill_missing = True, default_alphabet= 'ABCDEFGHIJK')
    for x in testing_windows
]

myarray = np.array([
    [ int(i) for i in elem.values()]
    for elem in foo
])

kdistances = KDTree(myarray, leaf_size = 20, metric = 'manhattan')
mydistances = kdistances.query(myarray, k =2, return_distance = False) # Returns the closest 2 neighbors to each

# TODO make this a function
G = nx.Graph()
G.add_nodes_from([ x for x, _ in enumerate(mydistances)])

edges_list = []
for x in mydistances:
    for y in x[1:]:
        edges_list.append((x[0], y))

G.add_edges_from(edges_list)

for x in G.nodes:
    G.nodes[x]['Sequence'] = testing_windows[x]

# END TODO make this a function

pos = nx.spring_layout(G, k=1/3.5, scale = 0.5)
mpl.figure(num=None, figsize=(5, 5), dpi=150, facecolor='w', edgecolor='k')
nx.draw(G, pos, labels = nx.get_node_attributes(G, 'Sequence'), font_size=12, with_labels=True)
mpl.savefig("mygraph.png")


"""
