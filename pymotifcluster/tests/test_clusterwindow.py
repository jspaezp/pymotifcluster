from pymotifcluster.clusterwindow import *
# import pytest
# import pytest_benchmark

testing_windows = [
"DSSASPEVVSDLPPSSPKGSPDRHDPSTSSP",
"VREQAVWALGNVAGDSPKCRDLVLSYGAMTP",
"LTSPIPRASVITNQNSPLSSRATRRTSISSG",
"VTPCKGSGDRSLGLTSPIPRASVITNQNSPL",
"AIKASSLSKEGLLFGSPKLSGGSSLYGKLQQ",
"GSFRKNLDTKDAIISSPKLMESGSGKLPVFE",
"SSIASVPITDTTHVKSETGSPRHSSSAKMDE",
"SVPITDTTHVKSETGSPRHSSSAKMDETNGR",
"GSLSKSPSFDWGEDVSPNIPLEKLLVYRLNV",
"DMSSIDGKETSRSGGSPNRAELRKRLSAAEE",
"FKSVFTEDLDPPETESESDSPKHSEEHEHPE",
"SVFTEDLDPPETESESDSPKHSEEHEHPEQE",
"FTEDLDPPETESESDSPKHSEEHEHPEQEHP",
"TGRLSPQTFTSSPSLSPSSSSGGSSFMARYA",
"SPQTFTSSPSLSPSSSSGGSSFMARYAMESS",
"PQTFTSSPSLSPSSSSGGSSFMARYAMESSK",
"NLPGNPDPEAEVIALSPKTLMATNRFLCEIC",
"SPRFSRQRTSPPSLHSPLRSLKEPKRQLIPQ",
"ESAASESGEKADEGRSQVDGSTEQSPKLESA",
"ESGEKADEGRSQVDGSTEQSPKLESASSTEP",
]


def test_SymetricDict():
    my_symdict = SymetricDict('PEPTIDES')
    # This tests that all values can be set appropiately
    for i, x  in enumerate(itertools.combinations_with_replacement(set('PEPTIDES'),2)):
        my_symdict[x] = i
        assert my_symdict[x] == i, "Value not being set appropiately"
        assert my_symdict[x[::-1]] == i, "Value not being set appropiately"

    for x in set('PEPTIDES'):
        for y in set('PEPTIDES'):
            assert my_symdict[x, y] == my_symdict[y, x], "Symdict is not symetrical"

def test_getPermutations():
    my_perm = getPermutations('123', 2)
    assert len(my_perm) == 9, 'Unexpected length of the permutations'


def test_getNmers():
    n_mers1 = getNmers("AAAAAAA", 2)
    assert len(n_mers1) == 1, 'More than 1 nmer for a sequence with only one unique character'
    n_mers2 = getNmers("AABAAAB", 3)
    assert len(n_mers2) == 4, 'Unexpected number of nmers, check manually'
    n_mers3  = getNmers("AABAAAB", 2, fill_missing=True, default_alphabet='ABC')
    print(str(len(n_mers3)))
    print(str(len(getPermutations('ABC', 2))))
    assert len(n_mers3) == len(getPermutations('ABC', 2)), 'Not filling nmers when asked'
    # TODO implement cropping border delimiters in getNmers
    #n_mers3 = getNmers("__AAAAA", 3)
    #assert len(n_mers1) == 1, 'More than 1 nmer for a sequence with only one unique character and border delimiters'

def test_getNmer_method_dispatching():
    assert len(getNmers('abc', 2)) == 2
    assert set(getNmers('abc', 2)) == set(getNmers('abc', 2))
    assert len(getNmers(['abc', 'cde'], 2)) == 2
    assert set(getNmers(['abc', 'cde'], 2)['abc']) == set(getNmers('abc', 2))
    assert set(getNmers(['abc', 'cde'], 2)['cde']) == set(getNmers('cde', 2))

def _test_trimer_clustering(testing_windows):
    foo = [getNmers(x, 3, fill_missing = True) for x in testing_windows]
    myarray = np.array([[ int(i) for i in elem.values()] for elem in foo])
    kdt = KDTree(myarray, leaf_size = 20, metric = 'manhattan')
    foo2 = kdt.query(myarray, k =2, return_distance = False)

def test_trimer_clustering_20(benchmark):
    benchmark(_test_trimer_clustering, testing_windows[0:10])

def test_trimer_clustering_20(benchmark):
    benchmark(_test_trimer_clustering, testing_windows)

def test_trimer_clustering_40(benchmark):
    benchmark(_test_trimer_clustering, testing_windows*2)

def test_trimer_clustering_60(benchmark):
    benchmark(_test_trimer_clustering, testing_windows*3)
