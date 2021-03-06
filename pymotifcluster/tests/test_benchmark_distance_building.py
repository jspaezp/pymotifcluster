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


def test_buildSequencesDistanceDict_10(benchmark):
    benchmark(Build_blosum_distance_dict, testing_windows[0:10])

def test_buildSequencesDistanceDict_20(benchmark):
    benchmark(Build_blosum_distance_dict, testing_windows)

def test_buildSequencesDistanceDict_40(benchmark):
    benchmark(Build_blosum_distance_dict, testing_windows*2)

def test_buildSequencesDistanceDict_60(benchmark):
    benchmark(Build_blosum_distance_dict, testing_windows*3)
