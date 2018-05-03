A Python package to cluster peptides
====================================

NOTE: Originally the intention was to have another way to find sequence similarities
among phosphorylation sites, nonetheless, in the current state it finds similarities 
between arbitrary sequences.

This program is fairly simple, given the set of peptides and the parameters, it
calls the clustering package internally and returns a json with the network of
the clustered peptides, which is then visualized using D3.js

This was more of a learning exercise for me (which took some sweet time ...) but there
is no reason to stop working on the project.

How to install
----------

.. code:: shell
  git clone https://github.com/jspaezp/pymotifcluster
  pip install -e .

How to run
----------

.. code:: python3

  from pymotifcluster.clusterwindow import shuffle_sequence, nmer_knn, knn_to_graph

  foregroundSequences = [
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


  # We generate decoy sequences by shuffling the foreground sequences (n times the size of the foreground)
  n=1
  decoys = [shuffle_sequence(x) for x in fg_sequence*n]

  # We add everything to a couple of lists
  sequences = []
  group = []

  sequences.extend(fg_sequence)
  sequences.extend(decoys)
  group.extend(['fg']*len(fg_sequence))
  group.extend(['decoys']*len(decoys))

  # Here we generate the clustering, connecting k nearest neighbors
  k = 2
  nmer_size = 3

  foo = nmer_knn(sequences, nmer_size=nmer_size, k=k)

  # And here we make it a graph (network) object
  G = knn_to_graph(foo, sequences)
  G.add_nodes_from([x for x, y in enumerate(group) if y == 'fg'], group=1)
  G.add_nodes_from([x for x, y in enumerate(group) if y == 'decoys'], group=2)

  # Here we can use networkx to plot it or ... anything else ...

TODO
====

- Implement an efficient way to provide weights to the nmers
- Implement a way to extract the clusters
- Implement a way to get whch nmers are the ones predominantly responsible for the grouping
- refactor the structure of the package
- clean out the dependencies (currently it depends in waaaay too many packages and uses only 1 or 2 functions from each)
- MANY MORE
