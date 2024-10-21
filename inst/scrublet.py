import matplotlib
matplotlib.use('agg')
import scrublet as scr
import scipy.io
import numpy as np
import os
from scipy.sparse import csc_matrix

def scrublet_py(i, j, val, dim, sim_doublet_ratio = 2.0, n_neighbors = 1,
                expected_doublet_rate = 0.1, stdev_doublet_rate = 0.002,
                random_state = 0, min_counts = 3, min_cells = 3,
                min_gene_variability_pctl = 85.0, n_prin_comps = 50):

    data = csc_matrix((val, (i, j)), shape = dim)
    scrub = scr.Scrublet(data, sim_doublet_ratio=int(sim_doublet_ratio), n_neighbors=int(n_neighbors), expected_doublet_rate=expected_doublet_rate, stdev_doublet_rate = stdev_doublet_rate, random_state = random_state)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=int(min_counts), min_cells=int(min_cells), min_gene_variability_pctl=min_gene_variability_pctl, n_prin_comps=int(n_prin_comps))
    return scrub

def call_scrublet_doublets(scrub_obj, threshold = 0.3):
    scrub_obj.call_doublets(threshold=threshold)
    return scrub_obj
