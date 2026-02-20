import numpy as np
import json
import math
from .nucleotide import Nucleotide, reset_nucleotide_index
from .strand import Strand, reset_strand_index


class System:
    def __init__(self, box):
        self._box = np.array(box, dtype=float)
        self._strands = []
        self._N = 0
        self._N_strands = 0
        self.isDNA = True
        reset_nucleotide_index()
        reset_strand_index()

    @property
    def N(self):
        return self._N

    @property
    def N_strands(self):
        return self._N_strands

    def add_strand(self, strand):
        self._strands.append(strand)
        self._N += strand.N
        self._N_strands += 1
        return True

    def _prepare(self):
        idx = 0
        for i, strand in enumerate(self._strands):
            idx = strand._prepare(i, idx)

    def calc_clusters(self, threshold=1.0018):
        """Compute connected clusters based on backbone distance."""
        breaks = set()
        prev_map = {}
        next_map = {}

        for strand in self._strands:
            for i in range(len(strand._nucleotides) - 1):
                n1 = strand._nucleotides[i]
                n2 = strand._nucleotides[i + 1]
                prev_map[id(n2)] = n1
                next_map[id(n1)] = n2
                diff = n1.distance(n2)
                if math.sqrt(np.dot(diff, diff)) > threshold:
                    breaks.add(id(n1))
                    breaks.add(id(n2))

        if len(breaks) == 0:
            # No breaks: everything is cluster 1
            for strand in self._strands:
                for nuc in strand._nucleotides:
                    nuc.cluster = 1
            return

        def neighbors(nuc):
            result = []
            if nuc.pair is not None:
                result.append(nuc.pair)
            p = prev_map.get(id(nuc))
            if p is not None:
                result.append(p)
            n = next_map.get(id(nuc))
            if n is not None:
                result.append(n)
            return result

        def merge_clusters(c1, c2):
            target = min(c1, c2)
            for strand in self._strands:
                for nuc in strand._nucleotides:
                    if nuc.cluster == c1 or nuc.cluster == c2:
                        nuc.cluster = target

        # Assign initial cluster IDs to break nucleotides
        break_nucs = []
        for strand in self._strands:
            for nuc in strand._nucleotides:
                if id(nuc) in breaks:
                    break_nucs.append(nuc)

        cluster_id = 1
        for nuc in break_nucs:
            nuc.cluster = cluster_id
            cluster_id += 1

        # BFS propagation
        for nuc in break_nucs:
            stack = [nuc]
            while stack:
                current = stack.pop()
                for nbr in neighbors(current):
                    if nbr.cluster != current.cluster:
                        if nbr.cluster is None:
                            nbr.cluster = current.cluster
                            stack.append(nbr)
                        elif not (id(current) in breaks and id(nbr) in breaks):
                            merge_clusters(current.cluster, nbr.cluster)

    def to_oxview_dict(self):
        """Generate oxView JSON as a Python dict."""
        self._prepare()
        box = np.round(self._box).astype(int).tolist()

        strands_out = []
        for strand in self._strands:
            monomers = []
            for i, nuc in enumerate(strand._nucleotides):
                m = {
                    'id': nuc.index,
                    'type': nuc.get_base(),
                    'class': 'DNA' if self.isDNA else 'RNA',
                    'p': nuc.cm_pos.tolist(),
                    'a1': nuc._a1.tolist(),
                    'a3': nuc._a3.tolist(),
                }
                # n3 neighbor
                if strand._circular:
                    n3 = (strand._nucleotides[-1].index if i == 0
                           else strand._nucleotides[i - 1].index)
                    n5 = (strand._nucleotides[0].index
                           if i == len(strand._nucleotides) - 1
                           else strand._nucleotides[i + 1].index)
                else:
                    n3 = -1 if i == 0 else strand._nucleotides[i - 1].index
                    n5 = (-1 if i == len(strand._nucleotides) - 1
                           else strand._nucleotides[i + 1].index)

                if n3 >= 0:
                    m['n3'] = n3
                if n5 >= 0:
                    m['n5'] = n5
                if nuc.pair is not None:
                    m['pair'] = nuc.pair.index
                if nuc.cluster is not None:
                    m['cluster'] = nuc.cluster
                if nuc.color is not None:
                    m['color'] = nuc.color
                monomers.append(m)

            strand_dict = {
                'id': strand.index,
                'end3': strand._nucleotides[0].index,
                'end5': strand._nucleotides[-1].index,
                'class': 'NucleicAcidStrand',
                'monomers': monomers,
            }
            strands_out.append(strand_dict)

        return {
            'box': box,
            'systems': [{'id': 0, 'strands': strands_out}]
        }

    def to_oxview_string(self):
        return json.dumps(self.to_oxview_dict())
