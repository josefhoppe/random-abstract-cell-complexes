# Adapted from https://github.com/josefhoppe/edge-flow-cell-complexes.
# Licensed under MIT License (https://github.com/josefhoppe/edge-flow-cell-complexes/blob/main/LICENSE)

# add code to import path
from pathlib import Path
import sys, resource
sys.path.append(str((Path(__file__).parent.parent.parent.parent).absolute()))

from cell_flower.detection import CellCandidateHeuristic
from edge_flow_cell_complexes.experiment import GROUND_TRUTH, RANDOM

method_map = {
    'max': CellCandidateHeuristic.MAX,
    'triangles': CellCandidateHeuristic.TRIANGLES,
    'ground_truth': GROUND_TRUTH,
    'similarity': CellCandidateHeuristic.SIMILARITY,
    'random': RANDOM
}
