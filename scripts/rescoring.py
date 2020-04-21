"""Rescore zDock poses 

Usage:
    rescoring.py --list <complex_list> --pdb <pdb_folder> --results <results_folder> --nb_frequencies <int> --nb_rescoring <int> --score <str>

Options:
    -h --help
    --list <complex_list> file with list of complexes names
    --pdb <pdb_folder> folder containing ligand and receptor pdbs for complexes
    --results <results_folder> folder containing zdock results for complexes
    --nb_frequencies <int> Number of poses to compute frequencies
    --nb_rescoring <int> Number of poses to rescore
    --score <str> Type of score. [default: all]

"""

from docopt import docopt

if __name__ == "__main__":
    ARGS = docopt(__doc__, version="1.0.0")
    print(ARGS)
