"""Rescore zDock poses 

Usage:
    rescoring.py --list <complex_list> --pdb <pdb_folder> --results <results_folder> --nb_frequencies <int> --nb_rescoring <int> --output <output_dir> [ --score <str> ] [ --thread <int> ]

Options:
    -h --help
    --list <complex_list> file with list of complexes names
    --pdb <pdb_folder> folder containing ligand and receptor pdbs for complexes
    --results <results_folder> folder containing zdock results for complexes
    --nb_frequencies <int> Number of poses to compute frequencies
    --nb_rescoring <int> Number of poses to rescore
    --output <output_dir> Output folder for write score
    --score <str> Type of score. [default: all available scores]
    --thread <int> Number of threads to use to compute contact map [default: 8]

"""

from docopt import docopt
import logging
import os
import DockingPP
import time
logging.basicConfig(level = logging.INFO, format='%(levelname)s\t%(filename)s:%(lineno)s\t%(message)s')

def args_gestion():
    args = docopt(__doc__, version="1.0.0")

    err = ""

    if not os.path.isfile(args["--list"]):
        err += f'{args["--list"]} doesn\'t exist\n'
    
    for folder in [args["--pdb"], args["--results"]]:
        if not os.path.isdir(folder):
            err += f"{folder} is not a directory\n"

    if not os.path.isdir(args["--output"]):
       os.mkdir(args["--output"]) 
    else:
        logging.warn(f'--output : {args["--output"]} already exists')
    
    args["--score"] = args["--score"] if args["--score"] else "all"
    args["--thread"] = args["--thread"] if args["--thread"] else 8

    for key in ["--thread", "--nb_frequencies", "--nb_rescoring"]:
        try:
            args[key] = int(args[key])
        except ValueError:
            err += f'{key}: {args[key]} is not int\n' 
    
    if err:
        raise DockingPP.error.InvalidArgument(err)

    return args
        

if __name__ == "__main__":
    start = time.time()
    ARGS = args_gestion()

    #THINK ABOUT THAT, WE NEED TO FIX SUFFIX AND WRITE IN DOC
    ligand_suffix = "_l_u.pdb" 
    receptor_suffix = "_r_u.pdb"
    results_suffix = ".zd3.0.2.fg.fixed.out"
    
    max_poses = max(ARGS["--nb_rescoring"], ARGS["--nb_frequencies"]) # Use the maximum number of poses to load zdock results and compute contact map

    with open(ARGS["--list"]) as complexes:
        for compl in complexes:
            compl = compl.rstrip()
            logging.info(f"== Complex {compl}")
            #Load zdock results
            DH = DockingPP.loadZdock( f'{ARGS["--results"]}/{compl}{results_suffix}', max_poses)
            #Set receptor and ligand
            DH.setReceptor (f'{ARGS["--pdb"]}/{compl}{receptor_suffix}')
            DH.setLigand (f'{ARGS["--pdb"]}/{compl}{ligand_suffix}')
            DH.computeContactMap(ARGS["--thread"], max_poses )
            DH.computeFrequencies(ARGS["--nb_frequencies"])
            DH.rescorePoses(ARGS["--nb_rescoring"], ARGS["--score"])
            DH.serializeRescoring(f'{ARGS["--output"]}/{compl}_scores.tsv', ARGS["--score"])
        
    logging.info(f"END in {time.time() - start} s")