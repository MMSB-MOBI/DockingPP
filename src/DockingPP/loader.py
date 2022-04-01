import re
import DockingPP.error as error
from DockingPP.rotation_utils import trans_matrix, eulerFromMatrix
import logging
from DockingPP.dockingHandlerZdock import DockingHandlerZdock
from DockingPP.dockingHandlerPDB import DockingHandlerPDB
import DockingPP.typecheck as typecheck
from typing import List
import pypstruct.coordinates as PDB
import ccmap


def loadZdock(zdock_results: str, nb_pose: int = -1) -> 'DockingPP.dockingHandler.DockingHandler':
    """Load zdock results into DockingHandler object

    Args:
        zdock_results (str): path to zdock results file
        nb_pose (int, optional): number of poses to load. If -1, all poses will be loaded. Defaults to -1.

    Raises:
        error.ZdockFormatError: Raise if error is detected in zdock format
        error.IncompatiblePoseNumber: Raise if you try to load more poses than possible.

    Returns:
        DockingPP.dockingHandler.DockingHandler: DockingHandler object for this zdock results.

    """
    typecheck.validFile(zdock_results)  # Check if the file exists
    logging.info(
        f"== Load zDock results ==\n path : {zdock_results}\n number of poses : {'All' if nb_pose == -1 else nb_pose}")
    reL1 = r'^([\d]+)[\s]+([\d\.]+)[\s]*$'
    reL2 = r'^[\s]*([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'
    reL3 = r'^[\s]*([\S]+)[\s]+([\.\d-]+)[\s]+([\d\.-]+)[\s]+([\.\d-]+)[\s]*$'

    reZPOSE = r'^([\d\.-]+)[\s]+([\d\.-]+)[\s]+([\d\.-]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+([\d]+)[\s]+.*'

    with open(zdock_results, 'r') as f:
        # Check header lines formatting
        re_line1 = re.match(reL1, f.readline())
        re_line2 = re.match(reL2, f.readline())
        re_line3 = re.match(reL3, f.readline())
        re_line4 = re.match(reL3, f.readline())

        if not re_line1:
            raise error.ZdockFormatError("Line 1 has wrong format")
        if not re_line2:
            raise error.ZdockFormatError("Line 2 has wrong format")
        if not re_line3:
            raise error.ZdockFormatError("Line 3 has wrong format")
        if not re_line4:
            raise error.ZdockFormatError("Line 4 has wrong format")

        # Set docking_collection attributes

        grid_dimension = int(re_line1.groups()[0])
        step = float(re_line1.groups()[1])
        initial_euler = (float(re_line2.groups()[0]), float(
            re_line2.groups()[1]), float(re_line2.groups()[2]))
        baryRec = (float(re_line3.groups()[1]), float(
            re_line3.groups()[2]), float(re_line3.groups()[3]))
        baryLig = (float(re_line4.groups()[1]), float(
            re_line4.groups()[2]), float(re_line4.groups()[3]))

        docking_collection = DockingHandlerZdock(
            grid_dimension, step, initial_euler, baryRec, baryLig)

        pose_index = 1
        # Parse poses lines
        for line in f:
            m = re.match(reZPOSE, line)
            if not m:
                raise error.ZdockFormatError("A pose line has wrong format")
            euler = (float(m.groups()[0]), float(
                m.groups()[1]), float(m.groups()[2]))

            if initial_euler != (0, 0, 0):  # Rotation has to be applied
                # Make rotation matrices
                rand_rot = trans_matrix(*initial_euler)
                pose_rot = trans_matrix(*euler)
                # Combine into one matrix
                double = pose_rot.dot(rand_rot)
                # Recover combined angles
                euler = eulerFromMatrix(double)

            _translation = [int(m.groups()[3]), int(
                m.groups()[4]), int(m.groups()[5])]
            translation = tuple(
                [t - grid_dimension if t > grid_dimension / 2 else t for t in _translation])
            # Copy Julia's calculations
            translation = tuple([-1 * t * step for t in translation])

            docking_collection.addPose(pose_index, euler, translation)

            if nb_pose == pose_index:
                return docking_collection

            pose_index += 1

        if nb_pose == -1:
            return docking_collection

        raise error.IncompatiblePoseNumber(
            f"You ask too much poses, only {pose_index} are present in the result file.")


def loadPDB(pairedList: str, directorypath: str = "") -> 'DockingPP.dockingHandler.DockingHandler':
    """Load Poses from a list of pdb files 

    Args:
        pairedList (str): path to a file with a list of ligand receptor
                            exemple : 
                                rec_refined_model_round1.BL00090001.pdb lig_refined_model_round1.BL00090001.pdb
                                rec_refined_model_round1.BL03860001.pdb lig_refined_model_round1.BL03860001.pdb
                                rec_refined_model_round1.BL02460001.pdb lig_refined_model_round1.BL02460001.pdb
                                rec_refined_model_round1.BL04120001.pdb lig_refined_model_round1.BL04120001.pdb
                                rec_refined_model_round1.BL02820001.pdb lig_refined_model_round1.BL02820001.pdb
        directorypath (str): path to paired files,  exemple : "Test_case/" , default value = ""

    Raises:
        error.ZdockFormatError: Raise if error is detected in zdock format
        error.IncompatiblePoseNumber: Raise if you try to load more poses than possible.

    Returns:
        DockingPP.dockingHandler.DockingHandler: DockingHandler object for this zdock results.

    """

    # Check if the file exists
    typecheck.validFile(pairedList)

    # docking_collection = DockingHandler( 0 , 0, (0, 0, 0), (0, 0, 0),(0, 0, 0) )
    nb_pose = len(open(pairedList, 'r').readlines())
    docking_collection = DockingHandlerPDB(nb_pose)

    pose_index = 1

    with open(pairedList, 'r') as f:
        for line in f:

            pair = re.split(' | \t', line.strip("\n"))
            typecheck.validFile(directorypath+pair[0])
            typecheck.validFile(directorypath+pair[1])

            # parse the first pdb "receptor"
            rec = directorypath+pair[0]
            docking_collection.addReceptor(rec)

            # parse the second pdb "ligand"
            lig = directorypath+pair[1]
            docking_collection.addLigand(lig)

            print(f"Pose #{pose_index} : [{lig.split('/')[-1]}, {rec.split('/')[-1]}]")
            docking_collection.addPose(pose_index)

            if nb_pose == pose_index:
                print("Done for", nb_pose, "poses")
                return docking_collection

            pose_index += 1
