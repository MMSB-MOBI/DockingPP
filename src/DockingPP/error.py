class ZdockFormatError(Exception):
    """Raised when zdock results have wrong format
    """
    pass

class IncompatiblePoseNumber(Exception):
    """Raised when the asked pose number is incompatible with data
    """
    pass

class InvalidRole(Exception):
    """Raised when invalid role (ligand|receptor) is provided
    """
    pass

class InvalidScore(Exception):
    """Raised when invalid score is provided
    """
    pass

class InvalidArgument(Exception):
    """Raised when invalid argument provided in command line
    """
    pass

class ClustersNotComputed(Exception):
    """Raised when clusters are needed and not computed
    """
    pass

class RescoringNotComputed(Exception):
    """Raised when rescoring is needed and not computed
    """
    pass

class PdbNotSet(Exception):
    """Raised when a pdb file is not set"""
    pass

class ContactMapNotComputed(Exception):
    """Raised when contact map is needed and not computed
    """
    pass

class FrequenciesNotComputed(Exception):
    """Raised when frequencies are needed and not computed
    """
    pass
