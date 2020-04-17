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