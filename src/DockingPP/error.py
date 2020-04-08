class ZdockFormatError(Exception):
    """Raised when zdock results have wrong format
    """
    pass

class IncompatiblePoseNumber(Exception):
    """Raised when the asked pose number is incompatible with data
    """
    pass