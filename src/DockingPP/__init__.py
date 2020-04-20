from .loader import loadZdock
import logging
import sys

logging.basicConfig(level = logging.INFO, format='%(levelname)s\t%(filename)s\t%(message)s', stream=sys.stdout)