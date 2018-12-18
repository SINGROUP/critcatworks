import logging

# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

logger=logging.getLogger("nomadcore")
logger.setLevel(logging.WARNING)
logger.addHandler(ch)

logger2=logging.getLogger("nomad")
logger2.setLevel(logging.WARNING)
logger2.addHandler(ch)

def debugToFile():
    "makes a full log to a file named detailed.log"
    fh = logging.FileHandler('detailed.log')
    fh.setLevel(logging.DEBUG)
    logger.setLevel(logging.DEBUG)
    logger2.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    logger2.addHandler(fh)
