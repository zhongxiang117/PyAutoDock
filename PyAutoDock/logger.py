import inspect
import logging


def mylogger():
    """Configure base logger"""
    frame = inspect.stack()[1]
    module_name = inspect.getmodule(frame[0]).__name__
    if module_name != '__main__':
        logger = logging.getLogger(module_name)
        if not logger.parent.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s [%(levelname)s] [%(filename)s:%(lineno)d]: %(message)s'
            )
            handler.setFormatter(formatter)
            logger.parent.addHandler(handler)
    else:
        logger = logging.getLogger('PyAutoDock')
    return logger


