import logging

log_level_conversion_table = {
    0: 'NOTSET',
    10: 'DEBUG',
    20: 'INFO',
    30: 'WARNING',
    40: 'ERROR',
    50: 'CRITICAL',
}


def configure_logging_from_args(args):
    """Checks args for --silent and --verbose mode and sets the logging level appropriately"""
    if getattr(args, 'verbose', None):
        log_level = logging.DEBUG
    elif getattr(args, 'silent', None):
        log_level = logging.WARNING
    else:
        log_level = logging.INFO
    logging.basicConfig(level=log_level)
    logger = logging.getLogger('cortexpy')
    logger.info('Log level is {}'.format(log_level_conversion_table[logger.getEffectiveLevel()]))


def configure_logging_from_args_and_get_logger(args, logger_name):
    configure_logging_from_args(args)
    return logging.getLogger(logger_name)
