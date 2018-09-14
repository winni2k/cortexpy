def get_shared_argparse():
    import argparse
    shared_parser = argparse.ArgumentParser(add_help=False)
    group = shared_parser.add_mutually_exclusive_group()
    group.add_argument('-v', '--verbose', help='Increase log level to debug', action='store_true')
    group.add_argument('-s', '--silent', help='Decrease log level to warnings and errors',
                       action='store_true')
    shared_parser.add_argument('-o', '--out', required=False, default='-',
                               help="Output cortexpy graph. '-' writes to stdout")

    return shared_parser
