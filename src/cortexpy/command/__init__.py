import os


def get_exit_code_yaml_path():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), 'exit_codes.yaml'))
