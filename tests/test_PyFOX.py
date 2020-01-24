import os
import pyfocs
import sh
import yaml
import pytest

path = os.path.dirname(os.path.abspath(__file__))
path_data = os.path.join(path, 'data')
path_data_external = os.path.join(path, 'data', 'multicore_demo', 'external')
example_config = os.path.join(path, 'data', 'example_configuration.yml')


def test_example_config_exists():
    '''
    Test that the example config exists.
    '''
    assert os.path.isfile(example_config)
    pass


def test_example_runs_on_pyfox():
    '''
    Run the actual PyFOX.py script on the example data.
    '''

    # Read the config file
    cfg = pyfocs.yamlDict(example_config)

    # Alter the two entries to work for an arbitrary system.
    cfg['directories']['local']['dir_pre'] = path_data
    cfg['directories']['external'] = path_data_external

    # Write to a temp file
    stream = open(os.path.join(path_data, 'temp.yml'), 'w')
    yaml.dump(cfg, stream)

    try:
        sh.python(['PyFOX.py', os.path.join(path_data, 'temp.yml')])
    except sh.ErrorReturnCode as e:
        print(e)
        pytest.fail(e)
    finally:
        os.remove(os.path.join(path_data, 'temp.yml'))
