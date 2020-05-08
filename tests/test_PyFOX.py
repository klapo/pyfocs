import os
import pyfocs
import sh
import yaml
import pytest

path = os.path.dirname(os.path.abspath(__file__))
path_data = os.path.join(path, 'data')
path_data_external = os.path.join(path, 'data', 'multicore_demo', 'external')

# Dictionary of configuration files
example_config_ss = os.path.join(path, 'data', 'example_configuration_steelfiber.yml')

example_dict = {'stainless steel': example_config_ss}

def test_example_config_exists():
    '''
    Test that the example config exists.
    '''
    assert os.path.isfile(example_config)
    pass


def test_pyfoc_LOVE19_examples():
    '''
    Run the actual PyFOX.py script on the example data from the LOVE19
    campaign. These are all single-ended examples.
    '''

    for ex_name, ex_config in example_dict.items():
        # Read the config file
        cfg = pyfocs.yamlDict(ex_config)

        # Alter the two entries to work for an arbitrary system.
        cfg['directories']['local']['dir_pre'] = path_data
        cfg['directories']['external'] = path_data_external

        # Write to a temp file
        stream = open(os.path.join(path_data, 'temp.yml'), 'w')
        yaml.dump(cfg, stream)

        # This block of code assumes that any problems will throw an error back to
        # python from the shell. This is not necessarily true. But, this also
        # explains why the code coverage report is so poor.
        try:
            sh.python(['PyFOX.py', os.path.join(path_data, 'temp.yml')])
        except sh.ErrorReturnCode as e:
            print(ex_name + ' had a problem.')
            print(e)
            pytest.fail(e)
        finally:
            os.remove(os.path.join(path_data, 'temp.yml'))
