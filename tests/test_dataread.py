import os
# import pyfocs
from pyfocs import xml_read


# Paths to example data
path = os.path.dirname(os.path.abspath(__file__))
path_data = os.path.join(path, 'data')
mc_demo = 'multicore_demo'
sf_demo = 'single_file_demo'
channel = 'channel 1'
sfn = 'channel 1_20190722000003996.xml'

path_data_singlefile = os.path.join(path_data, sf_demo, 'raw_xml', channel, sfn)
path_data_multifile = os.path.join(path_data, mc_demo)
yaml_file = os.path.join(path_data, 'example_configuration.yml')


def test_example_data_exists():
    '''
    Test that the example data is found along the hard coded paths.
    '''
    # The multicore data directory exists.
    assert os.path.isdir(path_data_multifile)

    # The single file we will be working with exists.
    assert os.path.isfile(path_data_singlefile)
    pass


def test_xmlread():
    '''
    Test the ability to read a single xml file and return the expected data.
    '''
    # Did we read the file without getting a CorruptedXMLError exception?
    assert xml_read(path_data_singlefile)
