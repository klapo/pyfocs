import os
# import pyfocs
from pyfocs import xml_read


# Multicore file example
path = os.path.dirname(os.path.abspath(__file__))
path_data_arch = os.path.join(path, 'data', 'multicore_demo')

# Single file
fn = 'archived/channel 1_20190722-0000/channel 1_20190722000003996.xml'
path_data_singlefile = os.path.join(path_data_arch, fn)

def test_example_data_exists():
    '''
    Test that the example data is found along the hard coded paths.
    '''
    # The data directory exists.
    assert os.path.isdir(path_data_arch)

    # The single file we will be working with exists.
    assert os.path.isfile(path_data_singlefile)
    pass


def test_xmlread():
    '''
    Test the ability to read a single xml file and return the expected data.
    '''
    # Did we read the file without getting a CorruptedXMLError exception?
    assert xml_read(path_data_singlefile)
