import os
import pandas as pd


def results(filenamePrefix, filenameSuffix):

    # Go to data directory
    os.chdir(dir_bmmflux_results_CSAT)

    # File names
    allResults = [filenamePrefix + '_results_' + filenameSuffix + '.csv']

    # Units (names + units = first two rows)
    results_cols = pd.read_csv(allResults[0], header=0, delimiter=',',
                               error_bad_lines=False, warn_bad_lines=False, nrows=2)
    results_units = results_cols.iloc[0]

    # Read all data
    resultsList = []

    # Loop through files
    for rFiles in allResults:
        df = pd.read_csv(rFiles, header=None, delimiter=',', skiprows=[0,1],
                         names=results_cols.columns.values, index_col=False)
        # Drop datetime_end column
        df.drop(['Datetime_start', 'Datetime_end', 'DOY'],
                axis=1, inplace=True)
        resultsList.append(df)

    # Concatenate into a single Pandas Dataframe
    bmmflux_df = pd.concat(resultsList)

    # Index using the datetime at the center of the interval (converted to a
    # pandas datetime object).
    bmmflux_df = bmmflux_df.reset_index(drop=True).set_index([matlabdn2datetime(bmmflux_df['Datetime_center'].values)])
    bmmflux_df.drop(['Datetime_center'], axis=1, inplace=True)

    # Files can get out of place, so resort the data
    bmmflux_df.sort_index(inplace=True)

     # Save units as separate dictionary
    for ru in results_units.keys():
        if '[1]' in results_units[ru]:
            results_units[ru] = 'unitless'

    return bmmflux_df, results_units
