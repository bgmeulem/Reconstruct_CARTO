"""File containing function to convert a pts file to add the point index as scalar data.
This makes paraview operations easier."""

import argparse
from tqdm import tqdm
import glob


def convert(filename):  # pts meshfile to paraview csv
    """
    Takes a pts file, adds point ID as extra data to each point, writes to csv file
    :param filename: the name of the .pts mesh file
    :return: 0 if successful
    """

    ofname = filename.split(".")[0]
    outfile = open(ofname + "_paraview.csv", "w+")
    df = open(filename)

    print('\n\tWriting {}_paraview.csv'.format(ofname))
    # csv header
    outfile.write("X,Y,Z,meshID\n")
    i = 0
    for line in tqdm(df.readlines()[1:], desc='        '):
        for e in line[:-2].split():  # x, y or z co-ordinate
            outfile.write(str(e) + ',')  # write co-ordinates
        outfile.write(str(i)+'\n')  # add point ID
        i += 1

    outfile.close()
    df.close()
    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # fetching the arguments
    parser.add_argument('--filename',
                        help="name of .pts file, to be converted to paraview .csv (including relative directory)",
                        type=str, default="")
    args = parser.parse_args()
    if not args.filename:
        fn = glob.glob("*.pts")[0]
    else:
        fn = args.filename
    convert(fn)
