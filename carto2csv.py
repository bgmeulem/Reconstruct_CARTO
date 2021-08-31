#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io


def cartoToCsv(inmesh, meshdir="", verbose=True):
    """Reads in a carto .mesh file and generates .csv files per header in real-time"""
    with io.open(inmesh, 'r', encoding='utf-8') as f:
        # variables used throughout loop:
        section = "Header"  # name of current section
        current_line = 0  # index of current line being read
        section_line = 0  # index of current section name
        section_ind = 0

        # Section 0 = Header, file format info
        # Section 1 = General Attributes, info about section data
        # Section 2 = VerticesSection, x y z coordinates plus data (normals and GroupID)
        # Section 3 = TrianglesSection,
        # section 4 = VerticesColorSection
        # Section 5 = VerticesAttributesSection

        for line in f:
            current_line += 1  # update line index

            if line[0] == "[":  # new section encountered
                if verbose:
                    print("\t" + str(section_ind) + '. ' + section)
                section = line[1:-2]  # section name
                section_line = current_line
                section_ind += 1

                # Opening or closing csv files
                if section_ind > 2:
                    of.close()  # past VerticesSection, an outfile has already been created (infra)
                if section_ind > 1:  # skip GeneralAttributes
                    of = open(meshdir + section + ".csv", "w+", encoding='utf-8')

                # New tqdm loops per section
                # if section_ind > 1:
                #     if section_ind > 2:
                #         time.sleep(1e-6)
                #         loop.close()  # close previous loop
                #     if section_ind == 3:
                #         time.sleep(1e-6)
                #         loop = tqdm(desc='        ' + section, total=n_triangles, position=0, leave=True)
                #     else:
                #         time.sleep(1e-6)
                #         loop = tqdm(desc='        ' + section, total=n_vertices, position=0, leave=True)

            # extract n_vertices and n_triangles to use for tqdm loops
            # elif section_ind == 1:
            #     if current_line == section_line + 3:  # n_verticess
            #         n_vertices = int(line.split("=")[1])
            #     elif current_line == section_line + 4:
            #         n_triangles = int(line.split('=')[1])

            elif section_ind > 1:  # useful data
                if section_ind == 4:
                    header_line = section_line + 2  # VerticesColorSection has extra line
                else:
                    header_line = section_line + 1

                if current_line == header_line:  # column names
                    column_names = line.split()[1:]  # first element is useless ";"
                    of.write("Index")
                    for name in column_names:
                        of.write("," + name)
                    of.write("\n")
                    # time.sleep(1e-8)

                elif len(line) > 1 and line[0] != ";":  # actual data, not empty line or column names
                    # time.sleep(1e-8)
                    # loop.update(1)
                    ind, data = line.split("=")  # line has shape "ind = x  y  z  nx  ny  nz  groupID
                    of.write(str(ind))
                    for el in data.split():  # ignore ind
                        of.write("," + str(el))
                    of.write("\n")
        # time.sleep(1e-8)
        # loop.close()
        # time.sleep(1e-8)
        of.close()


if __name__ == '__main__':
    infile = '1-AT1230CL.mesh'
    cartoToCsv(infile)

