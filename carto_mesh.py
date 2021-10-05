from mesh_tools import *
from tqdm import tqdm
import time
import glob
import matplotlib.pyplot as plt
import tetgen
import os
plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors


class CartoMesh:
    def __init__(self, name: str = "", directory: str = "./"):
        self.point_info = pd.DataFrame()
        self.triangle_info = pd.DataFrame()
        self.points = np.array([])
        self.triangles = np.array([])
        self.edges = np.array([])
        self.cells = []
        self.layers = None
        self.mesh = None
        self.dir = directory
        self.name = name
        self.thickness = 0
        self.ratio = None
        self.verbose = False  # TODO: implement in other functions
        self.switches = "-pYkAmNEFq2.5/20a2e+6"

    def __cartoToCsv(self, verbose: bool = True):
        """Reads in a carto .mesh file and generates .csv files per header in real-time in the same directory"""
        with io.open(self.name + '.mesh', 'r', encoding='utf-8') as f:
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
                        of = open(self.dir + section + ".csv", "w+", encoding='utf-8')

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

                    elif len(line) > 1 and line[0] != ";":  # actual data, not empty line or column names
                        ind, data = line.split("=")  # line has shape "ind = x  y  z  nx  ny  nz  groupID
                        of.write(str(ind))
                        for el in data.split():  # ignore ind
                            of.write("," + str(el))
                        of.write("\n")

            of.close()

    def __readVertTri(self):
        """Reads VerticesSection.csv and TrianglesSection.csv files and initializes class"""
        assert os.path.exists("{}/VerticesSection.csv".format(self.dir)) \
               and os.path.exists('{}/TrianglesSection.csv'.format(self.dir)), \
            "VerticesSection.csv or TrianglesSection.csv not found! Run cartoToCsv() first."
        vert = pd.read_csv("{}/VerticesSection.csv".format(self.dir), sep=',')
        tri = pd.read_csv('{}/TrianglesSection.csv'.format(self.dir), sep=',')
        return vert, tri

    def __getNeighboringFaces(self, index: int):
        """Finds all mesh facets containing some point.
        Args:
            index: the index referring to the mesh point.
        Returns:
            DataFrame: Pandas DataFrame containing the triangles that have the point at index. The DataFrame
            contains the following columns: 'Vertex0', 'Vertex1', 'Vertex2', 'NormalX', 'NormalY', 'NormalZ', 'GroupID'

        """
        facets_ = self.triangle_info.loc[
            (self.triangle_info[['Vertex0']].values == index) | (self.triangle_info[['Vertex1']].values == index) |
            (self.triangle_info[['Vertex2']].values == index)]
        return facets_

    def __update(self, mesh_):
        self.mesh = mesh_
        self.points = mesh_.points
        self.edges = mesh_.extract_all_edges()
        if type(mesh_) == pv.PolyData:
            self.triangles = pvToPmCells(mesh_.faces)
            self.cells = None
        elif type(mesh_) == pv.UnstructuredGrid:
            self.cells = mesh_.cells
            self.triangles = None

    def initialise(self, name: str = ""):
        """Reads in a .mesh file and writes out a .csv file per .mesh header. Reads in the VerticesSection and
        TrianglesSection and uses these to construct a mesh.
        Args:
            name: name of .mesh file to be initialised, including '.mesh' extension
        Returns:
            Nothing
        Todo:
            Add color to initial mesh
        """
        self.name = name if name else glob.glob(self.dir + '*.mesh')[0][:-5]  # first .mesh file without extension
        self.__cartoToCsv(verbose=False)
        vert, tri = self.__readVertTri()
        self.point_info = vert
        self.triangle_info = tri
        self.points = np.array([[1000.*e for e in (x, y, z)]
                                for x, y, z in self.point_info[['X', 'Y', 'Z']].values])
        self.triangles = self.triangle_info[['Vertex0', 'Vertex1', 'Vertex2']].values
        self.layers = 1
        self.mesh = pv.PolyData(self.points, pmToPvFaces(self.triangles))
        self.edges = self.mesh.extract_all_edges()

    def writeEndoEpi(self, thickness: float, ratio: float):
        """Reads in .csv files from carto2csv.py, calculates normals of mesh facets and adds two layers of points
        along these normals: the epi- and endocardium. Writes out these two layers as .txt files.
        Args:
            thickness: the distance between the inner (endo) and outer (epi) layer
            ratio: the ratio of distances from both layers to the original middle layer. A ratio of .8 will add an outer
            layer .8*thickness from the middle and an inner layer .2*thickness to the inside. Ratios > .5 are
            recommended to avoid self intersections.
        Returns:
            Nothing. Writes out two files: endo.txt and epi.txt
        """
        self.thickness = thickness
        self.ratio = ratio
        outside = ratio * thickness
        inside = (1 - ratio) * thickness
        endo = open(self.dir + 'endo.txt', 'w+')
        epi = open(self.dir + 'epi.txt', 'w+')
        for f in (endo, epi):
            f.write(self.dir + "Index,X,Y,Z\n")
        tqbar = tqdm(range(self.points.shape[0]), desc='        Adding points')
        for index, point in self.point_info.iterrows():
            tqbar.update(1)
            co = point[['X', 'Y', 'Z']]
            facets_ = self.__getNeighboringFaces(index)
            av_normal = calcAvNormal(facets_)
            newpoint = [e + outside * n for e, n in zip(co, av_normal)]
            newpoint2 = [e - inside * n for e, n in zip(co, av_normal)]
            for f, p in zip((endo, epi), (newpoint2, newpoint)):  # loop over files and corresponding points
                f.write(str(index))
                for comp in p:  # components of point
                    f.write(',' + str(comp * 1000.))  # Convert from mm to µm and write out
                f.write('\n')

        endo.close()
        epi.close()
        tqbar.close()

    def splitLayer(self, thickness: float = .5, ratio: float = .8):
        """Calls upon writeEndoEpi() to write 'endo.txt' and 'epi.txt': two files containing the point coordinates of
        the bounding mesh layers. These are read in again and used to construct two surface meshes. These two
        surface meshes are added together as one mesh. The class properties are updated to these points and faces.
        Args:
            thickness: the distance between the inner (endo) and outer (epi) layer
            ratio: the ratio of distances from both layers to the original middle layer. A ratio of .8 will add an outer
            layer .8*thickness from the middle and an inner layer .2*thickness to the inside. Ratios > .5 are
            recommended to avoid self intersections.
        """
        if not os.path.exists(self.dir + 'endo.txt') or not os.path.exists(self.dir + 'epi.txt'):
            self.writeEndoEpi(thickness, ratio)
        assert self.layers < 2, "Amount of layers is already >= 2. Splitting aborted."
        endo = pd.read_csv(self.dir + 'endo.txt', sep=',', index_col=0)
        epi = pd.read_csv(self.dir + 'epi.txt', sep=',', index_col=0)
        f = pmToPvFaces(self.triangles)
        fullmesh = pv.PolyData(endo.values, f) + pv.PolyData(epi.values, f)
        self.points = endo.append(epi).values
        self.triangles = np.concatenate((self.triangles, np.array([e + len(self.points) for e in self.triangles])))
        self.layers = 2
        self.mesh = fullmesh

    def cleanMesh(self, tol: float, max_iter: int = 10, print_si: bool = True):
        """Calls upon PyVista's built-in clean method to clean the mesh. Afterwards, calls upon PyMesh's
        built-in methods remove_degenerated_triangles() and detect_self_intersection() to remove
        duplicated faces and triangles, and resolve self-intersections.
        Args:
            tol: tolerance passed to PyVista's built-in clean() method (absolute tolerance)
            max_iter: max amount of iterations PyMesh is allowed to try and fix self-intersections,
            duplicate faces and duplicate vertices.
            print_si: print the progress of fixing self-intersections
        Returns:
            PyVista PolyData: The cleaned mesh"""
        mesh_ = self.mesh.clean(lines_to_points=True, polys_to_lines=True,
                                tolerance=tol, absolute=True)
        mesh_ = makePyMesh(mesh_)
        mesh_, info = pm.remove_degenerated_triangles(mesh_)
        si = len(pm.detect_self_intersection(mesh_))
        if print_si:
            print('\t' + str(si) + ' self-intersections detected')

        i = 0
        while i < max_iter and si != 0:  # there are still self-intersections, max 10 iterations
            mesh_, info = pm.remove_duplicated_faces(mesh_)
            mesh_, info = pm.remove_duplicated_vertices(mesh_, tol=tol)
            mesh_ = pm.resolve_self_intersection(mesh_)
            si = len(pm.detect_self_intersection(mesh_))
            if print_si:
                sys.stdout.write("\r" + '\tIteration {}: '.format(i + 1) + str(si) + ' self-intersection(s) left ')
                sys.stdout.flush()
            i += 1
        print("")
        mesh_, info = pm.remove_duplicated_vertices(mesh_, tol=tol)
        mesh_, info = pm.remove_duplicated_faces(mesh_)
        mesh_ = makePyVista(mesh_).clean(lines_to_points=True, polys_to_lines=True,
                                         tolerance=tol, absolute=True)
        self.__update(mesh_)
        self.point_info = pd.DataFrame()  # point info does no longer match
        self.triangle_info = pd.DataFrame()
        return mesh_

    def getEdgeLengths(self, mesh=None):
        """Gets all edge lengths from a PyMesh mesh (used in homogenizeMesh())"""
        if mesh:
            edges = mesh.extract_all_edges()
        else:
            edges = self.edges
        pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
        distances = []
        for pair in pmedges:
            co1, co2 = edges.points[pair]
            distances.append(dist(co1, co2))
        return np.array(distances)

    def homogenizeMesh(self, nsteps=10, boxplot=False, plot_mesh=False, verbose=False, return_dist=False,
                       min_edge=500., max_edge=1000.):
        """ Iteratively splits long edges and collapses short edges until they all have a length between
        min_edge and max_edge"""

        # Keep in mind that max allowed distance in split_long_edges should be about twice
        # the tolerance in order to get converging behaviour, since short edges become
        # at least twice as big as tolerance. Setting max_edge < 2 * min_edge wastes some
        # computational time, but can be useful in case the input mesh has self-intersections.

        # TODO: clean this

        def calcAlpha(n_, nsteps_, max_edge_, longest_edge_):
            """Calculate the longest allowed edge length (alpha) for a given iteration step.
            Alpha drops exponentially as the iterations progress. At step 0, alpha equals the
            longest edge of the input mesh. At the final step, alpha equals the maximum desired edge length."""
            return (1 - n_ / nsteps_) * longest_edge_ * np.exp(-3. * n_ / nsteps_) + max_edge_

        def calcTol(n_, nsteps_, min_edge_):
            """Calculate the shortest allowed edge for a given iteration step"""
            return 500. * (1 - n_ / nsteps_) + min_edge_  # min edge length drops linearly

        edge_lengths = self.getEdgeLengths()
        longest_edge = np.max(edge_lengths)
        mesh_ = makePyMesh(self.mesh)

        if boxplot:
            dist_range = [edge_lengths]  # initialize list of distances during refinement procedure
        for n in tqdm(range(1, nsteps + 1), desc='        Resizing edges'):  # don't use \t in tqdm descriptions
            # alpha can never be smaller than 2*tol. Keep some margin in-between
            # It will still run otherwise, but tetgen will probably detect self-intersections
            alpha = calcAlpha(n, nsteps, max_edge, longest_edge)
            tol = calcTol(n, nsteps, min_edge)
            splitmesh_, info = pm.split_long_edges(mesh_, max_edge_length=alpha)
            colmesh_, info = pm.collapse_short_edges(splitmesh_, tol, preserve_feature=True)
            mesh_ = colmesh_
            if boxplot:
                dist_range.append(self.getEdgeLengths(mesh_))

        if boxplot:
            def plotEdgelengths(dist_range, show=True, save=True):
                """Plots a boxplot of the edgelengths at each iteration step
                dist_range should be a 2D array, each array entry containing all the edge lengths at
                an iteration step"""
                plt.style.use('fivethirtyeight')
                medians = [np.median(d) for d in dist_range]
                means = [np.mean(d) for d in dist_range]
                sigmas = [np.sqrt(v) for v in [np.var(d) for d in dist_range]]

                fig, ax = plt.subplots(figsize=(8, 6))
                ax.set_ylabel("Edge length (µm)", size=22)
                ax.set_xlabel("Iteration step", size=22)
                plt.tight_layout(pad=1.8, h_pad=1.8)
                ax.fill((.5, nsteps + 1.5, nsteps + 1.5, .5), (max_edge, max_edge, min_edge, min_edge), color=colors[0],
                        alpha=.3, label='Desired edgelength\n({} - {} µm)'.format(min_edge, max_edge))
                # plt.axhline(max_edge, color='black', lw=2.5)
                # plt.axhline(min_edge, color='black', lw=2.5)
                polygonx = [.88] + list(range(2, nsteps + 2)) + list(reversed(range(2, nsteps + 2))) + [.88]
                polygony = [m + s for m, s in zip(means, sigmas)] + \
                           list(reversed([m - s for m, s in zip(means, sigmas)]))
                # ax.fill(polygonx, polygony, color=colors[0], alpha=.2, label="mean $\pm$ $\sigma$")
                ax.boxplot(dist_range, labels=range(len(dist_range)), whis=[0, 100],
                           boxprops=dict(color="gray", linewidth=2.5),
                           medianprops=dict(color="gray", linewidth=2.5),
                           whiskerprops=dict(color="gray", linewidth=2.5),
                           capprops=dict(color="gray", linewidth=2.5), zorder=1)
                ax.plot(range(1, nsteps + 2), [calcAlpha(step) for step in range(nsteps + 1)], color=colors[1],
                        )  # xaxis needs to be shifted by one for matplotlib bullshit reasons, hence range(1, nsteps+2)
                ax.plot(range(1, nsteps + 2), [calcTol(step) for step in range(nsteps + 1)], color=colors[1],
                        label=r'$\alpha$ and tolerance')
                ax.set_title("Mesh edge lengths during refinement", size=28)
                ax.legend(loc="upper right", prop={'size': 20})
                ax.tick_params(axis='x', labelsize=18)
                ax.tick_params(axis='y', labelsize=18)
                if show:
                    plt.show()
                if save:
                    fig.savefig("../Plots/EdgeLengthRefinement_{}-{}.png".format(min_edge, max_edge), dpi=300)

            plotEdgelengths(dist_range)
        print("\tCleaning mesh...")
        mesh_ = cleanMesh(makePyVista(mesh_), tol=min_edge / 2, iter=6)
        self.__update(mesh_)

        if return_dist:
            return mesh_, getEdgeLengths(makePyMesh(mesh_))
        self.name += '_{}-{}µm'.format(int(min_edge), int(max_edge))  # update name
        return mesh_  # return final version of mesh

    def tetrahedralise(self, switches: str, n_col=3):
        # TODO: what do i want to write out, what do i want to update in the class
        def runTetGen(name, switches: str):
            """Runs TetGen in terminal and captures output.
            Args:
                switches: Switches to be used together with TetGen bash command.
            Returns:
                """
            tetcommand = "tetgen {} {}.smesh".format(switches, name)
            command = ['stdbuf', '-o0', *tetcommand.split(" ")]
            print('\tNow running bash command {}\n'.format(" ".join(command)))
            print("\t-------------------- TetGen output --------------------")
            # adding stdbuf -o0 makes TetGen output real-time
            p = Popen(command, stdout=PIPE, stderr=PIPE, encoding='utf-8')
            colin_pair = []
            for line in p.stdout:
                l_ = str(line.rstrip())
                print("\t" + l_)
                # The following checks assume that, if an error is found, it's colinearity
                # This still needs to be expanded for overlapping lines and segments intersecting triangles
                if '1st' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                    # line is of form "1st: [24578,24581]."
                    # split at ":", drop first two (whitespace and "[") and last two ("]" and ".") chars
                    # split at comma to make list again
                    edge1 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                if '2nd' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                    # line is of form "2nd: [24578,24581]."
                    edge2 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                    colin_pair.append([edge1, edge2])
                p.stdout.flush()
            print("\t------------------ End TetGen output ------------------")
            return colin_pair

        colinear_pairs = runTetGen(name=self.name, switches=switches)
        i = 0
        mesh = self.mesh
        while i < n_col and len(colinear_pairs) > 0:  # check for collinear pairs
            print("\n\tCollinearity found! Adapting points...")
            mesh = makeNonCollinear(mesh, colinear_pairs)
            writeToSmesh(mesh, self.name)
            colinear_pairs = runTetGen(self.name, switches)
            i += 1

    def reconstruct(self, boxplot=False,
                    switches="-pYkAmNEFq2.5/20a2e+6", refine_steps=10,
                    keep_intmed=False, skip_reading=False,
                    min_edge=700., max_edge=1400., return_surfacemesh=True, ncv=1, n_col=1):
        """Reads in .mesh file and writes out a refined tetrahedron mesh in .vtk format and carp format.
        If 'speed.csv' exists in the cwd, also interpolates these speeds on the mesh. speed.csv should be a csv file with
        columns 'x', 'y', 'z' and 'speed', where the xyz coordinates refer to point coordinates. Can be calculated with
        DGM.


        @param refine_steps: <int> amount of refinement steps during the refinement process of the double-layered surface mesh
        @param boxplot: <>bool Make a boxplot of the edgelengths during the refinement procedure of the double-layered
        surface mesh
        @param switches: <str> switches to be used by TetGen
        @param keep_intmed: <bool> keep intermediate files generated by TetGen, carto2csv.py or applycv.py
        @param skip_reading: <bool> skip the conversion of CARTO .mesh files to .csv files
        (useful together with keep_intmed to save time)
        @param min_edge: <float> minimum allowed edge length for the mesh in µm
        @param max_edge: <float> maximum allowed edge length for the mesh in µm
        @param return_surfacemesh: <bool> stop the tetrahedralisation and return the refined double-layered surface mesh
        instead (for testing purposes)
        @param ncv: <int> amount of different conduction velocity distributions to be made. 1 corresponds with the
        original conduction velocities
        @param from_stl: <bool> Read in .stl file instead of .smesh file. Useful for preprocessing of meshes with e.g. MeshLab
        (not fully implemented yet)
        @param n_col: <int> amount of times TetGen needs to try and fix collinearities before giving up.
        """

        print("\n####### Creating 3D tetrahedralized {}\n".format(self.name))

        step = 1  # progress
        start = time.time()

        # Read in the meshes
        if not skip_reading:
            print("---- {}. Reading CARTO .mesh file and writing to csv\n".format(step))
            self.__cartoToCsv()
            step += 1

        # Add second layer of points and triangles
        print("\n---- {}. Adding endo and epi point layers\n".format(step))
        self.splitLayer()
        step += 1

        # Refine, make .smesh file
        self.homogenizeMesh(nsteps=refine_steps, boxplot=boxplot, min_edge=min_edge, max_edge=max_edge)
        writeToSmesh(self.mesh, self.name)
        # create_surface adapted filename to form: basefilename + "_<minedge>-<maxedge>µm"
        step += 1

        # Make tetrahedrons from double-layered refined surface mesh
        self.tetrahedralise(switches=switches)
        step += 1
        print('\n\tTetrahedralizing done.')

        print("\n---- {}. Converting to carp and paraview-friendly format\n".format(step))
        convertMesh_Meshtool(self.name)  # by default vtk to carp_txt
        # pts2paraview.convert(fn + '.pts')
        step += 1

        # clean created mesh
        print("\n---- {}. Cleaning with meshtool\n".format(step))
        cleanMesh_Meshtool(self.name, .2)
        print('\n\tAdapting .vtk ...')
        convertMesh_Meshtool(self.name, ifmt='carp_txt',
                             ofmt='vtk')  # Make sure .vtk file is same as carp file after cleaning
        step += 1

        # TODO
        # if ncv:
        #     print("\n---- {}. Applying conduction velocities\n".format(step))
        #     apply_cv.run(meshdir, fn + '.1.vtk', write_VTK_file=True, n_variations=ncv)
        #     step += 1

        if not keep_intmed:  # delete intermediate files
            print("\n----- {}. Deleting intermediate files\n".format(step))
            # remove intermediate files that weren't present before generating and aren't the .vtk file
            to_delete = ["mid.txt", "endo.txt", "epi.txt", "TrianglesSection.csv", "VerticesSection.csv",
                         self.name + ".smesh", self.name + '.mtr', self.name + '.1.mtr', self.name + '.1.p2t',
                         self.name + '.1.node', self.name + '.1.edge',
                         self.name + '.1.face', self.name + '.1.ele', 'VerticesAttributesSection.csv',
                         'VerticesColorsSection.csv', 'myo.csv', 'noncond.csv', 'Surface.stl']
            for trash in to_delete:
                if glob.glob(trash):
                    os.remove(trash)
                    print("\tDeleted ", trash)
            step += 1

        duration = time.time() - start
        print("\n\tMesh reconstructed in {0}m {1:.2f}s".format(int(duration // 60), duration % 60))
        self.__update(pv.read(self.name+'.1.vtk'))

    def plot(self):
        """Plots the mesh in its current state using PyVista's Plotter() method."""
        p = pv.Plotter()
        p.add_mesh(self.mesh)
        p.show()


if __name__ == '__main__':
    m = CartoMesh()
    m.initialise()
    m.reconstruct()
