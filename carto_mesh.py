from typing import Tuple
from mesh_tools import *
from tqdm import tqdm
import time
import glob
import matplotlib.pyplot as plt
from collections import OrderedDict
import os

plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors


class CartoMesh:
    """A class containing functions for initialising from file, plotting, reconstructing and
    applying conduction velocities."""

    def __init__(self, name: str = ""):
        self.name = ""
        self.dir = "./"
        self.point_info = pd.DataFrame()
        self.triangle_info = pd.DataFrame()
        self.points = np.array([])
        self.n_points = 0
        self.triangles = np.array([], dtype=int)
        self.n_triangles = 0
        self.edges = np.array([])
        self.cells = []
        self.n_cells = 0
        self.layers = None
        self.mesh = None
        self.myo = self.non_myo = None
        self.thickness = 0
        self.ratio = None
        self.reconstruction_parameters = {}
        self.cv_interpolation_parameters = {}
        self.__initialiseFromFile(name)
        self.__parseSettings()

    def __initialiseFromFile(self, name: str = ""):
        def parseName(fn):
            comp = fn.split(os.sep) if fn else ['.']
            if '.mesh' not in comp[-1] and '.vtk' not in comp[-1]:  # only directory given: look for file
                fns = glob.glob(os.path.join(*comp, '*.mesh'))  # take first .mesh file
                if len(fns):
                    comp = fns[0].split(os.sep)  # overwrite filename and dir components
                else:
                    raise FileNotFoundError("No .mesh or .vtk file found in directory \'{}\'"
                                            "".format(os.path.join(*comp)))
            d = os.path.join(*comp[:-1]) + os.sep
            n = comp[-1].split('.')[0]
            return n, d

        def initialiseFromMeshFile(self_, name_: str = "") -> None:
            """Initialize a carto mesh from .mesh file.
            Reads in a .mesh file and writes out a .csv file per .mesh header. Reads in the VerticesSection and
            TrianglesSection and uses these to construct a mesh.
            Args:
                self_: self, the CartoMesh() object, passed through from the overarching function.
                name_: name of .mesh file to be initialised, including '.mesh' extension
            Returns:
                None: Nothing. Updates mesh.
            """
            self_.name, self_.dir = parseName(name_)
            self_.__cartoToCsv(verbose=False)
            vert, tri = self_.__readVertTri()
            self_.point_info = vert
            self_.triangle_info = tri
            points = np.array([[1000. * e for e in (x, y, z)]
                               for x, y, z in self_.point_info[['X', 'Y', 'Z']].values])
            triangles = self_.triangle_info[['Vertex0', 'Vertex1', 'Vertex2']].values
            self_.layers = 1
            self_.__update(pv.PolyData(points, pmToPvFaces(triangles)))
            color = [e for e in self_.triangle_info['GroupID'].to_numpy()]
            self_.mesh['color'] = color
            self_.mesh = self_.mesh.ctp()  # triangle data to point data
            self_.myo, self_.non_myo = self_.extractMyoAndNonMyo()

        def initialiseFromVtkFile(self_, name_: str = "") -> None:
            """Initialize a carto mesh from .vtk file.
            Args:
                self_: self_: self, the CartoMesh() object, passed through from the overarching function.
                name_: name of .mesh file to be initialised, including '.vtk' extension
            Returns:
                None: Nothing. Updates mesh.
            """

            self_.name, self_.dir = parseName(name_)
            self_.layers = 1
            self_.__update(pv.read(self_.dir + self_.name + '.vtk'))

        if '.vtk' in name:
            initialiseFromVtkFile(self, name)
        else:
            initialiseFromMeshFile(self, name)


    def __cartoToCsv(self, verbose: bool = True):
        """Reads in a carto .mesh file and generates .csv files per header in real-time in the same directory"""
        with io.open(os.path.join(self.dir, self.name) + '.mesh', 'r', encoding='utf-8') as f:
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
                        of = open(os.path.join(self.dir, section + ".csv"), "w+", encoding='utf-8')

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
        self.n_points = len(self.points)
        self.edges = mesh_.extract_all_edges()
        if type(mesh_) == pv.PolyData:
            self.triangles = pvToPmCells(mesh_.faces)
            self.n_triangles = len(self.triangles)
            self.cells = None
            self.n_cells = 0
        elif type(mesh_) == pv.UnstructuredGrid:
            self.cells = pvToPmCells(mesh_.cells)
            self.n_cells = len(self.cells)
            self.triangles = None
            self.n_triangles = 0

    def writeEndoEpi(self, thickness: float, ratio: float) -> None:
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

    def splitLayer(self, thickness: float = .5, ratio: float = .8) -> None:
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
        endo, epi = pv.PolyData(endo.values, f), pv.PolyData(epi.values, f)
        endo['color'] = epi['color'] = self.triangle_info['GroupID']
        fullmesh = endo + epi
        fullmesh = fullmesh.ctp()
        self.points = fullmesh.points
        self.triangles = np.concatenate((self.triangles, np.array([e + len(self.points) for e in self.triangles])))
        self.layers = 2
        self.mesh = fullmesh

    def cleanMesh(self, tol: float, max_iter: int = 10, print_si: bool = True) -> pv.PolyData:
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

    def getEdgeLengths(self, mesh=None) -> np.ndarray:
        """Gets all edge lengths from a PyMesh mesh (used in homogenizeMesh())"""
        edges = mesh.extract_all_edges() if mesh else self.edges
        pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
        distances = []
        for pair in pmedges:
            co1, co2 = edges.points[pair]
            distances.append(dist(co1, co2))
        return np.array(distances)

    def homogenizeMesh(self, nsteps=10, boxplot=False, return_dist=False,
                       edge_range=(500., 1000.)) -> pv.PolyData:
        """ Iteratively splits long edges and collapses short edges until they all have a length between
        min_edge and max_edge"""

        # Keep in mind that max allowed distance in split_long_edges should be about twice
        # the tolerance in order to get converging behaviour, since short edges become
        # at least twice as big as tolerance. Setting max_edge < 2 * min_edge wastes some
        # computational time, but can be useful in case the input mesh has self-intersections.

        min_edge, max_edge = edge_range

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
            colmesh_.remove_attribute('face_sources')  # pymesh adds this attribute. Makes conversion difficult.
            mesh_ = colmesh_
            if boxplot:
                dist_range.append(self.getEdgeLengths(mesh_))

        if boxplot:
            def plotEdgelengths(dist_range, show=True, save=True):
                """Plots a boxplot of the edge lengths at each iteration step
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
        self.name += '_{}-{}µm'.format(int(min_edge), int(max_edge))  # update name

        if return_dist:
            return mesh_, getEdgeLengths(makePyMesh(mesh_))
        return mesh_  # return final version of mesh

    def tetrahedralise(self, switches: str, n_col=3) -> None:
        def runTetGen(name, switches: str) -> [[int, int]]:
            """Runs TetGen in terminal and captures output.
            Args:
                switches: Switches to be used together with TetGen bash command.
            Returns:
                """
            cwd = os.getcwd()
            os.chdir(self.dir)
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
                if '1st' in l_ and 'facet' not in l_ and 'seg' not in l_:
                    # line is of form "1st: [24578,24581]."
                    # split at ":", drop first two (whitespace and "[") and last two ("]" and ".") chars
                    # split at comma to make list again
                    edge1 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                if '2nd' in l_ and 'facet' not in l_ and 'seg' not in l_:
                    # line is of form "2nd: [24578,24581]."
                    edge2 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                    colin_pair.append([edge1, edge2])
                p.stdout.flush()
            print("\t------------------ End TetGen output ------------------")
            os.chdir(cwd)
            return colin_pair

        colinear_pairs = runTetGen(name=self.name, switches=switches)
        i = 0
        mesh = self.mesh
        while i < n_col and len(colinear_pairs) > 0:  # check for collinear pairs
            print("\n\tCollinearity found! Adapting points...")
            mesh = makeNonCollinear(mesh, colinear_pairs)
            writeToSmesh(mesh, self.dir + self.name)
            colinear_pairs = runTetGen(self.name, switches)
            i += 1

    def extractMyoAndNonMyo(self) -> Tuple[pv.PolyData, pv.PolyData]:
        """Checks scalar data of current mesh (clinical tags). Returns two point clouds: one of myocardium and one
        that's everything but myocardium"""
        myo = pv.PolyData()
        noncond = pv.PolyData()
        for tag in np.unique(self.mesh['color']):
            tagged_pointcloud = self.mesh.points[[self.mesh['color'] == tag]]
            if tag == 0:
                myo += pv.PolyData(tagged_pointcloud)
            else:
                noncond += pv.PolyData(tagged_pointcloud)
        return myo, noncond

    def writeTags(self) -> None:
        """Writes out all points in a mesh with a certain tag to a .csv file"""
        for tag in self.mesh['color'].unique():
            tagged_pointcloud = self.mesh.points[[self.mesh['color'] == tag]]
            with open(self.dir + "Tag_{}.csv".format(tag), 'w+') as of:
                csvWriter = csv.writer(of, delimiter=',')
                csvWriter.writerows(tagged_pointcloud)

    def getNonMyoIndices(self) -> [int]:
        if self.non_myo.npoints:
            self.non_myo["speed"] = self.non_myo.n_points * [0.]
            self.myo["speed"] = self.myo.n_points * [1.]
            c = self.non_myo + self.myo  # a mask that is 1 where the carto point is myocardium
            # for each meshpoint, find closest carto point
            tree = nb.KDTree(c.points)
            distances, indices = tree.query(self.mesh.points, k=1)  # these are c indices for each meshpoint
            nc_mesh_ind = [ind for ind in range(len(indices)) if c["speed"][indices[ind]] == 0.]
        else:
            nc_mesh_ind = []
        return nc_mesh_ind

    def writeAdjustNa2(self, tol: float = 1e-5, g_Na: float = 7.8, write_dat: bool = True) -> None:
        """Writes adjustment file to close off Na2+ channels in cells where CV ~ 0
        Args:
            tol: tolerance for when a point is considered to have a conduction velocity of 0. Points whose conduction velocity < tol are considered to be zero.
            g_Na: value for (maximum) Sodium channel conductance. Open Na-channels are set to this conductance. Default: 7.8 pS
            write_dat: write a .dat file for visual inspection in e.g. meshalyzer.
        Returns:
            None: Nothing. Writes out a .txt file and, if wanted, a .dat file"""

        def writeDat(d: str, mesh: pv.PolyData, name: str = "gNA2") -> None:
            datfile = open(d + name + ".dat", "w+")
            dat = np.zeros(mesh.n_points)
            dat[ptn == 0.] = 1
            for e in dat:
                datfile.write(str(e) + '\n')
            datfile.close()

        assert self.mesh.cells or self.mesh.triangles, "Mesh does not contain any cells or triangles."
        cells = pvToPmCells(self.mesh.cells) if self.cells else pvToPmCells(self.mesh.triangles)
        speed = self.mesh["speed"]

        ptn = np.ones(len(self.mesh.points)) * g_Na
        for i in range(len(speed)):
            if speed[i] < tol:  # if CV is ~ 0 -> close off Na channel
                vertices = cells[i]
                ptn[vertices] = 0.

        stimfile = open(self.dir + "gNA2.adj", "w+")
        stimfile.write(str(len(ptn)) + "\n" + "extra\n")
        for i in range(len(ptn)):
            stimfile.write(str(i) + " " + str(ptn[i]) + "\n")
        stimfile.close()

        if write_dat:
            writeDat(self.dir, self.mesh, "gNA2")

    def setNonCondByIndexFiles(self, region_dir='Regions', write_dat=False, index_col="meshID") -> None:
        """Opens directory region_dir and reads in .csv files there. These .csv files should contain the indices of
        points whose conduction velocity should be set to zero.
        Writes out a .dat file for each .csv file if wanted
        Args:
            region_dir: Name of directory containing the regions to be set to a conduction velocity of 0. These regions should be .csv files containing the indices of points to be set to CV=0
            index_col: Name of the column in the .csv files that contain the point indices. Default=\"meshID\" for easy workflow in correspondence with mesh_tools.ptsToParaview()"""

        def writeDat(data: pd.DataFrame, n_points, name):
            """Writes a .dat file for a given pd.DataFrame"""
            # write scar .dat file
            datfile = open(name + ".dat", "w+")
            dat = np.zeros(n_points)
            dat[data["meshID"]] = 1
            for e in dat:
                datfile.write(str(e) + '\n')
            datfile.close()

        def getNonCondRegions(meshdir, mesh, region_dir_=region_dir, write_dat_=write_dat):
            """"""
            non_cond_region_indices = []
            if os.path.isdir(meshdir + region_dir_):  # apply scars if directory exists
                for csvfile in glob.glob(meshdir + region_dir_ + "/*.csv"):
                    scar = pd.read_csv(csvfile)
                    if write_dat_:
                        writeDat(scar, csvfile.split('.')[0], len(mesh.points))
                    for index in scar[index_col].values:
                        non_cond_region_indices.append(index)
            return np.array(non_cond_region_indices)

        # Set manually selected scars to 0 velocity
        scars = getNonCondRegions(self.dir, self.mesh, 'Regions')
        self.mesh["speed"] = [0. if p in scars else self.mesh["speed"][p] for p in range(self.n_points)]

    def applyCV(self, cv_params=None, speed_file='speed.csv',
                write_csv=False, write_VTK_file=False, write_txt=True, write_dat=False, write_xyz=False,
                outdir='scale_factors/'):

        if cv_params is None:
            cv_params = {'radius': 4000,  # for interpolation
                         'sharpness': 1.5,  # for interpolation
                         'n_variations': 1,
                         'n_neighbors': 15,  # for randomization
                         'speed_limit': (0, 1.4)}

        def applySpeedLimit(data_: pd.DataFrame, limit: tuple, speed_col: str = 'speed') -> pd.DataFrame:
            speeds = data_[speed_col]
            for i in range(len(speeds)):
                s = speeds[i]
                if s < limit[0]:
                    speeds[i] = limit[0]
                if s > limit[1]:
                    speeds[i] = limit[1]
            data_[speed_col] = speeds
            return data_

        def randomizeCV(data_: pd.DataFrame, n_neighbors: int, speed_col: str = "speed") -> pd.DataFrame:
            """Given a mesh with conduction velocities, makes a variation on this conduction velocity distribution.
            A normal distribution is fitted to each point and their neighbors. From this distribution, a new
            conduction velocity is randomly sampled.
            Args:
                data_: the mesh in pd.DataFrame format. Must contain columns 'x', 'y' and 'z' with the coordinates in µm.
                n_neighbors: amount of neighbors to use for each normal distribution fitting.
                speed_col: name of the column containing the conduction velocity values for each mesh point.
            Returns:
                pd.DataFrame: DataFrame containing the original point coordinates and updated conduction velocities"""
            print("\n\t#### Variation ", n)
            # tweak conduction velocities of input CV file
            points = data_[["x", "y", "z"]].values
            tree = nb.KDTree(points)
            speeds = np.zeros(len(data_))
            for i in tqdm(range(len(points)), desc='        Calculating new velocities'):
                p = points[i]
                dist, neighbors = tree.query([p], k=n_neighbors)
                neighborCVs = input_data.loc[[int(e) for e in neighbors[0]]][speed_col]
                mean, sigma = np.mean(neighborCVs), np.std(neighborCVs)
                new_cv = np.random.normal(mean, sigma, 1)
                speeds[i] = np.abs(new_cv)
            data_[speed_col] = speeds
            return data_

        # read in data
        input_data = pd.read_csv(self.dir + speed_file,
                                 usecols=["speed", "x", "y", "z"])
        input_data = applySpeedLimit(input_data, cv_params['speed_limit'])
        med_speed = np.mean(input_data["speed"])
        # create new dataframe for mutable purposes
        calculated_data = input_data

        for n in range(cv_params['n_variations']):
            if n != 0:  # variations of conduction velocities
                calculated_data = randomizeCV(calculated_data, n_neighbors=cv_params['n_neighbors'])

            # Create PolyData to use in interpolation
            data = applySpeedLimit(calculated_data, cv_params['speed_limit'])
            pvdata = pv.PolyData(np.array([data["x"], data["y"], data["z"]]).T)
            pvdata["speed"] = data["speed"]

            # Interpolate on mesh
            print("\tInterpolating on mesh")
            mesh = self.mesh.interpolate(pvdata, radius=cv_params['radius'], sharpness=cv_params['sharpness'],
                                         strategy="null_value", null_value=med_speed,
                                         pass_point_arrays=False, pass_cell_arrays=False)

            pointdata = mesh["speed"]
            mesh = mesh.ptc()  # point data to cell data
            cell_data = mesh["speed"]
            sq_point_data = pd.DataFrame([e ** 2 for e in pointdata], columns=["squared speed"])
            sq_cell_data = pd.DataFrame([e ** 2 for e in cell_data], columns=["squared speed"])
            self.__update(mesh)

            if not os.path.exists(self.dir + outdir):
                os.mkdir(self.dir + outdir)
            # write to csv file
            if write_csv:
                print("\tWriting squared speed to csv")
                sq_point_data.to_csv(self.dir + outdir + "sq_CV_point_{}.csv".format(n), index_label="PointID")
                sq_cell_data.to_csv(self.dir + outdir + "sq_CV_cell_{}.csv".format(n), index_label="CellID")

            if write_xyz:
                print("\tWriting point cloud")
                of = open(self.dir + outdir + "input_points_{}.txt".format(n), "w+")
                of.write("X,Y,Z,speed\n")
                for i in tqdm(range(len(pvdata.points))):
                    p = pvdata.points[i]
                    d = pvdata["speed"][i]
                    for c in p:
                        of.write(str(c) + ",")
                    of.write(str(d) + '\n')
                of.close()

            # write text file to read in during simulation
            if write_txt:
                print("\tWriting txt file: {}scale_factor_{}.txt".format(outdir, n))
                of = open(self.dir + outdir + "scale_factor_{}.txt".format(n), "w+")
                for e in sq_cell_data["squared speed"]:
                    of.write(str(e) + '\n')
                of.close()

            # .dat file for visualisation purposes
            if write_dat:
                print("\tWriting dat file: {}scale_factor_{}.dat".format(outdir, n))
                of = open(self.dir + outdir + "scale_factor_{}.dat".format(n), "w+")
                for e in sq_point_data["squared speed"]:
                    of.write(str(e) + '\n')
                of.close()

            # write to vtk for inspection in paraview
            if write_VTK_file:
                print("\tWriting mesh to {}{}_CV{}.vtk".format(self.dir, self.name, n))
                pv.save_meshio(self.dir + "{}_CV{}.vtk".format(self.name, n), self.mesh)

    def reconstruct(self, boxplot=False,
                    switches="-pYkAmNEFq2.5/20a2e+6", refine_steps=10,
                    keep_intmed=False,
                    edge_range=(600., 1000.), ncv=1, n_col=1) -> None:
        """Reads in .mesh file and writes out a refined tetrahedron mesh in .vtk format and carp format.
        If 'speed.csv' exists in the cwd, also interpolates these speeds on the mesh. speed.csv should be a csv file with
        columns 'x', 'y', 'z' and 'speed', where the xyz coordinates refer to point coordinates. Can be calculated with
        DGM.

        Args:
          ncv: amount of conduction velocity distributions to make. The first one corresponds to the original distribution.
          edge_range: The desired values for the shortest and longest mesh edge length in µm
          refine_steps: Amount of times the mesh edges will be split and collapsed in homogenizeMesh()
          switches: The switches used in the TetGen command to tetrahedralize a mesh.
          boxplot: Keep track of the edge length distribution during homogenizeMesh() and plot a boxplot afterwards.
          keep_intmed: Keep files that were generated during the reconstruction process. These include .csv files from
          self.__cartoToCsv() and files generated by TetGen in tetrahedralise().
          n_col: Amount of times collinearities will be fixed during the tetrahedralisation process.

        Returns:
          Nothing.
          Writes out a .vtk file of the tetrahedron mesh. Writes out the same mesh in Carp text format.
          Updates the mesh to the tetrahedron mesh.

        """

        def getPipeline(self_, boxplot_, switches_, refine_steps_, edge_range_, ncv_, n_col_):
            pipeline_ = OrderedDict([
                (
                    "Adding endo and epi point layers",
                    {'function': [self_.splitLayer],
                     'args': [{}]}
                ),
                (
                    "Refining surface mesh",
                    {'function': [self_.homogenizeMesh],
                     'args': [{'nsteps': refine_steps_, 'boxplot': boxplot_, 'edge_range': edge_range_}]}

                ),
                (
                    "Writing to .smesh",
                    {'function': [writeToSmesh],
                     'args': [{'mesh': self_.mesh, 'name': self_.dir + self_.name}]}
                ),
                (
                    "Tetrahedralising with TetGen",
                    {'function': [self_.tetrahedralise],
                     'args': [{'switches': switches_, 'n_col': n_col_}]}
                ),
                (
                    "Converting to carp and paraview-friendly format",
                    {'function': [convertMesh_Meshtool, ptsToParaview],
                     'args': [{'meshname': self_.dir + self_.name},
                              {'filename': self_.dir + self_.name + '.pts', 'column_name': 'meshID'}]}
                ),
                (
                    "Cleaning with meshtool",
                    {'function': [cleanMesh_Meshtool, convertMesh_Meshtool],
                     'args': [{'meshname': self_.dir + self_.name, 'threshold': .2},
                              {'meshname': self_.dir + self_.name, 'ifmt': 'carp_txt', 'ofmt': 'vtk'}]}
                )
            ])
            if ncv_:
                cv_params = {'radius': 4000,  # for interpolation
                             'sharpness': 1.5,  # for interpolation
                             'n_variations': ncv_,
                             'n_neighbors': 15,  # for randomization
                             'speed_limit': (0, 1.4)}
                pipeline_.update({"Applying conduction velocities":
                                     {'function': [self.applyCV],
                                      'args': [{'cv_params': cv_params, 'write_VTK_file': True}]
                                      }
                                  })
            return pipeline_

        print("\n####### Creating 3D tetrahedralized {}\n".format(self.name))
        start = time.time()
        step = 1
        pipeline = getPipeline(self, boxplot, switches, refine_steps, edge_range, ncv, n_col)
        for action_name in pipeline:
            print("\n---- {}. {}\n".format(step, action_name))
            action = pipeline[action_name]
            for func, args in zip(action['function'], action['args']):
                func(**args)
            # Reinitialise pipeline dict to update function arguments
            pipeline = getPipeline(self, boxplot, switches, refine_steps, edge_range, ncv, n_col)
            step += 1

        if not keep_intmed:  # delete intermediate files
            print("\n----- {}. Deleting intermediate files\n".format(step))
            # remove intermediate files that weren't present before generating and aren't the .vtk file
            to_delete = ["mid.txt", "endo.txt", "epi.txt",
                         "TrianglesSection.csv", "VerticesSection.csv", 'VerticesAttributesSection.csv',
                         'VerticesColorsSection.csv',
                         self.name + ".smesh", self.name + '.1.mtr', self.name + '.1.mtr', self.name + '.1.p2t',
                         self.name + '.1.node', self.name + '.1.edge', self.name + '.1.face', self.name + '.1.ele',
                         'myo.csv', 'noncond.csv']
            for trash in to_delete:
                if os.path.exists(self.dir+trash):
                    os.remove(self.dir+trash)
                    print("\tDeleted ", trash)
            step += 1

        duration = time.time() - start
        print("\n\tMesh reconstructed in {0}m {1:.2f}s".format(int(duration // 60), duration % 60))
        self.__update(pv.read(self.dir+self.name + '.1.vtk'))

    def plot(self):
        """Plots the mesh in its current state using PyVista's Plotter() method."""
        p = pv.Plotter()
        p.add_mesh(self.mesh)
        p.show()


if __name__ == '__main__':
    m = CartoMesh('BlankMeshes/OC59')
    # m.reconstruct()
