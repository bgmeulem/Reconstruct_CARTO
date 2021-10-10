from typing import Tuple
from mesh_tools import *
from tqdm import tqdm
import time
import glob
import matplotlib.pyplot as plt
import os
plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors



class CartoMesh:
    """A class CartoMesh containing functions for initialising from file, plotting, reconstructing."""

    def __init__(self):
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
        self.dir = ""
        self.name = ""
        self.thickness = 0
        self.ratio = None
        self.verbose = False  # TODO: implement in other functions
        self.switches = "-pYkAmNEFq2.5/20a2e+6"

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
            self.n_cells = mesh_.n_cells
        elif type(mesh_) == pv.UnstructuredGrid:
            self.cells = pvToPmCells(mesh_.cells)
            self.n_cells = len(self.cells)
            self.triangles = None
            self.n_triangles = mesh_.n_faces

    def initialise(self, name: str = ""):
        """Initialize a carto mesh from .mesh file.
        Reads in a .mesh file and writes out a .csv file per .mesh header. Reads in the VerticesSection and
        TrianglesSection and uses these to construct a mesh.
        Args:
            name: name of .mesh file to be initialised, including '.mesh' extension
        Returns:
            Nothing
        Todo:
            Add color to initial mesh
        """

        def parseName(fn):
            comp = fn.split(os.sep) if fn else ['.']
            if '.mesh' not in comp[-1]:  # only directory given: look for .mesh file
                fns = glob.glob(os.path.join(*comp, '*.mesh'))
                if len(fns):
                    comp = fns[0].split(os.sep)  # overwrite filename and dir components
                else:
                    raise FileNotFoundError("No .mesh file found in directory \'{}\'".format(os.path.join(*comp)))
            d = os.path.join(*comp[:-1]) + os.sep
            n = comp[-1].split('.')[0]
            return n, d

        self.name, self.dir = parseName(name)
        self.__cartoToCsv(verbose=False)
        vert, tri = self.__readVertTri()
        self.point_info = vert
        self.triangle_info = tri
        self.points = np.array([[1000. * e for e in (x, y, z)]
                                for x, y, z in self.point_info[['X', 'Y', 'Z']].values])
        self.n_points = len(self.points)
        self.triangles = self.triangle_info[['Vertex0', 'Vertex1', 'Vertex2']].values
        self.layers = 1
        self.mesh = pv.PolyData(self.points, pmToPvFaces(self.triangles))
        color = [e for e in self.triangle_info['GroupID'].to_numpy()]
        self.mesh['color'] = color
        self.mesh = self.mesh.ctp()  # triangle data to point data
        self.myo, self.non_myo = self.extractMyoAndNonMyo()
        self.edges = self.mesh.extract_all_edges()

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

    def getEdgeLengths(self, mesh=None) -> np.ndarray[float]:
        """Gets all edge lengths from a PyMesh mesh (used in homogenizeMesh())"""
        edges = mesh.extract_all_edges() if mesh else self.edges
        pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
        distances = []
        for pair in pmedges:
            co1, co2 = edges.points[pair]
            distances.append(dist(co1, co2))
        return np.array(distances)

    def homogenizeMesh(self, nsteps=10, boxplot=False, plot_mesh=False, verbose=False, return_dist=False,
                       min_edge=500., max_edge=1000.) -> pv.PolyData:
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

        if return_dist:
            return mesh_, getEdgeLengths(makePyMesh(mesh_))
        self.name += '_{}-{}µm'.format(int(min_edge), int(max_edge))  # update name
        return mesh_  # return final version of mesh

    def tetrahedralise(self, switches: str, n_col=3) -> None:
        def runTetGen(name, switches: str) -> [[int, int]]:
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

    # This one works, but have to manually initialise mesh["speed'] first en set active scalar to "speed"
    def applyManualNonCondRegions(self, region_dir='Regions', write_dat=False, index_col="MeshID") -> None:
        """Opens directory 'Regions' and reads in the .csv files there.
        Sets the conduction velocity to zero for all points with an index that's present in these .csv files
        Writes out a .dat file for each .csv file if wanted"""
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

    # TODO: not tested. How to get a speed.csv? Extract from old vtk files?
    def applyCV(self, speed_limit=(0, 1.4), speed_file='speed.csv',
                radius=4000, sharpness=1.5,
                write_csv=False, write_VTK_file=False, write_txt=True, write_dat=False, write_xyz=False, write_adj=False,
                n_variations=10, n_neighbors=5,
                outdir='scale_factors/'):

        def convertMmtoMicron(data_: pd.DataFrame, co_cols=('x', 'y', 'z')) -> pd.DataFrame:
            for c in co_cols:
                data_["{}_um".format(c)] = [1000. * e for e in data_[c]]
            return data_

        def applySpeedLimit(data_: pd.DataFrame, limit: tuple, speed_col: str = 'speed') -> pd.DataFrame:
            speeds = data_[speed_col]
            for i in range(len(speeds)):
                s = speeds[i]
                if s < limit[0]:
                    speeds[i] = limit[0]
                if s > speed_limit[1]:
                    speeds[i] = limit[1]
            data_[speed_col] = speeds
            return data_

        def writeAdjust(meshdir, mesh) -> None:
            """Writes adjustment file to close off Na2+ channels in cells where CV ~ 0"""
            # mesh needs to be pyvista PolyData
            # meshname = glob.glob(meshdir+"*.1.vtk")[0].split('.')[0]
            cells = pvToPmCells(mesh.cells)  # mesh cells in PyMesh format
            speed = mesh["speed"]

            ptn = np.ones(len(mesh.points)) * 7.8  # default value for g_Na
            for i in range(len(speed)):
                if speed[i] < 1e-5:  # if CV is ~ 0 -> close off Na channel
                    vertices = cells[i]
                    ptn[vertices] = 0.

            stimfile = open(meshdir + "gNA2.adj", "w+")
            stimfile.write(str(len(ptn)) + "\n" + "extra\n")
            for i in range(len(ptn)):
                stimfile.write(str(i) + " " + str(ptn[i]) + "\n")
            stimfile.close()

            datfile = open(meshdir + "gNA2.dat", "w+")
            dat = np.zeros(len(mesh.points))
            dat[ptn == 0.] = 1
            for e in dat:
                datfile.write(str(e) + '\n')
            datfile.close()

        def randomizeCV(data: pd.DataFrame, speed_col: str = "speed") -> pd.DataFrame:
            print("\n\t#### Variation ", n)
            # tweak conduction velocities of input CV file
            points = data[["x_um", "y_um", "z_um"]].values
            tree = nb.KDTree(points)
            speeds = np.zeros(len(data))
            for i in tqdm(range(len(points)), desc='        Calculating new velocities'):
                p = points[i]
                dist, neighbors = tree.query([p], k=n_neighbors)
                neighborCVs = input_data.loc[[int(e) for e in neighbors[0]]][speed_col]
                mean, sigma = np.mean(neighborCVs), np.std(neighborCVs)
                new_cv = np.random.normal(mean, sigma, 1)
                speeds[i] = np.abs(new_cv)
            data[speed_col] = speeds
            return data

        # read in data
        input_data = pd.read_csv(self.dir + speed_file,
                                 usecols=["speed", "x", "y", "z"])  # TODO: needs to be present
        input_data = convertMmtoMicron(input_data)  # TODO: depends on input file
        input_data = applySpeedLimit(input_data, speed_limit)
        med_speed = np.mean(input_data["speed"])
        # create new dataframe for mutable purposes
        calculated_data = input_data

        ids = getGroupIds()  # Carto tags (excluding 0 and -10000): correspond to MV, LPV and RPV. At max 3
        print("\tDetected tags (non-conductive): ", ids)
        for n in range(n_variations):
            if n != 0:  # variations of conduction velocities
                calculated_data = randomizeCV(calculated_data)

            # Create PolyData to use in interpolation
            data = applySpeedLimit(calculated_data, speed_limit)
            pvdata = pv.PolyData(np.array([data["x_um"], data["y_um"], data["z_um"]]).T)
            pvdata["speed"] = data["speed"]

            # Interpolate on mesh
            print("\tInterpolating on mesh")
            mesh = self.mesh.interpolate(pvdata, radius=radius, sharpness=sharpness,
                                         strategy="null_value", null_value=med_speed,
                                         pass_point_arrays=False, pass_cell_arrays=False)

            # Set auto-detected non-conductive regions after interpolation
            nc_mesh_ind = self.getNonMyoIndices()
            mesh["speed"] = [0. if p in nc_mesh_ind else mesh["speed"][p] for p in range(mesh.n_points)]
            pointdata = mesh["speed"]
            mesh = mesh.ptc()  # point data to cell data
            cell_data = mesh["speed"]
            sq_point_data = pd.DataFrame([e ** 2 for e in pointdata], columns=["squared speed"])
            sq_cell_data = pd.DataFrame([e ** 2 for e in cell_data], columns=["squared speed"])
            self.__update(mesh)

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

            if write_adj:
                print("\tWriting file: gNA2.adj")
                writeAdjust(self.dir, mesh)

            # write to vtk for inspection in paraview
            if write_VTK_file:
                print("\tWriting mesh to {}_CV{}.vtk".format(self.dir.split('.')[0], n))
                polyd = pv.UnstructuredGrid(mesh.cells, np.array(len(pvToPmCells(mesh.cells)) * [10]), mesh.points)
                polyd["speed"] = cell_data
                pv.save_meshio(self.dir + "{}_CV{}.vtk".format(self.dir.split('.')[0], n), mesh)

            if len(ids) < 3:
                # don't calculate variations if you have to recalculate again with manual scars
                print("\n\t!!! Only {} out of 3 tag(s) found - Manual input needed !!!\n".format(len(ids)))
                print("\tVariations of CVs will not be calulated")
                break

    def reconstruct(self, boxplot=False,
                    switches="-pYkAmNEFq2.5/20a2e+6", refine_steps=10,
                    keep_intmed=False,
                    min_edge=700., max_edge=1400., ncv=1, n_col=1) -> None:
        """Reads in .mesh file and writes out a refined tetrahedron mesh in .vtk format and carp format.
        If 'speed.csv' exists in the cwd, also interpolates these speeds on the mesh. speed.csv should be a csv file with
        columns 'x', 'y', 'z' and 'speed', where the xyz coordinates refer to point coordinates. Can be calculated with
        DGM.

        Args:
          min_edge: The desired value for the shortest mesh edge length in µm
          max_edge: The desired value for the longest mesh edge length in µm
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

        print("\n####### Creating 3D tetrahedralized {}\n".format(self.name))

        step = 1  # progress
        start = time.time()

        # Read in the meshes
        if not all(
                [os.path.exists(os.path.join(self.dir, e)) for e in ('VerticesSection.csv', 'TrianglesSection.csv')]):
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
        self.tetrahedralise(switches=switches, n_col=n_col)
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

        # TODO:
        # if ncv:
        #     print("\n---- {}. Applying conduction velocities\n".format(step))
        #     applyCV(self.dir, self.name + '.1.vtk', write_VTK_file=True, n_variations=ncv)
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
        self.__update(pv.read(self.name + '.1.vtk'))

    def plot(self):
        """Plots the mesh in its current state using PyVista's Plotter() method."""
        p = pv.Plotter()
        p.add_mesh(self.mesh)
        p.show()


if __name__ == '__main__':
    m = CartoMesh()
    m.initialise('BlankMeshes/OC59')
    m.reconstruct()
