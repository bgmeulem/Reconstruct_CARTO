"""
This file contains the class CartoMesh with all its reconstruction-specific methods.
This class contains methods for initialising from a .mesh file, plotting, reconstructing and applying conduction velocities.
This file imports :mod:`mesh_tools` for all general mesh functions that are not specific to carto meshes.
:mod:`mesh_tools` also concerns itself with all python module imports.
"""
from mesh_tools import *


class CartoMesh:
    """
    A class containing functions for initialising from file, plotting, reconstructing and
    applying conduction velocities.
    """

    def __init__(self, name: str = ""):
        self.name = ""
        self.dir = "./"

        self.points = np.array([])
        self.n_points = 0
        self.point_info = pd.DataFrame()

        self.triangles = np.array([], dtype=int)
        self.n_triangles = 0
        self.triangle_info = pd.DataFrame()

        self.edges = np.array([])
        self.cells = []
        self.n_cells = 0

        self.layers = None
        self.mesh = None

        self.myo = self.non_myo = None

        self.thickness = 0
        self.ratio = None

        self.settings = {}
        self.__initialiseFromFile(name)
        self.__readSettings()

    def __initialiseFromFile(self, name: str = "") -> None:
        """
        Initialise the class from either a .mesh file or a .vtk file.

        Args:
            name: The filename (.mesh or .vtk), or the directory containing the .mesh file.

        Returns:
            None: Nothing
        """

        def parseName(fn) -> Tuple[str, str]:
            """
            Given a filename to initialise, this method parses the filename and extracts the filetype and
             parent directory.
            Args:
                fn: Filename, containing extension or not, or name of the parent directory

            Returns:
                Tuple: tuple containing the name of the file and the name of its parent directory
            """
            if '*' in fn:  # starred expression: file name not given
                fn_ = glob.glob(fn)[0]  # search for the file
            else:
                fn_ = fn
            comp = fn_.split(os.sep) if fn_ else ['.']
            if '.mesh' not in comp[-1] and '.vtk' not in comp[-1]:  # only directory given: look for file
                fns = glob.glob(os.path.join(*comp, '*.mesh'))  # take first .mesh file
                if len(fns):
                    comp = fns[0].split(os.sep)  # overwrite filename and dir components
                else:
                    raise FileNotFoundError("No .mesh or .vtk file found in directory \'{}\'"
                                            "".format(os.path.join(*comp)))
            d = os.path.join(*comp[:-1]) + os.sep
            n = '.'.join(comp[-1].split('.')[:-1])  # drop file extension
            return n, d

        def initialiseFromMeshFile(self_, name_: str = "") -> None:
            """
            Initialize a carto mesh from .mesh file.
            Reads in a .mesh file and writes out a .csv file per .mesh header. Reads in the VerticesSection and
            TrianglesSection and uses these to construct a mesh.

            Args:
                self_: self, the :class:`CartoMesh` object, passed through from the overarching function.

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
            """Initialize a mesh from .vtk file.

            Args:
                self_: self_: self, the :class:`CartoMesh` object, passed through from the overarching function.

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

    def __cartoToCsv(self, verbose: bool = True) -> None:
        """
        Reads in a carto .mesh file and generates .csv files per header in real-time in the same directory

        Args:
            verbose: Use verbose output. Default=True

        Returns:
            None: Nothing. Writes out .csv files.
        """
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

    def __readVertTri(self) -> None:
        """
        Reads VerticesSection.csv and TrianglesSection.csv files and initializes class

        Returns:
            None: Nothing. Updates mesh.
        """
        assert os.path.exists("{}/VerticesSection.csv".format(self.dir)) \
               and os.path.exists('{}/TrianglesSection.csv'.format(self.dir)), \
            "VerticesSection.csv or TrianglesSection.csv not found! Run cartoToCsv() first."
        vert = pd.read_csv("{}/VerticesSection.csv".format(self.dir), sep=',')
        tri = pd.read_csv('{}/TrianglesSection.csv'.format(self.dir), sep=',')
        return vert, tri

    def __getNeighboringFaces(self, index: int) -> pd.DataFrame:
        """
        Finds all mesh facets containing some point.

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

    def __update(self, mesh_: Union[pv.PolyData, pv.UnstructuredGrid]) -> None:
        """
        Updates the CartoMesh to match the given mesh.

        Args:
            mesh_: A PyVista mesh of the type PolyData or UnstructuredGrid .

        Returns:
            None: Nothing. Updates the mesh.
        """
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

    def __readSettings(self) -> None:
        """
        Reads the settings.ini file and updates the CartoMesh accordingly.

        Returns:
            None: Nothing. Updates the mesh settings attribute.
        """
        assert os.path.exists("settings.ini"), "No settings.ini file found"
        config = configparser.ConfigParser(inline_comment_prefixes='#')
        config.read('settings.ini')
        for section in config.sections():
            self.settings[section.lower()] = {key: eval(val) for key, val in config[section].items()}

    def writeEndoEpi(self, thickness: float, ratio: float) -> None:
        """
        Reads in .csv files from carto2csv.py, calculates normals of mesh facets and adds two layers of points
        along these normals: the epi- and endocardium. Writes out these two layers as .txt files.

        Args:
            thickness: the distance between the inner (endo) and outer (epi) layer in cm.

            ratio: the ratio of distances from both layers to the original middle layer. A ratio of .8 will add an outer
            layer .8*thickness from the middle and an inner layer .2*thickness to the inside. Ratios > .5 are
            recommended to avoid self intersections.

        Returns:
            None: Nothing. Writes out two files: endo.txt and epi.txt. These are two point-cloud meshes containing the
            mesh point coordinates of the epicardium and endocardium.
            Coordinates are written out in mm. Not in cm like the input file.
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
        """
        Calls upon :func:`writeEndoEpi` to write 'endo.txt' and 'epi.txt': two files containing the point coordinates of
        the bounding mesh layers. These are read in again and used to construct two surface meshes. These two
        surface meshes are added together as one mesh. The class properties are updated to these points and faces.

        Args:
            thickness: the distance between the inner (endo) and outer (epi) layer

            ratio: the ratio of distances from both layers to the original middle layer. A ratio of .8 will add an outer
            layer .8*thickness from the middle and an inner layer .2*thickness to the inside. Ratios > .5 are
            recommended to avoid self intersections.
        Returns:
            Nothing
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
        """
        Calls upon PyVista's built-in clean method to clean the mesh. Afterwards, calls upon PyMesh's
        built-in methods remove_degenerated_triangles() and detect_self_intersection() to remove
        duplicated faces and triangles, and resolve self-intersections.

        Args:
            tol: tolerance passed to PyVista's built-in clean() method (absolute tolerance)

            max_iter: max amount of iterations PyMesh is allowed to try and fix self-intersections,
            duplicate faces and duplicate vertices.

            print_si: print the progress of fixing self-intersections

        Returns:
            PyVista PolyData: The cleaned mesh
        """
        mesh_ = self.mesh.clean(lines_to_points=True, polys_to_lines=True, tolerance=tol, absolute=True)
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

    def getEdgeLengths(self, mesh: Union[pm.Mesh, bool] = None) -> np.ndarray:
        """
        Gets all edge lengths. If a mesh is given, the mesh must be a PyMesh Mesh object. If no mesh is given (default),
        the CartoMesh's mesh attribute is used (PyVista PolyData or UnstructuredGrid)

        Args:
            mesh: A PyMesh Mesh object. Default: None

        Returns:
            np.ndarray: An array containing the edge lengths of te mesh.
        """
        edges = mesh.extract_all_edges() if mesh else self.edges
        pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
        distances = []
        for pair in pmedges:
            co1, co2 = edges.points[pair]
            distances.append(dist(co1, co2))
        return np.array(distances)

    def homogenizeMesh(self, nsteps=10, boxplot=False, return_dist=False,
                       edge_range=(600., 1000.)) -> pv.PolyData:
        """
        Iteratively splits long edges and collapses short edges until they all have a length between
        min_edge and max_edge

        Args:
            nsteps: Amount of steps to take to iteratively adapt the mesh edge lengths

            boxplot: Make a boxplot of the mesh edge lengths for each iteration step after the refinement procedure.

            return_dist: Return a 2xnsteps array of the mesh edge lengths

            edge_range: Minimum and maximum allowed edge lengths

        Returns:
            Union(PolyData, Tuple[PolyData, List]): Either the refined surface mesh in PyVista's PolyData format, or both the refined mesh and a list containing all the edge lengths for each iteration step.
        """

        def calcAlpha(n_: int, nsteps_: int, max_edge_: float, longest_edge_: float) -> float:
            """
            Calculate the longest allowed edge length (alpha) for a given iteration step.
            Alpha drops exponentially as the iterations progress. At step 0, alpha equals the
            longest edge of the input mesh. At the final step, alpha equals the maximum desired edge length.

            Args:
                n_: Current iteration step

                nsteps_: Total amount of iteration steps

                max_edge_: Maximum allowed edge length at final iteration step

                longest_edge_: The longest edge length of the input mesh before any refinement

            Returns:
                float: The longest allowed edge length for the current iteration step
            """
            return (1 - n_ / nsteps_) * longest_edge_ * np.exp(-3. * n_ / nsteps_) + max_edge_

        def calcTol(n_: int, nsteps_: int, min_edge_: float) -> float:
            """
            Calculate the shortest allowed edge for a given iteration step.
            The shortest allowed edge augments linearly as the iteration steps progress.

            Args:
                n_: Current iteration step

                nsteps_: Total amount of iteration steps

                min_edge_: Minimum allowed edge length at the final iteration step.

            Returns:
                float: The minimum allowed edge length for the current iteration step.
            """
            return 500. * (1 - n_ / nsteps_) + min_edge_  # min edge length drops linearly

        def plotEdgelengths(dist_range: np.ndarray, show: bool = True, save: bool = True) -> None:
            """
            Plots a boxplot of the edge lengths at each iteration step
            dist_range should be a 2D array, each array entry containing all the edge lengths at
            an iteration step

            Args:
                dist_range: An NxM array containing the mesh edge lengths for each iteration step, where M is the amount
                of mesh edges at iteration step N.

                show: Show the plot. Default = True

                save: Save the plot. Default = True

            Returns:
                None: Nothing
            """
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
            plt.axhline(max_edge, color='black', lw=2.5)
            plt.axhline(min_edge, color='black', lw=2.5)
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

        min_edge, max_edge = edge_range
        edge_lengths = self.getEdgeLengths()
        longest_edge = np.max(edge_lengths)
        mesh_ = makePyMesh(self.mesh)

        if boxplot:
            dist_range = [edge_lengths]  # initialize list of distances during refinement procedure
        for n in tqdm(range(1, nsteps + 1), desc='        Resizing edges'):
            alpha = calcAlpha(n, nsteps, max_edge, longest_edge)
            tol = calcTol(n, nsteps, min_edge)
            splitmesh_, info = pm.split_long_edges(mesh_, max_edge_length=alpha)
            colmesh_, info = pm.collapse_short_edges(splitmesh_, tol, preserve_feature=True)
            colmesh_.remove_attribute('face_sources')  # pymesh adds this attribute. Makes conversion difficult.
            mesh_ = colmesh_
            if boxplot:
                dist_range.append(self.getEdgeLengths(mesh_))

        if boxplot:
            plotEdgelengths(dist_range)

        print("\tCleaning mesh...")
        mesh_ = cleanMesh(makePyVista(mesh_), tol=min_edge / 3., iterations=6)
        self.__update(mesh_)
        self.name += '_{}-{}µm'.format(int(min_edge), int(max_edge))  # update name

        if return_dist:
            return mesh_, getEdgeLengths(makePyMesh(mesh_))
        return mesh_  # return final version of mesh

    def tetrahedralise(self, switches, n_col_retry=3) -> None:
        """
        Runs TetGen as a bash command.

        Args:
            switches: Switches to use with the TetGen command. See https://www.wias-berlin.de/software/tetgen/switches.html

            n_col_retry: Amount of times to attempt to fix collinearities during the TetGen process.

        Returns:
            None: Nothing. Writes out a .vtk file.
        """
        def runTetGen(name: str, switches_: str) -> [[int, int]]:
            """
            Runs TetGen in terminal and captures output.

            Args:
                name: name of the base file

                switches_: Switches to be used together with TetGen bash command.

            Returns:
                [[int, int]]: A list containing pairs of indices. These pairs define a mesh segment that's colinear.
            """
            cwd = os.getcwd()
            os.chdir(self.dir)
            tetcommand = "tetgen {} {}.smesh".format(switches_, name)
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

        colinear_pairs = runTetGen(self.name, switches)
        i = 0
        mesh = self.mesh
        while i < n_col_retry and len(colinear_pairs) > 0:  # check for collinear pairs
            print("\n\tCollinearity found! Adapting points...")
            mesh = makeNonCollinear(mesh, colinear_pairs)
            writeToSmesh(mesh, self.dir + self.name)
            colinear_pairs = runTetGen(self.name, switches)
            i += 1

    def extractMyoAndNonMyo(self) -> Tuple[pv.PolyData, pv.PolyData]:
        """
        Checks scalar data of the mesh in its current state (clinical tags) and extracts the part of the mesh with
        tags that correspond to myocardium (tag == 0).

        Returns:
            Tuple[pv.PolyData, pv.PolyData]: A tuple containing two PyVista PolyData meshes: a pointcloud of points
            corresponding to the mesh myocardium and a pointcloud of the rest of the mesh.
        """
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
        """
        Writes out all points in a mesh with a certain tag to a .csv file

        Returns:
            None: Nothing. Writes out .csv files.
        """
        for tag in self.mesh['color'].unique():
            tagged_pointcloud = self.mesh.points[[self.mesh['color'] == tag]]
            with open(self.dir + "Tag_{}.csv".format(tag), 'w+') as of:
                csvWriter = csv.writer(of, delimiter=',')
                csvWriter.writerows(tagged_pointcloud)

    def writeAdjustNa2(self, tol: float = 1e-5, g_Na: float = 7.8, write_dat: bool = True) -> None:
        """
        Writes adjustment file to close off Na2+ channels in cells where CV ~ 0

        Args:
            tol: tolerance for when a point is considered to have a conduction velocity of 0. Points whose
            conduction velocity < tol are considered to be zero.

            g_Na: value for (maximum) Sodium channel conductance. Open Na-channels are set to this conductance.
            Default: 7.8 pS

            write_dat: write a .dat file for visual inspection in e.g. meshalyzer.

        Returns:
            None: Nothing. Writes out a .txt file and, if wanted, a .dat file
        """

        def writeDat(d: str, mesh: pv.PolyData, name: str = "gNA2") -> None:
            """
            Write a .dat file. This file is a column of 0s and 1s. The index column corresponds to the mesh point index.
            If the value equals to one on row i, then mesh.points[i] is a point whose Na2+ channel has been closed off.

            Args:
                d: Directory to write .dat file to

                mesh: A Pyvista PolyData mesh

                name: Name of the .dat file, without file extension. Default = 'gNA2'

            Returns:
                None: Nothing. Writes out .dat file <name>.dat
            """
            datfile = open(d + name + ".dat", "w+")
            dat = np.zeros(mesh.n_points)
            dat[ptn == 0.] = 1
            for e in dat:
                datfile.write(str(e) + '\n')
            datfile.close()

        assert self.mesh.n_cells or self.mesh.n_triangles, "Mesh does not contain any cells or triangles."
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
        """
        Opens directory region_dir and reads in .csv files there. These .csv files should contain the indices of
        points whose conduction velocity should be set to zero.
        Writes out a .dat file for each .csv file if wanted

        Args:
            region_dir: Name of directory containing the regions to be set to a conduction velocity of 0.
            These regions should be .csv files containing the indices of points to be set to CV=0

            index_col: Name of the column in the .csv files that contain the point indices.
            Default=\"meshID\" for easy workflow in correspondence with :func:`mesh_tools.ptsToParaview`

        Returns:
            None: Nothing
        """

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

    def applyCV(self, speed_file='speed.csv',
                write_csv=False, write_VTK_file=False, write_txt=True, write_dat=False, write_xyz=False,
                outdir='scale_factors/', region_dir='', ncv=None, speed_col="speed"):
        """
        Applies conduction velocities on a reconstructed tetrahedralised carto mesh.

        Args:
            speed_file: Name of the file containing the coordinates and conduction velocity values to interpolate.

            write_csv: Write out the coordinates and speeds to a .csv file

            write_VTK_file: Write out the resulting mesh to .vtk

            write_txt: Write out the conduction velocities as scale factors (readable by OpenCARP) to <outdir>

            write_dat: Write out the regions of the mesh with CV=0 to .dat file, to visualize in Meshalyzer
            write_xyz: Write point cloud of the conduction velocites, as interpolated on the mesh
            outdir: Name of the directory to contain the scale factors aka conduction velocities
            region_dir: Name of the directory containing .csv files with indices of mesh points whose conduction velocity will be set to 0.

        Returns:
            None: Nothing. Writes out one or more of the following files in <outdir>: .csv, .vtk, .txt, .dat or a
            pointcloud .csv file if <write_xyz> == True.
        """

        def applySpeedLimit(data_: pd.DataFrame, limit: tuple, speed_col: str = 'speed') -> pd.DataFrame:
            """
            Given a DataFrame, cuts off all conduction velocities whose conduction velocity in the column <speed_col>
            have a value outside the range of <limit>
            Args:
                data_: The DataFrame containing the conduction velocities
                limit: The minimum and maximum allowed conduction velocities as a tuple with length 2
                speed_col: Name of the DataFrame column containing the conduction velocities. Default = 'speed'

            Returns:
                DataFrame: DataFrame whose conduction velocities in the column <speed_col> have been cut off to match <speed_col>
            """
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
                DataFrame: DataFrame containing the original point coordinates and updated conduction velocities"""
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

        if not ncv:
            n_cv = self.settings['reconstruction_parameters']['n_cv']
        else:
            n_cv = ncv

        # read in data
        input_data = pd.read_csv(self.dir + speed_file,
                                 usecols=["speed", "x", "y", "z"])
        input_data = applySpeedLimit(input_data, self.settings['cv_interpolation_parameters']['speed_limit'])
        med_speed = np.mean(input_data["speed"])
        # create new dataframe for mutable purposes
        calculated_data = input_data

        for n in range(n_cv):
            if n != 0:  # variations of conduction velocities
                calculated_data = randomizeCV(calculated_data,
                                              n_neighbors=self.settings['cv_interpolation_parameters']['n_neighbors'])

            # Create PolyData to use in interpolation
            data = applySpeedLimit(calculated_data, self.settings['cv_interpolation_parameters']['speed_limit'])
            pvdata = pv.PolyData(np.array([data["x"], data["y"], data["z"]]).T)
            pvdata["speed"] = data["speed"]

            # Interpolate on mesh
            print("\tInterpolating on mesh")
            mesh = self.mesh.interpolate(pvdata, radius=self.settings['cv_interpolation_parameters']['radius'],
                                         sharpness=self.settings['cv_interpolation_parameters']['sharpness'],
                                         strategy="null_value", null_value=med_speed,
                                         pass_cell_arrays=False)
            self.__update(mesh)  # set mesh
            if region_dir:
                print("\tSetting regions to 0 conduction velocity")
                self.setNonCondByIndexFiles(region_dir, write_dat)

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

    def reconstruct(self) -> None:
        """
        Reads in .mesh file and writes out a refined tetrahedron mesh in .vtk format and carp format.
        If 'speed.csv' exists in the cwd, also interpolates these speeds on the mesh. speed.csv should be a csv file with
        columns 'x', 'y', 'z' and 'speed', where the xyz coordinates refer to point coordinates. Can be calculated with
        DGM.

        Args:

        Returns:
          Nothing.
          Writes out a .vtk file of the tetrahedron mesh. Writes out the same mesh in Carp text format.
          Updates the mesh to the tetrahedron mesh.
        """

        def getPipeline(self_, boxplot, switches, refine_steps, min_edge, max_edge, n_cv, n_col_retry, keep_intmed):
            """Given a set of parameters, generate the functions with corresponding parameters to reconstruct
            a carto mesh from zero to simulation-ready.
            Args:
                self_: self
                boxplot: bool: Whether or not to make a boxplot of the edgelengths during refinement in :func:`homogenizeMesh`
                switches: str: The switches used in the tetgen command
                refine_steps: int: amount of refinement steps during :func:`homogenizeMesh`
                min_edge: float: Minimum allowed edge length
                max_edge: float: Maximum allowed edge length
                n_cv: int: Amount of conduction velocity distribution files to make.
                n_col_retry: int: Amount of times you want to try to fix colinearities during TetGen command.
                keep_intmed: bool: Unused. For completeness."""
            pipeline_ = OrderedDict([
                (
                    "Adding endo and epi point layers",
                    {'function': [self_.splitLayer],
                     'args': [{}]}
                ),
                (
                    "Refining surface mesh",
                    {'function': [self_.homogenizeMesh],
                     'args': [{'nsteps': refine_steps, 'boxplot': boxplot, 'edge_range': (min_edge, max_edge)}]}

                ),
                (
                    "Writing to .smesh",
                    {'function': [writeToSmesh],
                     'args': [{'mesh': self_.mesh, 'name': self_.dir + self_.name}]}
                ),
                (
                    "Tetrahedralising with TetGen",
                    {'function': [self_.tetrahedralise],
                     'args': [{'switches': switches, 'n_col_retry': n_col_retry}]}
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
            if n_cv:
                pipeline_.update({"Applying conduction velocities":
                                      {'function': [self.applyCV],
                                       'args': [{'write_VTK_file': True}]
                                       }
                                  })
            return pipeline_

        print("\n####### Creating 3D tetrahedralized {}\n".format(self.name))
        start = time.time()
        step = 1
        pipeline = getPipeline(self, **self.settings['reconstruction_parameters'])
        for action_name in pipeline:
            print("\n---- {}. {}\n".format(step, action_name))
            action = pipeline[action_name]
            for func, args in zip(action['function'], action['args']):
                func(**args)
            # Reinitialise pipeline dict to update function arguments
            pipeline = getPipeline(self, **self.settings['reconstruction_parameters'])
            step += 1

        if not self.settings['reconstruction_parameters']['keep_intmed']:  # delete intermediate files
            print("\n----- {}. Deleting intermediate files\n".format(step))
            # remove intermediate files that weren't present before generating and aren't the .vtk file
            to_delete = ["mid.txt", "endo.txt", "epi.txt",
                         "TrianglesSection.csv", "VerticesSection.csv", 'VerticesAttributesSection.csv',
                         'VerticesColorsSection.csv',
                         self.name + ".smesh", self.name + '.1.mtr', self.name + '.1.mtr', self.name + '.1.p2t',
                         self.name + '.1.node', self.name + '.1.edge', self.name + '.1.face', self.name + '.1.ele',
                         'myo.csv', 'noncond.csv']
            for trash in to_delete:
                if os.path.exists(self.dir + trash):
                    os.remove(self.dir + trash)
                    print("\tDeleted ", trash)
            step += 1

        duration = time.time() - start
        print("\n\tMesh reconstructed in {0}m {1:.2f}s".format(int(duration // 60), duration % 60))
        self.__update(pv.read(self.dir + self.name + '.1.vtk'))

    def plot(self) -> None:
        """
        Plots the mesh in its current state using PyVista's Plotter() method.

        Returns:
            None: Nothing. Opens a VtkRenderer window to show the plot.
        """
        p = pv.Plotter()
        p.add_mesh(self.mesh)
        p.show()
