"""File containing some DGM functions"""

import os
import numpy as np
import logging

""" Collection of helper methods"""


# Set filename with path for DG with _p{length}.txt appendix
def filename(name, folder, length):
    f = "{}_p{}.txt".format(name, length)
    return os.path.join(folder, f)


# Set filename with path for DG with _t{plottime}_p{length}
def tp_filename(name, extension, folder, time, length):
    f = "{}_t{}_p{}.{}".format(name, time, length, extension)
    return os.path.join(folder, f)


def distance(r1, r2):
    """Calculate the euclidean distance between 2 vectors

		Args:
			r1 (tuple): x,y,z of the first coordinate
			r2 (tuple): x,y,z of the second coordinate

		Return:
			distance between r1 and r2
	"""
    (x1, y1, z1) = r1
    (x2, y2, z2) = r2
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


def check_data_array(at, time):
    """Calculate the first AT time after time t with aperiodical data

		Args:
			at: array of at values of the considered node
			time: reference time for which at such be larger

		Return:
			updated AT, -1 for scar tissue
	"""
    if type(at) == int or type(at) == np.float64:
        return int(at)
    elif len(at) == 1:
        return at[0]
    else:
        if at == []:
            return -1
        else:
            i = 0
            while at[i] < time and i < len(at):
                i += 1
            return at[i]


def create_color(visu_mode, color, var, intervar, range):
    switcher = {
        "var": 1000 - var * 100,
        "range": range * 100,
        "intervar": intervar * 100,
    }
    return switcher.get(visu_mode, color)


def collectKeysandArrows(keys_to_plot, color_to_plot, coordinates):
    """Create the relevant scalars, colors, xs, ys, zs, x_dirs, y_dirs, z_dirs to visualize in VtkWindow

		Args:
			keys_to_plot: keys which have to be plotted
			color_to_plot: colors corresponding with the keys which have to be plotted
			coordinateValues: coordinates of considered electrodes

		Return:
			scalars: array of the scalar values of the considered electrodes
			colors: array of the color values of the considered electrodes
			xs: array of the x coordinates of the considered electrodes
			ys: array of the y coordinates of the considered electrodes
			zs: array of the z coordinates of the considered electrodes
			xs: array of the x_dir direction along unit x of the considered electrodes (empty when a sphere object)
			ys: array of the y_dir direction along unit y of the considered electrodes (empty when a sphere object)
			zs: array of the z_dir direction along unit z of the considered electrodes (empty when a sphere object)
	"""
    scalars, colors, xs, ys, zs, x_dirs, y_dirs, z_dirs = [], [], [], [], [], [], [], []
    for key in keys_to_plot.keys():
        for j in range(len(keys_to_plot[key])):
            color = color_to_plot[key][j]
            key2 = keys_to_plot[key][j]
            x, y, z = coordinates[key]
            x1, y1, z1 = coordinates[key2]
            x_dir, y_dir, z_dir = x1 - x, y1 - y, z1 - z
            scalar = np.sqrt(x_dir * x_dir + y_dir * y_dir + z_dir * z_dir)
            scalars.append(scalar)
            colors.append(color)
            xs.append(x)
            ys.append(y)
            zs.append(z)
            x_dirs.append(x_dir)
            y_dirs.append(y_dir)
            z_dirs.append(z_dir)
    return scalars, colors, xs, ys, zs, x_dirs, y_dirs, z_dirs


def get_AT(AT, t, period):
    """Set valid AT

		Return:
			new AT
	"""
    ATnew = AT
    if AT < -2 * period:
        return -1  # scar
    while ATnew < t:
        ATnew += period
    return int(ATnew)


def select_spherical1(xyz, unip, time_sep, distance_sep, period):
    from scipy.stats import norm
    from scipy.spatial import distance
    from sklearn.cluster import KMeans

    D = distance.squareform(distance.pdist(xyz))
    closest = np.argsort(D, axis=1)
    kmeans = KMeans(n_clusters=1000, random_state=0).fit(xyz)

    xyz_new = kmeans.cluster_centers_
    unip_new = np.zeros(len(xyz_new))
    for i in range(len(unip)):
        unip[i] = get_AT(unip[i], 0, period)

    labels = kmeans.labels_
    g = k = gamma = 0

    xyz_final = []
    unip_final = []
    for g in range(len(xyz_new)):
        indices = []
        for k in range(len(xyz)):
            if labels[k] == g:
                indices.append(k)
        lats = [unip[k] for k in indices]
        lats_indexen = np.argsort(lats)
        lats = np.sort(lats)
        lats_derived = np.diff(lats)

        (mu, sigma) = norm.fit(lats)

        if len(lats) < 2:
            continue
        else:
            if sigma < 10:  # if np.amax(lats_derived) - np.amin(lats_derived ) < 20:
                unip_final.append(mu)
                xyz_final.append(xyz_new[g])

    return xyz_final, unip_final  # coords, at


def select_spherical2(xyz, unip, time_sep, distance_sep):
    xyz_new = [xyz[:, 0]]
    unip_new = [unip[0, 0]]
    for i in range(1, len(xyz[0])):
        r1 = xyz[:, i]
        LAT1 = unip[0, i]
        add = False
        for j in range(len(xyz_new)):
            r2 = xyz_new[j]
            LAT2 = unip_new[j]
            if not add:
                if np.abs(LAT2 - LAT1) < time_sep and df.distance(r1, r2) < distance_sep:
                    add = True
        if not add:
            xyz_new = np.append(xyz_new, [r1], axis=0)
            unip_new = np.append(unip_new, LAT1)
    return xyz_new, unip_new  # coords, at


def matlab_to_coord(threshold_mv, bipolar_ms, unipolar_ms, bipolar_mV, unipolar_mV, vertices, maps, period):
    """Convert matlab to coordinates

		Return:
			coordinates and bipolar list
	"""
    logging.debug("Converting matlab to coordinates")
    if maps[0][0]['cutoutMask'].any():  # array has items
        cutoutMask = maps[0][0]['cutoutMask']
        cutoutMask = cutoutMask[0]
        nup = len(cutoutMask) - np.count_nonzero(cutoutMask)
        # print(np.amax(bipolar_mV[0]))
        # print(np.amin(bipolar_mV[0]))
        mask = (bipolar_mV[0] > threshold_mv) & (cutoutMask == 0)
        dimensions = np.shape(vertices)
        if dimensions[0] < dimensions[1]:
            vertices = np.transpose(vertices)
        xyz = vertices[mask]
        bip = bipolar_ms[0][mask]
        unip = unipolar_ms[0][mask]
        xyz, bip = select_spherical1(xyz, bip, 5, 5, period)
    else:
        cutoutMask = maps[0][0]['cutoutMask']
        cutoutMask = cutoutMask[0]
        mask = (bipolar_mV[0] > threshold_mv) & (cutoutMask == 0)
        xyz = vertices[:, mask]
        bip = bipolar_ms[:, mask]
        xyz, bip = select_spherical2(xyz, bip, 5, 5)

    return xyz, bip
