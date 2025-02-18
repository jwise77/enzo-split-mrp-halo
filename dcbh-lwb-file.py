import ytree
import matplotlib.pyplot as plt
import h5py as h5
import glob
import unyt as u
import numpy as np
import yt.utilities.cosmology as co

# Load arbor and get unique DCBH host halo candidates
def load_halos(filename, dcbh_host=False, halo_id=None):
    arbor = ytree.load(filename)
    if dcbh_host:
        dcbh_halos = list(arbor.select_halos("tree['forest', 'dcbh_host'] == 1"))
        filtered_halos = []
        for h in dcbh_halos:
            nmarked = h["tree", "dcbh_host"].sum()
            if nmarked == 1:
                filtered_halos.append(h)
    elif halo_id != None:
        filtered_halos = list(arbor.select_halos(f"tree['forest', 'uid'] == {halo_id}"))
    return arbor, filtered_halos


def position_history(arbor, node):

    directory_prefix = "DD"
    prefix = "output_"

    if node["prog"] is None:
        line = []
    else:
        line = list(node["prog"])[::-1]
    while node.descendent is not None:
        line.append(node.descendent)
        node = node.descendent

    nnodes = len(line)
    pos = arbor.arr(np.zeros((nnodes, 3)), units="Mpc/h")
    vel = arbor.arr(np.zeros((nnodes, 3)), units="km/s")
    for i, halo in enumerate(line):
        pos[i, :] = halo["position"].to("Mpc/h")
        vel[i, :] = halo["velocity"].to("km/s")
    vavg = vel.mean(0)

    # Output redshifts
    znum = []
    zz = []
    prefixlen = len(directory_prefix)
    with open("redshifts.dat") as fp:
        for fline in fp:
            znum.append(int(fline[prefixlen:prefixlen+4]))
            zz.append(float(fline.split("=")[1]))
        zsum = np.array(znum)
        zz = np.array(zz)

    zmax = line[0]["redshift"]
    znum_max = np.argmin(np.abs(zmax - zz))

    # Extrapolation to higher redshifts if necessary, assuming constant (avg) velocity
    cos = co.Cosmology(
        hubble_constant=arbor.hubble_constant,
        omega_lambda=arbor.omega_lambda,
        omega_matter=arbor.omega_matter,
    )

    pos_ex = arbor.arr(np.zeros((znum_max, 3)), units="Mpc/h")
    current_pos = pos[0, :].copy()
    for num in range(znum_max, 0, -1):
        t0 = cos.t_from_z(zz[num])
        t1 = cos.t_from_z(zz[num - 1])
        dt = t1 - t0
        a = 1 / (zz[num] + 1)
        current_pos += (vavg / a) * dt
        pos_ex[num - 1, :] = current_pos

    allpos = arbor.arr(np.concatenate([pos_ex, pos]), "Mpc/h")
    redshift = {}
    center = {}
    time = {}
    for ii in range(allpos.shape[0]):
        ds = f"{prefix}{znum[ii]:04d}"
        redshift[ds] = zz[ii]
        center[ds] = allpos[ii]
        time[ds] = cos.t_from_z(zz[ii])
    return redshift, time, center


def load_star_data(stardir="star-data"):
    sdata = {}
    first_time = True
    attributes = dict(
        top=["time_unit", "boxsize"], output=["current_redshift", "current_time"]
    )
    for fn in sorted(glob.glob("%s/*.h5" % (stardir))):
        ds = fn.split("/")[-1].split("-stars")[0]
        sdata[ds] = {}
        fp = h5.File(fn, "r")
        if first_time:
            for attr in attributes["top"]:
                sdata[attr] = fp.attrs.get(attr)
            first_time = False
        for attr in attributes["output"]:
            sdata[ds][attr] = fp.attrs.get(attr)
        for k in fp.keys():
            sdata[ds][k] = fp[k][:]
        fp.close()

    for k, dd in sdata.items():
        if k in attributes["top"]:
            continue
        idsort = np.argsort(dd["particle_index"])
        for f in dd.keys():
            if f in attributes["output"]:
                continue
            dd[f] = dd[f][idsort]
        # Correct Pop III masses for SN explosions (1e-20 factor)
        SN3 = (dd["particle_mass"] < 1e-10) & (dd["type"] == 3)
        dd["particle_mass"][SN3] *= 1e20

    # If time_unit didn't exist in data file, we have to manually grab it from the parameter file
    if sdata["time_unit"] == None:
        sdata["time_unit"] = 6.5927815e14 / 3.1557e13  # Myr
    return sdata


def find_adj_ds(z, allz):
    nz = len(allz.keys())
    za = np.zeros(nz)
    ds = []
    _ds = []
    for i, k in enumerate(allz.keys()):
        za[i] = allz[k]
        _ds.append(k)
    zsort = np.argsort(za)
    za = za[zsort]
    for zs in zsort:
        ds.append(_ds[zs])
    iz = np.searchsorted(za, z)
    factor = 0.0
    dt = 0.0  # Myr
    if iz == 0:
        ds0 = ds[iz]
        ds1 = None
    elif iz == nz:
        ds0 = None
        ds1 = ds[iz - 1]
    else:
        ds1 = ds[iz - 1]
        ds0 = ds[iz]
        # Convert redshift to time (not normalized since we're dealing with fractional values)
        # Assumes an EdS universe
        t = (z + 1) ** (-2.0 / 3)
        t0 = (za[iz - 1] + 1) ** (-2.0 / 3)
        t1 = (za[iz] + 1) ** (-2.0 / 3)
        factor = (t - t0) / (t1 - t0)
    return ds0, ds1, factor


def p3life(mass):
    # Schaerer (2002)
    # Input: solar masses
    # Output: Lifetime in Myr
    x = np.log10(mass)
    return 10 ** (3.785 - 3.759 * x + 1.413 * x * x - 0.186 * x * x * x)


def p3lw(mass):
    # Schaerer (2002)
    # Input: solar masses
    # Output: LW photon luminosity (ph/s)
    x = np.log10(mass)
    return 10 ** (44.03 + 4.59 * x - 0.77 * x * x)


def interp_particles(data, ds0, ds1, dt, fint):
    nstars1 = data[ds1]["particle_index"].size
    luminosity_t = np.zeros(nstars1)
    position_t = np.zeros((nstars1, 3))

    #
    # Direct interpolation for star particles that exist in both datasets
    #
    exists = np.isin(
        data[ds0]["particle_index"], data[ds1]["particle_index"], invert=False
    )
    # idel = np.where(data[ds0]['particle_index'][deleted])
    match = np.isin(data[ds1]["particle_index"], data[ds0]["particle_index"][exists])
    dr = data[ds1]["particle_position"][match] - data[ds0]["particle_position"][exists]
    position_t[match] = data[ds0]["particle_position"][exists] + fint * dr
    # print position_t[match]

    # Search for massive star particles that are still radiating at the given time
    # Pop III stars
    p3_lifetime = p3life(data[ds1]["particle_mass"])
    p2_lifetime = 20.0
    p3 = (
        (data[ds1]["type"] == 3)
        & (data[ds1]["age"] - fint * dt < p3_lifetime)
        & (data[ds1]["age"] - fint * dt > 0)
    )
    # Metal-enriched stars
    p2 = (
        (data[ds1]["type"] == 2)
        & (data[ds1]["age"] - fint * dt < p2_lifetime)
        & (data[ds1]["age"] - fint * dt > 0)
    )

    # Calculate ionizing photon luminosity for the living star particles
    specific_LW = (
        1.288 * 1.12e46
    )  # 1.12e46 ionizing ph/s, 1.288x more LW photons (Schaerer 2003)
    luminosity_t[p3] = p3lw(data[ds1]["particle_mass"][p3])
    luminosity_t[p2] = specific_LW * data[ds1]["particle_mass"][p2]

    #
    # For stars that only exist in the latter dataset, we use the velocity to extrapolate
    # their positions backwards.
    #
    # Discard star particles that haven't been born yet
    new_stars = np.logical_not(match)
    position_t[new_stars] = (
        data[ds1]["particle_position"][new_stars]
        - fint * (dt / data["time_unit"]) * data[ds1]["particle_velocity"][new_stars]
    )
    # print position_t[new_stars]
    # if data[ds1]["age"].min() < 20:
    #     import pdb; pdb.set_trace()
    # Result dictionary with only radiating particle positions and LW luminosities
    nonz = np.nonzero(luminosity_t)[0]
    idata = dict(position=position_t[nonz], luminosity=luminosity_t[nonz])
    return idata


def scale(val):
    minv = 48.0
    maxv = 52.0
    srange = (1, 100)
    ds = (srange[1] - srange[0]) / (maxv - minv)
    return (np.log10(val) - minv) * ds


def calc_JLW_at_point(idata, point, cm_boxsize, z):
    N_to_J = 4.80814e-28  # ph/s/cm^2 -> erg/s/cm^2/Hz/sr
    if not isinstance(point, np.ndarray):
        point = np.array(point)
    if point.size != 3:
        raise RuntimeError("Input (point) must be 3D position")
    scale = cm_boxsize * 3.086e24 / (1 + z)
    r2 = ((point - idata["position"]) ** 2).sum(1) * scale**2

    # photon flux
    NLW = idata["luminosity"] / (4 * np.pi * r2)

    # Convert to intensity
    JLW = NLW * N_to_J

    return JLW.sum()


def plot_lwb(allz, JLW, zdcbh=None, hid=0):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.semilogy(allz, JLW * 1e21, lw=3)
    if zdcbh is not None:
        ax.axvline(zdcbh, lw=1, ls=":")
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"J$_{LW}$/J$_{21}$")
    plt.title(f'Halo {hid}')
    plt.savefig(f"dcbh-lwb-{hid}.png")
    plt.close(fig)

    # Background from Wise & Abel (2005), Wise et al. (2012b)
    zz = allz + 1
    log_LWB = (
        -23.56688
        + 4.5213e-1 * zz
        - 2.679832e-2 * zz**2
        + 5.88234e-4 * zz**3
        - 5.05576e-6 * zz**4
    )
    LWB = 10**log_LWB
    JLW_total = JLW + LWB

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.semilogy(allz, JLW * 1e21, lw=2, ls="--", label="Local Sources")
    ax.semilogy(allz, LWB * 1e21, lw=2, ls=":", label="BG (Wise+ 2012)")
    ax.semilogy(allz, JLW_total * 1e21, c="k", lw=4, ls="-", label="Combined")
    if zdcbh is not None:
        ax.axvline(zdcbh, lw=1, ls=":", label="DCBH Form Estimate")
    ax.set_xlabel("Redshift")
    ax.set_ylabel(r"J$_{LW}$/J$_{21}$")
    ax.legend(loc="best")
    plt.title(f'Halo {hid}')
    plt.savefig(f"dcbh-lwb-total-{hid}.png")
    plt.close(fig)


def calc_lwb(redshift, time, center, sdata, zrange, dz=0.05):
    allz = np.mgrid[zrange[0] : zrange[1] + 0.1 * dz : dz]
    boxsize = sdata["boxsize"]  # comoving Mpc
    JLW = np.zeros(allz.size)
    for i, z in enumerate(allz):
        ds0, ds1, factor = find_adj_ds(z, redshift)
        if ds0 == None or ds1 == None:
            print("Redshift (%f) out of range. No data. Skipping." % (z))
            continue
        # print z, redshift[ds0], redshift[ds1], ds0, ds1, factor
        dt = time[ds1] - time[ds0]
        dpoint = center[ds1] - center[ds0]
        point = center[ds0] + factor * dpoint
        print("Redshift = %f (%s -> %s; dt = %.1f Myr)" % (z, ds0, ds1, dt.to('Myr')))
        iparticles = interp_particles(sdata, ds0, ds1, dt.to("Myr").v, factor)
        JLW[i] = calc_JLW_at_point(iparticles, point.to("unitary").v, boxsize, z)
    return allz, JLW


def k31(J21):
    return 1.42e9 * (J21)


def write_grackle_file(hid, allz, JLW):
    outfile = f"dcbh-lwb-halo-{hid}.h5"
    info_string = f"dcbh-halo-{hid}"

    rate_names = [
        "k24",
        "k25",
        "k26",
        "k27",
        "k28",
        "k29",
        "k30",
        "k31",
        "piHI",
        "piHeII",
        "piHeI",
    ]

    outdata = dict(z=allz, k31=k31(JLW))
    rate0 = np.zeros(outdata["z"].size)

    fp = h5.File(outfile, "w")
    fp["UVBRates/Info"] = np.bytes_(info_string)
    fp["UVBRates/z"] = outdata["z"]

    for name in rate_names:
        if name == "z":
            fname = "/UVBRates/%s" % (name)
        elif name.startswith("k"):
            fname = "/UVBRates/Chemistry/%s" % (name)
        else:
            fname = "/UVBRates/Photoheating/%s" % (name)
        if name not in outdata.keys():
            fp[fname] = rate0
        else:
            fp[fname] = outdata[name]
    fp.close()


def plot_tree(node, hid=0):
    def my_node(halo):
        prog = list(halo.find_root()["prog", "uid"])
        color = "black"
        if "dcbh_host" in halo.arbor.field_list:
            if halo["dcbh_host"] == 1:
                color = "blue"
            elif halo["uid"] in prog:
                color = "red"
        else:
            if halo["uid"] in prog:
                color = "red"
        label = """
        id: %d // z = %.2f
        mass: %.2e Msun
        """ % (
            halo["uid"],
            halo["redshift"],
            halo["mass"].to("Msun"),
        )
        my_kwargs = {"label": label, "fontsize": 8, "shape": "square", "color": color}
        return my_kwargs

    p = ytree.TreePlot(node, dot_kwargs={"rankdir": "BT"}, node_function=my_node)
    p.save(f"dcbh-tree-{hid}.png")


def main():
    zrange = (11.18, 20)
    halo_id = 26813
    arbor_file = "full_arbor.h5"
    arbor, filtered_halos = load_halos(arbor_file, halo_id=halo_id)
    for tree in filtered_halos:
        hid = tree["uid"]
        zdcbh = tree["redshift"]
        redshift, time, center = position_history(arbor, tree)
        star_data = load_star_data()
        zz, lwb = calc_lwb(redshift, time, center, star_data, zrange)
        plot_lwb(zz, lwb, hid=hid, zdcbh=zdcbh)
        plot_tree(tree, hid)
        write_grackle_file(hid, zz, lwb)


if __name__ == "__main__":
    main()
