#!/usr/bin/env python3
"""Physics-aware checks for epw_plrn auxiliary files."""

import glob
import math
import re
import sys
from pathlib import Path

RYD2EV = 13.6057
RYD2MEV = 1000.0 * RYD2EV


def fail(message):
    raise RuntimeError(message)


def close(a, b, rtol=1.0e-7, atol=1.0e-10):
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def read_lines(path):
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        fail("{}: missing or empty".format(path))
    data = [line.strip() for line in p.read_text().splitlines() if line.strip()]
    if not data:
        fail("{}: empty".format(path))
    return p, data


def read_floats(fields, path):
    try:
        values = [float(x) for x in fields]
    except ValueError:
        fail("{}: non-numeric data row".format(path))
    if not all(math.isfinite(x) for x in values):
        fail("{}: NaN or Inf".format(path))
    return values


def parse_wf_header(data, path, uniform):
    head = data[0].split()
    if head[0] == "Scell":
        if len(head) != 4:
            fail("{}: bad Scell header".format(path))
        nktotf, nbndsub, nstate = [int(x) for x in head[1:]]
    else:
        if len(head) != 6:
            fail("{}: bad header".format(path))
        nkf1, nkf2, nkf3, nktotf, nbndsub, nstate = [int(x) for x in head]
        if uniform and nkf1 * nkf2 * nkf3 != nktotf:
            fail("{}: k-grid size mismatch".format(path))
    if min(nktotf, nbndsub, nstate) < 1:
        fail("{}: invalid dimensions".format(path))
    return nktotf, nbndsub, nstate


def check_amp(path="Amp.plrn"):
    p, data = read_lines(path)
    nktotf, nbndsub, nstate = parse_wf_header(data, path, True)
    rows = [read_floats(row.split(), p) for row in data[1:]]
    if len(rows) != nktotf * nbndsub * nstate or any(len(row) != 2 for row in rows):
        fail("{}: unexpected coefficient shape".format(path))

    norms = [0.0] * nstate
    for i, row in enumerate(rows):
        norms[i % nstate] += row[0] * row[0] + row[1] * row[1]
    # A_mp = (1/N_p) sum_nk A_nk exp(ikR_p) U^dagger_mnk, so by unitarity of U
    # sum_mp |A_mp|^2 = (1/N_p) sum_nk |A_nk|^2 = 1, since polaron_scf holds
    # sum_nk |A_nk|^2 at N_p through norm_plrn_wf(eigvec, REAL(nktotf)).
    expected = 1.0
    for i, norm in enumerate(norms, 1):
        if not close(norm, expected, rtol=3.0e-3, atol=3.0e-6):
            fail("{}: state {} norm {:.8g}, expected {:.8g}".format(path, i, norm, expected))
    return {"nktotf": nktotf, "nstate": nstate, "norms": norms}


def check_ank(path="Ank.plrn", require_norm=True, uniform=True):
    p, data = read_lines(path)
    nktotf, _nbndsub, nstate = parse_wf_header(data, path, uniform)
    rows = [read_floats(row.split(), p) for row in data[1:]]
    if not rows or any(len(row) != 6 for row in rows):
        fail("{}: expected k, band, energy, Re, Im, |A| rows".format(path))
    if len(rows) % (nktotf * nstate):
        fail("{}: row count is inconsistent with the header".format(path))

    nbnd = len(rows) // (nktotf * nstate)
    norms = [0.0] * nstate
    for i, row in enumerate(rows):
        ik = i // (nbnd * nstate) + 1
        ibnd = (i // nstate) % nbnd + 1
        if int(round(row[0])) != ik or int(round(row[1])) != ibnd:
            fail("{}: k/band ordering mismatch".format(path))
        magnitude = math.hypot(row[3], row[4])
        if not close(magnitude, row[5], rtol=3.0e-6, atol=3.0e-7):
            fail("{}: stored |A| disagrees with Re/Im".format(path))
        norms[i % nstate] += row[5] * row[5]

    if require_norm:
        # norm_plrn_wf(eigvec, REAL(nktotf, DP)) is applied every SCF iteration
        for i, norm in enumerate(norms, 1):
            if not close(norm, float(nktotf), rtol=3.0e-4, atol=3.0e-5):
                fail("{}: state {} norm {:.8g}, expected {}".format(path, i, norm, nktotf))
    return {"nktotf": nktotf, "nstate": nstate, "nstates": nktotf * nbnd, "norms": norms}


def check_bmat(path="Bmat.plrn", uniform=True):
    p, data = read_lines(path)
    head = data[0].split()
    if head[0] == "Scell":
        if len(head) != 3:
            fail("{}: bad Scell header".format(path))
        nqtotf, nmodes = [int(x) for x in head[1:]]
    else:
        if len(head) != 5:
            fail("{}: bad header".format(path))
        nqf1, nqf2, nqf3, nqtotf, nmodes = [int(x) for x in head]
        if uniform and nqf1 * nqf2 * nqf3 != nqtotf:
            fail("{}: q-grid size mismatch".format(path))

    rows = [read_floats(row.split(), p) for row in data[1:]]
    if len(rows) != nqtotf * nmodes or any(len(row) != 6 for row in rows):
        fail("{}: unexpected B(q,nu) shape".format(path))

    sum_b2 = 0.0
    sum_wb2 = 0.0
    for i, row in enumerate(rows):
        iq, mode, omega, re_part, im_part, magnitude = row
        if int(round(iq)) != i // nmodes + 1 or int(round(mode)) != i % nmodes + 1:
            fail("{}: q/mode ordering mismatch".format(path))
        if not close(math.hypot(re_part, im_part), magnitude, rtol=3.0e-8, atol=3.0e-10):
            fail("{}: stored |B| disagrees with Re/Im".format(path))
        sum_b2 += magnitude * magnitude
        sum_wb2 += omega / RYD2MEV * magnitude * magnitude
    if sum_b2 <= 0.0:
        fail("{}: zero B(q,nu) norm".format(path))
    return {"nrows": len(rows), "sum_b2": sum_b2, "sum_wb2": sum_wb2}


def check_dtau(path="dtau.plrn"):
    p, data = read_lines(path)
    head = data[0].split()
    if head[0] == "Scell":
        if len(head) != 3:
            fail("{}: bad Scell header".format(path))
        nqtotf, nmodes = [int(x) for x in head[1:]]
    else:
        if len(head) != 5:
            fail("{}: bad header".format(path))
        nqf1, nqf2, nqf3, nqtotf, nmodes = [int(x) for x in head]
        if nqf1 * nqf2 * nqf3 != nqtotf:
            fail("{}: q-grid size mismatch".format(path))
    rows = [read_floats(row.split(), p) for row in data[1:]]
    if len(rows) != nqtotf * nmodes or any(len(row) != 2 for row in rows):
        fail("{}: unexpected displacement shape".format(path))
    if not any(math.hypot(row[0], row[1]) > 1.0e-14 for row in rows):
        fail("{}: all displacements are zero".format(path))


def trapz(x, y):
    return sum(0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]) for i in range(1, len(x)))


def check_integral(name, value, expected, rtol=5.0e-2):
    atol = max(1.0e-5, 1.0e-5 * abs(expected))
    if not close(value, expected, rtol=rtol, atol=atol):
        fail("dos.plrn: {} integral {:.8g}, expected {:.8g}".format(name, value, expected))


def check_dos(ank, bmat, path="dos.plrn"):
    p, data = read_lines(path)
    rows = []
    for row in data:
        if row.startswith("#"):
            continue
        values = read_floats(row.split(), p)
        if len(values) != 7:
            fail("{}: expected seven columns".format(path))
        rows.append(values)
    if len(rows) < 100:
        fail("{}: too few DOS points".format(path))

    egrid = [row[0] for row in rows]
    edos = [row[1] for row in rows]
    edos_all = [row[2] for row in rows]
    pgrid = [row[3] for row in rows]
    pdos = [row[4] for row in rows]
    pdos_all = [row[5] for row in rows]
    sdos = [row[6] for row in rows]
    if any(b <= a for a, b in zip(egrid, egrid[1:])):
        fail("{}: electron grid is not increasing".format(path))
    if any(b <= a for a, b in zip(pgrid, pgrid[1:])):
        fail("{}: phonon grid is not increasing".format(path))
    for name, values in (("A^2", edos), ("edos", edos_all), ("B^2", pdos), ("pdos", pdos_all)):
        if min(values) < -1.0e-10:
            fail("{}: negative {}".format(path, name))

    check_integral("A^2", trapz(egrid, edos), ank["norms"][0])
    check_integral("electronic DOS", trapz(egrid, edos_all), ank["nstates"])
    check_integral("B^2", trapz(pgrid, pdos), bmat["sum_b2"])
    check_integral("phonon DOS", trapz(pgrid, pdos_all), bmat["nrows"])
    check_integral("S-DOS", trapz(pgrid, sdos), bmat["sum_wb2"], rtol=8.0e-2)


def check_psir():
    names = sorted(set(glob.glob("psir_plrn*.xsf") + glob.glob("psir_plrn*.csv")))
    if not names:
        fail("psir_plrn output is missing")
    bad = re.compile(r"(?i)(?<![A-Za-z])(?:nan|[+-]?inf(?:inity)?)(?![A-Za-z])")
    for name in names:
        p = Path(name)
        if p.stat().st_size == 0:
            fail("{}: empty".format(name))
        if bad.search(p.read_text(errors="replace")):
            fail("{}: NaN or Inf".format(name))


def check_scf():
    amp = check_amp()
    ank = check_ank()
    if amp["nktotf"] != ank["nktotf"]:
        fail("Amp.plrn and Ank.plrn use different k grids")
    for anorm, wnorm in zip(ank["norms"], amp["norms"]):
        if not close(wnorm * amp["nktotf"], anorm, rtol=3.0e-3, atol=3.0e-5):
            fail("Amp/Ank normalization conventions disagree")
    bmat = check_bmat()
    check_dtau()
    check_dos(ank, bmat)


def run(mode):
    if mode == "scf":
        check_scf()
    elif mode == "interp":
        check_ank("Ank.band.plrn", require_norm=False, uniform=False)
        check_bmat("Bmat.band.plrn", uniform=False)
        check_psir()
    elif mode == "psir":
        check_psir()
    else:
        fail("unknown check mode: {}".format(mode))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} MODE".format(Path(sys.argv[0]).name), file=sys.stderr)
        sys.exit(2)
    try:
        run(sys.argv[1])
    except (OSError, RuntimeError, ValueError) as exc:
        print("polaron file check failed: {}".format(exc), file=sys.stderr)
        sys.exit(1)
    print("polaron file check: {} ok".format(sys.argv[1]))
