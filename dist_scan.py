#!/usr/bin/env python3
"""
move.py -- generate a series of Gaussian input files in which the distance
between two atoms is stretched (or compressed) along the vector that joins
them.

    python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1

writes TEMPLATE_00.gjf ... TEMPLATE_10.gjf, in which d(3,7) is increased by
0.0, 0.1, ..., 1.0 Angstrom with respect to the template.

By default the whole molecular fragment that atom J belongs to is translated
rigidly, so an intermolecular scan pulls the complete molecule away instead
of tearing one atom out of it.

Everything else in the template (Link0, route, title, charge/multiplicity,
connectivity, Gen basis sets, ModRedundant sections, ...) is copied
verbatim; only the atom lines and the %chk name are touched.

Rewritten (Python 3, argparse, internuclear displacement vector, automatic
fragment detection) from Tharwan's original move.py.
"""

import argparse
import math
import os
import re
import subprocess
import sys

CHK_RE = re.compile(r"^\s*%chk\s*=\s*\S+\s*$", re.IGNORECASE)
MODRED_RE = re.compile(r"modredundant", re.IGNORECASE)
OPT_RE = re.compile(
    r"\b(?P<kw>opt)\b"
    r"(?:\s*=\s*\((?P<p1>[^)]*)\)"      # opt=(a,b)
    r"|\s*\((?P<p2>[^)]*)\)"             # opt(a,b)
    r"|\s*=\s*(?P<bare>[^\s,()]+))?",     # opt=tight
    re.IGNORECASE)

# Cordero et al., Dalton Trans. 2008, 2832 (Angstrom)
COVALENT_RADII = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76,
    "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58, "Na": 1.66, "Mg": 1.41,
    "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39,
    "Mn": 1.39, "Fe": 1.32, "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22,
    "Ga": 1.22, "Ge": 1.20, "As": 1.19, "Se": 1.20, "Br": 1.20, "Kr": 1.16,
    "Rb": 2.20, "Sr": 1.95, "Y": 1.90, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54,
    "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39, "Ag": 1.45, "Cd": 1.44,
    "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.40,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01,
    "Pm": 1.99, "Sm": 1.98, "Eu": 1.98, "Gd": 1.96, "Tb": 1.94, "Dy": 1.92,
    "Ho": 1.92, "Er": 1.89, "Tm": 1.90, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.70, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36,
    "Au": 1.36, "Hg": 1.32, "Tl": 1.45, "Pb": 1.46, "Bi": 1.48, "Po": 1.40,
    "At": 1.50, "Rn": 1.50,
}
SYMBOLS = list(COVALENT_RADII)          # index 0 -> H, so Z = index + 1


def die(msg):
    sys.stderr.write("move.py: error: %s\n" % msg)
    sys.exit(1)


def warn(msg):
    sys.stderr.write("move.py: warning: %s\n" % msg)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------
class Atom(object):
    __slots__ = ("label", "frozen", "x", "y", "z", "element")

    def __init__(self, label, frozen, x, y, z):
        self.label = label
        self.frozen = frozen          # Gaussian freeze code (-1/0), or None
        self.x, self.y, self.z = x, y, z
        self.element = element_of(label)

    def copy(self):
        return Atom(self.label, self.frozen, self.x, self.y, self.z)

    def shift(self, v):
        self.x += v[0]
        self.y += v[1]
        self.z += v[2]

    def line(self):
        if self.frozen is None:
            return "%-4s %14.8f %14.8f %14.8f" % (self.label, self.x, self.y, self.z)
        return "%-4s %4d %14.8f %14.8f %14.8f" % (
            self.label, self.frozen, self.x, self.y, self.z)


def element_of(label):
    """Element symbol from a Gaussian atom label ('C', 'C12', '6', 'C-CA')."""
    core = re.split(r"[(\-]", label, 1)[0]
    if core.isdigit():                              # atomic number given
        z = int(core)
        return SYMBOLS[z - 1] if 1 <= z <= len(SYMBOLS) else "X"
    m = re.match(r"[A-Za-z]{1,2}", core)
    if not m:
        return "X"
    tag = m.group(0)
    two = tag.capitalize()
    if len(tag) == 2 and two in COVALENT_RADII:
        # 'Cl' yes, 'C1' impossible here, but 'CA'-style force-field labels
        # should still resolve to carbon
        if tag[1].islower() or tag[0].upper() not in COVALENT_RADII:
            return two
    one = tag[0].upper()
    return one if one in COVALENT_RADII else "X"


def distance(a, b):
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def unit_vector(a, b):
    """Normalised vector pointing from atom a to atom b."""
    d = distance(a, b)
    if d < 1e-8:
        die("atoms %s and %s sit on top of each other; "
            "the displacement vector is undefined" % (a.label, b.label))
    return ((b.x - a.x) / d, (b.y - a.y) / d, (b.z - a.z) / d)


def scaled(v, f):
    return (v[0] * f, v[1] * f, v[2] * f)


# --------------------------------------------------------------------------
# Gaussian input parsing
# --------------------------------------------------------------------------
def parse_atom(line):
    """Parse a Cartesian atom line; return an Atom or None."""
    tok = line.split()
    if len(tok) not in (4, 5):
        return None
    try:
        x, y, z = (float(t) for t in tok[-3:])
    except ValueError:
        return None
    frozen = None
    if len(tok) == 5:
        try:
            frozen = int(tok[1])
        except ValueError:
            return None
    return Atom(tok[0], frozen, x, y, z)


def blocks_of(lines):
    """Index ranges [start, end) of the blank-line-separated blocks."""
    blocks = []
    start = None
    for n, line in enumerate(lines):
        if line.strip():
            if start is None:
                start = n
        elif start is not None:
            blocks.append((start, n))
            start = None
    if start is not None:
        blocks.append((start, len(lines)))
    return blocks


def locate_atoms(lines):
    """Return (atom line numbers, Atoms, blocks, index of the route block).

    A Gaussian input is: [Link0 +] route / blank / title / blank /
    charge-multiplicity + atoms / blank / optional extra sections.
    """
    blocks = blocks_of(lines)
    route = None
    for n, (s, e) in enumerate(blocks):
        if any(l.lstrip().startswith("#") for l in lines[s:e]):
            route = n
            break
    if route is None:
        die("no route section found (no line starting with '#')")
    if len(blocks) < route + 3:
        die("no molecule specification found after the title section")

    s, e = blocks[route + 2]
    idx, atoms = [], []
    for n in range(s + 1, e):          # first line is charge/multiplicity
        atom = parse_atom(lines[n])
        if atom is None:
            die("line %d is not a Cartesian atom specification:\n    %s\n"
                "(Z-matrix, ONIOM and fragment input are not supported; "
                "save the template with Cartesian coordinates)"
                % (n + 1, lines[n].rstrip()))
        idx.append(n)
        atoms.append(atom)
    if not atoms:
        die("the molecule specification contains no atoms")
    return idx, atoms, blocks, route


def parse_selection(spec, natoms, offset, what):
    """'5', '5-9', '1,4,7-9' -> sorted list of 0-based indices."""
    out = set()
    for chunk in spec.replace(" ", "").split(","):
        if not chunk:
            continue
        if "-" in chunk[1:]:
            lo, hi = chunk[1:].split("-", 1)
            lo = chunk[0] + lo
        else:
            lo = hi = chunk
        try:
            lo, hi = int(lo), int(hi)
        except ValueError:
            die("cannot understand the atom selection %r in %s" % (chunk, what))
        if lo > hi:
            lo, hi = hi, lo
        for k in range(lo, hi + 1):
            n = k - offset
            if not 0 <= n < natoms:
                die("atom %d in %s is out of range (the template has %d atoms)"
                    % (k, what, natoms))
            out.add(n)
    if not out:
        die("%s selects no atoms" % what)
    return sorted(out)


# --------------------------------------------------------------------------
# connectivity / fragments
# --------------------------------------------------------------------------
def bonds_from_block(lines, blocks, route, natoms):
    """Read the connectivity section written by geom=connectivity."""
    if len(blocks) < route + 4:
        return None
    s, e = blocks[route + 3]
    adj = [set() for _ in range(natoms)]
    seen = 0
    for n in range(s, e):
        tok = lines[n].split()
        if not tok:
            continue
        try:
            centre = int(tok[0]) - 1
            partners = [int(t) - 1 for t in tok[1::2]]
        except ValueError:
            return None
        if not 0 <= centre < natoms or any(not 0 <= p < natoms for p in partners):
            return None
        for p in partners:
            adj[centre].add(p)
            adj[p].add(centre)
        seen += 1
    return adj if seen == natoms else None


def bonds_from_distance(atoms, scale):
    """Simple covalent-radius criterion: d(A,B) < scale * (rA + rB)."""
    n = len(atoms)
    radii = []
    unknown = set()
    for a in atoms:
        if a.element in COVALENT_RADII:
            radii.append(COVALENT_RADII[a.element])
        else:
            unknown.add(a.label)
            radii.append(1.50)
    if unknown:
        warn("unknown element(s) %s; assuming a covalent radius of 1.50 A"
             % ", ".join(sorted(unknown)))
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if distance(atoms[a], atoms[b]) < scale * (radii[a] + radii[b]):
                adj[a].add(b)
                adj[b].add(a)
    return adj


def component(adj, start, forbidden_edge):
    """Atoms reachable from `start`, ignoring one edge."""
    cut = frozenset(forbidden_edge)
    seen = {start}
    stack = [start]
    while stack:
        cur = stack.pop()
        for nxt in adj[cur]:
            if frozenset((cur, nxt)) == cut or nxt in seen:
                continue
            seen.add(nxt)
            stack.append(nxt)
    return sorted(seen)


def formula(atoms, group):
    counts = {}
    for n in group:
        counts[atoms[n].element] = counts.get(atoms[n].element, 0) + 1
    order = ([e for e in ("C", "H") if e in counts]
             + sorted(e for e in counts if e not in ("C", "H")))
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


def compact(group, offset):
    """[3,4,5,9] -> '4-6,10' (with offset = 1)."""
    out, run = [], []
    for n in group:
        if run and n == run[-1] + 1:
            run.append(n)
            continue
        if run:
            out.append(run)
        run = [n]
    if run:
        out.append(run)
    return ",".join(str(r[0] + offset) if len(r) == 1
                    else "%d-%d" % (r[0] + offset, r[-1] + offset) for r in out)


def shortest_path(adj, start, goal):
    """List of atoms from start to goal, or None if disconnected."""
    prev = {start: None}
    queue = [start]
    while queue:
        cur = queue.pop(0)
        if cur == goal:
            break
        for nxt in adj[cur]:
            if nxt not in prev:
                prev[nxt] = cur
                queue.append(nxt)
    if goal not in prev:
        return None
    path, cur = [], goal
    while cur is not None:
        path.append(cur)
        cur = prev[cur]
    return path[::-1]


def cut_candidates(adj, i, j):
    """Bonds on the i->j path whose removal separates i from j."""
    path = shortest_path(adj, i, j)
    if path is None:
        return []
    out = []
    for a, b in zip(path, path[1:]):
        if i not in component(adj, j, (a, b)):
            out.append((a, b))
    return out


def fragments_of(adj):
    seen, frags = set(), []
    for n in range(len(adj)):
        if n in seen:
            continue
        frag = component(adj, n, (-1, -1))
        seen.update(frag)
        frags.append(frag)
    return frags


def print_info(atoms, adj, offset, source):
    print("bonds from : %s" % source)
    print("")
    print("%5s %-4s %11s %11s %11s   %s"
          % ("atom", "el", "x", "y", "z", "bonded to"))
    for n, a in enumerate(atoms):
        print("%5d %-4s %11.5f %11.5f %11.5f   %s"
              % (n + offset, a.label, a.x, a.y, a.z,
                 compact(sorted(adj[n]), offset) or "-"))
    print("")
    frags = fragments_of(adj)
    print("%d fragment(s):" % len(frags))
    for k, frag in enumerate(frags, 1):
        print("  %2d: %3d atoms [%-10s] %s"
              % (k, len(frag), formula(atoms, frag), compact(frag, offset)))


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------
def write_job(lines, atom_idx, atoms, chk, path, freeze_line=None,
              add_chk=False):
    out = list(lines)
    for n, atom in zip(atom_idx, atoms):
        out[n] = atom.line()
    if chk is not None:
        found = False
        for n, line in enumerate(out):
            if CHK_RE.match(line):
                out[n] = "%chk=" + chk
                found = True
        if not found and add_chk:
            for n, line in enumerate(out):
                if line.lstrip().startswith("#"):
                    out.insert(n, "%chk=" + chk)
                    break
    text = "\n".join(out).rstrip("\n")
    if freeze_line:
        text += "\n\n" + freeze_line
    text += "\n\n"                       # Gaussian wants a trailing blank line
    with open(path, "w") as fh:
        fh.write(text)


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="move.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Build a rigid scan along the vector joining two atoms. "
                    "By default the whole fragment containing atom J is "
                    "translated.",
        epilog="""examples:
  # pull the molecule containing atom 7 away from atom 3, 0.1 A per step
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1

  # move both fragments apart symmetrically
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1 --mode both

  # old behaviour: displace the single atom J only
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1 --rigid atom

  # spell out what moves
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 12 --size 0.2 --group-j 7-15

  # relaxed scan: freeze the distance and let the rest optimise
  # (route section must already contain Opt=ModRedundant)
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1 --freeze

  # write the files and queue them
  python move.py TEMPLATE.gjf --atoms 3 7 --steps 10 --size 0.1 \\
      --submit 'bsub -q Batch24 "cd {dir}; g16 {file}"'
""")
    p.add_argument("template", help="Gaussian input file used as a template")
    p.add_argument("--atoms", "-a", nargs=2, type=int,
                   metavar=("I", "J"),
                   help="the two atoms defining the distance (1-based, "
                        "same numbering as GaussView)")
    p.add_argument("--steps", "-n", type=int, metavar="N",
                   help="number of steps; N+1 files are written, "
                        "step 0 being the untouched template geometry")
    p.add_argument("--size", "-s", type=float, metavar="D",
                   help="increment of d(I,J) per step, in Angstrom "
                        "(negative to shorten)")
    p.add_argument("--start", type=float, default=0.0, metavar="D",
                   help="displacement of step 0 (default: 0.0, i.e. the "
                        "template geometry)")
    p.add_argument("--mode", choices=("second", "first", "both"),
                   default="second",
                   help="which side moves: 'second' moves J's side (default), "
                        "'first' moves I's side, 'both' moves each by half")
    p.add_argument("--rigid", choices=("fragment", "atom"), default="fragment",
                   help="what travels with I and J when no explicit group is "
                        "given: the whole connected fragment (default) or "
                        "just the atom itself")
    p.add_argument("--group-i", metavar="SPEC",
                   help="atoms that travel rigidly with I, e.g. '3,10-14' "
                        "(overrides --rigid)")
    p.add_argument("--group-j", metavar="SPEC",
                   help="atoms that travel rigidly with J (overrides --rigid)")
    p.add_argument("--cut", nargs=2, type=int, metavar=("A", "B"),
                   help="bond to break when splitting the system into two "
                        "rigid parts (default: the I-J pair itself). Use this "
                        "when I and J are not directly bonded.")
    p.add_argument("--info", action="store_true",
                   help="print the numbering, the detected bonds and the "
                        "fragments of the template, then exit")
    p.add_argument("--bonds", choices=("auto", "connectivity", "distance"),
                   default="auto",
                   help="how fragments are detected: from the connectivity "
                        "section if present, else from covalent radii "
                        "(default: auto)")
    p.add_argument("--bond-scale", type=float, default=1.2, metavar="F",
                   help="two atoms are bonded if d < F * (r1 + r2) "
                        "(default: 1.2)")
    p.add_argument("--zero-based", action="store_true",
                   help="atom numbering starts at 0 instead of 1")
    p.add_argument("--outdir", "-o", default=None, metavar="DIR",
                   help="where to write the files (default: next to the "
                        "template)")
    p.add_argument("--prefix", default=None,
                   help="base name of the generated files "
                        "(default: template name without extension)")
    p.add_argument("--ext", default=None,
                   help="extension of the generated files "
                        "(default: same as the template)")
    p.add_argument("--freeze", nargs="?", const="both",
                   choices=("both", "bond", "atoms"), metavar="WHAT",
                   help="constrain the scanned coordinate for a relaxed scan "
                        "and add Opt=ModRedundant to the route if missing. "
                        "'both' (the default) writes 'X I F', 'X J F' and "
                        "'B I J F'; 'bond' only the B line; 'atoms' only the "
                        "two X lines")
    p.add_argument("--no-chk", action="store_true",
                   help="do not rename the %%chk file of each job")
    p.add_argument("--add-chk", action="store_true",
                   help="insert a %%chk line if the template has none")
    p.add_argument("--submit", metavar="CMD",
                   help="shell command run for every generated file; "
                        "{file}, {name}, {dir} and {chk} are substituted")
    p.add_argument("--dry-run", action="store_true",
                   help="only print what would be done")
    return p


def route_lines(lines, blocks, route):
    """Line numbers of the route section (from the first '#' to end of block)."""
    s, e = blocks[route]
    for n in range(s, e):
        if lines[n].lstrip().startswith("#"):
            return list(range(n, e))
    return []


def ensure_modredundant(lines, blocks, route):
    """Make sure the route asks for Opt=ModRedundant. Returns (lines, note)."""
    rlines = route_lines(lines, blocks, route)
    if not rlines:
        return lines, None
    if MODRED_RE.search(" ".join(lines[n] for n in rlines)):
        return lines, None                      # already there

    out = list(lines)
    for n in rlines:
        m = OPT_RE.search(out[n])
        if not m:
            continue
        inner = m.group("p1") or m.group("p2") or m.group("bare")
        if inner:
            new = "%s=(%s,modredundant)" % (m.group("kw"), inner.strip())
        else:
            new = "%s=modredundant" % m.group("kw")
        out[n] = out[n][:m.start()] + new + out[n][m.end():]
        return out, "route: %s -> %s" % (m.group(0).strip(), new)

    out[rlines[-1]] = out[rlines[-1]].rstrip() + " opt=modredundant"
    return out, ("route: opt=modredundant appended (the template had no Opt "
                 "keyword -- check that this is really the job you want)")


def build_adjacency(args, atoms, lines, blocks, route):
    adj = None
    source = None
    if args.bonds in ("auto", "connectivity"):
        adj = bonds_from_block(lines, blocks, route, len(atoms))
        if adj is not None:
            source = "connectivity section"
        elif args.bonds == "connectivity":
            die("no usable connectivity section found; use --bonds distance")
    if adj is None:
        adj = bonds_from_distance(atoms, args.bond_scale)
        source = "covalent radii (scale %.2f)" % args.bond_scale
    return adj, source


def explain_connected(args, atoms, adj, i, j, offset):
    """Atoms I and J are in one piece: say why, and how to fix it."""
    msg = ["atoms %d and %d belong to the same connected fragment, so there "
           "is no separate molecule to translate." % (args.atoms[0],
                                                      args.atoms[1])]
    path = shortest_path(adj, i, j)
    if path:
        msg.append("They are linked through: %s."
                   % " - ".join(str(n + offset) for n in path))
    cands = cut_candidates(adj, i, j)
    if cands:
        msg.append("Choose the bond that should break with --cut:")
        for a, b in cands:
            grp = component(adj, j, (a, b))
            msg.append("    --cut %-3d %-3d  moves %2d atoms [%s] %s"
                       % (a + offset, b + offset, len(grp),
                          formula(atoms, grp), compact(grp, offset)))
    else:
        msg.append("No single bond separates them (they sit in a ring). "
                   "Use --group-i/--group-j to list what should travel, "
                   "or --rigid atom to displace the two atoms alone.")
    msg.append("Run with --info to see the numbering and the fragments.")
    die("\n  ".join(msg))


def resolve_groups(args, atoms, lines, blocks, route, i, j, offset):
    """Decide which atoms travel with I and which with J."""
    natoms = len(atoms)
    explicit_i = bool(args.group_i)
    explicit_j = bool(args.group_j)
    group_i = (parse_selection(args.group_i, natoms, offset, "--group-i")
               if explicit_i else None)
    group_j = (parse_selection(args.group_j, natoms, offset, "--group-j")
               if explicit_j else None)
    source = "explicit"

    if explicit_i or explicit_j:
        # one side spelled out: the other just anchors on its own atom
        if group_i is None:
            group_i = [i]
        if group_j is None:
            group_j = [j]
        return group_i, group_j, "explicit"

    if args.rigid == "fragment":
        adj, source = build_adjacency(args, atoms, lines, blocks, route)

        if args.cut:
            cut = tuple(k - offset for k in args.cut)
            for k, n in zip(args.cut, cut):
                if not 0 <= n < natoms:
                    die("atom %d in --cut is out of range" % k)
            if cut[1] not in adj[cut[0]]:
                warn("atoms %d and %d are not bonded according to the "
                     "detected connectivity; cutting there anyway"
                     % tuple(args.cut))
            source += ", cut %d-%d" % tuple(args.cut)
        else:
            cut = (i, j)

        frag_j = component(adj, j, cut)
        frag_i = component(adj, i, cut)
        if i in frag_j:
            explain_connected(args, atoms, adj, i, j, offset)
        if group_i is None:
            group_i = frag_i
        if group_j is None:
            group_j = frag_j
    else:
        group_i, group_j = [i], [j]
        source = "single atoms (--rigid atom)"

    return group_i, group_j, source


def describe(atoms, group, offset, share):
    if share == 0.0:
        return "fixed"
    return "%d atom%s [%s] %s%s" % (
        len(group), "s" if len(group) > 1 else "", formula(atoms, group),
        compact(group, offset), " at half speed" if share == 0.5 else "")


def main(argv=None):
    args = build_parser().parse_args(argv)

    with open(args.template) as fh:
        lines = fh.read().splitlines()
    atom_idx, atoms0, blocks, route = locate_atoms(lines)
    natoms = len(atoms0)
    offset = 0 if args.zero_based else 1

    if args.info:
        adj, source = build_adjacency(args, atoms0, lines, blocks, route)
        print("template   : %s (%d atoms)" % (args.template, natoms))
        print_info(atoms0, adj, offset, source)
        return 0

    missing = [n for n, v in (("--atoms", args.atoms), ("--steps", args.steps),
                              ("--size", args.size)) if v is None]
    if missing:
        die("the following arguments are required: %s" % ", ".join(missing))
    if args.steps < 0:
        die("--steps must be >= 0")

    i, j = (k - offset for k in args.atoms)
    for k, n in zip(args.atoms, (i, j)):
        if not 0 <= n < natoms:
            die("atom %d is out of range (the template has %d atoms)"
                % (k, natoms))
    if i == j:
        die("--atoms needs two different atoms")

    group_i, group_j, source = resolve_groups(
        args, atoms0, lines, blocks, route, i, j, offset)
    if set(group_i) & set(group_j):
        die("the two groups overlap; the distance would not change")
    if i not in group_i:
        die("the group of atom %d must contain atom %d"
            % (args.atoms[0], args.atoms[0]))
    if j not in group_j:
        die("the group of atom %d must contain atom %d"
            % (args.atoms[1], args.atoms[1]))

    # displacement direction: the actual I -> J vector of the template
    axis = unit_vector(atoms0[i], atoms0[j])
    d0 = distance(atoms0[i], atoms0[j])

    if args.mode == "second":
        share_i, share_j = 0.0, 1.0
    elif args.mode == "first":
        share_i, share_j = 1.0, 0.0
    else:
        share_i, share_j = 0.5, 0.5

    moving = set(group_j) if share_j else set()
    moving |= set(group_i) if share_i else set()
    spectators = sorted(set(range(natoms)) - set(group_i) - set(group_j))
    if spectators:
        warn("%d atom(s) belong to neither group and will stay put: %s"
             % (len(spectators), compact(spectators, offset)))

    # output naming
    directory = args.outdir or os.path.dirname(os.path.abspath(args.template))
    stem, ext = os.path.splitext(os.path.basename(args.template))
    prefix = args.prefix or stem
    ext = args.ext or ext or ".gjf"
    if not ext.startswith("."):
        ext = "." + ext
    width = max(2, len(str(args.steps)))
    if not args.dry_run and not os.path.isdir(directory):
        os.makedirs(directory)

    freeze_line = None
    route_note = None
    if args.freeze:
        parts = []
        if args.freeze in ("atoms", "both"):       # Gaussian numbering is 1-based
            parts.append("X %d F" % (i + 1))
            parts.append("X %d F" % (j + 1))
        if args.freeze in ("bond", "both"):
            parts.append("B %d %d F" % (i + 1, j + 1))
        freeze_line = "\n".join(parts)
        lines, route_note = ensure_modredundant(lines, blocks, route)

    print("template  : %s (%d atoms)" % (args.template, natoms))
    print("distance  : d(%d,%d) = %.6f A"
          % (args.atoms[0], args.atoms[1], d0))
    print("axis      : %.6f %.6f %.6f" % axis)
    print("fragments : %s" % source)
    if args.freeze:
        print("constraint: %s" % freeze_line.replace("\n", " | "))
        print("%s" % ("            " + route_note if route_note
                      else "            route already had Opt=ModRedundant"))
    print("  with %-4d: %s" % (args.atoms[0],
                               describe(atoms0, group_i, offset, share_i)))
    print("  with %-4d: %s" % (args.atoms[1],
                               describe(atoms0, group_j, offset, share_j)))
    print("")
    print("%-6s %12s %14s   %s" % ("step", "shift/A", "d(I,J)/A", "file"))

    for step in range(args.steps + 1):
        shift = args.start + step * args.size
        atoms = [a.copy() for a in atoms0]
        for n in group_i:
            atoms[n].shift(scaled(axis, -shift * share_i))
        for n in group_j:
            atoms[n].shift(scaled(axis, shift * share_j))

        name = "%s_%0*d" % (prefix, width, step)
        path = os.path.join(directory, name + ext)
        chk = name + ".chk"
        print("%-6d %12.4f %14.6f   %s"
              % (step, shift, distance(atoms[i], atoms[j]), path))
        if args.dry_run:
            continue

        write_job(lines, atom_idx, atoms,
                  None if args.no_chk else chk, path, freeze_line,
                  args.add_chk)

        if args.submit:
            cmd = args.submit.format(file=path, name=name,
                                     dir=directory, chk=chk)
            rc = subprocess.call(cmd, shell=True)
            if rc != 0:
                warn("submission of %s returned %d" % (path, rc))
    return 0


if __name__ == "__main__":
    sys.exit(main())
