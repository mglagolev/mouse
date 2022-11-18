#!/usr/bin/python3
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MOUSE
# Copyright (c) 2022 Mikhail Glagolev
#
# Released under the GNU Public Licence, v2 or any higher version

import argparse
import MDAnalysis as mda

"""
This utility reads LAMMPS data file, and writes out the configuration
in the PDB format

Possible options are:
    --no-pbc-bonds
    hide the bonds which are not between the nearest images
    of the particles, used for visualisation
"""

parser = argparse.ArgumentParser(
    description = 'Convert LAMMPS data file into PDB')

parser.add_argument(
    'input', metavar = 'LAMMPS_DATA', action = "store", help = "input file")

parser.add_argument(
    'output', metavar = 'PDB', action = "store", help = "output file")

parser.add_argument(
    "--hide-pbc-bonds", action = "store_true",
    help = "Hide the bonds transversing the periodic boundary conditions")

args = parser.parse_args()

u = mda.Universe(args.input, format = 'DATA')

if args.hide_pbc_bonds:
    minbox = min(u.dimensions) / 2.
    bonds_to_delete = [
        bond for bond in u.bonds if bond.length(pbc = False) > minbox]
    u.delete_bonds(bonds_to_delete)

u.add_TopologyAttr('name', ['']*u.atoms.n_atoms)

atomtypes = list(set([atom.type for atom in u.atoms]))

for atomtype in atomtypes:
    selection = u.select_atoms("type " + str(atomtype))
    selection.atoms.names = [str(atomtype)] * selection.atoms.n_atoms

u.atoms.write(args.output)
