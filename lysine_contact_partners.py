#!/usr/bin/env python3
"""
lysine_contact_partners.py
===========================
For a chosen set of lysines on one chain (default: FOXG1's 5 conserved
interface lysines), finds exactly which residue(s) on the partner chain
(HDAC1 / HDAC2 / HDAC6) are in direct atomic contact -- answering
"which HDAC residue touches K181?" rather than just "is K181 at the
interface?" (which is all dASA/interfaceResidues tells you).

Reports two things per lysine, since they answer slightly different
questions:
  - ANY heavy-atom contact within --cutoff  (who's nearby, structurally)
  - NZ-specific contacts within --nz-cutoff (who's near the reactive
    amine tip that HDAC's catalytic tunnel would need to grab onto --
    the mechanistically relevant one for a deacetylation hypothesis)

Also flags candidate salt bridges (NZ near a Asp/Glu carboxylate O,
within salt-bridge range) since that's a common and easy-to-spot
functional signature at a lysine-enzyme interface.

USAGE
-----
    python lysine_contact_partners.py \
        --cif fold_foxg1_hdac1_model_0.cif \
        --foxg1-chain A --hdac-chain B \
        --lysines 15,181,237,241,272

    # different partner, same lysines, different cutoff:
    python lysine_contact_partners.py --cif fold_foxg1_hdac6_model_0.cif \
        --foxg1-chain A --hdac-chain B --lysines 15,181,237,241,272 \
        --cutoff 5.0

DEPENDENCIES
------------
    PyMOL importable as pymol2 (conda install -c conda-forge pymol-open-source)
"""

import argparse
import sys

try:
    import pymol2
except ImportError:
    sys.exit(
        "ERROR: PyMOL (pymol2) not importable. Install with:\n"
        "  conda install -c conda-forge pymol-open-source"
    )


def dist(a, b):
    return ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2) ** 0.5


def analyze(cif, foxg1_chain, hdac_chain, lysines, cutoff, nz_cutoff, salt_cutoff):
    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(cif, "m")

        # sanity check the chain has these lysines
        space = {"present": set()}
        cmd.iterate(f"m and chain {foxg1_chain} and resn LYS and name CA",
                    "present.add(int(resi))", space=space)
        missing = [k for k in lysines if k not in space["present"]]
        if missing:
            print(f"WARNING: chain {foxg1_chain} has no LYS at residue(s) "
                  f"{missing} -- check --foxg1-chain / --lysines / numbering.")

        # partner-chain atom coordinates + identity, once
        pcoords = {}
        space = {"d": pcoords}
        cmd.iterate_state(1, f"m and chain {hdac_chain} and not elem H",
                          "d[(chain, int(resi), resn, name, elem)] = (x, y, z)",
                          space=space)

        for k in lysines:
            print(f"\n=== FOXG1 K{k} ===")
            lys_atoms = {}
            space = {"d": lys_atoms}
            cmd.iterate_state(
                1, f"m and chain {foxg1_chain} and resi {k} and resn LYS "
                   f"and not elem H",
                "d[name] = (x, y, z)", space=space)
            if not lys_atoms:
                print("  (residue not found / not LYS -- skipped)")
                continue

            any_hits = {}   # (chain,resi,resn) -> min dist among any atoms
            nz_hits = {}    # (chain,resi,resn) -> [(nz_dist, partner_atom)]
            salt_hits = []  # (partner residue, atom, dist)

            nz_xyz = lys_atoms.get("NZ")

            for lname, lxyz in lys_atoms.items():
                for (pchain, presi, presn, pname, pelem), pxyz in pcoords.items():
                    d = dist(lxyz, pxyz)
                    if d <= cutoff:
                        key = (pchain, presi, presn)
                        if key not in any_hits or d < any_hits[key]:
                            any_hits[key] = d
                    if lname == "NZ" and d <= nz_cutoff:
                        key = (pchain, presi, presn)
                        nz_hits.setdefault(key, []).append((round(d, 2), pname))
                        if presn in ("ASP", "GLU") and pname in (
                            "OD1", "OD2", "OE1", "OE2") and d <= salt_cutoff:
                            salt_hits.append((f"{presn}{presi}", pname, round(d, 2)))

            if not any_hits:
                print(f"  No {hdac_chain}-chain heavy atoms within {cutoff} A "
                      f"of any atom in K{k}.")
                continue

            print(f"  Any-atom contacts within {cutoff} A:")
            for (pchain, presi, presn), d in sorted(any_hits.items(),
                                                     key=lambda x: x[1]):
                print(f"    {presn}{presi} (chain {pchain})  {d:.2f} A")

            if nz_xyz is None:
                print("  (NZ atom missing on this lysine -- side chain may be "
                      "truncated/unresolved; skipping NZ-specific report)")
            elif not nz_hits:
                print(f"  No contacts within {nz_cutoff} A of the NZ tip "
                      f"specifically -- side chain likely points away from "
                      f"the partner despite being at the interface.")
            else:
                print(f"  NZ-tip-specific contacts within {nz_cutoff} A "
                      f"(mechanistically relevant for HDAC access):")
                for (pchain, presi, presn), hits in sorted(
                        nz_hits.items(), key=lambda x: min(h[0] for h in x[1])):
                    best = min(hits, key=lambda h: h[0])
                    print(f"    {presn}{presi}:{best[1]}  {best[0]:.2f} A")

            if salt_hits:
                print("  Candidate salt bridge(s) (NZ near Asp/Glu carboxylate):")
                for resname, atom, d in sorted(salt_hits, key=lambda x: x[2]):
                    print(f"    K{k}:NZ -- {resname}:{atom}  {d:.2f} A")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cif", required=True)
    ap.add_argument("--foxg1-chain", default="A")
    ap.add_argument("--hdac-chain", default="B")
    ap.add_argument("--lysines", default="15,181,237,241,272",
                    help="comma-separated FOXG1 residue numbers")
    ap.add_argument("--cutoff", type=float, default=4.5,
                    help="any-heavy-atom contact distance, Angstrom")
    ap.add_argument("--nz-cutoff", type=float, default=5.0,
                    help="NZ-tip-specific contact distance, Angstrom "
                         "(slightly looser: catalytic tunnels are deep, "
                         "and this flags residues lining the approach too)")
    ap.add_argument("--salt-cutoff", type=float, default=4.0,
                    help="salt-bridge distance for NZ...Asp/Glu O, Angstrom")
    args = ap.parse_args()

    lysines = [int(x) for x in args.lysines.split(",")]
    analyze(args.cif, args.foxg1_chain, args.hdac_chain, lysines,
            args.cutoff, args.nz_cutoff, args.salt_cutoff)


if __name__ == "__main__":
    main()
