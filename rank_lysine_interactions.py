#!/usr/bin/env python3
"""
rank_lysine_interactions.py
=============================
Ranks a set of lysines on one chain (default: FOXG1's 5 conserved interface
lysines) by structural "interaction strength" against a partner chain
(HDAC1 / HDAC2 / HDAC6), using one AF3 model_0.cif per complex.

IMPORTANT FRAMING
------------------
AlphaFold3 predicts a static structure, not a binding energy. There is no
"strength of binding" number in the AF3 output. What this script reports
instead are four *structural proxies* that structural biologists commonly
use to compare how tightly different interface residues engage a partner:

  1. Buried surface area (BSA, sq. Angstrom) -- how much of the lysine's
     surface becomes inaccessible upon complex formation. This is the
     standard, most defensible single metric for "how engaged is this
     residue in the interface" and is what the ranking is primarily
     sorted by.
  2. Contact atom count -- how many partner-chain heavy atoms sit within
     the contact cutoff of the lysine (packing density).
  3. Closest heavy-atom distance -- tightest point of contact, any atom.
  4. NZ-tip distance -- closest partner atom to the reactive epsilon-amine
     specifically (mechanistically relevant for an HDAC-substrate story).
  5. Salt bridge flag -- NZ near an Asp/Glu carboxylate O, a real
     electrostatic interaction and typically the "strongest" single
     contact type at a protein interface.

None of this is a substitute for an actual computed or measured binding
affinity (MM-GBSA rescoring, FEP, ITC, SPR, etc.) if the researchers need a
publication-grade strength claim -- this is a same-model structural ranking,
useful for prioritizing which lysines to test first.

USAGE
-----
    python rank_lysine_interactions.py \
        --cif fold_foxg1_hdac1_model_0.cif \
        --foxg1-chain A --hdac-chain B \
        --lysines 15,181,237,241,272 \
        --label HDAC1

Run once per complex (HDAC1 / HDAC2 / HDAC6) and compare the three ranked
tables side by side.

DEPENDENCIES
------------
    PyMOL importable as pymol2 (conda install -c conda-forge pymol-open-source)
"""

import argparse
import csv
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


def per_residue_bsa(cmd, foxg1_chain, hdac_chain, lysines):
    """Total per-residue buried surface area (isolated - complex SASA),
    summed over all atoms of each target lysine. Same physical quantity
    dASA/interfaceResidues checks per-atom, but totaled per residue here
    for a continuous, rankable number instead of a pass/fail flag."""
    cmd.create("tempC", f"m and polymer and (chain {foxg1_chain} or chain {hdac_chain})")
    cmd.get_area("tempC", load_b=1)
    cmd.alter("tempC", "q=b")
    cmd.extract("chA", f"tempC and chain {foxg1_chain}")
    cmd.extract("chB", f"tempC and chain {hdac_chain}")
    cmd.get_area("chA", load_b=1)
    cmd.get_area("chB", load_b=1)
    cmd.alter("chA or chB", "b=b-q")

    space = {"rows": []}
    cmd.iterate("chA or chB", "rows.append((chain, int(resi), resn, b))",
                space=space)

    totals = {}
    for chain, resi, resn, b in space["rows"]:
        if chain == foxg1_chain and resi in lysines:
            totals[resi] = totals.get(resi, 0.0) + b

    for tmp in ("tempC", "chA", "chB"):
        cmd.delete(tmp)
    return totals


def contact_metrics(cmd, foxg1_chain, hdac_chain, lysines, cutoff, salt_cutoff):
    pcoords = {}
    space = {"d": pcoords}
    cmd.iterate_state(1, f"m and chain {hdac_chain} and not elem H",
                      "d[(chain, int(resi), resn, name)] = (x, y, z)",
                      space=space)

    out = {}
    for k in lysines:
        lys_atoms = {}
        space = {"d": lys_atoms}
        cmd.iterate_state(
            1, f"m and chain {foxg1_chain} and resi {k} and resn LYS "
               f"and not elem H",
            "d[name] = (x, y, z)", space=space)

        n_contacts = 0
        min_any = None
        min_nz = None
        min_nz_partner = None
        salt_partner = None
        salt_dist = None

        for lname, lxyz in lys_atoms.items():
            for (pchain, presi, presn, pname), pxyz in pcoords.items():
                d = dist(lxyz, pxyz)
                if d <= cutoff:
                    n_contacts += 1
                    if min_any is None or d < min_any:
                        min_any = d
                if lname == "NZ":
                    if d <= cutoff and (min_nz is None or d < min_nz):
                        min_nz = d
                        min_nz_partner = f"{presn}{presi}:{pname}"
                    if (presn in ("ASP", "GLU")
                            and pname in ("OD1", "OD2", "OE1", "OE2")
                            and d <= salt_cutoff
                            and (salt_dist is None or d < salt_dist)):
                        salt_dist = d
                        salt_partner = f"{presn}{presi}"

        out[k] = {
            "n_contacts": n_contacts,
            "min_any": min_any,
            "min_nz": min_nz,
            "min_nz_partner": min_nz_partner,
            "salt_partner": salt_partner,
            "salt_dist": salt_dist,
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cif", required=True)
    ap.add_argument("--foxg1-chain", default="A")
    ap.add_argument("--hdac-chain", default="B")
    ap.add_argument("--lysines", default="15,181,237,241,272")
    ap.add_argument("--cutoff", type=float, default=4.5,
                    help="heavy-atom contact distance, Angstrom")
    ap.add_argument("--salt-cutoff", type=float, default=4.0)
    ap.add_argument("--label", default="", help="e.g. HDAC1 -- cosmetic, for CSV name")
    ap.add_argument("--csv", default=None, help="optional output CSV path")
    args = ap.parse_args()

    lysines = [int(x) for x in args.lysines.split(",")]

    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(args.cif, "m")

        bsa = per_residue_bsa(cmd, args.foxg1_chain, args.hdac_chain, lysines)
        contacts = contact_metrics(cmd, args.foxg1_chain, args.hdac_chain,
                                   lysines, args.cutoff, args.salt_cutoff)

    rows = []
    for k in lysines:
        c = contacts[k]
        rows.append({
            "lysine": f"K{k}",
            "buried_area_A2": round(bsa.get(k, 0.0), 1),
            "n_contact_atoms": c["n_contacts"],
            "min_any_dist": round(c["min_any"], 2) if c["min_any"] else None,
            "min_nz_dist": round(c["min_nz"], 2) if c["min_nz"] else None,
            "nz_partner": c["min_nz_partner"] or "-",
            "salt_bridge": f"{c['salt_partner']} ({c['salt_dist']:.2f} A)"
                          if c["salt_partner"] else "no",
        })

    # primary sort: buried area desc. ties: more contacts, then closer min dist.
    rows.sort(key=lambda r: (
        -r["buried_area_A2"],
        -r["n_contact_atoms"],
        r["min_any_dist"] if r["min_any_dist"] is not None else 999,
    ))

    label = args.label or "complex"
    print(f"\n=== Ranked lysine interaction strength vs {label} ===")
    print(f"(primary: buried surface area; ties broken by contact count, "
          f"then closest distance)\n")
    header = f"{'rank':<5}{'lysine':<8}{'BSA (A2)':<10}{'#contacts':<11}" \
            f"{'min dist':<10}{'NZ dist':<10}{'NZ partner':<16}{'salt bridge'}"
    print(header)
    for i, r in enumerate(rows, 1):
        print(f"{i:<5}{r['lysine']:<8}{r['buried_area_A2']:<10}"
              f"{r['n_contact_atoms']:<11}"
              f"{str(r['min_any_dist']):<10}{str(r['min_nz_dist']):<10}"
              f"{r['nz_partner']:<16}{r['salt_bridge']}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["rank"] + list(rows[0].keys()))
            w.writeheader()
            for i, r in enumerate(rows, 1):
                w.writerow({"rank": i, **r})
        print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
