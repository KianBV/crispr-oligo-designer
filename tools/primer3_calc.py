#!/usr/bin/env python3
"""Batch thermodynamic QC for oligos, using primer3-py.

Called by the R app (R/qc.R) as:
    python tools/primer3_calc.py <seqfile> <mv_mM> <dv_mM> <dntp_mM> <dna_nM> <temp_C>

<seqfile> has one oligo sequence (ACGT, uppercase) per line. Prints CSV to stdout:
    idx,tm,hp_dg,hp_tm,dimer_dg
where tm is the oligo melting temperature (primer3 calc_tm, deg C), hp_dg/hp_tm the
most stable hairpin (dG in kcal/mol, Tm in deg C), and dimer_dg the most stable
self/homo-dimer (dG in kcal/mol). A missing/unstable value is "NA".

If primer3 is not installed the import fails and the process exits non-zero, so the
R side transparently falls back to its built-in heuristic.
"""
import sys

try:
    import primer3
except Exception as exc:  # primer3-py not installed -> signal the caller
    sys.stderr.write("primer3 unavailable: %s\n" % exc)
    sys.exit(2)


def _fmt(x):
    try:
        if x != x:  # NaN
            return "NA"
    except Exception:
        return "NA"
    return "%.2f" % x


def main():
    seqfile = sys.argv[1]
    mv = float(sys.argv[2]); dv = float(sys.argv[3]); dntp = float(sys.argv[4])
    dna = float(sys.argv[5]); temp = float(sys.argv[6])

    with open(seqfile) as fh:
        seqs = [ln.strip().upper() for ln in fh]

    kw = dict(mv_conc=mv, dv_conc=dv, dntp_conc=dntp, dna_conc=dna)
    out = ["idx,tm,hp_dg,hp_tm,dimer_dg"]
    for i, s in enumerate(seqs):
        if not s or any(c not in "ACGT" for c in s):
            out.append("%d,NA,NA,NA,NA" % i); continue
        try:
            tm = primer3.calc_tm(s, **kw)
        except Exception:
            tm = float("nan")
        try:
            hp = primer3.calc_hairpin(s, temp_c=temp, **kw)
            hp_dg, hp_tm = hp.dg / 1000.0, hp.tm
        except Exception:
            hp_dg, hp_tm = 0.0, float("nan")
        try:
            dm = primer3.calc_homodimer(s, temp_c=temp, **kw)
            dm_dg = dm.dg / 1000.0
        except Exception:
            dm_dg = 0.0
        out.append("%d,%s,%s,%s,%s" % (i, _fmt(tm), _fmt(hp_dg), _fmt(hp_tm), _fmt(dm_dg)))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
