#!/usr/bin/env python3
"""Bitwise comparison of the EL_* diagnostic streams of two q2p1_el_pipeflow logs.
usage: el_stream_compare.py [--strict] <reference.log> <candidate.log>
       el_stream_compare.py --self-test

Per tag the lines are compared as text (bytes after the tag), not as parsed numbers:
  bitwise    every line identical, same line count
  trailing   every reference line is reproduced verbatim and the candidate appends labelled
             tokens ('key= value ...', a newer binary writing more fields); accepted unless
             --strict, and always named in the summary line
  DIFF       anything else: different line count, a missing, extra unlabelled, NaN/Inf or
             changed token (numeric or text); the worst relative numeric difference over the
             aligned tokens is printed for diagnosis
Streams absent from the reference are reported, not compared. Exit status 0 only when every
common stream is bitwise (or trailing without --strict), 1 otherwise, 2 on usage errors."""
import re, sys, tempfile, os

NUM = re.compile(r"^[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?$")
BAD = re.compile(r"nan|inf|\*{2,}", re.IGNORECASE)
TAGS = ("EL_MOMENTUM_ELEMINT", "EL_MEAN_SLIP t=", "EL_MEAN_SLIP_TAVG", "EL_NEWTON_PAIR", "EL_VOLUME_CONSERVATION",
        "EL_FORCE_BUDGET", "EL_FEEDBACK_CONSERVATION", "EL_TERMINAL_VEL", "EL_CONTACT_STATS", "EL_VA_SOURCE",
        "EL_SUSP_STRESS", "EL_CONT_STRESS", "EL_HALO_RECORDS", "EL_LUB_IMPULSE", "EL_VA_MOMENT",
        "EL_VA_TSRC", "EL_VA_RESIDUAL", "EL_VA_COMPAT", "EL_FLUID_SOURCE", "EL_VA_FLUX t=", "EL_VA_FLUXMOM")


def load(path):
    d = {t: [] for t in TAGS}
    with open(path, errors="replace") as f:
        for line in f:
            for t in TAGS:
                if line.startswith(t):
                    d[t].append(line[len(t):].rstrip("\r\n"))
                    break
    return d


def tofloat(tok):
    return float(tok.replace("D", "E").replace("d", "e"))


def compare_lines(ra, rb):
    """One reference line vs one candidate line -> (status, note, worst, example)."""
    if ra == rb:
        return "bitwise", "", 0.0, ""
    ta, tb = ra.split(), rb.split()
    if BAD.search(rb):
        return "DIFF", "NaN/Inf/overflow token in candidate", 0.0, ""
    if len(tb) > len(ta) and tb[:len(ta)] == ta and tb[len(ta)].endswith("="):
        return "trailing", f"+{len(tb) - len(ta)} trailing tokens", 0.0, ""
    if len(tb) < len(ta):
        return "DIFF", f"{len(ta) - len(tb)} tokens missing", 0.0, ""
    if len(tb) > len(ta):
        return "DIFF", f"{len(tb) - len(ta)} unlabelled extra tokens", 0.0, ""
    worst, ex = 0.0, ""
    for x, y in zip(ta, tb):
        if x == y:
            continue
        if NUM.match(x) and NUM.match(y):
            fx, fy = tofloat(x), tofloat(y)
            if fx == fy:
                return "DIFF", f"same value, different text: {x} vs {y}", 0.0, ""
            r = abs(fx - fy) / max(abs(fx), abs(fy), 1e-300)
            if r > worst:
                worst, ex = r, f"{x}->{y}"
        else:
            return "DIFF", f"text token changed: {x} vs {y}", 0.0, ""
    return "DIFF", "numeric difference", worst, ex


def compare(ref, cand, strict=False, out=print):
    a, b = load(ref), load(cand)
    allok, trailing_streams = True, 0
    for t in TAGS:
        la, lb = a[t], b[t]
        if not la and not lb:
            continue
        if not la:
            out(f"  {t:24s} lines {len(la)}/{len(lb)}  new in candidate, not compared")
            continue
        status, note, worst, ex = "bitwise", "", 0.0, ""
        if len(la) != len(lb):
            status, note = "DIFF", "line count differs"
        for ra, rb in zip(la, lb):
            s, n, w, e = compare_lines(ra, rb)
            if s == "bitwise":
                continue
            if s == "trailing":
                if status == "bitwise":
                    status, note = "trailing", n
                continue
            status = "DIFF"
            if w > worst:
                worst, ex = w, e
            if not note or note.startswith("+"):
                note = n
        if status == "trailing" and strict:
            status = "DIFF"
        if status == "trailing":
            trailing_streams += 1
        if status == "DIFF":
            allok = False
        detail = f"worst rel diff {worst:.2e} ({ex})" if ex else note
        out(f"  {t:24s} lines {len(la)}/{len(lb)}  {detail:48s} {status}")
    if allok and trailing_streams:
        out(f"ALL COMMON STREAMS BITWISE ({trailing_streams} stream(s) with trailing tokens added in the candidate)")
    elif allok:
        out("ALL COMMON STREAMS BITWISE")
    else:
        out("NOT BITWISE")
    return 0 if allok else 1


def self_test():
    ref = ("EL_MEAN_SLIP t=   1.000000E-03 step= 1 eps=  9.9E-01 uf_super=  1.0E+00  2.0E+00  3.0E+00\n"
           "EL_FLUID_SOURCE phase= pre step= 1 mode= explicit applied= yes rhs=  1.0E+00\n"
           "EL_VA_SOURCE t=   1.000000E-03 step= 1 gbar=  1.5E-04\n")
    cases = [  # (candidate text, strict, expected exit, expected summary prefix)
        (ref, False, 0, "ALL COMMON STREAMS BITWISE"),
        (ref.replace("3.0E+00", "NaN"), False, 1, "NOT BITWISE"),
        (ref.replace("applied= yes", "applied= no"), False, 1, "NOT BITWISE"),
        (ref.replace("  3.0E+00", ""), False, 1, "NOT BITWISE"),
        (ref.replace("3.0E+00\n", "3.0E+00 uz_mean=  2.5E+00\n"), False, 0, "ALL COMMON STREAMS BITWISE (1 stream"),
        (ref.replace("3.0E+00\n", "3.0E+00 uz_mean=  2.5E+00\n"), True, 1, "NOT BITWISE"),
        (ref.replace("3.0E+00\n", "3.0E+00  2.5E+00\n"), False, 1, "NOT BITWISE"),
        (ref.replace("1.5E-04", "1.50E-04"), False, 1, "NOT BITWISE"),
        (ref.replace("1.5E-04", "1.6E-04"), False, 1, "NOT BITWISE"),
        (ref + "EL_VA_SOURCE t=   2.000000E-03 step= 2 gbar=  1.5E-04\n", False, 1, "NOT BITWISE"),
        (ref + "EL_VA_FLUX t=   1.000000E-03 step= 1 S=  0.0E+00\n", False, 0, "ALL COMMON STREAMS BITWISE"),
    ]
    with tempfile.TemporaryDirectory() as d:
        rp = os.path.join(d, "ref.log")
        open(rp, "w").write(ref)
        failed = 0
        for i, (cand, strict, code, summary) in enumerate(cases):
            cp = os.path.join(d, f"cand{i}.log")
            open(cp, "w").write(cand)
            lines = []
            rc = compare(rp, cp, strict=strict, out=lines.append)
            ok = rc == code and lines[-1].startswith(summary)
            failed += not ok
            print(f"  case {i:2d}: {'ok  ' if ok else 'FAIL'} exit {rc} (want {code}) {lines[-1]}")
        print("self-test passed" if not failed else f"self-test FAILED ({failed} case(s))")
        return 0 if not failed else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["--self-test"]:
        sys.exit(self_test())
    strict = "--strict" in args
    args = [x for x in args if x != "--strict"]
    if len(args) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(compare(args[0], args[1], strict=strict))
