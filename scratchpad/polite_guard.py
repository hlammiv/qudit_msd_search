"""Pause our batch jobs while other people need the machine; resume when they are done.

python3 scratchpad/polite_guard.py SID[,SID...] OTHERS_MAX [RESUME_BELOW] [PERIOD_S]
  SID: session ids of our queue scripts (they are started with setsid; every process in those sessions is ours).
  Every PERIOD_S (default 60) the CPU used by all *other* processes is measured over 10 s from /proc.  If it exceeds
  OTHERS_MAX cores on two consecutive checks, every process of our sessions gets SIGSTOP; when it falls below
  RESUME_BELOW (default OTHERS_MAX / 2) on two consecutive checks they get SIGCONT.  Only explicit PIDs read from
  /proc are signalled.  Exits (resuming everything) when no process of the sessions is left.  Log on stdout.
"""
import os, signal, sys, time

HZ = os.sysconf("SC_CLK_TCK")


def snapshot():
    out = {}
    for p in os.listdir("/proc"):
        if not p.isdigit():
            continue
        try:
            f = open(f"/proc/{p}/stat").read()
            rest = f[f.rindex(")") + 2:].split()
            out[int(p)] = (int(rest[3]), int(rest[11]) + int(rest[12]))      # session id, utime + stime
        except (OSError, ValueError, IndexError):
            pass
    return out


def main():
    sids = {int(x) for x in sys.argv[1].split(",")}
    hi = float(sys.argv[2]); lo = float(sys.argv[3]) if len(sys.argv) > 3 else hi / 2
    period = float(sys.argv[4]) if len(sys.argv) > 4 else 60
    me = os.getpid(); paused = False; strikes = 0
    while True:
        a = snapshot(); time.sleep(10); b = snapshot()
        ours = [p for p, (s, _) in b.items() if s in sids and p != me]
        if not ours:
            print(time.strftime("%F %T"), "no jobs left, exiting", flush=True)
            return
        others = sum(t - a[p][1] for p, (s, t) in b.items() if s not in sids and p in a and p != me) / HZ / 10
        want = (others > hi) if not paused else not (others < lo)
        strikes = strikes + 1 if want != paused else 0
        if strikes >= 2:
            sig = signal.SIGCONT if paused else signal.SIGSTOP
            for p in ours:
                try:
                    os.kill(p, sig)
                except ProcessLookupError:
                    pass
            paused = not paused; strikes = 0
            print(time.strftime("%F %T"), "others %.1f cores ->" % others, "PAUSED" if paused else "RESUMED", len(ours), "procs", flush=True)
        time.sleep(period)


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *a: sys.exit(0))       # kill -TERM still resumes the jobs (finally below)
    try:
        main()
    finally:                                         # never leave our jobs stopped
        for p, (s, _) in snapshot().items():
            if s in {int(x) for x in sys.argv[1].split(",")}:
                try:
                    os.kill(p, signal.SIGCONT)
                except ProcessLookupError:
                    pass
