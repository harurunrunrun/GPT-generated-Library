#!/usr/bin/env python3
"""
interactive_runner.py

Usage:
  python3 interactive_runner.py [--timeout SEC] [--show-io] -- <judge cmd...> -- <solution cmd...>

Example:
  python3 interactive_runner.py --timeout 2 --show-io -- ./interactor input.txt tout.txt -- ./solution

Notes:
- The judge's stdout is piped to the solution's stdin.
- The solution's stdout is piped to the judge's stdin.
- Both stderr streams are forwarded to this runner's stderr.
- --show-io prints the interactive transcript (both directions) to stderr.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import threading
import time
from typing import List, Optional


def split_by_double_dash(argv: List[str]) -> List[List[str]]:
    parts: List[List[str]] = []
    cur: List[str] = []
    for a in argv:
        if a == "--":
            parts.append(cur)
            cur = []
        else:
            cur.append(a)
    parts.append(cur)
    return parts


def pump_stream(
    src,
    dst,
    *,
    tee_prefix: Optional[bytes] = None,
    tee_to_stderr: bool = False,
    close_dst: bool = True,
):
    try:
        while True:
            data = src.read(4096)
            if not data:
                break
            if dst:
                try:
                    dst.write(data)
                    dst.flush()
                except BrokenPipeError:
                    # Destination closed; keep draining src to avoid blocking.
                    dst = None
            if tee_to_stderr and tee_prefix is not None:
                try:
                    sys.stderr.buffer.write(tee_prefix + data)
                    sys.stderr.buffer.flush()
                except Exception:
                    pass
    finally:
        try:
            src.close()
        except Exception:
            pass
        if close_dst and dst:
            try:
                dst.close()
            except Exception:
                pass


def forward_stderr(prefix: str, proc: subprocess.Popen):
    # Forward proc.stderr to runner stderr with a text prefix per line.
    try:
        for line in iter(proc.stderr.readline, b""):
            if not line:
                break
            sys.stderr.buffer.write(prefix.encode("utf-8") + line)
            sys.stderr.buffer.flush()
    except Exception:
        pass
    finally:
        try:
            proc.stderr.close()
        except Exception:
            pass


def terminate_process(p: subprocess.Popen):
    if p.poll() is not None:
        return
    try:
        p.terminate()
    except Exception:
        return
    # Wait a bit, then kill if needed
    try:
        p.wait(timeout=0.5)
    except Exception:
        try:
            p.kill()
        except Exception:
            pass


def main() -> int:
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--timeout", type=float, default=None, help="timeout in seconds (wall clock)")
    parser.add_argument("--show-io", action="store_true", help="print IO transcript to stderr")
    parser.add_argument("rest", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    # Expect: -- <judge...> -- <solution...>
    parts = split_by_double_dash(args.rest)
    # parts[0] may be empty if user wrote: python ... -- judge ... -- sol ...
    # We accept either:
    #   rest starts with "--" => parts[0]==[] , parts[1]==judge, parts[2]==solution
    # or:
    #   rest without leading "--" => parts[0]==judge, parts[1]==solution
    if len(parts) == 3 and parts[0] == []:
        judge_cmd = parts[1]
        sol_cmd = parts[2]
    elif len(parts) == 2:
        judge_cmd = parts[0]
        sol_cmd = parts[1]
    else:
        sys.stderr.write(
            "Invalid arguments.\n"
            "Expected: [--timeout SEC] [--show-io] -- <judge cmd...> -- <solution cmd...>\n"
        )
        return 2

    if not judge_cmd or not sol_cmd:
        sys.stderr.write("Judge command and solution command must be non-empty.\n")
        return 2

    # Start both processes
    try:
        judge = subprocess.Popen(
            judge_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,  # unbuffered
        )
    except FileNotFoundError:
        sys.stderr.write(f"Judge executable not found: {judge_cmd[0]}\n")
        return 2

    try:
        sol = subprocess.Popen(
            sol_cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0,
        )
    except FileNotFoundError:
        terminate_process(judge)
        sys.stderr.write(f"Solution executable not found: {sol_cmd[0]}\n")
        return 2

    show_io = bool(args.show_io)

    # Pump stdout<->stdin
    t1 = threading.Thread(
        target=pump_stream,
        args=(sol.stdout, judge.stdin),
        kwargs={
            "tee_prefix": b"[sol->judge] " if show_io else None,
            "tee_to_stderr": show_io,
            "close_dst": True,
        },
        daemon=True,
    )
    t2 = threading.Thread(
        target=pump_stream,
        args=(judge.stdout, sol.stdin),
        kwargs={
            "tee_prefix": b"[judge->sol] " if show_io else None,
            "tee_to_stderr": show_io,
            "close_dst": True,
        },
        daemon=True,
    )

    # Forward stderrs
    t3 = threading.Thread(target=forward_stderr, args=("[judge stderr] ", judge), daemon=True)
    t4 = threading.Thread(target=forward_stderr, args=("[sol stderr] ", sol), daemon=True)

    for t in (t1, t2, t3, t4):
        t.start()

    start = time.time()
    timed_out = False

    try:
        while True:
            j_rc = judge.poll()
            s_rc = sol.poll()
            if j_rc is not None and s_rc is not None:
                break
            if args.timeout is not None and (time.time() - start) > args.timeout:
                timed_out = True
                break
            time.sleep(0.01)
    except KeyboardInterrupt:
        sys.stderr.write("\nInterrupted.\n")
        terminate_process(judge)
        terminate_process(sol)
        return 130

    if timed_out:
        sys.stderr.write(f"\nTimed out after {args.timeout} seconds.\n")
        terminate_process(judge)
        terminate_process(sol)
        return 124

    # Ensure processes are terminated (in case one finished early and the other is stuck)
    terminate_process(judge)
    terminate_process(sol)

    # Join threads briefly
    for t in (t1, t2, t3, t4):
        t.join(timeout=0.2)

    j_rc = judge.returncode if judge.returncode is not None else -1
    s_rc = sol.returncode if sol.returncode is not None else -1

    sys.stderr.write(f"\n[judge exit] {j_rc}\n")
    sys.stderr.write(f"[sol exit] {s_rc}\n")

    # Conventional: return judge's code if it failed; else solution's code; else 0.
    if j_rc != 0:
        return j_rc if j_rc > 0 else 1
    if s_rc != 0:
        return s_rc if s_rc > 0 else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
