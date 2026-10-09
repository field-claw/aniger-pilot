"""Upload a directory's tracked files to a GitHub repo via the contents API.

Why this exists (measured, 2026-10-09): the git protocol (github.com:443) was
unreachable from this machine -- "Empty reply from server", "Connection was
reset", "schannel: server closed abruptly" -- while the REST API worked. The
obvious replacement, POST /git/trees with a base64 `content` per entry, silently
DOUBLE-ENCODED the payload: `gh api --input file.json` base64-encodes the string
fields it sends, so already-base64 content arrived base64-of-base64 and every
Python file landed as unreadable text (SyntaxError on first line).

This script therefore uses the plain contents API, which takes the content
already base64-encoded via `-f content=...`, and it VERIFIES the round trip by
decoding what GitHub stored and comparing bytes with the local file. A mismatch
is a hard failure -- an upload that "succeeded" while shipping undecodable
bytes is exactly the failure mode this script exists to prevent.

Usage:
    python upload_via_contents_api.py <repo> <dir> [--dry-run]
"""
import base64
import hashlib
import json
import os
import subprocess
import sys


def gh(args, stdin=None):
    return subprocess.run(["gh", "api"] + args, input=stdin,
                          capture_output=True, text=True)


def remote_sha(repo, path):
    r = gh(["repos/%s/contents/%s" % (repo, path), "--jq", ".sha"])
    return r.stdout.strip() if r.returncode == 0 else None


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    repo, workdir = sys.argv[1], sys.argv[2]
    dry = "--dry-run" in sys.argv

    files = subprocess.run(["git", "-C", workdir, "ls-files"],
                           capture_output=True, text=True).stdout.split()
    if not files:
        print("no tracked files in %s" % workdir)
        return 1

    bad = []
    for f in files:
        local = open(os.path.join(workdir, f), "rb").read()
        b64 = base64.b64encode(local).decode()
        if dry:
            print("would upload %-52s %8d B" % (f, len(local)))
            continue

        sha = remote_sha(repo, f)
        args = ["--method", "PUT", "repos/%s/contents/%s" % (repo, f),
                "-f", "message=upload %s" % f, "-f", "content=" + b64]
        if sha:
            args += ["-f", "sha=" + sha]
        r = gh(args)
        if r.returncode != 0:
            print("FAIL %-50s %s" % (f, r.stderr.strip()[:90]))
            bad.append(f)
            continue

        # verify the round trip: what GitHub stored must decode to our bytes
        got = gh(["repos/%s/contents/%s" % (repo, f), "--jq", ".content"])
        stored = base64.b64decode("".join(got.stdout.split()))
        if hashlib.md5(stored).hexdigest() != hashlib.md5(local).hexdigest():
            print("MISMATCH %-46s local=%d stored=%d"
                  % (f, len(local), len(stored)))
            bad.append(f)
        else:
            print("ok   %-52s %8d B" % (f, len(local)))

    if bad:
        print("\n%d file(s) failed: %s" % (len(bad), ", ".join(bad)))
        return 1
    print("\nall %d files uploaded and byte-verified" % len(files))
    return 0


if __name__ == "__main__":
    sys.exit(main())