# AGENTS.md

## Ecosystem orientation

For cross-repository work, read the optional local
[agent overview](../gfactory/docs/agent-overview.md), then only the relevant
owner contract and current source. It gives ownership, reuse links and current
priorities; this repository retains its own interfaces and build instructions.
If private gfactory is absent, continue with local guidance; do not fetch it or
make builds, tests, runtime, public evidence or website rendering depend on it.
A cross-repository link does not authorize edits in another repository.

## Local responsibility

gtools owns the reviewed public documentation snapshot and separately approved
release artifacts. It does not own gsuite analysis source or numerical methods.
Read [README.md](README.md) for the publishing contract.

## Website changes and checks

Edit gsuite website content in its owning development checkout. When that
checkout is available, run `python tools/build/prepare_public_site.py --sync-to ../gtools`
from gsuite; it validates the allowlist and links before synchronizing. Do not
hand-edit generated site HTML or copy the private gsuite `_site` tree.
The updater manages `site/`, README.md and the Pages workflow, and preserves
this AGENTS.md. If the source checkout is unavailable, review the snapshot
and report needed source corrections; do not invent a render/build dependency.
Use `git diff --check` for local edits. No native build is needed.

Preserve unrelated changes. Do not commit or push unless requested. A push to
main affecting site/ or the Pages workflow deploys the website; it does not
release software. Never copy private source, reports or archive content into site/.
