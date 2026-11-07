---
name: batch2-ssh-needs-no-vpn
description: batch2 is reachable by plain ssh to batch2.ncifcrf.gov when the Mac is on the institute network (no VPN needed there); off-site (home DNS) the name does not resolve; no bare "batch2" alias
metadata:
  type: reference
---
`ssh -o BatchMode=yes batch2.ncifcrf.gov '...'` works from the Mac with no VPN (user 2026-10-07: "you dont need the vpn
to run it... just ssh"). The bare alias `batch2` is not in ~/.ssh/config and fails with "Host key verification failed";
the git remote named batch2 points at ssh://atifao@batch2.ncifcrf.gov/home/atifao/SIMPLE.
Login banners and module messages precede output: print a marker (`echo ==`) and filter with `sed -n '/^==/,$p'`.

**Why:** I told the user a week of work had to wait for the VPN; it did not.
**How to apply:** when batch2 seems unreachable, retry with the full hostname before saying it is down.

2026-10-07 23:30: from a home network (Comcast DNS 75.75.75.75) batch2.ncifcrf.gov is NXDOMAIN - the internal name only
resolves on the institute network or VPN. Check `host batch2.ncifcrf.gov` before telling the user the cluster is down.
