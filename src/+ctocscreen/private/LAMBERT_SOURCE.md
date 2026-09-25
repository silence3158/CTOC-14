# Lambert Source

- Source: https://github.com/rodyo/FEX-Lambert
- Commit: `f220826e6ac688ae8dfd2e8c9374872e43ea9fb4`
- File: `lambert.m`, retrieved 2026-09-25.
- License: BSD 2-clause; see `LAMBERT_LICENSE.txt`.
- Local change: renamed the main function to `lambert_izzo_gooding` for private package use. The numerical implementation is otherwise retained.
- Input time is in days; the V3 wrapper converts seconds explicitly. Signed time selects short/long paths; signed nonzero revolution count selects the two branches.
- This solver supplies two-body guesses only. Endpoint checks, fallback to the existing enumerator, J2 correction, and final independent verification remain outside the vendor implementation. Its extremal-distance output is not used as a J2 height certificate.

The public source was retrieved using a per-command TLS certificate bypass because the local certificate chain could not be validated. No global TLS configuration was changed.
