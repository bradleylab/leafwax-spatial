# Chordal-fit provenance

This directory preserves the exact calibration input and software record for
the 2026-07-28 chordal-distance calibration fit reported in the associated
manuscript.

`leafwax_d2h_c29_calibration_fit.csv` is the byte-for-byte fit input (MD5
`bb52649130a02b9b7d897325899b4f5a`). The current public calibration file was
subsequently corrected in metadata only: 27 Río Bermejo records received their
deposit-level source and DOI, and six corrupted location labels were repaired.
The coordinates, isotope measurements, archive classes, and every modeled
value are identical between the fit input and the public calibration file.

Run `Rscript scripts/audit_calibration_metadata_update.R` from the repository
root to verify the field-by-field differences. The command writes
`scripts/reference_outputs/calibration_metadata_update_audit.csv` and fails if
any numerical or modeling field differs.

`fit_environment.json` is generated from the archived fit-time record by
`Rscript scripts/export_fit_provenance.R`. It contains the R, CmdStan, and
package versions; the container checksum; code and input checksums; and hashes
for all 17 prepared model configurations. Cluster-specific paths and operator
metadata are intentionally omitted.
