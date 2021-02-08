[inputs]
method      | "apPhoto"
# method: "apPhoto' or "mostLikely"
plotout     | "PDF"

[wcs]
# MAP PLOTTING PARAMETERS FOR SMALER MAPS IUSED FOR APERTURE PHOTOMETRY AND
# SOURCE FITTING
crpix       | 120           120
cdelt       | -.2           .2
ctype       | GLON-CYP      GLAT-CYP

[sources]
# LIST OF SOURCES TO ANALYSE
{G21.5-0.9}
name        | "G21.5-0.9"
coords      | 021.500642    -00.885663

[maps]

{SOME GENERIC TERM}
name      | "NAME OF SURVEY"
file      | "DIRECTORY"
type      | 'healpix' or 'wcs'
freqGHz   | CENTRAL FREQUENCY IN GHz
