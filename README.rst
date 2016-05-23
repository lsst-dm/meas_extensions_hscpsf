======================
meas_extensions_hscpsf
======================

An experimental PSF determintion code written by Kendrick Smith in the context
of the HSC project. He writes:

  What's currently implemented on the master branch is a PCA-style PSF model,
  which was initially based on PcaPsf in the LSST pipeline but ended up
  deviating quite a bit.

  [...]

  Another thing which is nearly complete is a reimplementation of the psfex
  PSF model using low-level classes and libraries from the LSST/HSC pipeline.
  The idea is to have the functionality of psfex in a form which is familiar
  and easy to modify for anyone familiar with the HSC pipeline (the
  reimplementation will also much smaller by lines of code).

This code should be considered for incorporation into the LSST stack; see
`DM-6170 <https://jira.lsstcorp.org/browse/DM-6170>`_ for details.
