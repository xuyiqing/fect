## CI-only R profile. Activated via R_PROFILE_USER in GitHub Actions.
##
## pak/pkgdepends uses getOption("Ncpus") to decide worker count.
## Forcing single worker avoids rare cache/archive race conditions in CI.
options(Ncpus = 1)

