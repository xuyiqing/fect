## CI-only R profile. Activated via R_PROFILE_USER in GitHub Actions.
##
## pak/pkgdepends uses getOption("Ncpus") to decide worker count.
## Forcing single worker avoids rare cache/archive race conditions in CI.
options(Ncpus = 1)

# pkgcache (used by pak) download stability knobs
options(
  pkgcache_connecttimeout = 60,
  pkgcache_timeout = 7200,
  pkgcache_low_speed_limit = 1,
  pkgcache_low_speed_time = 60
)

# Allow workflow to override package type (e.g., force source on macOS).
pkg_type <- Sys.getenv("R_PKG_TYPE", unset = "")
if (nzchar(pkg_type)) {
  options(pkgType = pkg_type)
}

