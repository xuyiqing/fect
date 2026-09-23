## Package-local session state.
##
## Used for once-per-session notices instead of options(): CRAN's policy asks
## packages not to leave changes to the user's options() behind, and an
## environment inside the namespace is invisible to the user and reset on
## reload. Add a named flag here rather than a new options() key.
.fect_state <- new.env(parent = emptyenv())
