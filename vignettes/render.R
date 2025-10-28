library(bookdown)
preview_chapter("index.qmd")

preview_chapter("01-start.Rmd")
preview_chapter("02-fect.Rmd")
preview_chapter("03-plots.Rmd")
preview_chapter("04-gsynth.Rmd")
preview_chapter("05-panel.Rmd")
preview_chapter("06-sens.Rmd")
preview_chapter("aa-cheatsheet.Rmd")

# render the whole book
render_book("index.Rmd")