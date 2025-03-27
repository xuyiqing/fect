library(bookdown)
preview_chapter("index.qmd")

preview_chapter("01-Start.Rmd")
preview_chapter("02-FEct.Rmd")
preview_chapter("03-Gsynth.Rmd")
preview_chapter("04-Table.Rmd")
preview_chapter("05-references.Rmd")

# render the whole book
render_book("index.Rmd")