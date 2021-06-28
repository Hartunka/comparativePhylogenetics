#proj.dir = file.path('')
models.dir = file.path(proj.dir, 'R/comparativePhylogenetics/comparative-phylogenetics/models')
dataprep.dir = file.path(proj.dir, 'R/comparativePhylogenetics/comparative-phylogenetics/dataprep')
fit.dir = file.path(proj.dir, 'R/comparativePhylogenetics/comparative-phylogenetics/fit')
post.dir = file.path(proj.dir, 'R/comparativePhylogenetics/comparative-phylogenetics/posterior')
data.dir = file.path(proj.dir, 'data')

assertthat::assert_that(dir.exists(proj.dir),
                        dir.exists(fit.dir),
                        dir.exists(models.dir),
                        dir.exists(dataprep.dir),
                        dir.exists(data.dir),
                        dir.exists(post.dir)
                        )