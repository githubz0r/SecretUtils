# SecretUtils
Some utility functions for analyzing data with Pagoda2/Conos and simulating data with Splatter.

[Example analysis script](https://githubz0r.github.io/SecretUtils/docs/mouse_alzheimer_analysis.html) using subset of data from Keren-Shaul et al https://www.cell.com/cell/fulltext/S0092-8674(17)30578-0 i.e. some mice with alzheimer. Annotation was done based on markers in the paper but coarser (fewer clusters).

[Example simulation script](https://githubz0r.github.io/SecretUtils/docs/simulations_plots_intro.html). Note that Conos is used by these functions to integrate and normalize etc. the simulated datasets. Estimation of initial distribution parameters can be done from any count matrix. The plot grid is not displayed in the html since it often looks like garbage when displayed inside the script, however [here's a pdf showing a result for all the distances and factors I varied.](https://github.com/githubz0r/SecretUtils/blob/master/pdfs/grid_of_grids_final.pdf)


A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
