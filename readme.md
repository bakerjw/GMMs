# Ground Motion Models

This repository provides Matlab scripts to perform calculations from various ground motion models and correlation models.
Citations for each model are included within the scripts themselves as well as in the report cited below.

Matlab functions which run each ground motion model are included in the folder 'gmms'.
Matlab functions which run each correlation model are included in the folder 'correlations'.
They are named based on author's last initials, year of model publication, and the situation for which they are applicable.

Matlab scripts which call the models to plot figures are in the folder 'testing'.
These demonstrate how to set up the input variable objects as well as how to utilize the outputs to create figures.

Figures created by the testing scripts are in the folder 'figures'.


Further documentation can be found in the following report:

Mongold and Baker (2021). A software repository of Ground Motion Models.” Blume Center Report # 207. https://blume.stanford.edu/reports

Users desiring additional information, or looking for other GMMs, are referred to John Douglas's compilation at http://www.gmpe.org.uk/. Python implimentations of ground motion models are available at https://github.com/gem/oq-engine/tree/master/openquake/hazardlib/gsim and https://github.com/arkottke/pygmm/.

Note: Note these models were formerly referred to as "attenuation models," but the use of that name is being discouraged as it is misleading in what exactly these models provide. Users of these models should refer to them as "ground motion models." See http://daveboore.com/daves_notes/Thoughts%20on%20the%20acronyms%20GMPE,%20GMPM,%20and%20GMM.pdf
