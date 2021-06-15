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

Mongold and Baker (2021). A software repository of Ground Motion Models.” Blume Earthquake Engineering Center Technical Report 207. https://purl.stanford.edu/qy113my5899

The following Ground Motion Models are currently provided:

| Function               | Region             | Predicted metric               | Reference                                            |
|------------------------|--------------------|--------------------------------|------------------------------------------------------|
| a_2015_active.m        | Active crustal     | Spectral   acceleration        | (Atkinson, 2015)                                     |
| as_1997_active.m       | Active crustal     | Spectral acceleration          | (Abrahamson and Silva, 1997)                         |
| as_2008_active.m       | Active crustal     | Spectral acceleration          | (Abrahamson and Silva, 2008)                         |
| ask_2014_active.m      | Active crustal     | Spectral acceleration          | (Abrahamson et al., 2014)                            |
| ba_2008_active.m       | Active crustal     | Spectral acceleration          | (Boore and Atkinson, 2008)                           |
| bjf_1997_active.m      | Active crustal     | Spectral acceleration          | (Boore et al., 1997)                                 |
| bssa_2014_active.m     | Active crustal     | Spectral acceleration          | (Boore et al., 2014)                                 |
| c_1997_active.m        | Active crustal     | Spectral acceleration          | (Campbell, 1997)                                     |
| cb_2008_active.m       | Active crustal     | Spectral acceleration          | (Campbell and Bozorgnia, 2008)                       |
| cb_2014_active.m       | Active crustal     | Spectral acceleration          | (Campbell and Bozorgnia, 2014)                       |
| cy_2008_active.m       | Active crustal     | Spectral acceleration          | (Chiou and Youngs, 2008)                             |
| cy_2014_active.m       | Active crustal     | Spectral acceleration          | (Chiou and Youngs, 2014)                             |
| i_2008_active.m        | Active crustal     | Spectral acceleration          | (Idriss, 2008)                                       |
| i_2014_active.m        | Active crustal     | Spectral acceleration          | (Idriss, 2014)                                       |
| scemy_1997_active.m    | Active crustal     | Spectral acceleration          | (Sadigh et al., 1997)                                |
| z_2006_active.m        | Active crustal     | Spectral acceleration          | (Zhao et al., 2006)                                  |
| as_1996_duration.m     |                    | Duration                       | (Abrahamson and Silva, 1996;   Stewart et al., 2002) |
| as_2016_duration.m     |                    | Duration                       | (Afshari and Stewart, 2016)                          |
| bsa_2009_duration.m    |                    | Duration                       | (Bommer et al., 2009)                                |
| ab_2006_stable.m       | Stable continental | Spectral acceleration          | (Atkinson and Boore, 2006)                           |
| sp_2016_stable.m       | Stable continental | Spectral acceleration          | (Shahjouei and Pezeshk, 2016)                        |
| ab_2003_subduction.m   | Subduction         | Spectral acceleration          | (Atkinson and Boore, 2003)                           |
| aga_2016_subduction.m  | Subduction         | Spectral acceleration          | (Abrahamson et al., 2016)                            |
| gswy_2002_subduction.m | Subduction         | Spectral acceleration          | (Gregor et al., 2002)                                |
| ycsh_1997_subduction.m | Subduction         | Spectral acceleration          | (Youngs et al., 1997)                                |
| as_1997_vert.m         | Active crustal     | Vertical spectral acceleration | (Abrahamson and Silva, 1997)                         |
| c_1997_vert.m          | Active crustal     | Vertical spectral acceleration | (Campbell, 1997)                                     |

The following correlation models are currently provided:

| Function                | Cross-period | Spatial | Reference                  |
|-------------------------|:------------:|:-------:|----------------------------|
| a_2011_corr.m           |       y      |         | (Al Atik, 2011)            |
| asa_2014_corr.m         |       y      |         | (Akkar et   al., 2014)     |
| bb_2017_corr.m          |       y      |         | (Baker and Bradley, 2017)  |
| bc_2006_corr.m          |       y      |         | (Baker and Cornell, 2006)  |
| bj_2008_corr.m          |       y      |         | (Baker and Jayaram, 2008)  |
| ga_2009_corr.m          |       y      |         | (Goda and Atkinson, 2009)  |
| gh_2008_spatial_corr.m  |              |    y    | (Goda and Hong, 2008)      |
| hm_2019_spatial_corr.m  |              |    y    | (Heresi and Miranda, 2019) |
| jb_2009_spatial_corr.m  |              |    y    | (Jayaram and Baker, 2009)  |
| lb_2013_spatial_corr.m  |       y      |    y    | (Loth and Baker, 2013)     |
| mcb_2018_spatial_corr.m |       y      |    y    | (Markhvida et   al., 2018) |

Users desiring additional information, or looking for other GMMs, may find the following resources useful:
http://www.gmpe.org.uk/
https://www.risksciences.ucla.edu/nhr3/gmtools/
https://github.com/gem/oq-engine/tree/master/openquake/hazardlib/gsim
https://github.com/arkottke/pygmm/

Note: Note these models were formerly referred to as "attenuation models," but the use of that name is being discouraged as it is misleading in what exactly these models provide. Users of these models should refer to them as "ground motion models." See http://daveboore.com/daves_notes/Thoughts%20on%20the%20acronyms%20GMPE,%20GMPM,%20and%20GMM.pdf
