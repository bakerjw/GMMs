# Ground Motion Models

This repository provides Matlab scripts to perform calculations from various ground motion models. Citations for each
model are included within the scripts themselves as well as in the report cited below.

Matlab functions which run each model are included in the folder 'gmms'. They are named based on author's last initials,
year of model publication, and the event type to which they are applicable.

Matlab scripts which call the models to plot figures are in the folder 'testing'. These demonstrate how to set up the input
variable objects as well as how to utilize the outputs to create figures.

Figures created by the testing scripts are in the folder 'figures'.


Further documentation can be found in the following report:

Mongold and Baker... (2021). “XXX.” Blume Center Report # XXX. https://XXX

Users desiring additional information, or looking for other GMMs, are referred to John Douglas's excellent compilation at http://www.gmpe.org.uk/.

Note: Note these models were formerly referred to as "attenuation models," but the use of that name is being discouraged as it is misleading in what exactly these models provide. Users of these models should refer to them as "ground motion models." See http://daveboore.com/daves_notes/Thoughts%20on%20the%20acronyms%20GMPE,%20GMPM,%20and%20GMM.pdf 
