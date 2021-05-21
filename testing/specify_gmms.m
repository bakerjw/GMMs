% Script to specify a list of GMMs to be evaluated


% Vector of indices correspoding with the NGA 2 West 2014 models
nga2west = [4 7 9 11 13];

% specify names of GMMs to consider
gmm_name{1}  = 'a_2015';
gmm_name{2}  = 'as_1997';
gmm_name{3}  = 'as_2008';
gmm_name{4}  = 'ask_2014';
gmm_name{5}  = 'ba_2008';
gmm_name{6}  = 'bjf_1997';
gmm_name{7}  = 'bssa_2014';
gmm_name{8}  = 'cb_2008';
gmm_name{9}  = 'cb_2014';
gmm_name{10} = 'cy_2008';
gmm_name{11} = 'cy_2014';
gmm_name{12} = 'i_2008';
gmm_name{13} = 'i_2014';
gmm_name{14} = 'scemy_1997';
gmm_name{15} = 'z_2006';


% specify line styles for plotting
line_style{1}  = ':';
line_style{2}  = '-.';
line_style{3}  = '--';
line_style{4}  = '-';
line_style{5}  = '--';
line_style{6}  = '-.';
line_style{7}  = '-';
line_style{8}  = '--';
line_style{9}  = '-';
line_style{10} = '--';
line_style{11} = '-';
line_style{12} = '--';
line_style{13} = '-';
line_style{14} = '-.';
line_style{15} = ':';
 

% specify line colors for plotting

% general colors to use for particular authors
ascolor =   [0,     0.4470, 0.7410]; 
bacolor =   [0.85,  0.3250, 0.0980]; 
cbcolor =   [0.929, 0.694,  0.1250]; 
cycolor =   [0.494, 0.184,  0.556]; 
icolor =    [0.466, 0.674,  0.188]; 
excolor =   [0.301, 0.745,  0.933]; 
othercolor =[0.635, 0.078,  0.184];


line_color{1}  = bacolor;
line_color{2}  = ascolor;
line_color{3}  = ascolor;
line_color{4}  = ascolor;
line_color{5}  = bacolor;
line_color{6}  = bacolor;
line_color{7}  = bacolor;
line_color{8}  = cbcolor;
line_color{9}  = cbcolor;
line_color{10} = cycolor;
line_color{11} = cycolor;
line_color{12} = icolor;
line_color{13} = icolor;
line_color{14} = othercolor;
line_color{15} = excolor;
