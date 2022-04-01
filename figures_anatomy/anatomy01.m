%% show anatomical location of all units
% first show units & AP slices in 3D space
% generate flat images for each AP slice to make illutrator editing easier!

%% paths
addpath functions/
addpath functions/NIfTI_20140122/
addpath ../functions/
addpath ../preprocessing/functions/unit/ %% I think this should be the functions from figures_units

fname_sessions = '../info/recording_sessions.csv';

dir_data = '/media/StorageHDD/WallisLabData';
dir_bhv = 'bhv';
dir_pl2 = 'pl2_cut';
dir_info = 'info';
dir_out = '.';

dir_atlas = 'brainatlas';
fname_atlas = fullfile(dir_atlas,'atlas_rb.nii.gz');

which_value = 'pics_expval_bin'; % 'pics_expval' or 'pics_expval_bin';
dir_unitglms = ['singleunits/',which_value];
which_model = 'alltrials'; % 'alltrials' or 'free'


%% get session info
ses_info = readtable(fname_sessions,'Format','%s%s%s%s%d%d%d');

%% show units
show_region = 'CN'; % can be 'OFC' or 'CN'
anatomy_fig_setup;
[AP_slices,coords_disc,coords_size] = add_units(ses_info,dir_data,dir_info,show_region);

%% add atlas slices
[show_svgs, slice_images, corners_x, corners_z] = show_atlas(fname_atlas,show_region,AP_slices); 

%% make_illustrator templates
make_APslice_projections(AP_slices, coords_disc, coords_size,...
    show_svgs, slice_images, corners_x, corners_z, dir_out, show_region);

%% NEED TO CHECK THAT EVERYTHING BELOW WORKS WITH NEW FILE STRUCTURE: %%
%% for CN only: check out distribution of glm results by anatomical postion
% load glm outputs
glmdata = load_glmresults(ses_info,fullfile(dir_working,dir_data), ...
    dir_unitglms, which_model);

% show units
anatomy_fig_setup;
view([0 0]);
%%
close all
add_units_byglm2D(ses_info,dir_hdd,dir_info,'CN',glmdata,'DV','George');
%%
add_units_byglm3D(ses_info,dir_hdd,dir_info,'CN',glmdata);
%%
add_units_byglm1D(ses_info,dir_hdd,dir_info,'CN',glmdata,'ML');

