%{
This script uses MTEX to arrange the existing grain orientations such that the misorientation
distribution function (MDF) are statistically equivalent to the MDF found in the EBSD data.
Users will want download MTEX (https://mtex-toolbox.github.io/download) and to change
"path_mtex" to wherever it is unzipped in their system.
Additionally, if using a different material from Aluminum, the user may need to change
"crystal_symmetry".
If the algorithm has been running for quite a while without improvement, then the desired
KS statistic cannot be reached for some reason. This can be caused by an insufficient number
of grains or a bad ODF that does not allow for the desired MDF. It can also just be caused
by an overly ambitious ks_target. just wait until after the next full save and exit the
swapping algorithm to diagnose further.
%}

clear all;clc

%%% start mtex
path_mtex = "../../MTEX";
addpath(path_mtex)
startup_mtex
ppi = get(0,'ScreenPixelsPerInch');
fontSize = round(15 * ppi/100);
setMTEXpref('FontSize',fontSize);

%%% vars
ks_target        = 0.04;
iters_per_epoch  = 100;
iters_per_backup = 2000;
crystal_symmetry = crystalSymmetry('m-3m', [4.050 4.050 4.050], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]);
n_bins = 50;
save_histograms  = true;

%%% paths
% files
path_ebsd_file_input           = "./pipeline_output/3-feature_attributes_ebsd.dream3d";
path_synthetic_file_input      = "./pipeline_output/15-synthetic_grains.dream3d";
path_temp_file_eulerangles     = "./pipeline_output/16-eulerangles.txt";
path_temp_file_misorientations = "./pipeline_output/16-misorientations.txt";
path_fig_output                = "./pipeline_output/images/16-misorientations.png";
% groups
path_ebsd_hdf5_cellfeaturedata      = "/DataContainers/ImageDataContainer/CellFeatureData";
path_synthetic_hdf5_cellfeaturedata = "/DataContainers/ImageDataContainer/CellFeatureData";
% datasets
path_ebsd_hdf5_misorientationlist = path_ebsd_hdf5_cellfeaturedata     +"/"+"MisorientationList";
path_synthetic_hdf5_neighborlist  = path_synthetic_hdf5_cellfeaturedata+"/"+"NeighborList";
path_synthetic_hdf5_eulerangles   = path_synthetic_hdf5_cellfeaturedata+"/"+"AvgEulerAngles";

%%% find linked datasets
path_synthetic_hdf5_numneighbors = ...
    rsplit(path_synthetic_hdf5_neighborlist,"/")+"/"+ ...
        h5readatt(path_synthetic_file_input,path_synthetic_hdf5_neighborlist,"Linked NumNeighbors Dataset" ...
    );

% make histogram images folder
if save_histograms
    mkdir(rsplit(path_fig_output,'/'))
end

%%% read from hdf5
% grain/neighbor map from synthetic microstructure
neighborlist_synthetic          = transpose(read_dream3d_dataset(path_synthetic_file_input, path_synthetic_hdf5_neighborlist ));
numneighbors_synthetic          =           read_dream3d_dataset(path_synthetic_file_input, path_synthetic_hdf5_numneighbors );
% target misorientation distribution from EBSD
misorientations_target_unsorted = transpose(read_dream3d_dataset(path_ebsd_file_input     , path_ebsd_hdf5_misorientationlist));

%%% import initial/current orientation, misorientation
[orientations_previous, misorientations_previous] = initialize_orientations( ...
    path_temp_file_eulerangles, ...
    path_temp_file_misorientations, ...
    path_synthetic_file_input, ...
    path_synthetic_hdf5_eulerangles, ...
    neighborlist_synthetic, ...
    numneighbors_synthetic, ...
    crystal_symmetry ...
);

%%% initialize KS test
% [~, ~, ks_previous] = kstest2(misorientations_target, misorientations_previous);
[misorientations_target, YMisorientationReference] = initialize_RapidKS(misorientations_target_unsorted);
ks_initial = RapidKS(misorientations_previous',misorientations_target,YMisorientationReference);

%%% main loop
iter   = 0;
n_swaps_successful_epoch = 0;
orientations_current = orientations_previous;
misorientations_current = misorientations_previous;
ks_previous = ks_initial;
ks_current = ks_initial;
time_elapsed_seconds_run = 0.0;
timer = tic;
timer_total = tic;
while abs(ks_previous) > ks_target

    % swap two grains, calculate new misorientations, and calculate KS stat
    [orientations_current, indecies] = swap(orientations_current);
    misorientations_current = get_misorientations( ...
        transpose(indecies), ...
        orientations_current, ...
        misorientations_current, ...
        neighborlist_synthetic, ...
        numneighbors_synthetic ...
    );
    %[~, ~, ks_current] = kstest2(misorientations_target, misorientations_current);
    ks_current = RapidKS(misorientations_current',misorientations_target,YMisorientationReference);
    
    % if the KS statistic is better, keep the swap, else revert the swap
    if abs(ks_current) < abs(ks_previous)
        orientations_previous    = orientations_current;
        misorientations_previous = misorientations_current;
        ks_previous              = ks_current;
        n_swaps_successful_epoch = n_swaps_successful_epoch + 1;
    else
        orientations_current    = orientations_previous;
        misorientations_current = misorientations_previous;
        ks_current              = ks_previous;
    end

    % user feedback
    if mod(iter, iters_per_epoch) == 0 && iter ~= 0
        
        time_elapsed_seconds_epoch = toc(timer);
        time_elapsed_seconds_run   = time_elapsed_seconds_run+time_elapsed_seconds_epoch;
        ks_improvment_run          = ks_initial - ks_current;
        ks_improvment_per_second   = ks_improvment_run/time_elapsed_seconds_run;
        time_remaining        = seconds((ks_current - ks_target)/ks_improvment_per_second);
        time_remaining.Format = "hh:mm:ss";

        disp("Iteration: "+string(iter))
        disp("   "+"Time elapsed last epoch    : " + time_elapsed_seconds_epoch)
        disp("   "+"Successful Swaps last epoch: " + n_swaps_successful_epoch)
        disp("   "+"Current KS: " + ks_current)
        disp("   "+"KS improvment rate over run: " + ks_improvment_per_second)
        disp("   "+"ETA: " + string(time_remaining))

        timer = tic;
        n_swaps_successful_epoch = 0;

    end

    % save state periodically
    if mod(iter, iters_per_backup) == 0 % && iter ~= 0
		save_temp_files(orientations_previous, misorientations_previous, path_temp_file_eulerangles, path_temp_file_misorientations)
		if save_histograms
            save_histogram(iter, misorientations_target_unsorted, misorientations_current, n_bins, path_fig_output)
        end
    end

    iter = iter + 1;

end

%%% final write
save_temp_files(orientations_previous, misorientations_previous, path_temp_file_eulerangles, path_temp_file_misorientations)
if save_histograms
    save_histogram(iter, misorientations_target_unsorted, misorientations_current, n_bins, path_fig_output)
end

%%% display total elapsed time
time_elapsed_seconds_total = seconds(toc(timer_total));
time_elapsed_seconds_total.Format = "hh:mm:ss";
disp("Total Elapsed time: "+string(time_elapsed_seconds_total)) 

%%

function [orientations, misorientations] = initialize_orientations( ...
    path_temp_file_eulerangles, ...
    path_temp_file_misorientations, ...
    path_synthetic_file_input, ...
    path_synthetic_hdf5_eulerangles, ...
    neighborlist, ...
    numneighbors, ...
    crystal_symmetry ...
    )
    
    % continue from previous swap if possible
    if isfile(path_temp_file_eulerangles) && isfile(path_temp_file_misorientations)
        disp("Continuing previous run...")
        eulerangles     = readmatrix(path_temp_file_eulerangles    );
        orientations    = get_orientations(eulerangles, crystal_symmetry);
        misorientations = readmatrix(path_temp_file_misorientations);
        return
    end

    if isfile(path_temp_file_eulerangles)
        disp("Continuing previous run, but recalculating misorientations...")
        eulerangles     = readmatrix(path_temp_file_eulerangles    );
        orientations    = get_orientations(eulerangles, crystal_symmetry);
        misorientations = get_misorientations( ...
            transpose(2:size(numneighbors,1)), ...
            orientations, ...
            zeros(size(neighborlist,1),1), ...
            neighborlist, ...
            numneighbors ...
        );
        return
    end

    % otherwise start from scratch
    disp("Starting from scratch...")
    eulerangles     = read_dream3d_dataset(path_synthetic_file_input, path_synthetic_hdf5_eulerangles);
    orientations    = get_orientations(eulerangles, crystal_symmetry);
    misorientations = get_misorientations( ...
        transpose(2:size(numneighbors,1)), ...
        orientations, ...
        zeros(size(neighborlist,1),1), ...
        neighborlist, ...
        numneighbors ...
    );

end

function [path, name] = rsplit(string_, delimiter)
    string_ = transpose(split(string_, delimiter));
    path = join(string_(1:size(string_,2)-1), delimiter);
    name = string_(size(string_,2));
end

function [neighbor_ids, n_start, n_stop] = get_feature_neighbors(featureid, neighborlist, numneighbors)

    %sum all the featureids before the current
    n_start = sum(numneighbors(1:featureid-1,:))+1;
    n_stop  = n_start+numneighbors(featureid)-1;
    neighbor_ids = neighborlist(n_start:n_stop,:)+1;

end

function misorientations = get_misorientations(featureids, orientations, misorientations, neighborlist, numneighbors)
    
    for i = 1:size(featureids,1)

        featureid = featureids(i);
        [featureneighborids, start, stop] = get_feature_neighbors(featureid, neighborlist, numneighbors);
        
        for j = 1:size(featureneighborids,1)

            featureneighborid = featureneighborids(j);

            misorientations(start+j-1, 1) = get_misorientation( ...
                orientations(featureid), ...
                orientations(featureneighborid) ...
            );

        end
        
    end

end

function misorientation = get_misorientation(orientation1, orientation2)
    misorientation = angle(orientation1, orientation2)/degree;
end

function orientations = get_orientations(eulerangles, crystal_symmetry)
    orientations = orientation.byEuler(eulerangles(:,1), eulerangles(:,2), eulerangles(:,3), crystal_symmetry);
end

function eulerangles = get_eulers(orientations)
    eulerangles = [orientations(:).phi1, orientations(:).Phi, orientations(:).phi2];
end

function [misorientations_target, YMisorientationReference] = initialize_RapidKS(misorientations_target_Unsorted)
    %%% format misorientations, initialize misorientation reference
    %%% for Max's RapidKS test
    misorientations_target_nonUnique = sort([0;misorientations_target_Unsorted;63]); % Max Pinz
    YMisorientationReference_nonUnique = 0:1/(numel(misorientations_target_nonUnique)-1):1; % Max Pinz
    % okay so we have a unique issue here that needs to be managed  % Max Pinz
    [misorientations_target,Ia,Ic] = unique(misorientations_target_nonUnique); % Max Pinz
    YMisorientationReference = YMisorientationReference_nonUnique(Ia); % Max Pinz
end

function KS_Val = RapidKS(Xq,misorientations_target,YMisorientationReference)
    XqSorted = sort(Xq);
    Y_CDF = 0:1/(numel(XqSorted)-1):1;
    YQ = interp1(misorientations_target,YMisorientationReference,XqSorted);
    KS_Val = max(abs(YQ - Y_CDF));
end

function [orientations, indecies] = swap(orientations)
    % pick two uniformly distributed, non-repeating indecies
    % the -1, +1 is to shift to the right such that the first
    % index is not chosen
    indecies = randperm(size(orientations,1)-1,2)+1;
    % swap the two euler angles
    orientations(indecies,:) = orientations(indecies([2,1]),:);
end

function dataset = read_dream3d_dataset(name_file, path_dataset)
    dataset = h5read(name_file, path_dataset);
    dataset = permute(dataset, length(size(dataset)):-1:1);
end

function save_histogram_misorientations(path_output, misorientations_target, misorientations_current, n_bins)
    
    fig = figure('visible','off');
    hold on
    histogram(misorientations_target , n_bins, 'Normalization','probability')
    histogram(misorientations_current, n_bins, 'Normalization','probability')
    xlabel("Misorientation (Degrees)")
    ylabel("Probability")
    legend({"EBSD","Current Swap Iteration"})
    hold off
    
    [~, extension] = rsplit(path_output, '.');
    if any(contains({'pdf', 'eps'}, extension))
        exportgraphics(fig, path_output, 'ContentType', 'vector');
    elseif any(contains({'jpg', 'jpeg', 'png', 'tif'}, extension))
        saveas(fig, strcat(path,'_pole_figure_',label,'.png'));
    end

end

function save_temp_files(orientations_previous, misorientations_previous, path_temp_file_eulerangles, path_temp_file_misorientations)
	eulerangles = get_eulers(orientations_previous);
	disp("   "+"Writing to temp files...")
	writematrix(eulerangles             , path_temp_file_eulerangles    , 'Delimiter',',')
	writematrix(misorientations_previous, path_temp_file_misorientations, 'Delimiter',',')
end

function save_histogram(iter, misorientations_target_unsorted, misorientations_current, n_bins, path_fig_output)
	disp("   "+"Saving misorientation histograms...")
	[head, tail] = rsplit(path_fig_output, '.');
	path_fig_output_modified = head+"_"+string(iter)+'.'+tail;
	save_histogram_misorientations(path_fig_output_modified, misorientations_target_unsorted, misorientations_current, n_bins)
end