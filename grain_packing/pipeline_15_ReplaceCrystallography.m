clear all;clc

%%% start mtex
path_mtex = "../../MTEX";
addpath(path_mtex)
startup_mtex

%%% paths
% files
path_file_input_eulerangles = "./pipeline_output/5-feature_attributes_ebsd.dream3d";
path_file_input_dream3d     = "./pipeline_output/14-synthetic_grains.dream3d";
path_file_output_dream3d    = "./pipeline_output/15-synthetic_grains.dream3d";
% groups
path_hdf5_cellfeaturedata = "/DataContainers/ImageDataContainer/CellFeatureData";
% datasets
path_synthetic_hdf5_eulerangles   = path_hdf5_cellfeaturedata+"/"+"AvgEulerAngles";

% set crystallographic and planar reference data
crystal_symmetry = crystalSymmetry('m-3m', [4.050 4.050 4.050], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]);
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
psi = deLaValleePoussinKernel('halfwidth',10*degree);

% create new dream3d file
copyfile(path_file_input_dream3d, path_file_output_dream3d, "f")

% find the number of grains that must be overwritten
eulerangles_old = read_dream3d_dataset(path_file_input_dream3d, path_synthetic_hdf5_eulerangles);
n_samples = size(eulerangles_old, 1);

% find the orientation distribution function
eulerangles_population  = import_eulerangles(path_file_input_eulerangles, path_synthetic_hdf5_eulerangles);
orientations_population = get_orientations(eulerangles_population, crystal_symmetry);
odf = calcDensity(orientations_population);

% sample from the orientation distribution function
% to obtain replacement eulerangles
orientations_sample = odf.discreteSample(n_samples);
eulerangles = get_eulers(orientations_sample);

% write samples out to file
write_dream3d_dataset(path_file_output_dream3d, path_synthetic_hdf5_eulerangles, eulerangles)
% read back from file
% reading back from the file rather than reading directly from
% orientations_sample also ensures they wrote correctly
eulerangles_new = read_dream3d_dataset(path_file_output_dream3d, path_synthetic_hdf5_eulerangles);

% compare textures
orientations_old = get_orientations(eulerangles_old, crystal_symmetry);
orientations_new = get_orientations(eulerangles_new, crystal_symmetry);
plot_pole_figures(orientations_old       , 'old', path_file_output_dream3d)
plot_pole_figures(orientations_population, 'txt', path_file_output_dream3d)
plot_pole_figures(orientations_new       , 'new', path_file_output_dream3d)

function orientations = get_orientations(eulerangles, crystal_symmetry)
    orientations = orientation.byEuler(eulerangles(:,1), eulerangles(:,2), eulerangles(:,3), crystal_symmetry);
end

function eulerangles = get_eulers(orientations)
    eulerangles = [orientations(:).phi1, orientations(:).Phi, orientations(:).phi2];
end

function dataset = read_dream3d_dataset(name_file,path_dataset)
    dataset = double(h5read(name_file,path_dataset));
    dataset = permute(dataset,length(size(dataset)):-1:1);
end

function write_dream3d_dataset(name_file,path_dataset, data)
    data = permute(data,length(size(data)):-1:1);
    h5write(name_file,path_dataset, single(data));
end

function plot_pole_figures(orientations, label, path_output)

    odf = calcDensity(orientations);
    miller_indecies = [Miller(1,0,0,odf.CS),Miller(1,1,0,odf.CS),Miller(1,1,1,odf.CS)];

    pole_figure = figure();
    plotPDF(odf, miller_indecies, 'antipodal', 'silent')
    mtexColorMap WhiteJet
    mtexColorbar
	
	path = rsplit(path_output, '.');
    saveas(pole_figure, strcat(path,'_pole_figure_',label,'.png'))

end

function eulerangles = import_eulerangles(path_file, path_avgeulerangles)
	[~, ext] = rsplit(path_file, '.');
	if strcmp(ext, "txt")
		eulerangles = readmatrix(path_file);
	elseif strcmp(ext, "dream3d")
		eulerangles = read_dream3d_dataset(path_file, path_avgeulerangles);
	end
end

function [path, name] = rsplit(string_, delimiter)
    string_ = transpose(split(string_, delimiter));
    path = join(string_(1:size(string_,2)-1), delimiter);
    name = string_(size(string_,2));
end