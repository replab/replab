% use https://github.com/bastibe/Violinplot-Matlab
addpath('~/software/violinplot');
addpath('~/software/export_fig');
data = load('unitaryError.mat');
h = figure;
violinplot(log10(data.errors), data.values_n)
xlabel('Matrix row/column size');
ylabel('Error on individual coefficients (log10)');
export_fig unitaryError.pdf
