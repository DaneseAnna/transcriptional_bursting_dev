clear;clc;
addpath(genpath(pwd));

% load('results\results_MEF\result_gene_26.mat', 'result');
% load('results\results_MEF\result_gene_26.mat', 'data');
tic;
load('results\results_test\result_gene_Atf4.mat', 'result');
load('results\results_test\result_gene_Atf4.mat', 'data');

figureResult2(data,result(:,end));
toc;


