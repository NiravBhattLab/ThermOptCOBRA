% This code gets the results on list of stochiometrically blocked and thermodynamically blocked
% reactions in bigg database. Also it gets the list of TICs and the respective direactions
% Results of 99 models are given in the folder 'BiggBlkdResults'. For
% remaining models that ran for more than 24 hrs only partially identifed
% TICs are provided (These models have the name prefix 'Incom')
clear
A = importdata('BiggModels.txt');
for i=1:numel(A)
    load(A{i})
    [Tcc,TICs,Dir] = ThermOptCC(model,1e-5);
    fcc = fastcc(model,1e-5,0);
    save(['./BiggBlkdResults/',A{i}],'Tcc','TICs','Dir','fcc')
end
