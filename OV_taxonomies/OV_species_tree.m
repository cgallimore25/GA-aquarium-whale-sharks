%% OV species phy tree

% .phy file generated from:
% https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
% by adding the species list 'OV_species_sci_names.txt' file and exporting
% as 'phylip tree'

file= 'phyliptree.phy';

tree= phytreeread(file);

[sp_matrix, id, dist]= getmatrix(tree);

matrix= full(sp_matrix); 

% figure;
% imagesc(matrix)

% force symmetry of sparse matrix
% from: https://stackoverflow.com/questions/38170193/adjacency-matrix-must-be-symmetric
sm= sp_matrix; 
sm= sm.' | sm; 

G= graph(sm, id);

figure; 
pG= plot(G);

[leafclusts, nodeclusts]= cluster(tree, [], 'criterion', 'gain', 'maxclust', 3);

% h= plot(tree, 'type', 'radial', 'leaflabels', true, 'branchlabels', true);
h= plot(tree, 'leaflabels', true, 'branchlabels', true);
set(h.BranchLines(nodeclusts == 1), 'Color', 'b')
set(h.BranchLines(nodeclusts == 2), 'Color', 'r')
set(h.BranchLines(nodeclusts == 3), 'Color', 'g')

ind= getbyname(tree, 'Branch');  % index generic branches
