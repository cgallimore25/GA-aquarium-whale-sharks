%SUMMARY
% Author: Connor Gallimore
% 02/05/2023

% input:
% matrices 1 and 2 are sparse matrices (though the inputs are 'full') of 
% some measurement organized as trials x timepoints x subjects, with 
% non-measured timepoints represented as NaNs.
% nreps: number of bootstrap repetitions
% nsubj: number of subjects

% for each subject, each iteration, and each time point, bootstrap over 
% non-nan indices of matrix 1's trial dimension, and take the time window 
% average. then average over these same indices for matrix 2 
% (which will contain nans). repeat the same procedure by reversing the 
% places of matrix 1 and matrix 2. seed rng for reproducibility.

% output:
% a subjects x nreps x 2 matrix where the 3rd dim corresponds to the 2 runs
% required to get non-nan data for pre- and post- periods. 
%--------------------------------------------------------------------------


function [boot_fin1, boot_fin2] = bootstrap_same_indices(matrix1, matrix2, nreps, nsubj)

nruns= 2;
tpts= size(matrix1, 2); 

npts_mat1= min( min(sum(~isnan(matrix1), 1), [], 3) ); 
npts_mat2= min( min(sum(~isnan(matrix2), 1), [], 3) ); 

npts_boot= min([npts_mat1 npts_mat2]);

% pre-allocate matrices
t_means= zeros(1, tpts);
boot_m1= zeros(nsubj, nreps, nruns); 
boot_m2= boot_m1; 
boot_fin1= boot_m1;
boot_fin2= boot_m1;
nn_ix= cell(1, tpts, nreps); % non-nan indices
ran_ix= zeros(npts_boot, tpts, nreps);  % get random index

% pre-allocate random seeds for reproducibility 
rs= reshape(1:prod([nruns nsubj tpts nreps]), nruns, nsubj, tpts, nreps);

for j= 1:nruns 

    if j == 1
        m1= matrix1;
        m2= matrix2;
    else
        m1= matrix2;
        m2= matrix1; 
    end

    for s= 1:nsubj
        subj_mat= m1(:, :, s);
        for i= 1:nreps
            for t= 1:tpts
                nn_ix{1, t, i}= find(~isnan(subj_mat(:, t)));
                rng(rs(j, s, t, i), 'twister');
                ran_ix(:, t, i)= randi(length(nn_ix{1, t, i}), 1, npts_boot);
                t_means(t)= mean(subj_mat(nn_ix{1, t, i}(ran_ix(:, t, i)), t));
            end
            boot_m1(s, i, j)= mean(t_means);
        end
    end
    
    for s= 1:nsubj
        subj_mat= m2(:, :, s);
        for i= 1:nreps
            for t= 1:tpts
                t_means(t)= mean(subj_mat(nn_ix{1, t, i}(ran_ix(:, t, i)), t), 'omitnan');
            end
            boot_m2(s, i, j)= mean(t_means, 'omitnan');
        end
    end

    if j == 1
        boot_fin1(:, :, j)= boot_m1(:, :, j);
        boot_fin2(:, :, j)= boot_m2(:, :, j);
    else
        boot_fin1(:, :, j)= boot_m2(:, :, j);
        boot_fin2(:, :, j)= boot_m1(:, :, j);
    end

end

end