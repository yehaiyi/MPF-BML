function [col_start, col_end] = res2col(amino_single_combine_array, res_idx)

col_start=0;
col_end=0;

for i = 1: res_idx-1
    col_start=length(amino_single_combine_array{i,1})-1 + col_start;
end
col_end = col_start+length(amino_single_combine_array{res_idx}) -1;
col_start = col_start+1;