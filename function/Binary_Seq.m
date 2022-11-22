function seq_bin = Binary_Seq(seq, amino_single_combine_array, ind_conserve)
% input: sequences with char
% ouput:  sequences with binary elements to represent mutations
switch nargin
  case 2
    ind_conserve = [];
end
if ~isempty(ind_conserve)
    seq(:,ind_conserve)=[];
end

seq_bin=[];

for i = 1:size(seq,1)
    bin=[];
    for j=1:length(amino_single_combine_array)
        b = char(amino_single_combine_array{j})==seq(i,j);
        if sum(b)==0   
            b(end)=1;
        end
        b = b(2:end);
        b = flip(b);
        bin = [bin b'];
    end
    seq_bin = [seq_bin;bin];
end

end