function mut_seq = mutate(seq,res,aa)

if res>479 && res<579
    res=res+1;
else if res>579
    res=res+6;
end
end

mut_seq = seq;
mut_seq(res-383)=aa;

end