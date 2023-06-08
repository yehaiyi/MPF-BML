%% load 1a model
load('results\Model_1a.mat')

%% load 1b model
load('results\Model_1b.mat')

%% load inferred 3a E2 model 
load('results/workspace/test_1226.mat')

%% 
entropy=[];

for kk = 0:9

    %%%% MCMC
    stream = RandStream("mcg16807","Seed",kk);
    rng(kk)
    tic;
    samples=[];
    N=1e7+1e4;
    Burn_in=1e4;
    Thin_in=1e3;
    USE_H=0;

    num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
    num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];

    positions = size(num_mutants_combine_array,2);    
    
    msa_bin_unique = unique(msa_bin,'row');
    randvalue = randi(stream,[1 size(msa_bin_unique,1)]);
    x = msa_bin_unique(randvalue,:);      % randomly pick a sequence from unique MSA sequence
    n=0;

    if(USE_H)
        energy = x*H';
    else
        energy = x*J_MPF_BML*x';   % engery of initial sequence 
    end

    while n<N
        xnew =x;
        energy_new = energy;
        pos = randi(stream,positions);
        index = randi(stream,[0,num_mutants_combine_array(1,pos)]); % index=0 i.e. consensus aa
        Start_site = num_mutants_combine_array_acc_all(pos)+1;
        End_site = num_mutants_combine_array_acc_all(pos+1);
        one = find(xnew(Start_site:End_site)>0);     % 目前seq x 在pos上mutant所在的位置
        if(~isempty(one) && ~logical(USE_H)) % 如果x在pos上不是consensus aa
            energy_new = energy_new-2*xnew*J_MPF_BML(:,Start_site+one-1)+J_MPF_BML(Start_site+one-1,Start_site+one-1);
        end
        if(index==0) 
            if(USE_H)
                energy_new = energy_new-xnew(Start_site:End_site)*H(Start_site:End_site)';
                xnew(Start_site:End_site)=0;                        
            else            
                xnew(Start_site:End_site)=0; 
            end
              
        else
            if(USE_H)
                energy_new = energy_new-xnew(Start_site:End_site)*H(Start_site:End_site)';
            end
            xnew(Start_site:End_site)=0;
            xnew(Start_site+index-1)=1;
            if(USE_H)
                energy_new = energy_new+xnew(Start_site+index-1)*H(Start_site+index-1)';
            else
                energy_new = energy_new+2*xnew*J_MPF_BML(:,Start_site+index-1)-J_MPF_BML(Start_site+index-1,Start_site+index-1); 
            end
        end
        trans_rate = exp(energy-energy_new);
        %trans_rate = 1/(1+exp((-energy+energy_new)));
        u = rand(stream);
        if energy_new<=energy
            x=xnew;
            energy=energy_new;
        elseif u < trans_rate
            x=xnew;
            energy=energy_new;
        end

        % sample
        if(n+1>Burn_in)
            if(mod((n+1),Thin_in)==0)
                samples=[samples;x];
            end
        end

        n=n+1;
        %     if rand < 1/(1+exp(xnew*triu(J_MPF_BML)*xnew'-x*triu(J_MPF_BML)*x'))

    end
    toc;

    %%%%%%%%%%% entropy
    [samples_unique, ind1, ind2]= unique(samples,'rows');


    weight_seq = ones(1,size(samples,1))/size(samples,1);
    weight_seq_unique=[];
    for indi_bin = 1:length(ind1)
        num_term = ind2(ind1(indi_bin));
        ind_values = find(ind2==num_term);
        weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
    end

    Z=0;
    E=0;
    for i =1:size(samples_unique,1)
        Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
        E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
    end

    entropy_maxent = E+log(Z);
    entropy= [ entropy; entropy_maxent];

end
mean(entropy)

%%  MY MCMC CODE
entropy=[];

for kk = 0:9
    stream = RandStream("mcg16807","Seed",kk);
    rng(kk)
    tic;
    samples=[];
    acp_rate=[];
    N=1000; %1e7
    Burn_in=0; %1e4
    Thin_in=1; %1e3
    USE_H=0;

    current_e=[];
    new_e=[];

    num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
    num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];

    positions = size(num_mutants_combine_array,2);    

    msa_bin_unique = unique(msa_bin,'row');
    randvalue = randi(stream,[1 size(msa_bin_unique,1)]);
    x = msa_bin_unique(randvalue,:);      % randomly pick a sequence from unique MSA sequence
    n=0;

    if(USE_H)
        energy = Calculate_Energy(x, J_MPF_BML,H);
    else 
        energy = Calculate_Energy(x, J_MPF_BML);
    end

    while n<N
        xnew =x;
        current_e = [current_e; energy];
        pos = randi(stream,positions);
        index = randi(stream,[0,num_mutants_combine_array(1,pos)]);  % new mutant index=0 i.e. consensus aa
        Start_site = num_mutants_combine_array_acc_all(pos)+1;
        End_site = num_mutants_combine_array_acc_all(pos+1);
        xnew(Start_site:End_site)=0;           
        if index~=0
            xnew(Start_site+index-1)=1;
        end

        if(USE_H)
        energy_new = Calculate_Energy(xnew, J_MPF_BML,H);
        else 
        energy_new = Calculate_Energy(xnew, J_MPF_BML);
        end

        new_e =[new_e; energy_new];
   
        trans_rate = exp(energy-energy_new);
        acp_rate =[ acp_rate; trans_rate] ;
        %trans_rate = 1/(1+exp((-energy+energy_new)));
        
        u=rand(stream);
        if energy_new<=energy
            x=xnew;
            energy=energy_new;
        elseif u < trans_rate
            x=xnew;
            energy=energy_new;
        end

        % sample
        if(n+1>Burn_in)
            if(mod((n+1),Thin_in)==0)
                samples=[samples;x];
            end
        end

        n=n+1;
       

    end
    toc;

    %%%%%%%%% entropy
    [samples_unique, ind1, ind2]= unique(samples,'rows');


    weight_seq = ones(1,size(samples,1))/size(samples,1);
    weight_seq_unique=[];
    for indi_bin = 1:length(ind1)
        num_term = ind2(ind1(indi_bin));
        ind_values = find(ind2==num_term);
        weight_seq_unique(indi_bin) = sum(weight_seq(ind_values));
    end

    Z=0;
    E=0;
    for i =1:size(samples_unique,1)
        Z = Z+exp(-samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)');
    
        E = E+samples_unique(i,:)*J_MPF_BML*samples_unique(i,:)'*weight_seq_unique(i);
    end

    entropy_maxent = E+log(Z);
    entropy= [ entropy; entropy_maxent];
end

mean(entropy)

%% 
std(entropy,0)
std(entropy,1)

%% Energy Trace 
rng(0)
stream = RandStream("mcg16807","Seed",0);
tic;
samples=[];
N=1e6;
Burn_in=0; %1e4
Thin_in=1; %1e3
USE_H=0;

num_mutants_combine_array_acc = cumsum(num_mutants_combine_array);
num_mutants_combine_array_acc_all = [0 num_mutants_combine_array_acc];
positions = size(num_mutants_combine_array,2);    
msa_bin_unique = unique(msa_bin,'row');
randvalue = randi(stream,[1 size(msa_bin_unique,1)]);
x = msa_bin_unique(randvalue,:);      % randomly pick a sequence from unique MSA sequence
n=0;
energy =  x * J_MPF_BML * x';
Energy_trace = [energy];

while n<N
    xnew =x;
    pos = randi(stream,positions);
    index = randi(stream,[0,num_mutants_combine_array(1,pos)]);  % new mutant index=0 i.e. consensus aa
    Start_site = num_mutants_combine_array_acc_all(pos)+1;
    End_site = num_mutants_combine_array_acc_all(pos+1);
    xnew(Start_site:End_site)=0;           
    if index~=0
        xnew(Start_site+index-1)=1;
    end

    energy_new = xnew * J_MPF_BML * xnew';

    trans_rate = exp(energy-energy_new);
    %trans_rate = 1/(1+exp((-energy+energy_new)));
    u= rand(stream);
    if energy_new<=energy
        x=xnew;
        energy=energy_new;
    elseif u < trans_rate
        x=xnew;
        energy=energy_new;
    end

    % sample
    %if(n+1>Burn_in)
        %if(mod((n+1),Thin_in)==0)
            %samples=[samples;x];
        %end
    %end
    Energy_trace = [Energy_trace; energy];
    n=n+1;
%     if rand < 1/(1+exp(xnew*triu(J_MPF_BML)*xnew'-x*triu(J_MPF_BML)*x'))

end
toc;
%% 
 plot(0:N,Energy_trace)