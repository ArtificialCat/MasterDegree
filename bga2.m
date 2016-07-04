clear;
clc;

warning('off','all');
load('141026_data.mat');
Org_Data=M;

win_size = 3;
Date = unique(Org_Data(:,end));


for date_idx=38:size(Date,1)-win_size

    trn_idx = ismember(Org_Data(:,end),Date(date_idx:date_idx+win_size,1));
%     tst_idx = ismember(M(:,end),Date(date_idx+win_size+1:date_idx+2*win_size,1));
    
    InputDATA = Org_Data(trn_idx(:,1)==1,1:end-2);
    CLAIM_INJ = Org_Data(trn_idx(:,1)==1,end-1);
    

%% GA parameter settings
N = size(InputDATA,2); % number of bits in a chromosome
M = 20; % number of chromosomes must be even
last = 20; % number of generations
sel = 0.5; % selection rate
M2 = 2*ceil(sel*M/2); % number of chromosomes kept
mutrate = 0.01; % mutation rate
nmuts = mutrate*N*(M-1); % number of mutations

%% creates M random chromosomes with N bits
pop = round(rand(M,N)); % initial population

pop_set = cell(last+1,1);
fitness_set = cell(last,1);
parameter_set = cell(last,1);

pop_set{1,1} = pop;

%% Iteration
    for ib = 1:last
        fprintf('Generation: %d\n', ib);
        [fitness, parameter] = fitness_auc(pop, InputDATA, CLAIM_INJ); % Fitness function

        fitness_set{ib,1} = fitness;
        parameter_set{ib,1} = parameter;

        % ranks results and chromosomes
        cost = 1-fitness;
        [cost,ind] = sort(cost);
        pop_parents = pop(ind(1:M2),:);
    %     [ib cost(1)]
        if ib == 1
            fitness_1gen = fitness;
            para_1gen = parameter;
        else
            fitness_2gen = fitness;
            para_2gen = parameter;
        end


        % mate
        cross=ceil((N-1)*rand(M2,1));

        % pairs chromosomes and performs crossover
        pop_children = zeros(M2,N);

        for ic=1:2:M2
    %         pop(ceil(M2*rand),1:cross)=pop(ic,1:cross);
    %         pop(ceil(M2*rand),cross+1:N)=pop(ic+1,cross+1:N);
    %         pop(ceil(M2*rand),1:cross)=pop(ic+1,1:cross);
    %         pop(ceil(M2*rand),cross+1:N)=pop(ic,cross+1:N);

            pop_children(ceil(M2*rand),1:cross)=pop_parents(ic,1:cross);
            pop_children(ceil(M2*rand),cross+1:N)=pop_parents(ic+1,cross+1:N);
            pop_children(ceil(M2*rand),1:cross)=pop_parents(ic+1,1:cross);
            pop_children(ceil(M2*rand),cross+1:N)=pop_parents(ic,cross+1:N);
        end

        % If chromosomes do not have 1s, make random chromosomes 
        for j = 1:size(pop_children,1)
            if sum(pop_children(j,:)) == 0
                pop_children(j,:) = round(rand(1,N));
            end
        end %j

        pop = [pop_parents; pop_children];

        % mutate
        for ic=1:floor(nmuts)
            ix=ceil(M*rand);
            iy=ceil(N*rand);
            pop(ix,iy)=1-pop(ix,iy);
        end %ic

        pop_set{ib+1,1} = pop;
    end %ib
    save(sprintf('%s_mog_result_%d.mat',datestr(date,'yymmdd'),date_idx));
end   


