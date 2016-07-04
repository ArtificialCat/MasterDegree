function [y, parameter] = fitness_auc(x, input, target)

warning('off','all');
M = size(x,1); % number of chromosomes
numPar = 50;
y = zeros(M,1);
parameter = y;
for chr = 1:M
    fprintf('\tchromosome %d...\n', chr);
    idx = x(chr,:)==1;
    temp = dataset(input(:,idx), target);
%     data = oc_set(temp,'0');
    auc_mean = zeros(numPar,1);  
    parfor k = 1:numPar
        try
            theta = k * 0.01;
            fprintf('\t\ttheta = %.2f\n',theta);
            numBags = 5;
            I = numBags;
            auc = zeros(1,numBags);        
            for i = 1:numBags
                [trn, tst, I] = dd_crossval(data,I);
                w = mog_dd(trn,theta);
                auc(1,i) = dd_auc(tst*w*dd_roc);
            end            
            auc_mean(k,1) = mean(auc);
        catch MExc
            auc_mean(k,1) = 0;
            fprintf('\t\t\tException...\n');
        end
    end    
    y(chr,1) = max(auc_mean(:,1));
    parameter(chr,1) = find(auc_mean(:,1) == y(chr,1), 1) * 0.01;    
end


