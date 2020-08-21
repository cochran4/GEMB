function DopamineSensitivity

% Initialization
Vda0   = velocityDopamine;
y0_da  = [41; 319; 126; 0.36; 2.65; 81; 0.002; 7.6780; 942.8183];

% Indices for modulated parameters (pathway,parameter #)
ix  = [1,5; 2,3; 2,6; 3,2; 4,2; 5,2; 6,2; 14,1; 7,1; 8,2; 1,6];

% Initialize parameters
np    = size(ix,1); 
nruns = 1000;
T     = [0,100]; % 200 simulated hours


% Build design matrix
rng(0)
X     = 1 + lhsnorm(zeros(np,1),eye(np)*0.05^2,nruns);

% Run through design matrix
for i=1:nruns
    % Display counter
    if mod(i-1,100)==0, disp(i), end
    
    % Rescale parameters
    Vda    = Vda0;
    for j=1:np
        Vda(ix(j,1)).param(ix(j,2)) = Vda(ix(j,1)).param(ix(j,2))*X(i,j);
    end
    
    % Run simulation
    hdl    = @(t,y)Monoamine(t,y,Vda);
    [t,Y]=ode23t(hdl,T,y0_da);
    
    % Collect output of extracellular dopamine
    eda(i,1) = Y(end,7);
end

% Calculate PRCC
rho           = partialcorr( [X,eda],'Type','Spearman');
dopamine_prcc = rho(1:end-1,end);

% Save weights
save('dopamine_prcc','dopamine_prcc')


% Turn into weights
weights     = zeros(12,1);
weights(1)  = abs(dopamine_prcc(1));                 % TH
weights(2)  = mean( abs(dopamine_prcc([2,3])) );     % QDPR
weights(3)  = abs(dopamine_prcc(4));                 % SLC7A5
weights(4)  = abs(dopamine_prcc(5));                 % DDC
weights(5)  = abs(dopamine_prcc(6))/2;               % SLC18A1
weights(6)  = abs(dopamine_prcc(6))/2;               % SLC18A1
weights(7)  = mean( abs(dopamine_prcc([7,8])) );     % SLC6A3
weights(8)  = mean( abs(dopamine_prcc([9,10])) )/3;  % COMT
weights(9)  = mean( abs(dopamine_prcc([9,10])) )/3;  % MAOA
weights(10) = mean( abs(dopamine_prcc([9,10])) )/3;  % MAOB
weights(11) = abs(dopamine_prcc(11))/2;              % DRD2
weights(12) = abs(dopamine_prcc(11))/2;              % DRD3
tmp.Gene   = {'TH', 'QDPR', 'SLC7A5', 'DDC',  ...
              'SLC18A1', 'SLC18A2', 'SLC6A3', ...
              'COMT', 'MAOA', 'MAOB', 'DRD2', 'DRD3'}';
tmp.Weight = weights/sum(weights);

% Save weights
writetable( struct2table(tmp), 'Weights_Dopamine.xlsx');