function v = velocityDopamine
%--------------------------------------------------------------------------
% Velocity information for Dopamine  
%--------------------------------------------------------------------------

% Tryosine hydroxylase
v(1).idx   = [1,1;2,-1; 3,-1; 4,1]; % Substance index and positive/negative contribution
v(1).fname = 'v_TH(y(2),y(3),y(5),y(7),v(1).param)';
v(1).param = [46,60,110,160,125,4.5];

% Dihydropteridine reductase
v(2).idx   = [1, -1;2,1];
v(2).fname = 'v_4Var(y(1),125,y(2),10,v(2).param)';
v(2).param = [100 75 200 10 75 80];

% Neutral amino acid transporter
v(3).idx   = [3,1];
v(3).fname = 'v_1Var(97,v(3).param)';
v(3).param = [64, 400];

% Aromatic amino acid decarboyxlase
v(4).idx   = [4, -1; 5, 1];
v(4).fname = 'v_1Var(y(4),v(4).param)';
v(4).param = [130 10000];

% Vesicular monoamine transporter 
v(5).idx   = [5, -1; 6, 1];
v(5).fname = 'v_1Var(y(5),v(5).param(1:2))-v(5).param(3)*y(6)';
v(5).param = [3 7082  40];%115.8841]; %[40 0.198 3500];

% Dopamine transporter
v(6).idx   = [5, 1; 7, -1];
v(6).fname = 'v_1Var(y(7),v(6).param)';
v(6).param = [0.2 8000];

% Catabolism of c5ht
v(7).idx   = [5,-1; 8, 1];
v(7).fname  = 'v(7).param(1)*y(5)'; %'v_1Var(y(5),v(7).param)';
v(7).param = 10; %[95,1000];

% Catabolism of e5ht
v(8).idx   = [7,-1; 8, 1];
v(8).fname = 'v_1Var(y(7),v(8).param)';
v(8).param = [3,30];

% Release of vda out of vescicles
v(9).idx   = [6,-1;7,1];
v(9).fname = 'fire(t)*y(6)';
v(9).param = [];

% Catabolism of hva
v(10).idx   = [8,-1];
v(10).fname = 'v(10).param*y(8)';
v(10).param = 3.45;

% Catabolism of tyr-pool
v(11).idx   = [9,-1];
v(11).fname = 'v(11).param*y(9)';
v(11).param = 0.2;

% Catabolism of typ
v(12).idx   = [3,-1];
v(12).fname = 'v(12).param*y(3)';
v(12).param = 0.2;

% Flow in and out of trp-pool
v(13).idx   = [3,-1; 9,1];
v(13).fname = 'v(13).param(1)*y(3) - v(13).param(2)*y(9)';
v(13).param = [6,0.6];

% Removal of external Serotonin
v(14).idx   = [7,-1];
v(14).fname = 'v(14).param*y(7)';
v(14).param = 400;
