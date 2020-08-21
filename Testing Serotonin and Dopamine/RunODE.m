function [T1,Y1,T2,Y2]=RunODE
%--------------------------------------------------------------------------
% Set up and run ODE model of Serotonin and Dopamine
% Based on models from J. Best
% Substance ID
% 1 - bh2
% 2 - bh4
% 3 - trp/tyr
% 4 - 5htp/l-dopa
% 5 - c5ht/cda (Monoamine in cell outside of vescicles)
% 6 - v5ht/vda (Monoamine in vescicles within cell)
% 7 - e5ht/eda (Monoamine outside of cell)
% 8 - 5hiaa/hva
% 9 - trp-pool/tyr-pool
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set-up ODE
%--------------------------------------------------------------------------

V5ht = velocitySerotonin;
Vda  = velocityDopamine;

% Supposed steady-states for serotonin model
y0_5ht  = [0.14; 0.86; 20.6; 2.26; 0.5; 21.45; 0.000768; 5.26; 144.9];

% Steady-states for dopamine model (hva and tyr-pool not mentioned)
% To get hva and tyr-pool, started with 0 and 0, ran 30s to steady-state,
% plugged into steady-state hva and tyr-pool values, and repeated process 
% until convergence
y0_da  = [41; 319; 126; 0.36; 2.65; 81; 0.002; 7.6780; 942.8183];

%--------------------------------------------------------------------------
% Run ODE
%--------------------------------------------------------------------------

% Serotonin
hdl = @(t,y)Monoamine(t,y,V5ht);
T   = [0,30];
[T1,Y1]=ode23t(hdl,T,y0_5ht);
%[T,Y]=ode15s(hdl,T,y0);
%[T,Y]=ode23(hdl,T,y0);

T0 = 10;
subplot(3,1,1)
plot(T1(T1<T0),Y1(T1<T0,5))
ylabel('Cystolic Serotonin')
subplot(3,1,2)
plot(T1(T1<T0),Y1(T1<T0,6))
ylabel('Vescicular Serotonin')
subplot(3,1,3)
plot(T1(T1<T0),Y1(T1<T0,7))
ylabel('Extracellular Serotonin')


% Dopamine
hdl = @(t,y)Monoamine(t,y,Vda);
T   = [0,10];
[T2,Y2]=ode23t(hdl,T,y0_da);

f = [Monoamine(1,y0_da,Vda), y0_da, Y2(end,:)'];

figure
subplot(3,1,1)
plot(T2(T2<T0),Y2(T2<T0,5))
ylabel('Cystolic Dopamine')
subplot(3,1,2)
plot(T2(T2<T0),Y2(T2<T0,6))
ylabel('Vescicular Dopamine')
subplot(3,1,3)
plot(T2(T2<T0),Y2(T2<T0,7))
ylabel('Extracellular Dopamine')

disp([y0_da(:)'; Y2(end,:)])
disp([])
