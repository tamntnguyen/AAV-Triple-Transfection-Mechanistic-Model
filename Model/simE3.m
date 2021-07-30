function [t,x] = simE3(x0,k,tEnd,outputTimes)
% Initial conditions - number of plasmids 
tEx = 6;
options = odeset('RelTol',1e-2,'AbsTol',1e-3,'NonNegative',1);

ic = zeros(21,1); 
ic(1)  = x0(1); % pRC
ic(5)  = x0(2); % pVector
ic(9)  = x0(3); % pHelper

odefun  = @(t,x) ode_viralProd(t,x,k);

[t1,x1]  = ode15s(odefun,[0 tEx],ic,options);

% Media exchange
ic2     = x1(end,:)'; 
ic2(1)  = 0; % Packaging
ic2(5)  = 0; % Vector
ic2(9)  = 0; % Helper




if isempty(outputTimes)
    [t2, x2]    = ode15s(odefun,[tEx tEnd],ic2,options);
    t           = [t1; t2];
    x           = [x1; x2];
else
    [t2, x2]    = ode15s(odefun,[tEx; outputTimes],ic2,options);
    t = outputTimes;
    x = x2(2:end,:);
end
 

end
