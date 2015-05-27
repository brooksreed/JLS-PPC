% generate and save an ocean process run
% save Ap,Cp,Wp,ApLoners,CpLoners,WpLoners,(Wn,Z),ApInt,CpInt,WpInt
% xp(Nx,Ns) ,wpt(Nw,Ns)
% pb(2,Nv,Ns) (reference origins)
% gb(2,Nv,Ns) (normalized ref gradients -- unit vector)

% BR, 9/16/2014

Nv = 3;
Ns = 220;


% controlled integrator chasing chained mass
%kiCross = .0405; % for Wn1 .3487
% %wp = 0.5;  % for 1.6*8 sec
% %wp = .1;   % for 0.8*7 sec 
%wp = .35;    % for 1.4*7 sec = 9.8 m


% for testing -- slower
kiCross = .02; %
%wp = 0.5;  % for 1.6*8 sec
%wp = .1;   % for 0.8*7 sec 
%wp = .3;    % for 1.2*7 sec = 8.4 m

saveStr = sprintf('expRefs_test2_Ns%d_kc%s,wp%s_%s',Ns,printNumFile(kiCross,4),printNumFile(wp,3),dateString('DHM'));


%% REFERENCES
% local MOOS x,y  grid 

pb = zeros(2,Nv,Ns);
gb = zeros(2,Nv,Ns);

% to start, using fixed references
% (if change to moving, add in adjustments to speed constraints)

% spread out in front of dock
%pb1 = [-30,-200];
%pb2 = [85,-160];
%pb3 = [190,-120];

% (for 9/23 first part of day)
%pb1 = [210,-140];
%pb2 = [240,-130];
%pb3 = [270,-120];

% for 9/23 end of day
% pb1 = [160,-180];
% pb2 = [180,-170];
% pb3 = [200,-160];
% gbAng = 120;

% parallel to dock close in
pb1 = [30,-60];
pb2 = [40,-80];
pb3 = [50,-100];
gbAng = 22;



% gb1 = [0,1];
% gb2 = [0,1];
% gb3 = [0,1];

gb1 = [cos(deg2rad(gbAng)),sin(deg2rad(gbAng))];
gb2 = [cos(deg2rad(gbAng)),sin(deg2rad(gbAng))];
gb3 = [cos(deg2rad(gbAng)),sin(deg2rad(gbAng))];

pb(:,1,:) = repmat(pb1',[1,Ns]);
pb(:,2,:) = repmat(pb2',[1,Ns]);
pb(:,3,:) = repmat(pb3',[1,Ns]);
gb(:,1,:) = repmat(gb1',[1,Ns]);
gb(:,2,:) = repmat(gb2',[1,Ns]);
gb(:,3,:) = repmat(gb3',[1,Ns]);




%% PERTURBATIONS

[Ap,Cp,Wp,ApLoners,CpLoners,WpLoners,Wn,Z,ApInt,CpInt,WpInt] = ...
    generateChainedMasses(Nv,kiCross,wp);

% deterministic forcing -- drive middle mass at Wn(1)
forcing = zeros(size(Wp,1),Ns);
forcing(Nv+2,:) = 0.1*sin(Wn(1)*(0:(Ns-1)));
disp('with deterministic forcing')

% generate and check specific perturbations
% save the wp used to generate this

x0p = -1+10*rand(Nv,1);

OKtrial = 0;
while(~OKtrial)
    wpt = real(sqrtm(Wp))*randn(2*Nv,Ns);
    wpt(:,1) = -abs(wpt(:,1)*2);    % first noise negative (offset sine)
    wpt = wpt+forcing;

    xp = zeros(2*Nv,Ns);    xp(:,1) = Cp'*x0p;
    for t = 1:Ns-1
        xp(:,t+1) = Ap*xp(:,t) + wpt(:,t);
    end

    % check speeds
    speeds=zeros(Nv,Ns-1);
    for i = 1:Nv
        speeds(i,:) = diff(Cp(1,:)*xp);
    end
    speeds = reshape(speeds,[1,Nv*(Ns-1)]);

    hf = figure;
    subplot(3,1,1)
    hist(speeds)
    title('speeds')
    xlabel('horz units/time step')
    subplot(3,1,[2 3])
    stairs((Cp*xp)')
    xlabel('time steps')

    okt = input('OK trial (y) or retry (n):  ','s');
    if( ~isempty(strfind(okt,'y')) || ~isempty(strfind(okt,'Y')))
        OKtrial=1;
    else
        OKtrial=0;
    end
    close(hf)

end

%save

save([saveStr,'.mat'],'Nv','Ns','Wn','Z','kiCross','wp','Ap','Cp','Wp',...
    'ApLoners','CpLoners','WpLoners','ApInt','CpInt','WpInt',...
    'xp','wpt','pb','gb')
% save Ap,Cp,Wp,ApLoners,CpLoners,WpLoners,(Wn,Z),ApInt,CpInt,WpInt
% xp(Nx,Ns) ,wpt(Nw,Ns)
% pb(2,Nv,Ns) (reference origins)
% gb(2,Nv,Ns) (normalized ref gradients -- unit vector)

