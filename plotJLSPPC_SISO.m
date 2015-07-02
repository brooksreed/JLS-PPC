function [hf] = plotJLSPPC_SISO(r)
% plots (for SISO systems)
% r structure of results
% also r.sys: system settings (only TAU_C required right now)

% in r struct, P is a proxy for size of A
if(size(r.P,1)==2)
    CPlot = [1 0];  % output for plotting
else
    CPlot = 1;
end

plot_XhMPC = 1;  % plots prediction used for computing control
plot_losses = 1; % plots packet losses for c, m (no a right now)

% (could be used with a saved r struct too)
NX_SYS = size(r.P,1);    % underlying system states (no buffer)
SIM_LENGTH = size(r.X,2);

hf = figure;

subplot(3,1,[1 2])
hx = plot(0:SIM_LENGTH-1,CPlot*r.X(1:NX_SYS,:));
hold on
hxh = plot(0:SIM_LENGTH-1,CPlot*r.Xh(1:NX_SYS,:),'g.--');
%title(sprintf('integrator sys, alphaBar = %0.2f, W=%.1f, V=%.1f',alpha_cBar,W,V))

hb=[];hc=[];
if(plot_losses)
    for k = 1:SIM_LENGTH
        if( (r.PI_C(k) - r.alpha_c(k)*r.PI_C(k)) == 1 )
            %hc = plot([k-1+tc k-1+tc], [-0.25 0],'r');
            hc = plot(k-1+r.sys.TAU_C,0,'ro');
        end
        if( (r.PI_M(k) - r.alpha_m(k)*r.PI_M(k)) ==1 )
            %hb = plot([k-1 k-1],[-.75 -.5],'m');
            hb = plot(k-1,-2,'mo');
        end
%         if( (r.Pi_a(k) - r.alpha_a(k)*r.Pi_a(k)) ==1 )
%             % need to save/load tac, do on per-veh. basis
%             ha = plot([k-1-tc-tac k-1-tc-tac],[-.75 -.5],'c');
%         end
    end
end

if(plot_XhMPC)
    hxhm = plot(0:SIM_LENGTH-1,CPlot*squeeze(r.XhMPC(1:NX_SYS,end,1:SIM_LENGTH)),'k.');
    legend([hx hxh hxhm hc hb],'X','XHat','XHatMPC','c loss','m loss')
else
    legend([hx hxh hc hb],'X','XHat','c loss','m loss')
end

subplot(3,1,3)
stairs(0:SIM_LENGTH-1,r.u,'k')
hold on
stairs(0:SIM_LENGTH-1,r.u_no_loss,'b:')
legend('u true','u planned')
%stairs(0:Ns-1,r.w,'b:')
%legend('u','w')
xlabel('time step')

if(size(r.P,1)==1)
    figure
    plot(sqrt(squeeze(r.P)))
    xlabel('time step')
elseif(size(r.P,1)==2)
    P11 = sqrt(squeeze(r.P(1,1,:)));
    P22 = sqrt(squeeze(r.P(2,2,:)));
    P12 = sqrt(squeeze(r.P(1,2,:)));
    figure
    plot([P11,P22,P12])
    xlabel('time step')
    legend('P_{11}^{(1/2)}','P_{22}^{(1/2)}','P_{12}^{(1/2)}')
end
    

