function [vehicleReport,vehRepOut,trueNav,cmdAck] = grabVehicleReport
% waits for postings from mexmoos -- breaks on reception of vehicleReport
% tries to grab pNav and cmdAck as well 
% returns raw vehicle report string, and parsed VehRepOut cell
% [vehicleReport,vehRepOut,trueNav,cmdAck] = grabVehicleReport
%
% Vehicle report - posted at end of each slot by MOOS timing program:
% id:acomms_status:command_ack:cycle_count:slot,x,y:slot,x,y
% id = -1 if this slot is for sending commands (no RX of measurements)
% trueNav is true vehicle nav at end of slot
% trueNav should be posted BEFORE vehicleReport


% BR,9/19/2014

vehicleNames = {'NOSTROMO','SILVANA','KESTREL'};
Nv = length(vehicleNames);

% (preallocate larger than expected)
msgSave = cell(1,12);

% always waiting on report
gotRep = 0;
count = 1;
while(~gotRep)
    
    % grab mail from mexmoos
    % (most recent messages are lower indices in the buffer)
    msgs=mexmoos('FETCH');
    
    % save msgs
    if(~isempty(msgs))
        msgSave{count} = msgs;
        count = count+1;
    end
    
    % msgs SHOULD first have nav (beginning of slot)
    % and then vehicle report (with the same nav sent over acomms)
    for k = 1:length(msgs)
        msgStr = msgs(k).KEY;
        if(strcmp(msgStr,'PURSUIT_VEHICLE_REPORT'))
            vehicleReport = msgs(k).STR;
            gotRep = 1;
            fprintf('\nGot Vehicle Report')
            vehRepOut = textscan(vehicleReport,...
                '%d %d %d %d %s %s %s %s','Delimiter',':');
        end
    end
end

% now fill up trueNav and cmdAck...
% should be in beginning of msgSave (before vehicle report)
trueNav = zeros(2,Nv);
cmdAck = zeros(1,Nv);
gotAll = zeros(1,3*Nv);

for j = 1:length(msgSave)
    tmpmsg = msgSave{j};
    for k = 1:length(tmpmsg)
        msgStr = tmpmsg(k).KEY;
        % grab nav and ACKs
        for i = 1:Nv
            if(strfind(msgStr,vehicleNames{i}))
                % NAV_X
                if(strfind(tmpmsg(k).KEY,'X'))
                    trueNav(1,i) = tmpmsg(k).DBL;
                    gotAll(i) = 1;
                % NAV_Y
                elseif(strfind(tmpmsg(k).KEY,'Y'))
                    trueNav(2,i) = tmpmsg(k).DBL;
                    gotAll(Nv+i) = 1;
                end
                % COMMAND_RECEIVE
                if(strfind(tmpmsg(k).KEY,'COMMAND'))
                    cmdAck(i) = tmpmsg(k).DBL;
                    gotAll(2*Nv+i) = 1;
                end
            end
        end
    end
end

disp('nav and cmdAck')
disp(gotAll)




