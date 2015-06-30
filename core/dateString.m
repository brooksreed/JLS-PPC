function datestring=dateString(size)
% creates a string based on the date and time 
% (useful for unique filenames when automating save)
% for example, for creating file names to save
% datestring=dateString(size)
% 2 options for size input
% 'DH' outputs time with just date and hour
% 'DHM' outputs time with just date hour min
% 'DHMS' outputs time with date hour min sec

% BR 
% changed 6/21/2011 to put year - month - day (for better sorting by name)
% 3/18/2013: added 'DHM' option

today=date;
out = textscan(today,'%s %s %s','delimiter','-');
day=char(out{1});
month=char(out{2});
year=char(out{3});
%
clockvec=clock;
hour=clockvec(4);
min=clockvec(5);
sec=round(clockvec(6));
%
switch size
    case 'DH'
        datestring=sprintf('%s-%s-%s-%d',year,month,day,hour);
    case 'DHM'
        datestring=sprintf('%s-%s-%s-%d-%d',year,month,day,hour,min);
    case 'DHMS'
        datestring=sprintf('%s-%s-%s-%d-%d-%d',year,month,day,hour,min,sec);
end