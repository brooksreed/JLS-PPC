function E_termsSquared = computePstars(nStars)

% BR 6/22/2015
% (first add after initial local dev)

% confirm that terms are similar
% just run this for largest nStars (not loop)

G = zeros(nStars+1,nStars);
G(1,1:nStars) = ones(1,nStars);
G(2:nStars+1,1:nStars) = diag(-ones(1,nStars));

syms ab
Eat = ab;
syms at

singleCoeffs{1} = [(at-ab);(1-at)-(1-ab)];

for star = 1:nStars
    % looping through whole thing (should be fast)
    % but maybe unnecessary?  
    if(star>1)
        str = sprintf('syms atm%d',star-1);
        eval(str)
        dingus = ['(1-at)*'];
        dangus = ['(1-at)*'];
        for i = 1:(star-2)
            dingus = [dingus sprintf('(1-atm%d)*',i)];
            dangus = [dangus sprintf('(1-atm%d)*',i)];
        end
        dingus = [dingus sprintf('atm%d-(1-ab)^%d*ab',star-1,star-1)];
        dangus = [dangus sprintf('(1-atm%d)-(1-ab)^%d',star-1,star)];

        singleCoeffs{star} = sym('tmp',[star+1,1]);
        singleCoeffs{star}(1:star-1,:) = singleCoeffs{star-1}(1:star-1,:);
        singleCoeffs{star}(star,:) = dingus;
        singleCoeffs{star}(star+1,:) = dangus;
    end
end



%{
% FOR NOW - hardcoded in needed syms
% and definitions of expected values
% (will there be any individual E__ that don't equal ab?)
% P**
%syms ab at atm1 
Eat = ab;
Eatpow2 = ab;
Eatm1 = ab;
Eatm1pow2 = ab;
% P***
%syms atm2
Eatm2 = ab;
Eatm2pow2 = ab;
%}


% (note -- atm* are unique as input here, so don't end up with ^power)
diffCoeffs = G\singleCoeffs{nStars};
for i = 1:nStars
    term{i} = diffCoeffs(i);
end
    
for i = 1:nStars
    
    termsSquared{i} = expand(term{i}^2);    

    % PREP FOR EXPECTATIONS BY SUBSTITUTING VARIABLE NAMES
    % (^2 after something that is NOT ab --> pow2
    expr = '(?<!ab)\^2';
    tmp = regexprep(char(termsSquared{i}),expr,'pow2');

    % check to make sure no higher powers
    exprTest = '(?<!ab)\^\d';
    tmpTest = regexprep(tmp,exprTest,'TEST');
    if(~isempty(strfind(tmpTest,'TEST')))
        disp('WARNING:  higher powers than 2 found')
    end
    
    % turn any term with "at" to Eat
    % (keeping * in between terms)
    %   leveraging the fact that Eatm* = Eat = ab (defined above)
    % and the fact that powers of a single term already converted
    %   since E[at^2] = E[at^*]... = E[at] = ab
    expr2 = '([\w]*at[\w]*)';
    termStr{i} = regexprep(tmp,expr2,'Eat');
    
    % OLD:     
    % add E to beginning of any term with "at"
    %termStr{i} = regexprep(tmp,expr2,'E$0');
    
end

for i = 1:nStars
    E_termsSquared{i} = eval(termStr{i});
    fprintf('%d:   ',i);disp(factor(E_termsSquared{i}))
    disp(E_termsSquared{i})
end

