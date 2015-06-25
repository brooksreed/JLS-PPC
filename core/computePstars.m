function E_termsSquared = computePstars(nStars,alphaBar)
% Numerically evaluate Pstar coefficients for a given alphaBar

% BR 6/22/2015
% (much faster than the symbolic version)

% G describes the rearranging of terms
% u(t,t), u(t,t-1), etc. to (u(t,t)-u(t,t-1)) etc.  
G = zeros(nStars+1,nStars);
G(1,1:nStars) = ones(1,nStars);
G(2:nStars+1,1:nStars) = diag(-ones(1,nStars));

% Define some symbolic vars
% ab = alphaBar, at = alpha(t)
syms ab
Eat = ab;
syms at

singleCoeffs = cell(nStars,1);
% coeffs for u(t), u(t,t-1) before rearranging
singleCoeffs{1} = [(at-ab);(1-at)-(1-ab)];
% loop through recursively creating the estimator error expressions
for star = 2:nStars
    
    % new alpha(t-(star-1)) var
    str = sprintf('syms atm%d',star-1);
    eval(str)
    
    % construct new terms (2nd-to-last and last) 
    dingus = ['(1-at)*'];
    for i = 1:(star-2)
        dingus = [dingus sprintf('(1-atm%d)*',i)];
    end
    dangus = [dingus sprintf('atm%d-(1-ab)^%d*ab',star-1,star-1)];
    dungus = [dingus sprintf('(1-atm%d)-(1-ab)^%d',star-1,star)];
    
    % construct coefficient vector
    singleCoeffs{star} = sym('tmp',[star+1,1]);
    singleCoeffs{star}(1:star-1,:) = singleCoeffs{star-1}(1:star-1,:);
    singleCoeffs{star}(star,:) = dangus;
    singleCoeffs{star}(star+1,:) = dungus;
    
end

% solve linear system to rearrange terms
% (note -- atm* are unique as input here, so don't end up with ^power)
diffCoeffs = G\singleCoeffs{nStars};
term = cell(nStars,1);
for i = 1:nStars
    term{i} = diffCoeffs(i);
end

% square terms and store in string form
termsSquared = cell(nStars,1);
termStr = cell(nStars,1);
for i = 1:nStars
    
    termsSquared{i} = expand(term{i}^2);    

    % PREP FOR EXPECTATIONS BY SUBSTITUTING VARIABLE NAMES
    % (^2 after something that is NOT ab --> pow2
    expr = '(?<!ab)\^2';
    tmp = regexprep(char(termsSquared{i}),expr,'pow2');

    % check to make sure no higher powers - shouldn't be
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
    
end

% plug in and evaluate strings 
ab = alphaBar;
Eat = alphaBar;
E_termsSquared = zeros(nStars,1);
for i = 1:nStars
    E_termsSquared(i)= eval(termStr{i});
    %fprintf('%d:   ',i);disp(factor(E_termsSquared{i}))
    %disp(E_termsSquared{i})
end

%{
% plug in and evaluate strings 
E_termsSquared = cell(nStars,1);
for i = 1:nStars
    E_termsSquared{i} = eval(termStr{i});
    fprintf('%d:   ',i);disp(factor(E_termsSquared{i}))
    disp(E_termsSquared{i})
end
%}
