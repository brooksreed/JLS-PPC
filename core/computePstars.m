function E_terms_squared = computePstars(n_stars,ALPHAC_BAR,MODE)
% compute Pstar coefficients for a given alphaBar
% E_termsSquared = computePstars(n_stars,ALPHAC_BAR,MODE)
% MODE= 'NUMERIC';
% MODE = 'ANALYTIC'; (alphaBar not used - any input is fine)
%
% Symbolic version used to generate fast analytic lookup
%   (experienced problems with symbolic when nStars>11)
% Numeric version with a given alphaBar
%   (much faster than the symbolic version, but lookup is faster)


% G describes the rearranging of terms
% u(t,t), u(t,t-1), etc. to (u(t,t)-u(t,t-1)) etc.  
G = zeros(n_stars+1,n_stars);
G(1,1:n_stars) = ones(1,n_stars);
G(2:n_stars+1,1:n_stars) = diag(-ones(1,n_stars));

% Define some symbolic vars
% ab = alphaBar, at = alpha(t)
syms ab
Eat = ab;
syms at

single_coeffs = cell(n_stars,1);
% coeffs for u(t), u(t,t-1) before rearranging
single_coeffs{1} = [(at-ab);(1-at)-(1-ab)];
% loop through recursively creating the estimator error expressions
for star = 2:n_stars
    
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
    single_coeffs{star} = sym('tmp',[star+1,1]);
    single_coeffs{star}(1:star-1,:) = single_coeffs{star-1}(1:star-1,:);
    single_coeffs{star}(star,:) = dangus;
    single_coeffs{star}(star+1,:) = dungus;
    
end

% solve linear system to rearrange terms
% (note -- atm* are unique as input here, so don't end up with ^power)
diff_coeffs = G\single_coeffs{n_stars};
term = cell(n_stars,1);
for i = 1:n_stars
    term{i} = diff_coeffs(i);
end

% square terms and store in string form
terms_squared = cell(n_stars,1);
terms_str = cell(n_stars,1);

if(strcmp(MODE,'ANALYTIC'))
    E_terms_squared = cell(n_stars,1);
end

for i = 1:n_stars
    
    terms_squared{i} = expand(term{i}^2);    

    % PREP FOR EXPECTATIONS BY SUBSTITUTING VARIABLE NAMES
    % (^2 after something that is NOT ab --> pow2
    expr = '(?<!ab)\^2';
    tmp = regexprep(char(terms_squared{i}),expr,'pow2');

    % check to make sure no higher powers - shouldn't be
    expr_test = '(?<!ab)\^\d';
    tmp_test = regexprep(tmp,expr_test,'TEST');
    if(~isempty(strfind(tmp_test,'TEST')))
        disp('WARNING:  higher powers than 2 found')
    end
    
    % turn any term with "at" to Eat
    % (keeping * in between terms)
    %   leveraging the fact that Eatm* = Eat = ab (defined above)
    % and the fact that powers of a single term already converted
    %   since E[at^2] = E[at^*]... = E[at] = ab
    expr2 = '([\w]*at[\w]*)';
    terms_str{i} = regexprep(tmp,expr2,'Eat');
    
    if(strcmp(MODE,'ANALYTIC'))
        
        % plug in and evaluate strings
        E_terms_squared{i} = eval(terms_str{i});
        fprintf('%d:   ',i);%disp(factor(E_termsSquared{i}))
        disp(E_terms_squared{i})
        
        %
        diary off
        diary on
        %
        
    end
    
end

if(strcmp(MODE,'NUMERIC'))
    
    % plug in and evaluate strings 
    ab = ALPHAC_BAR;
    Eat = ALPHAC_BAR;
    E_terms_squared = zeros(n_stars,1);
    for i = 1:n_stars
        E_terms_squared(i)= eval(terms_str{i});
    end

end

