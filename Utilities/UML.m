classdef UML <handle
    
% A Matlab toolbox for the updated-maximum-likelihood (UML) procedure
% This software is developed by 
% Wei Dai and Yi Shen,
% The Hearing Lab,
% The Department of Cognitive Sciences,
% The University of California, Irvine
%
% The UML procedure is an adaptive procedure for the efficient estimation
% of the psychometric function.
%
% A detail description on the UML procedure can be found in:
% Shen, Y. and Richards, V.M. (2012). "A maximum-likelihood procedure for
% estimating psychometric functions: Thresholds, slopes, and lapses of
% attention," J. Acoust. Soc. Am.,132(2), 957-967.
% 
% Shen, Y., Dai, W., & Richards, V. M. (2015). "A MATLAB toolbox for the 
% efficient estimation of the psychometric function using the updated 
% maximum-likelihood adaptive procedure. Behavior research methods, 47(1), 
% 13-26.
%
% The UML toolbox is freely available and freely redistributable, according
% to the conditions of the GNU General Public License
% (http://www.gnu.org/licenses/gpl-3.0.txt). You may not distribute the
% software, in whole or in part, in conjunction with proprietary code. That
% means you only have my permission to distribute a program that uses my
% code if you also make freely available (under the terms of the GNU GPL)
% the source code for your whole project. You may not pass on the software
% to another party in its current form or any altered, embellished or
% reduced form, without acknowledging the author and including a copy of
% this license. The software does not come with ANY WARRANTY.
%
% Any question or suggestion regarding this software should be directed to:
% Yi Shen
% shen2@indiana.edu
%
% Class definition for the UML class.
% To construct an object under the UML class:
%
% uml = UML(par)
% 
% The input argument "par" takes a specific form, here is an example of how 
% to set up the data structure "par":
% 
% ---------------------------------------------------------------------------
% % par.model = 'logit';    % the model for the psychometric function, choose 
% %                         % between 'logit' for the logistic psychometric
% %                         % function, 'weibull' for the Weibull 
% %                         % psychometric function, or 'gaussian' for the
% %                         % cumulative gausian function.
% %                         
% % The logistic psychometric function is given by:
% % 
% % p = gamma+(1-gamma-lambda).*(1+exp(-(x-alpha).*beta)).^(-1),
% % 
% % where p is proportion of positive responses, x is the signal strength, 
% % gamma indicates the lower asymptote of the psychometric function (e.g., 
% % gamma = 0.5 for 2AFC experiments), lambda is the lapse rate (difference 
% % between the upper asymptote of the psychometric function and 1, alpha and
% % beta are the threshold and beta parameters, respectively.
% % 
% % The Weibull psychometric fucntion is given by:
% % 
% % p = gamma+(1-gamma-lambda).*(1-exp(-(x./alpha).^beta)).
% % 
% % Note that, for Weibull psychometric function, the alpha parameter has to
% % be strictly greater than zero, and the beta parameter has to be greater 
% % or equal to one.
% % 
% % par.ndown = 2;  % the parameter for the up-down sweetpoint selection rule
% % par.mup = 1;  % the parameter for the up-down sweetpoint selection rule
% % par.x0 = 30;    % the initial signal strength
% % 
% % par.alpha = struct(...
% %     'limits',[-10 30],...       %range of the parameter space for alpha
% %     'N',61,...                %number of alpha values. If this value is set to 1, then the first element of alpha.limits would be the assumed alpha and the alpha parameter is not estimated.
% %     'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
% %     'dist','flat',...         %prior distribution of the alpha parameter. Choose between 'norm' and 'flat'.
% %     'mu',0,...                %mean of the prior distribution.
% %     'std',20 ...              %standard deviation of the prior distribution.  
% %     );
% % 
% % par.beta = struct(...
% %     'limits',[.1 10],...      %range of the parameter space for beta
% %     'N',41,...                %number of beta values. If this value is set to 1, then the first element of beta.limits would be the assumed beta and the beta parameter is not estimated.
% %     'scale','log',...         %the linear or log spacing. Choose between 'lin' and 'log'.
% %     'dist','flat',...         %prior distribution of the beta parameter. Choose between 'norm' and 'flat'.
% %     'mu',1,...                %mean of the prior distribution.
% %     'std',5 ...               %standard deviation of the prior distribution.
% %     );
% % 
% % par.gamma = 0.5;
% % 
% % par.lambda = struct(...
% %     'limits',[0 0.2],...      %range of the parameter space for lambda
% %     'N',1,...                 %number of lambda values. If this value is set to 1, then the first element of lambda.limits would be the assumed lambda and the lambda parameter is not estimated.
% %     'scale','lin',...         %the linear or log spacing. Choose between 'lin' and 'log'.
% %     'dist','flat',...         %prior distribution of the lambda parameter. Choose between 'norm' and 'flat'.
% %     'mu',0,...                %mean of the prior distribution.
% %     'std',0.1 ...             %standard deviation of the prior distribution.  
% %     );
% % 
% % Note that, when 'scale' is set to 'lin' and 'dist' is set to 'norm', the 
% % prior distribution is constructed based on:
% % p0 ~ N(mu,std)
% % On the other hand, when 'scale' is set to 'lin' and 'dist' is set to 
% % 'norm', the prior distribution is constructed based on:
% % log10(p0) ~ N(mu,std);
% ------------------------------------------------------------------------------------
%
% A UML object holds the following variables:
%     uml.x         %the signal strength
%     uml.xnext     %the next stimulus level based on the collected data
%     uml.p         %the posterior parameter distribution
%     uml.phi0      %the parameters to be estimated in the weibull psychometric function
%     uml.phi       %the space of all possible parameter combinations
%     uml.swpts     %the sweet points
%     uml.r         %the listener's responses
%     uml.PC        %the current proportion of positive responses
%     uml.par       %the data structure containing the configurations of the parameter space   
%     uml.n         %the trial number
%     uml.psycfun   %the function handle for the psychometric-function model
%     uml.alpha     %the array of potential alpha values in the parameter space
%     uml.beta      %the array of potential beta values in the parameter space
%     uml.gamma     %the value of gamma
%     uml.lambda    %the array of potential lambda values in the parameter space
%
% The following methods are available for the UML class:
%     uml.reset(par);           % clear the data and reset the parameter distributions
%     uml.setX0(x0);            % set the initial signal strength
%     uml.update(r);            % update the signal strength and other variables based on a new response r
%     conf = uml.getConf(prct); % return the Bayesian confidence limits for paramters (e.g., conf = uml.getConf([0.025 0.975])).
% 
%     uml.setPhi0(phi0);        % the true psychometric function parameters for a virtual observal, used for simulations
%     r = uml.simulateResponse(x); % simulate a response from the virtual observer
    
    properties
        %basic properties
        x   %the signal strength
        xnext % the next signal strength based on the previous data;
        p   %the posterior parameter distribution
        phi0 = [0.2,1,0.5,0.1]; %the parameters to be estimated in the weibull psychometric function
        phi %the space of all possible parameter combinations
        swpts % the sweet points
        r   %the listener's responses [1 - positive (higher, louder, correct), 0 - negative (lower, weaker, incorrect)]
        PC
        par %the data structure containing the configurations of the parameter space   
        n=0;   %the trial number
        psycfun %the function handle for the psychometric-function model
        
        alpha
        beta
        gamma
        lambda
        
        % hidden parameters used by the algorithm (flags, counters, and indices)
        swpts_idx %the indexes of the valid sweet points
        down_step_flag %flag is used to count the previous positive reponses, before making a negative reponse  
        up_step_flag %flag is used to count the previous negative reponses, before making a positive reponse 
        track_direction % whether the track is going up or down
        current_step %the index of the current sweet point
        rev_flag % flag for the reversals in the track
        nrev % the total number of reversal points
        a
        b
        g
        l
        
        % variable names reserved for user defined data
        userdata01
        userdata02
        userdata03
        userdata04
    end
    
    methods
        
        %constructor
        function uml = UML(par)
            reset(uml, par);
        end  
        
        % INITIALIZATION
        
        %set parameter ndown used in the updown algorithm
        function reset(uml, par)
            uml.par = par;
            uml.setP0();
            uml.x = [];
            uml.r = [];
            uml.phi = [];
            uml.n = 0;
            
            %initlize the parameters in updown algorithm
            if strcmp(par.model,'logit')
                uml.psycfun = @(x,a,b,g,l) g+(1-g-l).*(1+exp(-(x-a).*b)).^(-1);
            elseif strcmp(par.model,'weibull')
                if min(uml.par.alpha.limits)<=0 || min(uml.par.beta.limits)<1
                    error('The lower bound of alpha has to be larger than zero, and the lower bound of beta has to be larger than one.');
                end
                uml.psycfun = @(x,a,b,g,l) g+(1-g-l).*(1-exp(-(x./a).^b));
            elseif strcmp(par.model,'gaussian')
                uml.psycfun = @(x,a,b,g,l) g+(1-g-l).*(1+erf((x-a)./sqrt(2*b.^2)))/2;
            else
                error('Please choose among "logit", "weibull", and "gaussian" for "par.model".');
            end
                
            uml.swpts_idx = [];
            if uml.par.alpha.N > 1
                uml.swpts_idx(end+1) = 3;
            end
            if uml.par.beta.N > 1
                uml.swpts_idx(end+1:end+2) = [2 4];
            end
            if uml.par.gamma.N > 1
                uml.swpts_idx(end+1) = 1;
            end
            if uml.par.lambda.N > 1
                uml.swpts_idx(end+1) = 5;
            end
            uml.swpts_idx = sort(uml.swpts_idx);
            uml.down_step_flag = 0;
            uml.up_step_flag = 0;
            uml.track_direction = 0;
            uml.rev_flag = [];
            uml.nrev = 0;
            
            % set the first sweet point location
            if uml.par.x0 < mean(uml.par.x_lim)
                uml.current_step = min(uml.swpts_idx);
            elseif uml.par.x0 >= mean(uml.par.x_lim)
                uml.current_step = max(uml.swpts_idx);
            end
            uml.xnext = uml.par.x0;
                
        end
        
        %set the hypothesis made on the parameters to be estimated
        function setP0(uml)
            %set the space for phi
            uml.alpha = setParSpace(uml.par.alpha);
            uml.beta = setParSpace(uml.par.beta);
            uml.gamma = setParSpace(uml.par.gamma);
            uml.lambda = setParSpace(uml.par.lambda);
            
            [uml.a,uml.b,uml.g, uml.l] = ndgrid(uml.alpha, uml.beta, uml.gamma, uml.lambda);
            
            %set the prior value and the space for hypo phi
            A = setPrior(uml.a, uml.par.alpha);
            B = setPrior(uml.b, uml.par.beta);
            G = setPrior(uml.g, uml.par.gamma);
            L = setPrior(uml.l, uml.par.lambda);

            uml.p = log(prepare_prob(A.*B.*G.*L));
        end
        
        %set the initial stimulus level used in the algorithm
        function setX0(uml, x0)
            uml.par.x0 = x0;
            % set the first sweet point location
            if uml.par.x0 < mean(uml.par.x_lim)
                uml.current_step = min(uml.swpts_idx);
            elseif uml.par.x0 >= mean(uml.par.x_lim)
                uml.current_step = max(uml.swpts_idx);
            end
            uml.xnext = uml.par.x0;
            uml.down_step_flag = 0;
            uml.up_step_flag = 0;
            uml.track_direction = 0;
            uml.rev_flag = [];
            uml.nrev = 0;
        end
        
        % ITERATION
        
        %Update likelihood based on the previous likelihood and the
        %new response r. Find the according new estimated k, beta, lambda based on the new
        %likelihood. Return the next stimulus level.
        function update(uml, r) 
            
            nresps = length(r);
            uml.n = uml.n + 1;
            uml.x(end+1,:) = uml.xnext;
            uml.r(end+1,:) = sum(r)+1j*nresps;
            for iresp = 1:nresps
                %Update likelihood based on the previous likelihood and the
                %new response r. Find new alpha, beta, lambda based on the new
                %likelihood.
                update_posterior(uml, uml.xnext, r(iresp));
            end
            uml.phi(end+1,:) = get_phi(uml);
            
            %Find the next signal strength at the appropriate sweet point
            if strcmp(uml.par.model,'logit')
                swpt = logit_sweetpoints(uml.phi(end,:));  % note that the lambda sweet point is always at the maximum of the stimulus space
            elseif strcmp(uml.par.model,'weibull')
                swpt = weibull_sweetpoints(uml.phi(end,:));  % note that the lambda sweet point is always at the maximum of the stimulus space
            elseif strcmp(uml.par.model,'gaussian')
                swpt = gaussian_sweetpoints(uml.phi(end,:));  % note that the lambda sweet point is always at the maximum of the stimulus space
            else
                error('Choose among "logit", "weibull", and "gaussian" for "uml.par.model".');
            end 
            
            swpt_1 = 2*swpt(1)-swpt(2);     % defining the gamma sweet point so that the gamma, lower-beta, and alpha sweet points are evenly spaced along the signal strength axis
            swpt_5 = 2*swpt(3)-swpt(2);     % defining the lambda sweet point so that the alpha, upper-beta, and lambda sweet points are evenly spaced along the signal strength axis
            swpt = [swpt_1, swpt, swpt_5];
            swpt = max(min(swpt,uml.par.x_lim(2)),uml.par.x_lim(1)); % limit the sweet points to be within the stimulus space

            %Standard n-down, 1-up algorithm
            uml.rev_flag(end+1,:) = 0;
            if mean(r)>0.5   % positive
                if uml.down_step_flag == uml.par.ndown-1 % reached ndown
                    uml.current_step = max(uml.current_step-1,min(uml.swpts_idx));
                    uml.xnext = swpt(uml.current_step);
                    uml.down_step_flag = 0;
                    if uml.track_direction == 1
                        uml.rev_flag(end,:) = 1;
                        uml.nrev = uml.nrev + 1;
                    end
                    uml.track_direction = -1;
                elseif uml.down_step_flag <uml.par.ndown-1 % has not reached ndown
                    uml.xnext = swpt(uml.current_step);
                    uml.down_step_flag = uml.down_step_flag+1;
                end
                uml.up_step_flag = 0;
            elseif mean(r)<=0.5    % negative
                if uml.up_step_flag == uml.par.mup-1 % reached mup
                    uml.current_step = min(uml.current_step+1,max(uml.swpts_idx));
                    uml.xnext = swpt(uml.current_step);
                    uml.up_step_flag = 0;
                    if uml.track_direction == -1
                        uml.rev_flag(end,:) = 1;
                        uml.nrev = uml.nrev + 1;
                    end
                    uml.track_direction = 1;
                elseif uml.up_step_flag < uml.par.mup-1 % has not reached mup
                    uml.xnext = swpt(uml.current_step);
                    uml.up_step_flag = uml.up_step_flag+1;
                end
                uml.down_step_flag = 0;
            end
            uml.swpts(end+1,:) = swpt;
            
            % update percent positive
            uml.PC = sum(real(uml.r))/sum(imag(uml.r));
            
        end
        
        function update_posterior(uml, x, r)
            %Update likelihood based on the previous likelihood and the
            %new response r. Find new alpha, beta, lambda based on the new
            %likelihood.
            
            uml.p = uml.p +...
                log(prepare_prob(uml.psycfun(x,uml.a, uml.b, uml.g, uml.l)).^r) + ...
                log(prepare_prob(1-uml.psycfun(x,uml.a, uml.b, uml.g, uml.l)).^(1-r));

            uml.p = uml.p-max(max(max(max(uml.p))));
            
        end
        
        function phi = get_phi(uml)
            if strcmp(uml.par.method, 'mode')
                [~,idx] = max(reshape(uml.p,numel(uml.p),1));
                phi = [uml.a(idx), uml.b(idx), uml.g(idx), uml.l(idx)];
            elseif strcmp(uml.par.method, 'mean')
                pdf_tmp = exp(uml.p);
                pdf_tmp = pdf_tmp/sum(sum(sum(sum(pdf_tmp))));
                % alpha
                alpha_est_tmp = sum(sum(sum(sum(pdf_tmp.*uml.a))));
                % beta
                beta_est_tmp = sum(sum(sum(sum(pdf_tmp.*uml.b))));
                % gamma
                gamma_est_tmp = sum(sum(sum(sum(pdf_tmp.*uml.g))));
                % lambda
                lambda_est_tmp = sum(sum(sum(sum(pdf_tmp.*uml.l))));
                % combine together
                phi = [alpha_est_tmp, beta_est_tmp, gamma_est_tmp, lambda_est_tmp];
            else
                error('Choose between "mean" and "mode" for "uml.par.method".');
            end
        end
        
        % SIMULATION
        
        %set the parameters in the psychometric function of the virtual
        % observer for simulation
        function setPhi0(uml, phi0)
            uml.phi0 = phi0;
        end
        
        %given a stimulus, get a simulated virtual response
        function r = simulateResponse(uml, x)
            %given a stimulus, get a response
            r = rand<uml.psycfun(x,uml.phi0(1),uml.phi0(2),uml.phi0(3),uml.phi0(4));
        end
        
        % RETRIEVAL AND PLOTTING
        
        % get Bayesian confidence limits (credible limits) for each
        % parameter
        function conf = getConf(uml, prct)
            
            ptmp = exp(uml.p)/sum(sum(sum(sum(exp(uml.p)))));
            conf = zeros(length(prct),4);
            for i = 1:length(prct)
                if length(uml.beta) == 1 && length(uml.lambda) == 1 && length(uml.gamma) == 1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = uml.beta;
                    gamma_conf = uml.gamma;
                    lambda_conf = uml.lambda;
                elseif length(uml.beta) == 1 && length(uml.gamma)==1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = uml.beta;
                    gamma_conf = uml.gamma;
                    lambda_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),3))),uml.lambda,prct(i),'linear','extrap');
                elseif length(uml.beta) == 1 && length(uml.lambda)==1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = uml.beta;
                    gamma_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),4))),uml.gamma,prct(i),'linear','extrap');
                    lambda_conf = uml.lambda;
                elseif length(uml.gamma) == 1 && length(uml.lambda)==1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = interp1(cumtrapz(sum(sum(sum(ptmp,1),3),4)),uml.beta,prct(i),'linear','extrap');
                    gamma_conf = uml.gamma;
                    lambda_conf = uml.lambda;
                elseif length(uml.beta) == 1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = uml.beta;
                    gamma_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),4))),uml.gamma,prct(i),'linear','extrap');
                    lambda_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),3))),uml.lambda,prct(i),'linear','extrap');
                elseif length(uml.gamma) == 1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = interp1(cumtrapz(sum(sum(sum(ptmp,1),3),4)),uml.beta,prct(i),'linear','extrap');
                    gamma_conf = uml.gamma;
                    lambda_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),3))),uml.lambda,prct(i),'linear','extrap');
                elseif length(uml.lambda) == 1
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = interp1(cumtrapz(sum(sum(sum(ptmp,1),3),4)),uml.beta,prct(i),'linear','extrap');
                    gamma_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),4))),uml.gamma,prct(i),'linear','extrap');
                    lambda_conf = uml.lambda;
                else
                    alpha_conf = interp1(cumtrapz(sum(sum(sum(ptmp,2),3),4)),uml.alpha,prct(i),'linear','extrap');
                    beta_conf = interp1(cumtrapz(sum(sum(sum(ptmp,1),3),4)),uml.beta,prct(i),'linear','extrap');
                    gamma_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),4))),uml.gamma,prct(i),'linear','extrap');
                    lambda_conf = interp1(squeeze(cumtrapz(sum(sum(sum(ptmp,1),2),3))),uml.lambda,prct(i),'linear','extrap');
                end
                conf(i,:) = [alpha_conf, beta_conf, gamma_conf, lambda_conf];
            end
        end       
    end
end

% Some service functions

function space = setParSpace(s)

%set parameter space on a linear scale, or a logarithmic scale

    if s.N == 1
        space = s.limits(1);
    else
        if strcmp(s.scale, 'lin')
            space = linspace(s.limits(1),s.limits(2),s.N)';
        elseif strcmp(s.scale, 'log')
            space = logspace(log10(s.limits(1)),log10(s.limits(2)),s.N)';
        end          
    end
end

function p0 = setPrior(phi, s)           

%set the prior distribution

    if strcmp(s.dist, 'norm')
        if strcmp(s.scale, 'lin')
            p0 = normpdf(phi,s.mu,s.std);
        elseif strcmp(s.scale, 'log')
            p0 = normpdf(log10(phi),s.mu,s.std);
        end
    elseif strcmp(s.dist, 'flat')
        p0 = ones(size(phi));
    end
end

% eliminate the 1 or 0 probability to prevent the nummerical errors due to
% taking log of zeros.
function p = prepare_prob(p)
    p = p*(1-1e-10);
    p = p/sum(sum(sum(sum(p))));
end

% Find the sweet points for logistic psychometric function
function swpts = logit_sweetpoints(phi)

%calculate the sweet points for a logit psychometric function, for
%parameter alpha, beta and lambda.

alpha = phi(1);
beta = phi(2);
gamma = phi(3);
lambda = phi(4);

swpts(1) = fminsearch(@(x) betavar_est(x,alpha,beta,gamma,lambda)+(x>=alpha)*1e10,alpha-10);
swpts(3) = fminsearch(@(x) betavar_est(x,alpha,beta,gamma,lambda)+(x<=alpha)*1e10,alpha+10);
swpts(2) = fminsearch(@(x) alphavar_est(x,alpha,beta,gamma,lambda),alpha);
swpts = sort(swpts);

    function sigmaalphasq = alphavar_est(x,alpha,beta,gamma,lambda)

    term1 = exp(2*beta*(alpha-x));
    term2 = (1+exp(beta*(x-alpha))).^2;
    term3 = -gamma+(lambda-1)*exp(beta*(x-alpha));
    term4 = 1-gamma+lambda*exp(beta*(x-alpha));
    term5 = beta^2*(-1+gamma+lambda)^2;

    sigmaalphasq = -term1.*term2.*term3.*term4./term5;
    end

    function sigmabetasq = betavar_est(x,alpha,beta,gamma,lambda)

    term1 = exp(2*beta*(alpha-x));
    term2 = (1+exp(beta*(x-alpha))).^2;
    term3 = -gamma+(lambda-1)*exp(beta*(x-alpha));
    term4 = 1-gamma+lambda*exp(beta*(x-alpha));
    term5 = (x-alpha).^2*(-1+gamma+lambda)^2;

    sigmabetasq = -term1.*term2.*term3.*term4./term5;
    end

end


% Find the sweet points for weibull psychometric function
function swpts = weibull_sweetpoints(phi)

%calculate the sweet points for a weibull psychometric function, for
%parameter k, beta and lambda.

k = phi(1);
beta = phi(2);
gamma = phi(3);
lambda = phi(4);

opt = optimset('Display','off');
swpts(1) = fminsearch(@(x) betavar_est(x,k,beta,gamma,lambda)+(x>=k)*1e10,k/2,opt);
swpts(3) = fminsearch(@(x) betavar_est(x,k,beta,gamma,lambda)+(x<=k)*1e10,k*2,opt);
swpts(2) = fminsearch(@(x) kvar_est(x,k,beta,gamma,lambda),k,opt);
swpts = max(sort(swpts),0);

    function sigmaksq = kvar_est(x,k,beta,gamma,lambda)

    term1 = k.^2.*(x./k).^(-2*beta);
    term2 = -1+gamma-exp((x./k).^beta).*(-1+lambda)+lambda;
    term3 = -1+gamma+lambda-exp((x./k).^beta).*lambda;
    term4 = beta.^2.*(-1+gamma+lambda).^2;

    sigmaksq = -term1.*term2.*term3./term4;

    end

    function sigmabetasq = betavar_est(x,k,beta,gamma,lambda)

    term1 = (x./k).^(-2*beta);
    term2 = -1+gamma-exp((x./k).^beta).*(-1+lambda)+lambda;
    term3 = -1+gamma+lambda-exp((x./k).^beta).*lambda;
    term4 = (-1+gamma+lambda).^2.*log(x./k).^2;

    sigmabetasq = -term1.*term2.*term3./term4;

    end

end


% Find the sweet points for gaussian psychometric function
function swpts = gaussian_sweetpoints(phi)

    % Return the sweet points for the mu and sigma parameters for a cumulative
    % gaussian psychometric function. The sweet points are locations on the
    % psychometric function corresponding to minimum expected variablity in the
    % parameter estimates.
    % Numerical search for the minimum is implemented here using fminsearch.
    
    sigma2_mu = @(x) psycfunc(x,phi).*(1-psycfunc(x,phi))./...
        (psycfunc_derivative_mu(x,phi)).^2;

    sigma2_sigma = @(x) psycfunc(x,phi).*(1-psycfunc(x,phi))./...
        (psycfunc_derivative_sigma(x,phi)).^2;

    swpt_mu = fminsearch(sigma2_mu, phi(1));
    swpt_sigma_L = fminsearch(sigma2_sigma, phi(1)-10);
    swpt_sigma_H = fminsearch(sigma2_sigma, phi(1)+10);

    swpts = [swpt_sigma_L, swpt_mu, swpt_sigma_H];

    function p = psycfunc(x,phi)
        mu = phi(1);
        sigma = phi(2);
        gamma = phi(3);
        lambda = phi(4);
        p = gamma+((1-gamma-lambda)/2).*(1+erf((x-mu)./sqrt(2*sigma.^2)));
    end
    function dpdm = psycfunc_derivative_mu(x,phi)
        mu = phi(1);
        sigma = phi(2);
        gamma = phi(3);
        lambda = phi(4);
        dpdm = -(1-gamma-lambda)./(sqrt(2*pi)*sigma).*exp(-(x-mu).^2./(2*sigma^2));
    end
    function dpds = psycfunc_derivative_sigma(x,phi)
        mu = phi(1);
        sigma = phi(2);
        gamma = phi(3);
        lambda = phi(4);
        dpds = -(1-gamma-lambda).*(x-mu)./(sqrt(2*pi)*sigma.^2).*exp(-(x-mu).^2./(2*sigma^2));
    end

end

%eof