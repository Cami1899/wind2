function varargout = autocorr(varargin)
%AUTOCORR Sample autocorrelation
%
% Syntax:
%
%   [acf,lags] = autocorr(y)
%   ACFTbl = autocorr(Tbl)
%   [...,bounds] = autocorr(...)
%   [...,bounds,h] = autocorr(...)
%   [...] = autocorr(...,param,val,...)
%   [...] = autocorr(ax,...)
%   autocorr(...)
%
% Description:
%
%   Compute the sample autocorrelation function (ACF) of a univariate 
%   time series. AUTOCORR optionally plots the ACF with confidence bounds.
%
% Input Arguments:
%
%   y - Univariate time series data, specified as a numeric vector.
%
%   Tbl - Time series data, specified as a table or timetable. Specify a
%       single series for y using the 'DataVariable' parameter.
%
%   ax - Axes object in which to plot. If unspecified, AUTOCORR plots to
%       the current axes (gca).
%
% Optional Input Parameter Name/Value Arguments:
%
%  'NumLags' Positive integer that determines the number of lags at which the
%            ACF is computed. The lags used to compute the ACF are 0:NumLags.
%            The default is min[20,N-1], where N is the effective sample size
%            of y.
%
%  'NumMA'   For computing confidence bounds, a nonnegative integer less
%            than NumLags specifying the number of lags in a theoretical
%            MA(NumMA) model of y. For lags > NumMA, AUTOCORR uses
%            Bartlett's approximation [1] to compute the standard error
%            under the model assumption. The default is 0, in which case
%            the standard error is 1/sqrt(N), for Gaussian white noise.
%
%  'NumSTD'  For computing confidence bounds, a nonnegative scalar multiple
%            specifying an interval of +/-(NumSTD) times the computed
%            standard error. The default is 2 (approximate 95% confidence).
%
%  'DataVariable' Variable in Tbl to use for y, specified as a name in
%            Tbl.Properties.VariableNames. Variable names are character
%            vectors, string scalars, integers or logical vectors. The
%            default is the last variable in Tbl.
%
% Output Arguments:
%
%   acf - Sample ACF. Vector of length NumLags+1 of values computed at lags
%       0,1,2,...,NumLags. For all y, acf(1) = 1 at lag 0.
%
%   lags - Vector of lag numbers of length NumLags+1 used to compute acf.
%
%   ACFTbl - When input is Tbl, outputs lags and acf are returned in
%       table ACFTbl.
%
%   bounds - Two-element vector of approximate upper and lower confidence
%       bounds, assuming that y is an MA(NumMA) process.
%
%   h - Vector of handles to plotted graphics objects. AUTOCORR plots the
%       ACF when the number of output arguments is 0 or 4.
%
% Notes:
%
%  o Specify missing observations of y using NaN. AUTOCORR treats these
%    values as "missing completely at random."
%
%  o If y is fully observed, without NaNs, AUTOCORR uses a Fourier
%    transform to compute the ACF in the frequency domain, then converts
%    back to the time domain using an inverse Fourier transform.
%
%  o In the presence of NaNs, AUTOCORR computes the ACF at lag k in the
%    time domain, including in the sample average only those terms for
%    which the cross product y(t)*y(t+k) exists, so that the effective
%    sample size at any lag is a random variable.
%
% Example:
%
%   % Create an MA(2) process from a sequence of 1000 Gaussian deviates,
%   % and assess whether the ACF is effectively zero for lags > 2:
%
%   x = randn(1000,1);         % 1000 Gaussian deviates ~ N(0,1)
%   y = filter([1 -1 1],1,x);  % Create an MA(2) process
%   autocorr(y,'NumMA',2)      % Inspect the ACF with 95% confidence
%
% References:
%
%   [1] Box, G. E. P., G. M. Jenkins, and G. C. Reinsel. Time Series
%       Analysis: Forecasting and Control. 3rd edition. Upper Saddle River,
%       NJ: Prentice-Hall, 1994.
%
%   [2] Hamilton, J.D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
% See also PARCORR, CROSSCORR, FILTER.

% Copyright 2022 The MathWorks, Inc.   

% Preprocess varargin for target axes:

try
    
    [ax,args] = internal.econ.axesparser(varargin{:});
    
catch ME

	throwAsCaller(ME)

end

% This function produces a single plot:

if ~isempty(ax) && ~isscalar(ax)
    
    error(message('econ:internal:econ:axesparser:InvalidParent'));
    
end

% Parse inputs and set defaults:

Data = args{1};
args = args(2:end);
isTabular = istable(Data) || istimetable(Data);

if ~isTabular && (numel(args) > 0) && ...
   isnumeric(args{1}) % Deprecated positional syntax

    % Deprecated positional syntax:
    % [...] = autocorr(y,numLags,numMA,numSTD)

    % Parse:
    
    try

        parser = checkPositionalInputs(Data,args{:});

    catch ME

        throwAsCaller(ME)

    end

    y = parser.Results.Data;
    rowSeries = (size(y,1) == 1);
	numLags = parser.Results.NumLags;
   
else % Documented syntax
        
    % Construct parser:

    parser = inputParser;
    parser.addRequired('Data',...
                       @(x)validateattributes(x,{'double','table','timetable'},{'nonempty','2d'}));
    parser.addParameter('NumLags',20,...
                        @(x)validateattributes(x,{'double'},{'scalar','integer','>',0}));
    parser.addParameter('NumMA',0,...
                        @(x)validateattributes(x,{'double'},{'scalar','integer','>=',0}));
    parser.addParameter('NumSTD',2,...
                        @(x)validateattributes(x,{'double'},{'scalar','>=',0}));
    parser.addParameter('DataVariable',[],...
                        @(x)validateattributes(x,{'double','logical','char','string'},{'vector'}));

    % Parse:

    try
    
        parser.parse(Data,args{:});
  
    catch ME
    
        throwAsCaller(ME)
  
    end
  
    varSpec = parser.Results.DataVariable;

    % Select y with 'DataVariable':

    if isnumeric(Data)

        y = Data;
        rowSeries = (size(y,1) == 1);
        
        if isvector(y) && ~isempty(varSpec)

            warning(message('econ:autocorr:DataVariableUnused'))
        
        end

    else % Tabular data

        if ~isempty(varSpec)

            try

                y = Data(:,varSpec);

            catch ME

                throwAsCaller(ME)

            end

        else

            y = Data(:,end); % Default

        end
        
        try

            internal.econ.TableAndTimeTableUtilities.isTabularFormatValid(y,'y')
            internal.econ.TableAndTimeTableUtilities.isTabularDataSinglePath(y,'y')

       	catch ME

            throwAsCaller(ME)

        end
        
     	y = table2array(y);
        y = double(y);

    end
    
    if ~isvector(y)

    	error(message('econ:autocorr:NonVectorInput'))

   	end
    
   	% Set numLags:

    numLags = parser.Results.NumLags;

    if any(strcmpi('NumLags',parser.UsingDefaults)) % Default lags
        
        numLags = min(numLags,sum(~isnan(y))-1);
        
    else % User-specified lags
        
        if numLags > (sum(~isnan(y))-1)
          
            error(message('econ:autocorr:LagsTooLarge'))
        
        end
      
    end
    
end

% Preprocess validated inputs:

y = y(:);
y = y - mean(y,"omitnan");
N = sum(~isnan(y)); % Effective sample size

% Set numMA, numSTD:

numMA = parser.Results.NumMA;

if numMA >= numLags
    
	error(message('econ:autocorr:NumMATooLarge'))
    
end

numSTD = parser.Results.NumSTD;

% Compute ACF:

if N < length(y) % Missing data

% Compute ACF in the presence of missing data.
%
%   To compute the jth autocovariance, compute the cross product of y(t)
%	and y(t+j) with sample size (T) adjusted for missing data. The 
%	autocovariances computed in the presence of missing data follow
%	expression 2.1.10 on page 31 of [1].
%
%	The procedure assumes data is missing at random. Asymptotically, if P
%	is the probability that any element of y(t) is not missing, then for
%	all j > 0:
%
%	sum(~isnan(cross))/length(y) --> P^2,
%
%   or
%
%	length(y)/(T-sum(~isnan(y(1:j)))) --> 1/P^2.
%
%	If the input series y(t) has no missing data, time-domain
%	autocovariances are identical to those in the FFT approach.

	acf = nan(numLags+1,1);

    for j = 0:numLags
        
        cross   = y(1:end-j).*y(j+1:end);
        iNonNaN = ~isnan(cross);

        if any(iNonNaN)
            
            T        = sum(iNonNaN)+sum(~isnan(y(1:j)));
            acf(j+1) = sum(cross,"omitnan")/T;
            
        end
        
    end

else % No missing data

% Compute ACF using expression 2.1.10 on page 31 of [1], but perform the
% computation in the frequency domain using an FFT:

    nFFT = 2^(nextpow2(length(y))+1);
    F = fft(y,nFFT);
    F = F.*conj(F);
    acf = ifft(F);
    acf = acf(1:(numLags+1)); % Retain nonnegative lags
    acf = real(acf);
   
end

acf = acf./acf(1); % Normalize
   
% Compute approximate confidence bounds using the approach in [1],
% equations 2.1.13 and 6.2.2, pp. 33 and 188, respectively:

sigmaNMA = sqrt((1+2*(acf(2:numMA+1)'*acf(2:numMA+1)))/N);  
bounds = sigmaNMA*[numSTD;-numSTD];
lags = (0:numLags)';

% Perform nargout-dependent operations:

if isTabular
    
    nargoutchk(0,3)
    
else

    nargoutchk(0,4)
    
end

% Create plot:

needPlot = (nargout == 0) || ...
           (isTabular && (nargout == 3)) || ...
           (~isTabular && (nargout == 4));

if needPlot
    
    % Plot to gca if no parent axes is specified:

    if isempty(ax)
    
        ax = gca;
    
    end

    % Store NextPlot flag (and restore on cleanup):

    next = get(ax,'NextPlot');
    cleanupObj = onCleanup(@()set(ax,'NextPlot',next));

    % Plot the sample ACF:

    hPlot = stem(ax,lags,acf,'filled','r-o','MarkerSize',4,'Tag','ACF');

    % Plot confidence bounds under the hypothesis that y is an MA(numMA)
    % process. Bartlett's approximation gives an indication of whether the
    % ACF is effectively zero beyond lag numMA. Confidence bounds appear
    % over the ACF only for lags greater than numMA, where the null
    % hypothesis is assumed to hold.
    
    set(ax,'NextPlot','add');
    hBounds = plot(ax,[numMA+0.5 numMA+0.5; numLags numLags],...
                      [bounds([1 1]) bounds([2 2])],'-b',...
                      'Tag','Confidence Bound');
    hXAxis = plot(ax,[0 numLags],[0 0],'-k','Tag','X Axis');
    
    % Return "plot object":

    h = [hPlot;hBounds;hXAxis];
    
    % Modify axes properties conditional on NextPlot flag:
    
    ax.Tag = 'ACFPlot';

    switch next
    
        case {'replace','replaceall'}
    
            grid(ax,'on')
            xlabel(ax,'Lag')
            ylabel(ax,'Sample Autocorrelation')
            title(ax,'Sample Autocorrelation Function')
            ax.YLim(2) = 1;
        
        case {'replacechildren','add'}
        
            % Do not modify axes properties
    
    end

end

% Re-format outputs to conform to row/column orientation of y:

if isTabular
    
    % Outputs are rows of Tbl
    
else
    
    if nargout > 0

        if rowSeries

            acf = acf';
            lags = lags';
            bounds = bounds';

       end

    end

end

% Create output table:

if isTabular
    
    ACFTbl = table(lags,acf,...
                   'VariableNames',{'Lags','ACF'});
    
end   

% Suppress assignment to ans:

if nargout > 0
    
    if isTabular

        if needPlot

            varargout = {ACFTbl,bounds,h};

        else

            varargout = {ACFTbl,bounds};

        end 

    else

        if needPlot

            varargout = {acf,lags,bounds,h};

        else

            varargout = {acf,lags,bounds};

        end

    end

end

% -------------------------------------------------------------------------
function parser = checkPositionalInputs(Data,numLags,numMA,numSTD)

% Ensure the sample data is a vector:

[rows,columns] = size(Data);

if ((rows ~= 1) && (columns ~= 1))
    
	error(message('econ:autocorr:NonVectorInput'))
    
end

N = sum(~isnan(Data)); % Effective sample size

% Ensure numLags is a positive integer or set default:

if (nargin >= 2) && ~isempty(numLags)
    
	if numel(numLags) > 1
       
        error(message('econ:autocorr:NonScalarLags'))
      
	end
   
    if (round(numLags) ~= numLags) || (numLags <= 0)
       
        error(message('econ:autocorr:NonPositiveInteger'))
      
    end
   
    if numLags > (N-1)
       
        error(message('econ:autocorr:LagsTooLarge'))
      
    end
   
else
    
    numLags = min(20,N-1); % Default
   
end

% Ensure numMA is a nonnegative integer or set default:

if (nargin >= 3) && ~isempty(numMA)
    
	if numel(numMA) > 1
       
        error(message('econ:autocorr:NonScalarNumMA'))
        
	end
   
    if (round(numMA) ~= numMA) || (numMA < 0)
       
        error(message('econ:autocorr:NegativeIntegerNumMA'))
      
    end
   
    if numMA >= numLags
       
        error(message('econ:autocorr:NumMATooLarge'))
      
    end
   
else
    
    numMA = 0; % Default
   
end

% Ensure numSTD is a positive scalar or set default:

if (nargin >= 4) && ~isempty(numSTD)
    
	if numel(numSTD) > 1
       
        error(message('econ:autocorr:NonScalarSTDs'))
      
	end
   
    if numSTD < 0
       
        error(message('econ:autocorr:NegativeSTDs'))
      
    end
   
else
    
    numSTD = 2; % Default
   
end

parser.Results.Data    = Data;
parser.Results.NumLags = numLags;
parser.Results.NumMA   = numMA;
parser.Results.NumSTD  = numSTD;