function TE = detect_spikes_wavelet(...
    Signal, SFr, Wid, Ns, option, L, wname, PltFlg, CmtFlg)

% DETECT_SPIKES_WAVELET wavelet based algorithm for detection of transients
% from neural data.
%
%   TE=DETECT_SPIKES_WAVELET(Signal,SFr,Wid,Ns,option,L,wname,PltFlg,CmtFlg)
%
%   Signal - extracellular potential data to be analyzed 1 x Nt;
%
%   SFr - sampling frequency [kHz];
%
%   Wid - 1 x 2 vector of expected minimum and maximum width [msec] of transient 
%   to be detected Wid=[Wmin Wmax]. For most practical purposes Wid=[0.5 1.0];
%
%   Ns - (scalar): the number of scales to use in detection (Ns >= 2);
%
%   option - (string): the action taken when no coefficients survive hard 
%   thresholding 
%   'c' means conservative and returns no spikes if P(S) is found to be 0
%   'l' means assume P(S) as a vague prior (see the original reference)
%
%   L is the factor that multiplies [cost of comission]/[cost of omission].
%   For most practical purposes -0.2 <= L <= 0.2. Larger L --> omissions
%   likely, smaller L --> false positives likely. For unsupervised
%   detection, the suggested value of L is close to 0.  
%
%   wname - (string): the name of wavelet family in use
%           'bior1.5' - biorthogonal
%           'bior1.3' - biorthogonal
%           'db2'     - Daubechies
%           'sym2'    - symmlet
%           'haar'    - Haar function
%   Note: sym2 and db2 differ only by sign --> they produce the same
%   result;
%
%   PltFlg - (integer) is the plot flag: 
%   PltFlg = 1 --> generate figures, otherwise do not;
%  
%   CmtFlg - (integer) is the comment flag, 
%   CmtFlg = 1 --> display comments, otherwise do not;
%
%   TE is the vector of event occurrence times;
%
%   Reference: Z. Nenadic and J.W. Burdick, Spike detection using the 
%   continuous wavelet transform, IEEE T. Bio-med. Eng., vol. 52,
%   pp. 74-87, 2005.

%   Originally developed by:
%   Zoran Nenadic
%   California Institute of Technology
%   May 2003
%
%   Modified by:
%   Zoran Nenadic
%   University fo California, Irvine
%   February 2008


%admissible wavelet families (more wavelets could be added)
wfam = {'bior1.5','bior1.3','sym2','db2','haar'};

if sum(strcmp(wname,wfam)) == 0
    error('unknown wavelet family')
elseif CmtFlg == 1
    disp(['wavelet family: ' wname])
    to = clock;
end

%make sure signal is zero-mean
Signal = Signal - mean(Signal);

Nt = length(Signal);      %# of time samples

%define relevant scales for detection
W = determine_scales(wname,Wid,SFr,Ns);

%initialize the matrix of thresholded coefficients
ct = zeros(Ns,Nt);

%get all coefficients 
c = cwt(Signal,W,wname);  

%define detection parameter
Lmax = 36.7368;       %log(Lcom/Lom), where the ratio is the maximum 
                    %allowed by the current machine precision
L = L * Lmax;

%initialize the vector of spike indicators, 0-no spike, 1-spike
Io = zeros(1,Nt);

%loop over scales
for i = 1:Ns
    
    %take only coefficients that are independent (W(i) apart) for median
    %standard deviation
    
    Sigmaj = median(abs(c(i,1:round(W(i)):end) - mean(c(i,:))))/0.6745;
    Thj = Sigmaj * sqrt(2 * log(Nt));     %hard threshold
    index = find(abs(c(i,:)) > Thj);
    if isempty(index) & strcmp(num2str(option),'c')
        %do nothing ct=[0];
    elseif isempty(index) & strcmp(num2str(option),'l')
        Mj = Thj;
        %assume at least one spike
        PS = 1/Nt;
        PN = 1 - PS;
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];    %decision threshold
        DTh = abs(DTh) * (DTh >= 0);                 %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        if isempty(ind)
            %do nothing ct=[0];
        else
            ct(i,ind) = c(i,ind);
        end
    else
        Mj = mean(abs(c(i,index)));       %mean of the signal coefficients
        PS = length(index)/Nt;            %prior of spikes
        PN = 1 - PS;                        %prior of noise
        DTh = Mj/2 + Sigmaj^2/Mj * [L + log(PN/PS)];   %decision threshold
        DTh = abs(DTh) * (DTh >= 0);         %make DTh>=0
        ind = find(abs(c(i,:)) > DTh);
        ct(i,ind) = c(i,ind);
    end
    
    %find which coefficients are non-zero
    Index = ct(i,:) ~= 0;
    
    %make a union with coefficients from previous scales
    Index = or(Io,Index);
    Io = Index;
end

TE = parse(Index,SFr,Wid);

if PltFlg == 1
    close all
    figure(1)
    scale = 64./[max(abs(c),[],2) * ones(1,Nt)];
    temp = zeros(1,Nt);
    temp(TE) = 1;
    image(flipud(abs(c)) .* scale)
    colormap pink
    ylabel('Scales')
    Wt = [fliplr(W)];
    set(gca,'YTick',1:Ns,'YTickLabel',Wt,'Position',[0.1 0.2 0.8 0.6], ...
        'XTick',[])
    title(['|C| across scales: ' num2str(W)])
    ah2 = axes;
    set(ah2,'Position',[0.1 0.1 0.8 0.1])
    plot(temp,'o-m','MarkerSize',4,'MarkerFaceColor','m')
    set(gca,'YTick',[],'XLim',[1 Nt])
    xlabel('Time (samples)')
    ylabel('Spikes')
    
    figure(2)
    plot(Signal,'Color',[0.7 0.7 0.7],'LineWidth',2)
    hold on
    plot(ct','-o','LineWidth',1,'MarkerFaceColor','k', ...
        'MarkerSize',4)
    xlabel('Time (samples)')
    ylabel('Coefficients')
    set(gca,'XLim',[1 Nt])
end

if CmtFlg == 1
    disp([num2str(length(TE)) ' spikes found'])
    disp(['elapsed time: ' num2str(etime(clock,to))])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Scale = determine_scales(wname,Wid,SFr,Ns)
  
%Ns - # of scales  

dt = 1/SFr;  %[msec]

%signal sampled @ 1 KHz  
Signal = zeros(1,1000);
%create Dirac function
Signal(500) = 1;
  
Width = linspace(Wid(1),Wid(2),Ns);

%infinitesimally small number
Eps = 10^(-15);

ScaleMax = 5;
ScaleMax = ScaleMax*SFr;

switch num2str(wname)
  
 case 'haar'
  for i = 1:Ns
    Scale(i) = Width(i)/dt - 1; 
  end
 case 'db2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'sym2'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of positive slope zero crossings
    IndZeroCross = find(IndDer == 1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.3'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
   for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 case 'bior1.5'
  Scales = 2:ScaleMax;
  c = cwt(Signal,Scales,wname);
  for i = 1:length(Scales)
    %indicators of positive coefficients
    IndPos = (c(i,:) > 0);
    %indicators of derivative
    IndDer = diff(IndPos);
    %indices of negative slope zero crossings
    IndZeroCross = find(IndDer == -1);
    IndMax = IndZeroCross > 500;
    Ind(2) = min(IndZeroCross(IndMax))+1;
    IndMin = IndZeroCross < 500;
    Ind(1) = max(IndZeroCross(IndMin));
    WidthTable(i) = diff(Ind) * dt;
  end
  WidthTable = WidthTable + [1:length(Scales)] * Eps;
  %look-up table
  Scale = round(interp1(WidthTable,Scales,Width,'linear'));
 otherwise
  error('unknown wavelet family')
end

NaNInd = isnan(Scale);

if sum(NaNInd) > 0
  warning(['Your choice of Wid is not valid given' ...
        ' the sampling rate and wavelet family'])
  if NaNInd(1) == 1
    disp(['Most likely Wid(1) is too small'])
  elseif NaNInd(Ns) == 1
    disp(['Most likely Wid(2) is too large'])
    disp(['Change the value on line: ''ScaleMax = 2'' to something larger'])
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn = parse(Index,SFr,Wid);

%This is a special function, it takes the vector Index which has 
%the structure [0 0 0 1 1 1 0 ... 0 1 0 ... 0]. This vector was obtained
%by coincidence detection of certain events (lower and upper threshold
%crossing for threshold detection, and the appearance of coefficients at
%different scales for wavelet detection). 
%The real challenge here is to merge multiple 1's that belong to the same
%spike into one event and to locate that event

Refract = 1.5 * Wid(2);    %[ms] the refractory period -- can't resolve spikes 
                           %that are closer than Refract;
Refract = round(Refract * SFr);

Merge = mean(Wid);      %[ms] merge spikes that are closer than Merge, since 
                        %it is likely they belong to the same spike

Merge = round(Merge * SFr);   


Index([1 end]) = 0;   %discard spikes located at the first and last samples

ind_ones = find(Index == 1);    %find where the ones are

if isempty(ind_ones)
    TE = [];
else
    temp = diff(Index);  %there will be 1 followed by -1 for each spike
    N_sp = sum(temp == 1); %nominal number of spikes
    
    lead_t = find(temp == 1);  %index of the beginning of a spike
    lag_t = find(temp == -1);  %index of the end of the spike
    
    for i = 1:N_sp
        tE(i) = ceil(mean([lead_t(i) lag_t(i)]));
    end
   
    i = 1;        %initialize counter
    while 0 < 1
        if i > (length(tE) - 1)
            break;
        else
            Diff = tE(i+1) - tE(i);
            if Diff < Refract & Diff > Merge
                tE(i+1) = [];      %discard spike too close to its predecessor
            elseif Diff <= Merge
                tE(i) = ceil(mean([tE(i) tE(i+1)]));  %merge
                tE(i+1) = [];                         %discard
            else
                i = i+1;
            end
        end
    end 
    TE = tE;
end

fcn = TE;