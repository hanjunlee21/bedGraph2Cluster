function bedGraph2Cluster(bedGraphs_Target, bedGraphs_Nontarget, bedGraphs_Control, Outdir, BED_Bin, k, QNorm, Workingdir)
%% bedGraph2Cluster
% e.g., bedGraph2Cluster("bam/RB.WT.filtered.bedgraph,bam/RB.dCDK.filtered.bedgraph", "bam/E2F1.filtered.bedgraph,bam/CTCF.shSCR.filtered.bedgraph,bam/c-Jun.shSCR.filtered.bedgraph", "bam/INPUT.WT.filtered.bedgraph,bam/INPUT.dCDK.filtered.bedgraph", "test_output", "bed/hg19.200bp.bed", "8", "true", "../")
% 
% Required arguments
%     bedGraphs_Target (string): comma-delimited list of bedGraph files to be included during k-means clustering
%     bedGraphs_Nontarget (string): comma-delimited list of bedGraph files to be excluded during k-means clustering
%     bedGraphs_Control (string): comma-delimited list of bedGraph files to be used as controls for peak calling
%     Outdir (string): path to the output directory
%     BED_Bin (string): path to the BED file used for binned bedGraph generation
%     k (string): number of clusters during k-means clustering
%     QNorm (string): whether to perform QNorm normalization ("true": QNorm, "false": CPM)
%
% Optional arguments 
%     Workingdir (string): path to the output directory
% 
%% MIT License
%
% Copyright (c) 2022 Hanjun Lee (MIT/Broad/MGH), Michael S. Lawrence (Broad/MGH)
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%% Normalization
if strcmp(convertCharsToStrings(QNorm),"true")
    QNorm = true;
elseif strcmp(convertCharsToStrings(QNorm),"false")
    QNorm = false;
else
    error('QNorm has an inappropriate value')
end

if ~exist('Workingdir','var')
    Workingdir = strcat(convertCharsToStrings(pwd),"/");
else
    Workingdir = strcat(tostringmatrix(Workingdir),"/");
end

%% Reading input files
% Validating paths
target = strcat(Workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Target),',').'));
nontarget = strcat(Workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Nontarget),',').'));
control = strcat(Workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Control),',').'));
bin = strcat(Workingdir, tostringmatrix(BED_Bin));
outputdir = strcat(Workingdir, tostringmatrix(Outdir));
k = str2double(k);
if k <= 0
    error('k has to be a positive integer')
end
areallpathsvalid(target);
areallpathsvalid(nontarget);
areallpathsvalid(control);
areallpathsvalid(bin);
[~,~] = mkdir(outputdir);
msg = "All paths are valid";
disp(msg)

%% Creating structure of binned coverage data
% Creation of structure
[X, rtarget, rnontarget, rcontrol] = structify(target, nontarget, control, bin);

% QNorm normalization
X.samp = struct;
samp_names = [target; nontarget; control];
X.samp.name = cell(size(samp_names,1),1);
for i = 1:size(samp_names,1), X.samp.name{i,1} = convertStringsToChars(samp_names(i,1));, end
X.samp.totct = sum(X.bin.ct_raw,1)';
samps_for_median_totct = [rtarget, rnontarget, rcontrol];
X.bin.ct_norm = bsxfun(@rdivide,X.bin.ct_raw,X.samp.totct'/median(X.samp.totct(samps_for_median_totct)));
X.hist.bin = [0;unique(round(geometric_series(1,2000,200)))];
X.samp.hist_norm = nan(slength(X.samp),slength(X.hist));
for i = 1:slength(X.samp)
    X.samp.hist_norm(i,:) = histc(X.bin.ct_norm(:,i),X.hist.bin);
end
X.samp.cf_norm = cumsum(X.samp.hist_norm,2)/slength(X.bin);
X.hist.cf_avg_norm = X.samp.cf_norm(1,:)';

% Find fits for QNorm normalization
for i=1:slength(X.samp)
    fprintf('%d/%d ',i,slength(X.samp));
    x0 = log10(1.5+X.hist.bin); y0 = log10(1-X.hist.cf_avg_norm); y0(isinf(y0))=nan;
    yi = log10(1-X.samp.cf_norm(i,:)');
    lins = geometric_series(0.1,10,100); exps = geometric_series(0.1,10,100);
    score = nan(length(lins),length(exps));
    for lini=1:length(lins),lin=lins(lini);
        for expi=1:length(exps),exp=exps(expi);
            adjbins = lin*(X.hist.bin.^exp); x2 = log10(1.5+adjbins);
            d = bsxfun(@minus,x2,x0'); [~,map] = min(abs(d),[],1);
            score(lini,expi) = nanmean((y0-yi(map)).^2);
        end
    end
    [~,ord] = min(score(:)); [lini,expi] = ind2sub(size(score),ord); lin=lins(lini); exp=exps(expi);
    adjbins = lin*(X.hist.bin.^exp); x2 = log10(1.5+adjbins);
    X.samp.lin(i,1)=lin; X.samp.exp(i,1)=exp;
end
fprintf('\n');
X.samp.exp(rcontrol)=1; X.samp.lin(rcontrol)=1;

% Apply transformation for QNorm normalization
X.bin.ct = bsxfun(@times,bsxfun(@power,X.bin.ct_norm,X.samp.exp'),X.samp.lin');
X.samp.hist = nan(slength(X.samp),slength(X.hist));
for i=1:slength(X.samp)
    X.samp.hist(i,:) = histc(X.bin.ct(:,i),X.hist.bin);
end
X.samp.cf = cumsum(X.samp.hist,2)/slength(X.bin);
save(strcat(outputdir,"/tiles_200_data.mat"),'X','-v7.3');

% Peak selection for target
damp=16; thresh=1; maxgap=3;
X.bin.avgct_rb = mean(X.bin.ct(:,rtarget),2);
X.bin.maxct_rb = max(X.bin.ct(:,rtarget),[],2);
X.bin.maxct_ctl = max(X.bin.ct(:,rcontrol),[],2);
X.bin.log2fc = log2((damp+X.bin.maxct_rb)./(damp+X.bin.maxct_ctl));
bidx = find(X.bin.log2fc>=thresh);
dat = reorder_struct(keep_fields(X.bin,{'chr','pos'}),bidx); dat.bidx = bidx;
dat.diff = difff(dat.bidx); dat.samechr = (difff(dat.chr)==0); dat.diff(~dat.samechr)=inf;
dat.peakst = (dat.diff>maxgap | isnan(dat.diff));
z = nan(slength(dat),1); dat.st=z; dat.en=z; dat.bidx_first=z; dat.bidx_last=z; dat.orig_npeakbins=z;
for i=1:slength(dat), if ~dat.peakst(i), continue; end
    n = find(dat.peakst(i+1:end),1,'first'); if isempty(n), n = slength(dat)-i+1; end
    j = i+n-1; dat.orig_npeakbins(i) = n;
    dat.st(i)=dat.pos(i)-99; dat.en(i)=dat.pos(j)+100; dat.bidx_first(i)=dat.bidx(i); dat.bidx_last(i) = dat.bidx(j);
end
dat.len = dat.en-dat.st+1; dat.orig_bidx_last = dat.bidx_first+dat.orig_npeakbins-1;
dat = reorder_struct(dat,dat.peakst); dat=rmfield(dat,{'bidx','diff','samechr','peakst'});
for i=slength(dat):-1:1
    a = X.bin.avgct_rb(dat.bidx_first(i):dat.bidx_last(i));
    [mx,mxi] = max(a);
    dat.bidx_max(i,1) = mxi+dat.bidx_first(i)-1;
    halfmax = find(a>=0.5*mx);
    dat.bidx_halfmax{i,1} = halfmax+dat.bidx_first(i)-1;
    dat.nbin_halfmax(i,1) = length(halfmax);
    pct80max = find(a>=0.8*mx);
    dat.bidx_pct80max{i,1} = pct80max+dat.bidx_first(i)-1;
    dat.nbin_pct80max(i,1) = length(pct80max);
end
align_to_max = true;
width = 10000/200;
peak = rename_fields(dat,{'bidx_first','bidx_last'},{'bidx_st','bidx_en'});
if ~align_to_max
    peak.bidx_mid = round((peak.bidx_st + peak.orig_bidx_last)/2);
else
    peak.bidx_mid = peak.bidx_max;
end
peak.pos = X.bin.pos(peak.bidx_mid);
peak.bidx_first = peak.bidx_mid-(width/2); peak.bidx_last = peak.bidx_first + width-1;
peak.chr_first = X.bin.chr(peak.bidx_first); peak.chr_last = X.bin.chr(peak.bidx_last); all(peak.chr_first==peak.chr_last) % 1  --> good, no peaks cross a chr boundary
peak = rmfield(peak,{'chr_first','chr_last'});
height = slength(X.samp);
peak.raw_scalar = zeros(slength(peak),slength(X.samp));
peak.dat_scalar = zeros(slength(peak),slength(X.samp));
peak.raw = zeros(slength(peak),width*height);
peak.dat = zeros(slength(peak),width*height);
for i=1:slength(peak)
    peak.raw_scalar(i,:) = sum(X.bin.ct_norm(peak.bidx_st(i):peak.bidx_en(i),:),1);
    peak.dat_scalar(i,:) = sum(X.bin.ct(peak.bidx_st(i):peak.bidx_en(i),:),1);
    peak.raw(i,:) = reshape(X.bin.ct_norm(peak.bidx_first(i):peak.bidx_last(i),:),height*width,1);
    peak.dat(i,:) = reshape(X.bin.ct(peak.bidx_first(i):peak.bidx_last(i),:),height*width,1);
end
X.peak=peak; X.peak.avgct_rb = nan(slength(X.peak),1);for i=1:slength(X.peak),X.peak.avgct_rb(i)=mean(X.bin.avgct_rb(X.peak.bidx_first(i):X.peak.bidx_last(i)));end
X.pixel = []; for samp=1:slength(X.samp), X.pixel.samp((samp-1)*50+[1:50],1)=samp; end; X.pixel.dist = repmat([-4900:200:4900]',slength(X.samp),1);
X.peak = mf2a(X.peak,'pos','chr');

% Clustering
samps_to_cluster = [rtarget, rnontarget]; pixels_to_cluster = find(ismember(X.pixel.samp,samps_to_cluster));
randinit(1234);
X.peak.(['clust',num2str(k)]) = kmeansd(double(1e-5+X.peak.dat(:,pixels_to_cluster)),k,'distance','cosine','maxiter',1000);
X = rmfield(X,'bin'); X.peak.dat=single(X.peak.dat); X.peak.raw=single(X.peak.raw);
save(strcat(outputdir,"/tiles_200_data_peak.mat"),'X','-v7.3');
samps_to_show = [rtarget, rnontarget, rcontrol]';
neighborhood_to_show = 10000; clustfld = ['clust',num2str(k)];
X.peak = sort_struct(X.peak,{clustfld,'avgct_rb'},[1 -1]);
figure(1),clf,hold on,ff,viscap=30;colorscheme=2;
pixels_to_show={}; for i=1:length(samps_to_show), pixels_to_show{i} = find(X.pixel.samp==samps_to_show(i)&abs(X.pixel.dist)<neighborhood_to_show/2); end
pixels_to_show = cat(1,pixels_to_show{:}); dat = X.peak.dat(:,pixels_to_show)';
if colorscheme==1      % original colorscheme
    img=nan(size(dat,1),size(dat,2),3);
    for row=1:size(dat,1),for rgb=1:3,c=min(1,dat(row,:)/viscap);img(row,:,rgb)=(1-c)*X.samp.clr_bkgd(X.pixel.samp(pixels_to_show(row)),rgb);end,end
elseif colorscheme==2  % blue-yellow
    img=convert_1d_colors(dat,parula,0,viscap,[0.8 0.8 0.8]);
end
image(img);set(gca,'ydir','rev','position',[0.135 0.025 0.85 0.97]); xlim(0.5+[0 size(dat,2)]);ylim(0.5+[0 size(dat,1)]);
xlabels_by_group(X.peak.(clustfld));ylabels_by_group(X.samp.name(X.pixel.samp(pixels_to_show)));
w=12;h=7;set(gcf,'papersize',[w h],'paperposition',[0.2 0.2 w-0.4 h-0.4]);
print_to_file(strcat(outputdir,"/clustering_heatmap.pdf"));

end

function output = tostringmatrix(input)
if ischarorstring(input)
    stringmatrix = convertCharsToStrings(input);
    if size(stringmatrix,2) > 1
        if size(stringmatrix,1) > 1
            disp(input)
            error('The above variable is not in the appropriate string or matrix of strings format');
        else
            output = stringmatrix.';
        end
    else
        output = stringmatrix;
    end
else
    disp(input)
    error('The above variable is not in the appropriate string or matrix of strings format');
end
end

function output = ischarorstring(input)
output = ischar(input) | isstring(input);
end

function areallpathsvalid(input)
for i = 1:size(input,1)
    ispathvalid(input(i,1));
end
end

function ispathvalid(input)
if ~exist(input, 'file')
    disp(input);
    error('The above file does not exist');
end
end

function [X, rtarget, rnontarget, rcontrol] = structify(target, nontarget, control, bin)
% Read bins
X = struct;
fileID = fopen(bin,'r');
dataArray = textscan(fileID, '%s%s%s%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
bed = [dataArray{1:end-1}];
X.bin = struct;
X.bin.chr = zeros(size(bed,1),1);
for i = 1:size(X.bin.chr,1)
    X.bin.chr(i,1) = removechr(bed(i,1));
end
X.bin.pos = str2double(bed(:,2))+100;
clearvars fileID dataArray bed ans;

% Define ranges
rtarget = 1:size(target,1);
rnontarget = (size(target,1)+1):(size(target,1)+size(nontarget,1));
rcontrol = (size(target,1)+size(nontarget,1)+1):(size(target,1)+size(nontarget,1)+size(control,1));

% Read bedGraphs
X.bin.ct_raw = zeros(size(X.bin.chr,1),(size(target,1)+size(nontarget,1)+size(control,1)));
for i = 1:size(target,1)
    fileID = fopen(target(i,1),'r');
    dataArray = textscan(fileID, '%*s%*s%*s%f%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'EmptyValue', 0, 'ReturnOnError', false);
    fclose(fileID);
    X.bin.ct_raw(:,rtarget(1,i)) = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
end
for i = 1:size(nontarget,1)
    fileID = fopen(nontarget(i,1),'r');
    dataArray = textscan(fileID, '%*s%*s%*s%f%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'EmptyValue', 0, 'ReturnOnError', false);
    fclose(fileID);
    X.bin.ct_raw(:,rnontarget(1,i)) = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
end
for i = 1:size(control,1)
    fileID = fopen(control(i,1),'r');
    dataArray = textscan(fileID, '%*s%*s%*s%f%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'EmptyValue', 0, 'ReturnOnError', false);
    fclose(fileID);
    X.bin.ct_raw(:,rcontrol(1,i)) = [dataArray{1:end-1}];
    clearvars fileID dataArray ans;
end

msg = "All bedGraph files are imported";
disp(msg)
end

function r = geometric_series(st,en,num)
f = en/st;
step = f.^(1/(num-1));
r = as_column(st*step.^(0:num-1));
r(1) = st;
r(end) = en;
end

function d = difff(x,n,dim)
if ~exist('n','var'), n = 1; end
if ~exist('dim','var')
    if size(x,2)>1 && size(x,1)==1
        dim=2;
    else
        dim=1;
    end
end
d = diff(x,n,dim);
if dim==1
    d = cat(1,nan(n,size(d,2)),d);
elseif dim==2
    d = cat(2,nan(size(d,1),n),d);
elseif dim==3
    fprintf('WARNING: difff with dim>2 is untested, please double-check!\n');
    d = cat(3,nan(size(d,1),size(d,2),n),d);
else
    error('diff is untested with dim>2');
end
if length(size(d))~=length(size(x)) || ~all(size(d)==size(x))
    error('difff failed');
end
end

function withoutchr = removechr(withchr)
chrs = ["chr1";"chr2";"chr3";"chr4";"chr5";"chr6";"chr7";"chr8";"chr9";"chr10";"chr11";"chr12";"chr13";"chr14";"chr15";"chr16";"chr17";"chr18";"chr19";"chr20";"chr21";"chr22";"chrX";"chrY"];
withoutchr = find(strcmpi(withchr,chrs));
end

function [s2,ord]=sort_struct(s1,keyfield,order)
if slength(s1)==0
  s2 = s1;
  ord = [];
  return
end
if length(keyfield)==0, return; end
if ~iscell(keyfield)
  keyfield = {keyfield};
end
if ~exist('order','var')
  order = repmat(1,length(keyfield),1);
end
if ischar(order) && strcmpi(order,'descend')
  order = [-1];
end
if length(order) ~= length(keyfield)
  error('order and keyfield must have same number of elements');
end
if any(order~=1 & order~=-1) error('unknown order type'); end
orig_len = slength(s1);
ord=(1:orig_len)';
fields = fieldnames(s1);
nf = length(fields);
rank = nan(orig_len,nf);
for k=1:length(keyfield)
  f = getfield(s1,keyfield{k});
  if length(f)<orig_len, error('Attempted to sort on truncated field "%s"',keyfield{k}); end
  if islogical(f), f=1*f; end
  if isnumeric(f)
    [u ui uj] = unique(f,'rows');
    [tmp ordi] = sortrows(u);   
  else
    [u ui uj] = unique(f);
    [tmp ordi] = sort(u);
  end
  if order(k)==-1, ordi=ordi(end:-1:1); end
  rank(:,k) = ordi(uj);
end
[tmp ord] = sortrows(rank);
s2 = reorder_struct(s1,ord);
end

function [s,order]=reorder_struct(s,order)
if nargin~=2, error('reorder_struct(s,order)'); end
if islogical(order), order = find(order); end
if ischar(order)
    if strcmpi(order,'end')
        order = slength(s);
    else
        error('invalid index parameter');
    end
end
order = as_column(order);
nanflag = any(isnan(order));
fields = fieldnames(s);
nf = length(fields);
for i=1:nf
    f = getfield(s,fields{i});
    if nanflag
        f = nansub(f,order);
    else
        f = f(order,:,:,:,:,:,:,:,:,:);
    end
    s = setfield(s,fields{i},f);
end
end

function Y = nansub(X,idx,filler)
if length(size(X))==2 && size(X,1)==1 && size(X,2)>1
%   fprintf('note: converting first argument to column vector\n');
  X = X';
end
if iscellstr(X) && size(X,1)==1 && size(X,2)>1
  X=X';
end
if islogical(X)
  type = 0;
elseif isnumeric(X)
  type = 1;
elseif iscell(X)
  type = 2;
else
  error('Unsupported array type');
end
if ~exist('filler','var')
  if type==0
    filler = false;
  elseif type==1
    filler = nan;
  elseif type==2
    filler = {''};
  else
    error('Inconsistent behavior with "type"');
  end
end
if type==0
  if ~islogical(filler)
    error('Inappropriate filler for logical array');
  end
elseif type==1
  if ~isnumeric(filler)
    error('Inappropriate filler for numeric array');
  end
elseif type==2
  if ischar(filler)
    filler = {filler};
  end
  if ~iscell(filler)
    error('Inappropriate filler for cell array');
  end
else
  error('Inconsistent behavior with "type"');
end
sz = size(X); sz(1) = length(idx);
Y = repmat(filler,sz);
idx2 = find(~isnan(idx) & idx>=1 & idx<=length(X));
Y(idx2,:,:,:,:,:,:,:,:) = X(idx(idx2),:,:,:,:,:,:,:,:);
end

function varargout = mf2a(varargin)
if nargout==0
    move_field_to_after(varargin{:});
elseif nargout>1
    varargout = cell(nargout,1);
    [varargout{:}] = move_field_to_after(varargin{:});
else
    [varargout{1}] = move_field_to_after(varargin{:});
end
end

function Y = move_field_to_after(X,fld1,fld2)
demand_fields(X,{fld1,fld2});
f = fieldnames(X);
for i=1:length(f)
    if strcmp(f{i},fld1), continue; end
    Y.(f{i}) = X.(f{i});
    if strcmp(f{i},fld2)
        Y.(fld1) = X.(fld1);
    end
end
end

function demand_fields(varargin)
require_fields(varargin{:});
end

function require_fields(T,fields)
if ~iscell(fields)
    fields = {fields};
end
for i=1:length(fields)
    if ~isfield(T,fields{i})
        error(['Structure is missing required field "' fields{i} '"']);
    end
end
end

function [idx,c,sumd,d] = kmeansd(x,k,varargin)
initcen = rand(k,size(x,2));
[idx,c,sumd,d] = kmeans(x,k,'start',initcen,varargin{:});
end

function varargout = kmeans(X, k, varargin)
if nargin > 0
    X = convertStringsToChars(X);
end
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end
if nargin < 2
    error(message('stats:kmeans:TooFewInputs'));
end
if ~isreal(X)
    error(message('stats:kmeans:ComplexData'));
end
wasnan = any(isnan(X),2);
hadNaNs = any(wasnan);
if hadNaNs
    warning(message('stats:kmeans:MissingDataRemoved'));
    X = X(~wasnan,:);
end
[n, p] = size(X);
pnames = {   'distance'  'start' 'replicates' 'emptyaction' 'onlinephase' 'options' 'maxiter' 'display'};
dflts =  {'sqeuclidean' 'plus'          []  'singleton'         'off'        []        []        []};
[distance,start,reps,emptyact,online,options,maxit,display] ...
    = internal.stats.parseArgs(pnames, dflts, varargin{:});
distNames = {'sqeuclidean','cityblock','cosine','correlation','hamming'};
distance = internal.stats.getParamVal(distance,distNames,'''Distance''');
switch distance
    case 'cosine'
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error(message('stats:kmeans:ZeroDataForCos'));
        end
        X =  X./Xnorm;
    case 'correlation'
        X = X - mean(X,2);
        Xnorm = sqrt(sum(X.^2, 2));
        if any(min(Xnorm) <= eps(max(Xnorm)))
            error(message('stats:kmeans:ConstantDataForCorr'));
        end
        X =  X./Xnorm;
    case 'hamming'
        if  ~all( X(:) ==0 | X(:)==1)
            error(message('stats:kmeans:NonbinaryDataForHamm'));
        end
end
Xmins = [];
Xmaxs = [];
CC = [];
if ischar(start)
    startNames = {'uniform','sample','cluster','plus','kmeans++'};
    j = find(strncmpi(start,startNames,length(start)));
    if length(j) > 1
        error(message('stats:kmeans:AmbiguousStart', start));
    elseif isempty(j)
        error(message('stats:kmeans:UnknownStart', start));
    elseif isempty(k)
        error(message('stats:kmeans:MissingK'));
    end
    start = startNames{j};
    if strcmp(start, 'uniform')
        if strcmp(distance, 'hamming')
            error(message('stats:kmeans:UniformStartForHamm'));
        end
        Xmins = min(X,[],1);
        Xmaxs = max(X,[],1);
    end
elseif isnumeric(start)
    CC = start;
    start = 'numeric';
    if isempty(k)
        k = size(CC,1);
    elseif k ~= size(CC,1)
        error(message('stats:kmeans:StartBadRowSize'));
    end
    if size(CC,2) ~= p
        error(message('stats:kmeans:StartBadColumnSize'));
    end
    if isempty(reps)
        reps = size(CC,3);
    elseif reps ~= size(CC,3)
        error(message('stats:kmeans:StartBadThirdDimSize'));
    end
    if isequal(distance, 'correlation')
        CC = CC - mean(CC,2);
    end
else
    error(message('stats:kmeans:InvalidStart'));
end
emptyactNames = {'error','drop','singleton'};
emptyact = internal.stats.getParamVal(emptyact,emptyactNames,'''EmptyAction''');
[~,online] = internal.stats.getParamVal(online,{'on','off'},'''OnlinePhase''');
online = (online==1);
if ~isempty(display)
    options = statset(options,'Display',display);
end
if ~isempty(maxit)
    options = statset(options,'MaxIter',maxit);
end
options = statset(statset('kmeans'), options);
display = find(strncmpi(options.Display, {'off','notify','final','iter'},...
    length(options.Display))) - 1;
maxit = options.MaxIter;
if ~(isscalar(k) && isnumeric(k) && isreal(k) && k > 0 && (round(k)==k))
    error(message('stats:kmeans:InvalidK'));
elseif n < k
    error(message('stats:kmeans:TooManyClusters'));
end
if isempty(reps)
    reps = 1;
elseif ~internal.stats.isScalarInt(reps,1)
    error(message('stats:kmeans:BadReps'));
end
[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(options,true);
usePool = useParallel && poolsz>0;
if display > 1
    if usePool
        internal.stats.parallel.distributeToPool( ...
            'workerID', num2cell(1:poolsz) );
        if display == 3
            warning(message('stats:kmeans:displayParallel2'));
            fprintf('    worker\t  iter\t phase\t     num\t         sum\n' );
        end
    else
        if useParallel
            warning(message('stats:kmeans:displayParallel'));
        end
        if display == 3
            fprintf('  iter\t phase\t     num\t         sum\n');
        end
    end
end
if issparse(X) || ~isfloat(X) || strcmp(distance,'cityblock') || ...
        strcmp(distance,'hamming')
    [varargout{1:nargout}] = kmeans2(X,k, distance, emptyact,reps,start,...
        Xmins,Xmaxs,CC,online,display, maxit,useParallel, RNGscheme,usePool,...
        wasnan,hadNaNs,varargin{:});
    return;
end
emptyErrCnt = 0;
loopbody = @loopBody;
totsumD = 0;
iter = 0;
X = X';
Xmins = Xmins';
Xmaxs = Xmaxs';
ClusterBest = internal.stats.parallel.smartForReduce(...
    reps, loopbody, useParallel, RNGscheme, 'argmin');
varargout{1} = ClusterBest{5};
varargout{2} = ClusterBest{6}';
varargout{3} = ClusterBest{3};
totsumDbest = ClusterBest{1};
if nargout > 3
    varargout{4} = ClusterBest{7};
end
if display > 1
    fprintf('%s\n',getString(message('stats:kmeans:FinalSumOfDistances',sprintf('%g',totsumDbest))));
end
if hadNaNs
    varargout{1} = statinsertnan(wasnan, varargout{1});
    if nargout > 3
        varargout{4} = statinsertnan(wasnan, varargout{4});
    end
end
    function cellout = loopBody(rep,S)
        if isempty(S)
            S = RandStream.getGlobalStream;
        end
        if display > 1 % 'iter'
            if usePool
                dispfmt = '%8d\t%6d\t%6d\t%8d\t%12g\n';
                labindx = internal.stats.parallel.workerGetValue('workerID');
            else
                dispfmt = '%6d\t%6d\t%8d\t%12g\n';
            end
        end
        cellout = cell(7,1);
        cellout{1} = Inf;
        cellout{2} = rep;
        switch start
            case 'uniform'
                C = Xmins(:,ones(1,k)) + rand(S,[k,p])'.*(Xmaxs(:,ones(1,k))-Xmins(:,ones(1,k)));
                if isequal(distance, 'correlation')
                    C = C - mean(C,1);
                end
                if isa(X,'single')
                    C = single(C);
                end
            case 'sample'
                C = X(:,randsample(S,n,k));
            case 'cluster'
                k0 = floor(max(k, 0.1*n));
                Xsubset = X(:,randsample(S,n,k0));
                if k<k0
                    optIndex = find(strcmpi('options',varargin));
                    if isempty(optIndex)
                        opts = statset('Display','off');
                        varargin = [varargin,'options',opts];
                    else
                        varargin{optIndex+1}.Display = 'off';
                    end
                    [~, C] = kmeans(Xsubset', k, varargin{:}, 'start','sample', 'replicates',1);
                    C = C';
                else
                    C = Xsubset;
                end
            case 'numeric'
                C = CC(:,:,rep)';
                if isa(X,'single')
                    C = single(C);
                end
            case {'plus','kmeans++'}
                index = zeros(1,k);
                [C(:,1), index(1)] = datasample(S,X,1,2);
                minDist = inf(n,1);
                for ii = 2:k
                    minDist = min(minDist,distfun(X,C(:,ii-1),distance));
                    denominator = sum(minDist);
                    if denominator==0 || isinf(denominator) || isnan(denominator)
                        C(:,ii:k) = datasample(S,X,k-ii+1,2,'Replace',false);
                        break;
                    end
                    sampleProbability = minDist/denominator;
                    [C(:,ii), index(ii)] = datasample(S,X,1,2,'Replace',false,...
                        'Weights',sampleProbability);
                end
        end
        if ~isfloat(C)
            C = double(C);
        end
        D = distfun(X, C, distance, 0, rep, reps);
        [d, idx] = min(D, [], 2);
        m = accumarray(idx,1,[k,1])';
        try
            converged = batchUpdate();
            if online
                converged = onlineUpdate();
            end
            if display == 2
                fprintf('%s\n',getString(message('stats:kmeans:IterationsSumOfDistances',rep,iter,sprintf('%g',totsumD) )));
            end
            if ~converged
                if reps==1
                    warning(message('stats:kmeans:FailedToConverge', maxit));
                else
                    warning(message('stats:kmeans:FailedToConvergeRep', maxit, rep));
                end
            end
            nonempties = find(m>0);
            D(:,nonempties) = distfun(X, C(:,nonempties), distance, iter, rep, reps);
            d = D((idx-1)*n + (1:n)');
            sumD = accumarray(idx,d,[k,1]);
            totsumD = sum(sumD(nonempties));
            cellout = {totsumD,rep,sumD,iter,idx,C,D}';
        catch ME
            if reps == 1 || (~isequal(ME.identifier,'stats:kmeans:EmptyCluster')  && ...
                    ~isequal(ME.identifier,'stats:kmeans:EmptyClusterRep'))
                rethrow(ME);
            else
                emptyErrCnt = emptyErrCnt + 1;
                warning(message('stats:kmeans:EmptyClusterInBatchUpdate', rep, iter));
                if emptyErrCnt == reps
                    error(message('stats:kmeans:EmptyClusterAllReps'));
                end
            end
        end
        function converged = batchUpdate()
            moved = 1:n;
            changed = 1:k;
            previdx = zeros(n,1);
            prevtotsumD = Inf;
            iter = 0;
            converged = false;
            while true
                iter = iter + 1;
                [C(:,changed), m(changed)] = gcentroids(X, idx, changed, distance);
                D(:,changed) = distfun(X, C(:,changed), distance, iter, rep, reps);
                empties = changed(m(changed) == 0);
                if ~isempty(empties)
                    if strcmp(emptyact,'error')
                        if reps==1
                            error(message('stats:kmeans:EmptyCluster', iter));
                        else
                            error(message('stats:kmeans:EmptyClusterRep', iter, rep));
                        end
                    end
                    switch emptyact
                        case 'drop'
                            if reps==1
                                warning(message('stats:kmeans:EmptyCluster', iter));
                            else
                                warning(message('stats:kmeans:EmptyClusterRep', iter, rep));
                            end
                            D(:,empties) = NaN;
                            changed = changed(m(changed) > 0);
                        case 'singleton'
                            for i = empties
                                d = D((idx-1)*n + (1:n)');
                                [~, lonely] = max(d);
                                from = idx(lonely);
                                if m(from) < 2
                                    from = find(m>1,1,'first');
                                    lonely = find(idx==from,1,'first');
                                end
                                C(:,i) = X(:,lonely);
                                m(i) = 1;
                                idx(lonely) = i;
                                D(:,i) = distfun(X, C(:,i), distance, iter, rep, reps);
                                [C(:,from), m(from)] = gcentroids(X, idx, from, distance);
                                D(:,from) = distfun(X, C(:,from), distance, iter, rep, reps);
                                changed = unique([changed from]);
                            end
                    end
                end
                totsumD = sum(D((idx-1)*n + (1:n)'));
                if prevtotsumD <= totsumD
                    idx = previdx;
                    [C(:,changed), m(changed)] = gcentroids(X, idx, changed, distance);
                    iter = iter - 1;
                    break;
                end
                if display > 2
                    if usePool
                        fprintf(dispfmt,labindx,iter,1,length(moved),totsumD);
                    else
                        fprintf(dispfmt,iter,1,length(moved),totsumD);
                    end
                end
                if iter >= maxit
                    break;
                end
                previdx = idx;
                prevtotsumD = totsumD;
                [d, nidx] = min(D, [], 2);
                moved = find(nidx ~= previdx);
                if ~isempty(moved)
                    moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
                end
                if isempty(moved)
                    converged = true;
                    break;
                end
                idx(moved) = nidx(moved);
                changed = unique([idx(moved); previdx(moved)])';
            end
        end
        function converged = onlineUpdate()
            changed = find(m > 0);
            lastmoved = 0;
            nummoved = 0;
            iter1 = iter;
            converged = false;
            Del = NaN(n,k);
            while iter < maxit
                switch distance
                    case 'sqeuclidean'
                        for i = changed
                            mbrs = (idx == i)';
                            sgn = 1 - 2*mbrs;
                            if m(i) == 1
                                sgn(mbrs) = 0;
                            end
                            Del(:,i) = (m(i) ./ (m(i) + sgn)) .* sum((X - C(:,i)).^2, 1);
                        end
                    case {'cosine','correlation'}
                        normC = sqrt(sum(C.^2, 1));
                        if any(normC < eps(class(normC)))
                            if reps==1
                                error(message('stats:kmeans:ZeroCentroid', iter));
                            else
                                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
                            end
                        end
                        for i = changed
                            XCi =  C(:,i)'*X;
                            mbrs = (idx == i)';
                            sgn = 1 - 2*mbrs;
                            Del(:,i) = 1 + sgn .*...
                                (m(i).*normC(i) - sqrt((m(i).*normC(i)).^2 + 2.*sgn.*m(i).*XCi + 1));
                        end
                end
                previdx = idx;
                [minDel, nidx] = min(Del, [], 2);
                moved = find(previdx ~= nidx);
                moved(m(previdx(moved))==1)=[];
                if ~isempty(moved)
                    moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
                end
                if isempty(moved)
                    if (iter == iter1) || nummoved > 0
                        iter = iter + 1;
                        if display > 2 % 'iter'
                            if usePool
                                fprintf(dispfmt,labindx,iter,2,length(moved),totsumD);
                            else
                                fprintf(dispfmt,iter,2,length(moved),totsumD);
                            end
                        end
                    end
                    converged = true;
                    break;
                end
                moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;
                if moved <= lastmoved
                    iter = iter + 1;
                    if display > 2 % 'iter'
                        if usePool
                            fprintf(dispfmt,labindx,iter,2,length(moved),totsumD);
                        else
                            fprintf(dispfmt,iter,2,length(moved),totsumD);
                        end
                    end
                    if iter >= maxit, break; end
                    nummoved = 0;
                end
                nummoved = nummoved + 1;
                lastmoved = moved;
                oidx = idx(moved);
                nidx = nidx(moved);
                totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);
                idx(moved) = nidx;
                m(nidx) = m(nidx) + 1;
                m(oidx) = m(oidx) - 1;
                switch distance
                    case {'sqeuclidean','cosine','correlation'}
                        C(:,nidx) = C(:,nidx) + (X(:,moved) - C(:,nidx)) / m(nidx);
                        C(:,oidx) = C(:,oidx) - (X(:,moved) - C(:,oidx)) / m(oidx);
                end
                changed = sort([oidx nidx]);
            end
        end
    end
end

function D = distfun(X, C, dist, iter,rep, reps)
switch dist
    case 'sqeuclidean'
        D = internal.stats.pdist2mex(X,C,'sqe',[],[],[],[]);
    case {'cosine','correlation'}
        normC = sqrt(sum(C.^2, 1));
        if any(normC < eps(class(normC)))
            if reps==1
                error(message('stats:kmeans:ZeroCentroid', iter));
            else
                error(message('stats:kmeans:ZeroCentroidRep', iter, rep));
            end
        end
        C = C./normC;
        D = internal.stats.pdist2mex(X,C,'cos',[],[],[],[]);
end
end

function [centroids, counts] = gcentroids(X, index, clusts, dist)
p = size(X,1);
num = length(clusts);
centroids = NaN(p,num,'like',X);
counts = zeros(1,num,'like',X);
for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);
        switch dist
            case {'sqeuclidean','cosine','correlation'}
                centroids(:,i) = sum(X(:,members),2) / counts(i);
        end
    end
end
end

function print_to_file(filename,varargin)
if nargin == 1
    [dev res] = interpret_print_filename(filename);
else
    dev = interpret_print_filename(filename);
    res = varargin{1};
end
if nargin == 1
    [dev res] = interpret_print_filename(filename);
else
    dev = interpret_print_filename(filename);
    res = varargin{1};
end
fprintf('Outputting figure to %s\n', filename);
print(['-d' dev], ['-r' num2str(res)], filename);
end

function [dev,res] = interpret_print_filename(filename)
tmp = regexp(filename, '(\.[^\.]*)$', 'tokens');
if isempty(tmp) || isempty(tmp{1})
    error('Please specify output file with extension .png, .jpg, .eps, .pdf, or .tif'); ...
end
ext = tmp{1}{1}(2:end);
if strcmpi(ext,'jpeg') || strcmpi(ext,'jpg')
    dev = 'jpeg';
    res = 300;
elseif strcmpi(ext,'eps')
    dev = 'epsc';
    res = 180;
elseif strcmpi(ext,'tif') || strcmpi(ext,'tiff')
    dev = 'tiff';
    res = 180;
elseif strcmpi(ext,'png')
    dev = 'png';
    res = 180;
elseif strcmpi(ext,'pdf')
    dev = 'pdf';
    res = 1200;
else
    error('Unknown output format: please specify .png, .jpg, .eps, .pdf, or .tif');
end
end

function l = slength(S,quiet)
if ~exist('quiet','var'), quiet=false; end
l=NaN;
if isstruct(S)
    l = 0;
    if ~isempty(S) && ~isempty(fieldnames(S))
        f = fields(S);
        nf = length(f);
        len = nan(nf,1);
        for i=1:nf
            f1 = getfield(S,f{i});
            len(i) = size(f1,1);
        end
        ulen = unique(len);
        if length(ulen)==1, l = ulen;
        else
            if ~quiet, fprintf('Warning: deprecated use of slength for structure with fields of nonuniform length\n'); end
            l = len(1);
        end
    end
end
end

function S2 = keep_fields(S,flds)
if ischar(flds), flds = {flds}; end
S2=[];
for i=1:length(flds)
  if isempty(S)
    f = [];
  else
    f = getfield(S,flds{i});
  end
  S2=setfield(S2,flds{i},f);
end
end

function S = rename_fields(S, oldname, newname)
S = rename_field(S,oldname,newname);
end

function S = rename_field(S, oldname, newname)
if iscell(oldname) && iscell(newname)
  if ~iscell(newname) || length(oldname)~=length(newname), error('lists must be same length'); end
elseif ~iscell(oldname) && ~iscell(newname)
  oldname = {oldname};
  newname = {newname};
else
  error('improper parameters');
end
flds = fieldnames(S);
for i=1:length(oldname)
  f = getfield(S, oldname{i});
  S = setfield(S, newname{i}, f);
  if ~strcmp(oldname{i},newname{i})
    S = rmfield(S, oldname{i});
  end
  idx = find(strcmp(flds,oldname{i}));
  if length(idx)~=1, error('unexpected behavior'); end
  flds{idx} = newname{i};
end
S = order_fields_first(S,unique_keepord(flds));
end

function [u ui uj] = unique_keepord(x,varargin);
if exist('varargin','var') && length(varargin)>=1 && ischar(varargin{1}) && (strcmpi(varargin{1},'first')|strcmpi(varargin{1},'last'))
  error('please do not specify "first" or "last" with this function.  (default is "first")');
end
[u1 ui1 uj1] = unique(x,'first',varargin{:});
[ui ord] = sort(ui1);
u = x(ui1(ord));
[tmp ord2] = sort(ord);
uj = ord2(uj1);
return
if iscell(x)
  if any(~strcmp(x,u(uj))) || any(~strcmp(u,x(ui))), error('unique_keepord not working properly!!!'); end
else
  if any(x~=u(uj)) || any(u~=x(ui)), error('unique_keepord not working properly!!!'); end
end
end

function X = order_fields_first(varargin)
X = orderfields_first(varargin{:});
end

function S = orderfields_first(S,first_flds)
if ischar(first_flds), first_flds = {first_flds}; end
all_flds = fieldnames(S);
if ~isempty(setdiff(first_flds,all_flds)), error('Some of those fields don''t exist'); end
rest_flds = all_flds;
rest_flds(ismember(rest_flds,first_flds)) = [];
S = orderfields(S,[as_column(first_flds);as_column(rest_flds)]);
end

function rgb = convert_1d_colors(orig,cmap,cmin,cmax,nancolor)
if ~exist('cmap','var'), cmap=colormap(); end
if length(size(orig))>2, error('only works on 1D or 2D input'); end
nrow=size(orig,1);
ncol=size(orig,2);
if ~exist('cmin','var'), cmin=min(orig(:)); end
if ~exist('cmax','var'), cmax=max(orig(:)); end
crange=cmax-cmin;
orig(orig<cmin)=cmin;
orig(orig>cmax)=cmax;
if ~exist('nancolor','var'), nancolor = [0.7 0.7 0.7]; end
nc=size(cmap,1);
ci = 1+floor((nc-1)*(orig-cmin)/crange);
cmap(nc+1,:) = nancolor;
ci(isnan(ci)) = nc+1;
rgb = nan(nrow,ncol,3);
for i=1:3, rgb(:,:,i) = reshape(cmap(ci(:),i),nrow,ncol); end
end

function xlabels_by_group(labels,lines_color,varargin)
rot=90; halign='right'; valign='middle';
if exist('lines_color','var') && ischar(lines_color) && strcmpi('hor',lines_color)
  rot=0; halign='center'; valign='top';
  clear lines_color
end
if ~exist('lines_color','var'), lines_color = [0 0 0]; end
if isnumeric(labels), labels = num2cellstr(labels); end
labels = as_column(labels);
l1 = [labels(1:end)];
l2 = [labels(2:end);'***'];
b = find(~strcmp(l1,l2));
set(gca,'xtick',[0;b+0.5],'xticklabel',{});
yl = ylim();
if strcmp(get(gca,'ydir'),'reverse')
  ypos = yl(2);
else
  ypos = yl(1);
end
b = [0;b];
for i=1:length(b)
  if i>1
    xpos = (b(i-1) + b(i)) / 2;
text(xpos,ypos,labels{b(i)},'rotation',rot,'clipping','off','horizontalalignment',halign,'verticalalignment',valign,'interpreter','none',varargin{:});
  end
end
if ~isempty(lines_color)
  for i=0:length(b)
    if i==0
      line([0 0]+0.5,ylim,'color',lines_color,'clipping','off');
    else
      line([b(i) b(i)]+0.5,ylim,'color',lines_color);
    end
  end
end
end

function ylabels_by_group(labels,lines_color,varargin)
if ~exist('lines_color','var'), lines_color = [0 0 0]; end
if isnumeric(labels), labels = num2cellstr(labels); end
labels = as_column(labels);
l1 = [labels(1:end)];
l2 = [labels(2:end);'***'];
b = find(~strcmp(l1,l2));
set(gca,'ytick',[0;b+0.5],'yticklabel',{});
xl = xlim();
if strcmp(get(gca,'xdir'),'reverse')
  xpos = xl(2);
else
  xpos = xl(1) - 0.01*(xl(2)-xl(1));
end
for i=1:length(b)
  if i==1
    y1=0;
  else
    y1=b(i-1);
  end
  y2 = b(i);
  ypos = y1 + 0.5*(y2-y1);
text(xpos,ypos,labels{b(i)},'clipping','off','verticalalign','middle',varargin{:},'horizontalalign','right','interpreter','none');
end
if ~isempty(lines_color)
  for i=0:length(b)
    if i==0
      line(xlim,[0 0]+0.5,'color',lines_color,'clipping','off');
    else
      line(xlim,[b(i) b(i)]+0.5,'color',lines_color);
    end
  end
end
end

function x=as_column(x)
if size(x,2)>1
    x=x';
end
end

function randinit(randseed)
if ~exist('randseed','var'), randseed=1234; end
rand('twister',randseed);
randn('seed',randseed);
end

function ff
% finish figure
set(gca,'tickdir','out','linewidth',1.5,'xcolor',[0 0 0],'ycolor',[0 0 0],'zcolor',[0 0 0]);
set(gcf,'color',[1 1 1]);
end

function A = num2cellstr(a);
A = cell(length(a),1);
for i=1:length(a)
  A{i} = num2str(a(i));
end
end
