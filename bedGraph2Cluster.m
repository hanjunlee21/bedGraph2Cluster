function bedGraph2Cluster(bedGraphs_Signal, bedGraphs_Control, bedGraphs_Cluster, bedGraphs_Heatmap, outdir, bed_bin, fold_change, normalization_method, k, distance_method, clustering_method, workingdir)
%% bedGraph2Cluster
% e.g., bedGraph2Cluster("bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph", "bedgraph/INPUT.WT.bedgraph,bedgraph/INPUT.dCDK.bedgraph", "bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph", "bedgraph/RB.WT.bedgraph,bedgraph/RB.dCDK.bedgraph,bedgraph/H3K4me3.WT.bedgraph,bedgraph/H3K4me3.dCDK.bedgraph,bedgraph/H3K4me.WT.bedgraph,bedgraph/H3K4me.dCDK.bedgraph,bedgraph/H3K27ac.WT.bedgraph,bedgraph/H3K27ac.dCDK.bedgraph,bedgraph/E2F1.bedgraph,bedgraph/CTCF.shSCR.bedgraph,bedgraph/c-Jun.shSCR.bedgraph", "output", "bed/hg19.200bp.bed", "2", "true", "8", "cosine", "1", "../")
% 
% Required arguments
%     bedGraphs_Signal (string): comma-delimited list of bedGraph files to be included during peak calling
%     bedGraphs_Control (string): comma-delimited list of bedGraph files to be used as controls for peak calling  
%     bedGraphs_Cluster (string): comma-delimited list of bedGraph files to be included during k-means clustering
%     bedGraphs_Heatmap (string): comma-delimited list of bedGraph files to be included in heatmap
%     outdir (string): path to the output directory
%     bed_bin (string): path to the BED file used for binned bedGraph generation
%     fold_change (string): threshold for the fold change of signal over control during peak calling
%     normalization_method (string): normalization method to utilize ("QNorm": QNorm, "CPM": CPM)
%     k (string): number of clusters during k-means clustering
%     distance_method (string): distance metric for k-means clustering ("sqeuclidean", "cityblock", "cosine", "correlation")
%     clustering_method (string): clustering method to utilize ("1" = profile, "2" = profile+scalar, "3" = symmetry_collapsed_profile+scalar)
%
% Optional arguments 
%     workingdir (string): path to the working directory
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
warning('off','all')
if strcmpi(convertCharsToStrings(normalization_method),"true") || strcmpi(convertCharsToStrings(normalization_method),"yes") || strcmpi(convertCharsToStrings(normalization_method),"QNorm")
    normalization_method = true;
elseif strcmpi(convertCharsToStrings(normalization_method),"false") || strcmpi(convertCharsToStrings(normalization_method),"no") || strcmpi(convertCharsToStrings(normalization_method),"CPM")
    normalization_method = false;
else
    error('QNorm has an inappropriate value')
end

if ~strcmp(convertCharsToStrings(distance_method),"sqeuclidean") && ~strcmp(convertCharsToStrings(distance_method),"cityblock") && ~strcmp(convertCharsToStrings(distance_method),"cosine") ...
  && ~strcmp(convertCharsToStrings(distance_method),"correlation") && ~strcmp(convertCharsToStrings(distance_method),"hamming")
    error('distance has to be one of following: sqeuclidean, cityblock, cosine, correlation, hamming')
end

if ~strcmp(convertCharsToStrings(clustering_method),"1") && ~strcmp(convertCharsToStrings(clustering_method),"2") && ~strcmp(convertCharsToStrings(clustering_method),"3")
    error('clustering_method has to be one of following: ("1" = profile, "2" = profile+scalar, "3" = symmetry_collapsed_profile+scalar)') 
end

if ~exist('Workingdir','var')
    workingdir = strcat(convertCharsToStrings(pwd),"/");
else
    workingdir = strcat(tostringmatrix(workingdir),"/");
end

%% Reading input files
% Validating paths
signal = strcat(workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Signal),',').'));
control = strcat(workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Control),',').'));
cluster = strcat(workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Cluster),',').'));
heatmap = strcat(workingdir, tostringmatrix(strsplit(tostringmatrix(bedGraphs_Heatmap),',').'));
bin = strcat(workingdir, tostringmatrix(bed_bin));
outputdir = strcat(workingdir, tostringmatrix(outdir));
fold_change = str2double(fold_change);
if fold_change <= 1
    error('fold_change has to be greater than one')
end
k = str2double(k);
if k <= 0
    error('k has to be a positive integer')
end
areallpathsvalid(signal);
areallpathsvalid(control);
areallpathsvalid(cluster);
areallpathsvalid(heatmap);
areallpathsvalid(bin);
[~,~] = mkdir(outputdir);
msg = "All paths are valid";
disp(msg)

%% Creating structure of binned coverage data
% Creation of structure
[X, rsignal, rcontrol, rcluster, rheatmap] = structify(signal, control, cluster, heatmap, bin);

% Normalization to median total counts
X.samp.totct = sum(X.bin.ct_raw,1)';
samps_for_median_totct = unique([rsignal;rcontrol;rcluster]);
median_totct = median(X.samp.totct(samps_for_median_totct));
X.bin.ct_norm = bsxfun(@rdivide,X.bin.ct_raw,X.samp.totct'/median_totct);

% Normalization method. QNorm is a two parameter normalization method that was inspired by the S3norm method of Xiang et al. https://doi.org/10.1093/nar/gkaa105
if normalization_method == true
  fprintf('QNorm: ');
  qnorm_ref_sample = rsignal(1);  % first sample in the list of "signal" bedgraphs will be used as the reference to normalize other samples to
  X.hist.bin = [0;unique(round(geometric_series(1,2000,200)))];
  X.samp.hist_norm = nan(slength(X.samp),slength(X.hist));
  for i = 1:slength(X.samp)
    X.samp.hist_norm(i,:) = histc(X.bin.ct_norm(:,i),X.hist.bin);
  end
  X.samp.cf_norm = cumsum(X.samp.hist_norm,2)/slength(X.bin);
  X.hist.cf_ref_norm = X.samp.cf_norm(qnorm_ref_sample,:)';
  X.samp.lin = nan(slength(X.samp),1); X.samp.exp = nan(slength(X.samp),1);
  for i=1:slength(X.samp)
    fprintf('%d/%d ',i,slength(X.samp));
    x0 = log10(1.5+X.hist.bin); y0 = log10(1-X.hist.cf_ref_norm); y0(isinf(y0))=nan;
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
  X.bin.ct = bsxfun(@times,bsxfun(@power,X.bin.ct_norm,X.samp.exp'),X.samp.lin');
else
  X.bin.ct = X.bin.ct_norm;
end

% Save file with all data
save(strcat(outputdir,"/tiles_200_data.mat"),'X','-v7.3');

% Peak selection
thresh=log(fold_change)/log(2); maxgap=3;
X.bin.avgct_rb = mean(X.bin.ct(:,rsignal),2);
pseudocount = round(quantile(X.bin.avgct_rb,1-0.005)); % The algorithm assumes the top 0.5 percent of counts per bin as signals and the bottom 99.5% of counts per bin as potential noise
X.bin.maxct_rb = max(X.bin.ct(:,rsignal),[],2);
X.bin.maxct_ctl = max(X.bin.ct(:,rcontrol),[],2);
X.bin.log2fc = log2((pseudocount+X.bin.maxct_rb)./(pseudocount+X.bin.maxct_ctl));
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

% Extract data
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
samps_to_cluster = [rcluster];
pixels_to_cluster = []; for i=1:length(samps_to_cluster), pixels_to_cluster=[pixels_to_cluster;find(X.pixel.samp==samps_to_cluster(i))]; end
if strcmp(convertCharsToStrings(clustering_method),"1")
    % "1" = profile (basic)
    kmean_input = double(1e-5+X.peak.dat(:,pixels_to_cluster));
elseif strcmp(convertCharsToStrings(clustering_method),"2")
    % "2" = profile + scalar
    scalar = double(1e-5+X.peak.dat_scalar(:,samps_to_cluster));
    profile = double(1e-5+X.peak.dat(:,pixels_to_cluster));
    kmean_input = [scalar profile];
elseif strcmp(convertCharsToStrings(clustering_method),"3")
    % "3" = symmetry_collapsed_profile + scalar
    scalar = double(1e-5+X.peak.dat_scalar(:,samps_to_cluster));
    profile = double(1e-5+X.peak.dat(:,pixels_to_cluster));
    for i=1:50:size(profile,2),profile(:,i:i+24)=profile(:,i:i+24)+profile(:,i+49:-1:i+25);profile(:,i+49:-1:i+25)=nan;end;profile(:,all(isnan(profile),1))=[];
    kmean_input = [scalar profile];
end
randinit(1234);
X.peak.(['clust',num2str(k)]) = kmeansd(kmean_input,k,'distance',convertStringsToChars(distance_method),'maxiter',1000);

% Save file with peaks only
X = rmfield(X,'bin'); X.peak.dat=single(X.peak.dat); X.peak.raw=single(X.peak.raw);
save(strcat(outputdir,"/tiles_200_data_peak.mat"),'X','-v7.3');

% Heatmap
samps_to_show = [rheatmap]';
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
print_to_file(strcat(outputdir,"/clustering_heatmap.pdf"),300);
close all;

% Export peakset to BED files
chrs = ["chr1";"chr2";"chr3";"chr4";"chr5";"chr6";"chr7";"chr8";"chr9";"chr10";"chr11";"chr12";"chr13";"chr14";"chr15";"chr16";"chr17";"chr18";"chr19";"chr20";"chr21";...
        "chr22";"chrX";"chrY"];
fileID = fopen(convertStringsToChars(strcat(outputdir,"/peaks.bed")), 'w');
for i = 1:size(X.peak.chr,1)
    fprintf(fileID, '%s\t',chrs(X.peak.chr(i,1),1));
    fprintf(fileID, '%s\t',convertCharsToStrings(num2str(X.peak.st(i,1)+1)));
    fprintf(fileID, '%s\t',convertCharsToStrings(num2str(X.peak.en(i,1))));
    fprintf(fileID, '%s\t',convertCharsToStrings(num2str(i)));
    fprintf(fileID, '%s\t',"1000");
    fprintf(fileID, '%s\n',".");
end
fclose(fileID);
for clusterID = 1:k
    fileID = fopen(convertStringsToChars(strcat(outputdir,"/peaks.clust",convertCharsToStrings(num2str(k)),".",convertCharsToStrings(num2str(clusterID)),".bed")), 'w');
    for i = 1:size(X.peak.chr,1)
        if X.peak.(['clust',num2str(k)])(i,1) == clusterID
            fprintf(fileID, '%s\t',chrs(X.peak.chr(i,1),1));
            fprintf(fileID, '%s\t',convertCharsToStrings(num2str(X.peak.st(i,1)+1)));
            fprintf(fileID, '%s\t',convertCharsToStrings(num2str(X.peak.en(i,1))));
            fprintf(fileID, '%s\t',convertCharsToStrings(num2str(i)));
            fprintf(fileID, '%s\t',"1000");
            fprintf(fileID, '%s\n',".");
        end
    end
    fclose(fileID);
end
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

function [X, rsignal, rcontrol, rcluster, rheatmap] = structify(signal, control, cluster, heatmap, bin)

X = struct;
  
% Make list of unique bedgraphs (so we load each only once)
X.samp = struct;
X.samp.bedgraph = unique([signal; control; cluster; heatmap]);
for i = 1:slength(X.samp)
  [~,samp_name,~] = fileparts(X.samp.bedgraph(i,1));
  X.samp.name{i,1} = convertStringsToChars(samp_name);
end

% Define sample ranges
rsignal = listmap(signal,X.samp.bedgraph);
rcontrol = listmap(control,X.samp.bedgraph);
rcluster = listmap(cluster,X.samp.bedgraph);
rheatmap = listmap(heatmap,X.samp.bedgraph);

% Read bins
fprintf('Loading bins file %s\n',bin);
fileID = fopen(bin,'r');
bed = textscan(fileID, '%s%f%f%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string',  'ReturnOnError', false);
fclose(fileID);
X.bin = struct;
X.bin.chr = convert_chr(bed{1});
X.bin.pos = bed{2}+100;
X.bin = sort_struct(X.bin,{'chr','pos'});

% Read bedGraphs
X.bin.ct_raw = zeros(slength(X.bin),slength(X.samp));
for i = 1:slength(X.samp)
    fname = X.samp.bedgraph(i);
    fprintf('Importing bedGraph %d/%d   %s\n',i,slength(X.samp),fname);
    fileID = fopen(fname,'r');
    bg = textscan(fileID, '%s%f%f%f%[^\n\r]', 'Delimiter', '\t', 'TextType', 'string', 'EmptyValue', 0, 'ReturnOnError', false);
    fclose(fileID);
    B = struct;
    B.chr = convert_chr(bg{1});
    B.pos = bg{2}+100;
    B.ct = bg{4};
    B = sort_struct(B,{'chr','pos'});
    if ~all(B.chr==X.bin.chr & B.pos==X.bin.pos), error('Bins don''t match'); end
    if ~all(B.ct==round(B.ct) & B.ct>=0), warning('Expecting integer counts data'); end
    X.bin.ct_raw(:,i) = B.ct;
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

function withoutchr = convert_chr(withchr)
chrs = ["chr1";"chr2";"chr3";"chr4";"chr5";"chr6";"chr7";"chr8";"chr9";"chr10";"chr11";"chr12";"chr13";"chr14";"chr15";"chr16";"chr17";"chr18";"chr19";"chr20";"chr21";...
        "chr22";"chrX";"chrY"];
withoutchr = listmap(withchr,chrs);
if ~all(isnan(withoutchr)), return; end

chrs = ["1";"2";"3";"4";"5";"6";"7";"8";"9";"10";"11";"12";"13";"14";"15";"16";"17";"18";"19";"20";"21";"22";"X";"Y"];
withoutchr = listmap(withchr,chrs);
if ~all(isnan(withoutchr)), return; end

error('Failed to match chromosome names to those expected for human 1-22,XY');
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

function [idx c sumd d] = kmeansds(x,k,varargin)
[idx c sumd d] = kmeansd(x,k,varargin{:});
% cluster centroids by similarity
ord = clust1d(c,'cosine');
[tmp map] = sort(ord);
map = as_column(map);
% reorder clusters by similarity
idx = map(idx);
c = c(ord,:);
sumd = sumd(ord);
d = d(:,ord);
end

function xord = clust1d(data,dist,link)
% get dendrogram order from 1-D hierarchical clustering
% dist = distance metric
% link = linkage method
if ~exist('dist','var'), dist='euclidean'; end
if ~exist('link','var'), link='complete'; end
xord = dendrogram_get_perm(linkage(pdist(data,dist),link),0);
end

function perm = dendrogram_get_perm(Z,varargin)
m = size(Z,1)+1;
if nargin < 2
    p = 0;
%    p = 30;
%    fprintf('Cutting at maximum 30 leafs.\n');
end
if nargin == 2
    p = varargin{1};
end
orientation = 't';
horz = false;
color = false;
obslabels = [];
threshold = 0.7 * max(Z(:,3));
leafOrder = [];
if nargin > 2
    if isnumeric(varargin{1})
        p = varargin{1};
        offset = 1;
    else
        p = 30;
        offset = 0;
    end
    if rem(nargin - offset,2)== 0
        error('stats:dendrogram:BadNumArgs',...
              'Incorrect number of arguments to DENDROGRAM.');
    end
    okargs = {'orientation' 'colorthreshold' 'labels','reorder'};
    for j=(1 + offset):2:nargin-2
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('stats:dendrogram:BadParameter',...
                  'Unknown parameter name:  %s.',pname);
        elseif length(k)>1
            error('stats:dendrogram:BadParameter',...
                  'Ambiguous parameter name:  %s.',pname);
        else
            switch(k)
                case 1  % orientation
                    if ~isempty(pval)
                        if ischar(pval)
                            orientation = lower(pval(1));
                        else
                            orientation = 0;    % bad value
                        end
                    end
                    if ~ismember(orientation,{'t','b','r','l'})
                        orientation = 't';
                        warning('stats:dendrogram:BadOrientation',...
                            'Unknown orientation specified, using ''top''.');
                    end
                    if ismember(orientation,{'r','l'})
                        horz = true;
                    else
                        horz = false;
                    end
                case 2  % colorthreshold
                    color = true;
                    if ischar(pval)
                        if ~strncmpi(pval,'default',length(pval))
                            warning('stats:dendrogram:BadThreshold',...
                                'Unknown threshold specified, using default');
                        end
                    end
                    if isnumeric(pval)
                        threshold = pval;
                    end
                case 3 % labels
                    if ischar(pval)
                        pval = cellstr(pval);
                    end
                    if ~iscellstr(pval)
                        error('stats:dendrogram:BadLabels',...
                            'Value of ''labels'' parameter is invalid');
                    end
                    if ~isvector(pval) || numel(pval)~=m
                        error('stats:dendrogram:InputSizeMismatch',...
                            'Must supply a label for each observation.');
                    end
                    obslabels = pval(:);
                case 4 % leaf order
                    if ~isvector(pval) || numel(pval)~=m
                        error('stats:dendrogram:InputSizeMismatch',...
                            'Leaforder is not a valid permutation.');
                    end
                    leafOrder = pval;
            end
        end
    end
end
% For each node currently labeled m+k, replace its index by
% min(i,j) where i and j are the nodes under node m+k.
Z = transz(Z);
T = (1:m)';
% If there are more than p nodes, the dendrogram looks crowded.
% The following code will make the last p link nodes into leaf nodes,
% and only these p nodes will be visible.
if (m > p) && (p ~= 0)
    Y = Z((m-p+1):end,:);         % get the last nodes
    R = unique(Y(:,1:2));
    Rlp = R(R<=p);
    Rgp = R(R>p);
    W(Rlp) = Rlp;                 % use current node number if <=p
    W(Rgp) = setdiff(1:p, Rlp);   % otherwise get unused numbers <=p
    W = W(:);
    T(R) = W(R);
    % Assign each leaf in the original tree to one of the new node numbers
    for i = 1:p
        c = R(i);
        T = clusternum(Z,T,W(c),c,m-p+1,0); % assign to its leaves.
    end
    % Create new, smaller tree Z with new node numbering
    Y(:,1) = W(Y(:,1));
    Y(:,2) = W(Y(:,2));
    Z = Y;
    m = p; % reset the number of node to be 30 (row number = 29).
end
A = zeros(4,m-1);
B = A;
n = m;
X = 1:n;
Y = zeros(n,1);
r = Y;
% arrange Z into W so that there will be no crossing in the dendrogram.
W = zeros(size(Z));
W(1,:) = Z(1,:);
nsw = zeros(n,1); rsw = nsw;
nsw(Z(1,1:2)) = 1; rsw(1) = 1;
k = 2; s = 2;
while (k < n)
    i = s;
    while rsw(i) || ~any(nsw(Z(i,1:2)))
        if rsw(i) && i == s
            s = s+1;
        end
        i = i+1;
    end

    W(k,:) = Z(i,:);
    nsw(Z(i,1:2)) = 1;
    rsw(i) = 1;
    if s == i
        s = s+1;
    end
    k = k+1;
end
g = 1;
for k = 1:m-1 % initialize X
    i = W(k,1);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
    i = W(k,2);
    if ~r(i),
        X(i) = g;
        g = g+1;
        r(i) = 1;
    end
end
% if a leaf order is specified use
if ~isempty(leafOrder)
   [dummy, X] = sort(leafOrder); %#ok
end
[u,perm]=sort(X);   %#ok perm is the third output value
return
label = num2str(perm');
if ~isempty(obslabels)
   label = cellstr(label);
   % label(:) = {''};   % to make sure non-singletons get an empty label
   singletons = find(histc(T,1:m)==1);
   for j=1:length(singletons)
      sj = singletons(j);
      label(perm==sj) = obslabels(T==sj);
   end
end
% set up the color
theGroups = 1;
groups = 0;
cmap = [0 0 1];
if color
    groups = sum(Z(:,3)< threshold);
    if groups > 1 && groups < (m-1)
        theGroups = zeros(m-1,1);
        numColors = 0;
        for count = groups:-1:1
            if (theGroups(count) == 0)
                P = zeros(m-1,1);
                P(count) = 1;
                P = colorcluster(Z,P,Z(count,1),count);
                P = colorcluster(Z,P,Z(count,2),count);
                numColors = numColors + 1;
                theGroups(logical(P)) = numColors;
            end
        end
        cmap = hsv(numColors);
        cmap(end+1,:) = [0 0 0];
    else
        groups = 1;
    end

end
if  isempty(get(0,'CurrentFigure')) || ishold
    figure('Position', [50, 50, 800, 500]);
else
    newplot;
end
col = zeros(m-1,3);
h = zeros(m-1,1);
for n = 1:(m-1)
    i = Z(n,1); j = Z(n,2); w = Z(n,3);
    A(:,n) = [X(i) X(i) X(j) X(j)]';
    B(:,n) = [Y(i) w w Y(j)]';
    X(i) = (X(i)+X(j))/2; Y(i)  = w;
    if n <= groups
        col(n,:) = cmap(theGroups(n),:);
    else
        col(n,:) = cmap(end,:);
    end
end
ymin = min(Z(:,3));
ymax = max(Z(:,3));
margin = (ymax - ymin) * 0.05;
n = size(label,1);
if(~horz)
    for count = 1:(m-1)
        h(count) = line(A(:,count),B(:,count),'color',col(count,:));
    end
    lims = [0 m+1 max(0,ymin-margin) (ymax+margin)];
    set(gca, 'Xlim', [.5 ,(n +.5)], 'XTick', 1:n, 'XTickLabel', label, ...
        'Box', 'off');
    mask = logical([0 0 1 1]);
    if strcmp(orientation,'b')
        set(gca,'XAxisLocation','top','Ydir','reverse');
    end
else
    for count = 1:(m-1)
        h(count) = line(B(:,count),A(:,count),'color',col(count,:));
    end
    lims = [max(0,ymin-margin) (ymax+margin) 0 m+1 ];
    set(gca, 'Ylim', [.5 ,(n +.5)], 'YTick', 1:n, 'YTickLabel', label, ...
        'Box', 'off');
    mask = logical([1 1 0 0]);
    if strcmp(orientation, 'l')
        set(gca,'YAxisLocation','right','Xdir','reverse');
    end
end
if margin==0
    if ymax~=0
        lims(mask) = ymax * [0 1.25];
    else
        lims(mask) = [0 1];
    end
end
axis(lims);
if nargout>0
    hout = h;
end
end

function T = clusternum(X, T, c, k, m, d)
% assign leaves under cluster c to c.
d = d+1;
n = m; flag = 0;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = clusternum(X, T, c, k, n,d); % trace back left subtree
        T = clusternum(X, T, c, X(n,2), n,d);
        flag = 1; break;
    end
end
if flag == 0 && d ~= 1 % row m is leaf node.
    T(X(m,1)) = c;
    T(X(m,2)) = c;
end
end

function T = colorcluster(X, T, k, m)
% find local clustering
n = m;
while n > 1
    n = n-1;
    if X(n,1) == k % node k is not a leave, it has subtrees
        T = colorcluster(X, T, k, n); % trace back left subtree
        T = colorcluster(X, T, X(n,2), n);
        break;
    end
end
T(m) = 1;
end

function Z = transz(Z)
m = size(Z,1)+1;
for i = 1:(m-1)
    if Z(i,1) > m
        Z(i,1) = traceback(Z,Z(i,1));
    end
    if Z(i,2) > m
        Z(i,2) = traceback(Z,Z(i,2));
    end
    if Z(i,1) > Z(i,2)
        Z(i,1:2) = Z(i,[2 1]);
    end
end
end

function a = traceback(Z,b)
m = size(Z,1)+1;
if Z(b-m,1) > m
    a = traceback(Z,Z(b-m,1));
else
    a = Z(b-m,1);
end
if Z(b-m,2) > m
    c = traceback(Z,Z(b-m,2));
else
    c = Z(b-m,2);
end
a = min(a,c);
end

function [idx,c,sumd,d] = kmeansd(x,k,varargin)
% kmeans deterministic
initcen = rand(k,size(x,2));
[idx,c,sumd,d] = kmeans(x,k,'start',initcen,varargin{:});
end

function varargout = kmeans(X, k, varargin)
% by Mathworks
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

function A = num2cellstr(a)
A = cell(length(a),1);
for i=1:length(a)
  A{i} = num2str(a(i));
end
end

function m = listmap(a,b)
na = length(a);
nb = length(b);
x = [];
x.val = [as_column(a);as_column(b)];
x.a = [true(na,1);false(nb,1)];
x.idx = as_column([1:na 1:nb]);
[tmp ui x.uj] = unique(x.val);
x = sort_struct(x,{'uj','a'});
m = nan(na,1);
aidx = nan;
bidx = nan;
for i=1:slength(x)
  if i>1 && x.uj(i)~=x.uj(i-1), aidx=nan; bidx=nan; end
  if x.a(i), aidx=x.idx(i); else bidx=x.idx(i); end
  if ~isnan(aidx), m(aidx)=bidx; end
end
end

