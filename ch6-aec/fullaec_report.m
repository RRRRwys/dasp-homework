% Partitioned block frequency domain adaptive filtering NLMS and 
% standard time-domain sample-based NLMS 
% 分块频域自适应滤波器
  fid=fopen('aecFar.pcm', 'rb'); % Load far end
  rrin=fread(fid,inf,'int16');
  fclose(fid); 
  fid=fopen('aecNear.pcm', 'rb'); % Load near end
  ssin=fread(fid,inf,'int16');
  fclose(fid);
  rand('state',13);
  fs=16000;
  mult=fs/8000;
  if fs == 8000
  	cohRange = 2:3;
  elseif fs==16000
  	cohRange = 2;
  end
% Flags
  NLPon=1;  % NLP 
  CNon=1; % Comfort noise
  M = 16; % Number of partitions 块数
  N = 64; % Partition length 块长度
  L = M*N; % Filter length 
  if fs == 8000
      mufb = 0.6;
  else
      mufb = 0.5;  
  end

% init
  %mufb=1;  
  VADtd=48;
  alp = 0.1; % Power estimation factor alc = 0.1; % Coherence estimation factor
  beta = 0.9; % Plotting factor 
  %% Changed a little %%
  step = 0.3;%0.1875; % Downward step size 
  %%
  if fs == 8000
      threshold=2e-6;  % DTrob threshold
  else
      threshold=1.5e-6; end
  if fs == 8000
      echoBandRange = ceil(300*2/fs*N):floor(1800*2/fs*N);
  else
      echoBandRange = ceil(300*2/fs*N):floor(1800*2/fs*N);
  end
  %echoBandRange = ceil(1600*2/fs*N):floor(1900*2/fs*N);
  %echoBandRange = ceil(2000*2/fs*N):floor(4000*2/fs*N);
  suppState = 1;
  transCtr = 0;
  Nt=1;
  vt=1;
  ramp = 1.0003; % Upward ramp
  rampd = 0.999; % Downward ramp
  cvt = 20; % Subband VAD threshold;
  nnthres = 20; % Noise threshold 
  shh=logspace(-1.3,-2.2,N+1)';
  sh=[shh;flipud(shh(2:end-1))]; % Suppression profile
  len=length(ssin);
  w=zeros(L,1); % Sample-based TD NLMS 
  WFb=zeros(N+1,M); % Block-based FD NLMS
  WFbOld=zeros(N+1,M); % Block-based FD NLMS
  YFb=zeros(N+1,M);
  erfb=zeros(len,1);
  erfb3=zeros(len,1);
  ercn=zeros(len,1);
  zm=zeros(N,1);
  XFm=zeros(N+1,M);
  YFm=zeros(N+1,M);
  pn0=10*ones(N+1,1);
  pn=zeros(N+1,1);
  NN=len;
  Nb=floor(NN/N)-M;
  erifb=zeros(Nb+1,1)+0.1;
  erifb3=zeros(Nb+1,1)+0.1;
  ericn=zeros(Nb+1,1)+0.1;
  dri=zeros(Nb+1,1)+0.1;
  start=1;
  xo=zeros(N,1);
  do=xo;
  eo=xo;
  echoBands=zeros(Nb+1,1);
  cohxdAvg=zeros(Nb+1,1);
  cohxdSlow=zeros(Nb+1,N+1);
  cohedSlow=zeros(Nb+1,N+1);
  %overdriveM=zeros(Nb+1,N+1);
  cohxdFastAvg=zeros(Nb+1,1);
  cohxdAvgBad=zeros(Nb+1,1);
  cohedAvg=zeros(Nb+1,1);
  cohedFastAvg=zeros(Nb+1,1);
  hnledAvg=zeros(Nb+1,1);
  hnlxdAvg=zeros(Nb+1,1);
  ovrdV=zeros(Nb+1,1);
  dIdxV=zeros(Nb+1,1);
  SLxV=zeros(Nb+1,1);
  hnlSortQV=zeros(Nb+1,1);
  hnlPrefAvgV=zeros(Nb+1,1);
  mutInfAvg=zeros(Nb+1,1);
  %overdrive=zeros(Nb+1,1);
  hnled = zeros(N+1, 1);
  weight=zeros(N+1,1);
  hnlMax = zeros(N+1, 1);
  hnl = zeros(N+1, 1);
  overdrive = ones(1, N+1);
  xfwm=zeros(N+1,M);
  dfm=zeros(N+1,M);
  WFbD=ones(N+1,1);
  fbSupp = 0;
  hnlLocalMin = 1;
  cohxdLocalMin = 1;
  hnlLocalMinV=zeros(Nb+1,1);
  cohxdLocalMinV=zeros(Nb+1,1);
  hnlMinV=zeros(Nb+1,1);
  dkEnV=zeros(Nb+1,1);
  ekEnV=zeros(Nb+1,1);
  ovrd = 2;
  ovrdPos = floor((N+1)/4);
  ovrdSm = 2;
  hnlMin = 1;
  minCtr = 0;
  SeMin = 0;
  SdMin = 0;
  SeLocalAvg = 0;
  SeMinSm = 0;
  divergeFact = 1;
  dIdx = 1;
  hnlMinCtr = 0;
  hnlNewMin = 0;
  divergeState = 0;
  Sy=ones(N+1,1);
  Sym=1e7*ones(N+1,1);
  wins=[0;sqrt(hanning(2*N-1))];
  ubufn=zeros(2*N,1);
  ebuf=zeros(2*N,1);
  ebuf2=zeros(2*N,1);
  ebuf4=zeros(2*N,1);
  mbuf=zeros(2*N,1);
  cohedFast = zeros(N+1,1);
  cohxdFast = zeros(N+1,1);
  cohxd = zeros(N+1,1);
  Se = zeros(N+1,1);
  Sd = zeros(N+1,1);
  Sx = zeros(N+1,1);
  SxBad = zeros(N+1,1);
  Sed = zeros(N+1,1);
  Sxd = zeros(N+1,1);
  SxdBad = zeros(N+1,1);
  hnledp=[];
  cohxdMax = 0;
  progressbar(0);

% Nb 是数据总块数减 16
for kk=1:Nb
  pos = N * (kk-1) + start; % 数据块首元素的位置
  
  % FD block method
  % ----------------------   Organize data
  % xk 和 dk 是读到的远端和近端 64 个时域样点
  xk = rrin(pos:pos+N-1);
  dk = ssin(pos:pos+N-1);

  xx = [xo;xk]; % 和上一轮的64点数据拼起来，用于重叠相加法
  xo = xk;
  tmp = fft(xx);
  XX = tmp(1:N+1); % fft 取 65 个频点

  dd = [do;dk];  % Overlap
	do = dk;
  tmp = fft(dd); % Frequency domain 
  DD = tmp(1:N+1); % DD 用于参与背景噪声估计
  
  % ------------------------  Power estimation 远端功率谱
  % 之前功率谱占85%
  pn0 = (1 - alp) * pn0 + alp * real(XX.* conj(XX));
  pn = pn0;
  % 背景噪声估计
	if (CNon)
		Yp =  real(conj(DD).*DD); % Instantaneous power 
		Sy =  (1 - alp) * Sy + alp * Yp; % Averaged power    

		mm = min(Sy,Sym);  
		diff = Sym - mm;
		if (kk>50)
			Sym = (mm + step*diff) * ramp; % Estimated background noise power   
		end
	end
  
  % ----------------------   Filtering   
  XFm(:,1) = XX; % XFm 65x16 是最近的16块远端频谱
  for mm=0:(M-1)
      m=mm+1; 
      YFb(:,m) = XFm(:,m) .* WFb(:,m); % YFb 是滤波器的滤波结果，WFb 是自适应滤波器的频域表示
  end
  % 将估计的频谱按列求和，包含了最近16块频谱估计信息，是回声的估计信号
  % 容许延迟 16x(64x16000)x1000 = 64 ms
  yfk = sum(YFb,2); 
  tmp = [yfk ; flipud(conj(yfk(2:N)))];
  ykt = real(ifft(tmp));
  ykfb = ykt(end-N+1:end); 
  
  % ----------------------   Error estimation 
  % 近端信号减估计回声，得到误差信号
  ekfb = dk - ykfb; 
  erfb(pos:pos+N-1) = ekfb; % for plot
  tmp = fft([zm;ekfb]);      % FD version for cancelling part (overlap-save)
  Ek = tmp(1:N+1);
  % ------------------------  Adaptation 
  Ek2 = Ek ./(M*pn + 0.001); % Normalized error

  % Limt Ek2 Amplitude 
  absEf = max(abs(Ek2), threshold);
  absEf = ones(N+1,1)*threshold./absEf;
  Ek2 = Ek2.*absEf;

  % see formula
	mEk = mufb.*Ek2;
  PP = conj(XFm).*(ones(M,1) * mEk')'; 
  tmp = [PP ; flipud(conj(PP(2:N,:)))];
	IFPP = real(ifft(tmp));
  PH = IFPP(1:N,:);
  tmp = fft([PH;zeros(N,M)]);
  FPH = tmp(1:N+1,:);
  WFb = WFb + FPH;

  % 10 * mult = 20, kk % 20 == 0
  % 20 块更新一次，权重矩阵WFb计算每个频点的功率，按列求和，找到16个块中累加和最大的索引
  if mod(kk, 10*mult) == 0
      WFbEn = sum(real(WFb.*conj(WFb)));
      [tmp, dIdx] = max(WFbEn);
      % WFbD for plot
      WFbD = sum(abs(WFb(:, dIdx)),2);
      WFbD = min(max(WFbD, 0.5), 4);
  end
  dIdxV(kk) = dIdx;
  % NLP
  if (NLPon)
    % 误差信号与之前一块拼起来
    ee = [eo;ekfb]; 
    eo = ekfb;
	  window = wins; % hanning
    if fs == 8000
      gamma = 0.9;
    else
	    gamma = 0.93;
    end
    % xx 远端 dd 近端 ee 误差
    % xf df ef 对应频谱（加窗）
  	tmp = fft(xx.*window);
  	xf = tmp(1:N+1);
  	tmp = fft(dd.*window);
  	df = tmp(1:N+1);
  	tmp = fft(ee.*window);
  	ef = tmp(1:N+1);
  
    % 维护最近16块的信号
    xfwm(:,1) = xf;
    % xf 改为dIdx所在的那一列，dIdx就是之前计算的权重矩阵能量最大的那列
    xf = xfwm(:,dIdx);
    dfm(:,1) = df;
  
    % 递归平均法计算ef,df,xf的功率谱，都是65x1
    SxOld = Sx;
  	Se = gamma*Se + (1-gamma)*real(ef.*conj(ef));
  	Sd = gamma*Sd + (1-gamma)*real(df.*conj(df));
  	Sx = gamma*Sx + (1-gamma)*real(xf.*conj(xf));
  
  	% coherence
    % Sxd 远端近端互功率谱
    Sxd = gamma*Sxd + (1 - gamma)*xf.*conj(df);
    % Sed 误差信号近端互功率谱
  	Sed = gamma*Sed + (1-gamma)*ef.*conj(df);
    % 相关性
  	cohed = real(Sed.*conj(Sed))./(Se.*Sd + 1e-10);
  	cohxd = real(Sxd.*conj(Sxd))./(Sx.*Sd + 1e-10);
  
    % 不相关性 hnled，用于衡量回声大小，回声越大 hnled 越小
    % cohxd 远端信号与近端信号相关性越大，表示包含的回声越大
    % cohed 误差信号与近端信号的相关性越小，表示包含的回声越大，考虑消干净的情况
    hnled = min(1 - cohxd, cohed);
  
    % for plot
    if kk > 1
        cohxdSlow(kk,:) = 0.99*cohxdSlow(kk-1,:) + 0.01*cohxd';
        cohedSlow(kk,:) = 0.99*cohedSlow(kk-1,:) + 0.01*(1-cohed)';
    end
  
    % echoBandRange 是回声可能处于的频率区间
    % cohedMean 平均相干性
    cohedMean = mean(cohed(echoBandRange));
    % 对 1-cohxd(echoBandRange) 排序，越靠前回声越大
    [hnlSort, hnlSortIdx] = sort(1-cohxd(echoBandRange));
    % 对远端信号功率谱 Sx 排序
    [xSort, xSortIdx] = sort(Sx);
    % 近端远端不相关性的怕平均值
    hnlSortQ = mean(1 - cohxd(echoBandRange));
    % 对 hnled 排序
    [hnlSort2, hnlSortIdx2] = sort(hnled(echoBandRange));
    % 取 hnled 1/2 和 3/4 处2个值
    hnlQuant = 0.75;
    hnlQuantLow = 0.5;
    qIdx = floor(hnlQuant*length(hnlSort2));
    qIdxLow = floor(hnlQuantLow*length(hnlSort2));
    hnlPrefAvg = hnlSort2(qIdx);
    hnlPrefAvgLow = hnlSort2(qIdxLow);
  
    % 不相关性较大时，关闭噪声抑制
    if cohedMean > 0.98 & hnlSortQ > 0.9
        suppState = 0;
    elseif cohedMean < 0.95 | hnlSortQ < 0.8
        suppState = 1;
    end
  
    % 计算一些参数
    if hnlSortQ < cohxdLocalMin & hnlSortQ < 0.75
        cohxdLocalMin = hnlSortQ;
    end
    if cohxdLocalMin == 1
        ovrd = 3;
        hnled = 1-cohxd;
        hnlPrefAvg = hnlSortQ;
        hnlPrefAvgLow = hnlSortQ;
    end
    if suppState == 0
        hnled = cohed;
        hnlPrefAvg = cohedMean;
        hnlPrefAvgLow = cohedMean;
    end
    if hnlPrefAvgLow < hnlLocalMin & hnlPrefAvgLow < 0.6
        hnlLocalMin = hnlPrefAvgLow;
        hnlMin = hnlPrefAvgLow;
        hnlNewMin = 1;
        hnlMinCtr = 0;
    end
    if hnlNewMin == 1
        hnlMinCtr = hnlMinCtr + 1;
    end
    if hnlMinCtr == 2 
        hnlNewMin = 0;
        hnlMinCtr = 0;
        ovrd = max(log(0.00001)/(log(hnlMin + 1e-10) + 1e-10), 3);
    end
    hnlLocalMin = min(hnlLocalMin + 0.0008/mult, 1);
    cohxdLocalMin = min(cohxdLocalMin + 0.0004/mult, 1);
    if ovrd < ovrdSm
        ovrdSm = 0.99*ovrdSm + 0.01*ovrd;
    else
        ovrdSm = 0.9*ovrdSm + 0.1*ovrd;
    end
    
    % 计算误差能量ekEn，近端信号能量dkEn
    ekEn = sum(Se);
    dkEn = sum(Sd);
    % 发散处理，根据ekEn 和 dkEn 大小决定是否发散和当前输出
    if divergeState == 0
        if ekEn > dkEn
            ef = df;
            divergeState = 1;
        end
    else
        if ekEn*1.05 < dkEn  
            divergeState = 0;
        else
            ef = df;
        end
    end
    
    % > 13dB认为发散，滤波器置零
    if ekEn > dkEn*19.95
        WFb=zeros(N+1,M); % Block-based FD NLMS
    end
  
    % for plot
    ekEnV(kk) = ekEn;
    dkEnV(kk) = dkEn;
    hnlLocalMinV(kk) = hnlLocalMin;
    cohxdLocalMinV(kk) = cohxdLocalMin;
    hnlMinV(kk) = hnlMin;
  
    % 利用之前的计算的参数，对hnled做尺度变换
    % 拉大近端和远端的差距
    aggrFact = 0.3;
    wCurve = [0; aggrFact*sqrt(linspace(0,1,N))' + 0.1];
    weight = wCurve;
    hnled = weight.*min(hnlPrefAvg, hnled) + (1 - weight).*hnled;
    od = ovrdSm*(sqrt(linspace(0,1,N+1))' + 1);
    sshift = ones(N+1,1);
    hnled = hnled.^(od.*sshift);
    hnl = hnled;
    ef = ef.*(hnl);
    
    % save
    ovrdV(kk) = ovrdSm;
    hnledAvg(kk) = 1-mean(1-cohed(echoBandRange));
    hnlxdAvg(kk) = 1-mean(cohxd(echoBandRange));
    hnlSortQV(kk) = hnlPrefAvgLow;
    hnlPrefAvgV(kk) = hnlPrefAvg;
  
  	% Comfort noise
  	if (CNon)
      % Sym 是之前计算的背景噪声功率谱
  		snn=sqrt(Sym);
  		snn(1)=0; % Reject LF noise
  		Un=snn.*exp(j*2*pi.*[0;rand(N-1,1);0]);
  		% Weight comfort noise by suppression
      % 误差信号抑制的越厉害，噪声越大
  		Un = sqrt(1-hnled.^2).*Un;
  		Fmix = ef + Un;
  	else
  		Fmix = ef;
  	end
  
  	% Overlap and add in time domain for smoothness 
    % 重叠相加法获得时域信号
  	tmp = [Fmix ; flipud(conj(Fmix(2:N)))];
  	mixw = wins.*real(ifft(tmp));
  	mola  = mbuf(end-N+1:end) + mixw(1:N);
  	mbuf = mixw;
  	ercn(pos:pos+N-1) = mola; 
  end % NLPon

  XFm(:,2:end) = XFm(:,1:end-1);
  YFm(:,2:end) = YFm(:,1:end-1);
  xfwm(:,2:end) = xfwm(:,1:end-1);
  dfm(:,2:end) = dfm(:,1:end-1);
end

