clear 
rng(5)
%% Introduction
% This function set the parameters in radar ROC curve simulation
%% RF setup
param.c = 3e8; % speed of the light
param.fc = 77*1e9; % central frequency, 77GHz
param.lambda= param.c/param.fc; % wavelength
%% Target properties
% Adapted from: radar/FMCWExample and phased/BasicMonostaticRadarExample
% Webs: https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
% https://www.mathworks.com/help/phased/ug/designing-a-basic-monostatic-pulse-radar.html
param.dTarget = [20+rand*70,20+rand*70]; % target distance, m
param.vTarget = [-5+rand*10,-5+rand*10]; % target velocity, m/s
param.azTarget = [-50+rand*100,-50+rand*100]; % target azimuth angle, degree 
param.elTarget = [0,0]; % target elevation angle, degree 
param.RCS = db2pow(min(10*log10(param.dTarget)+5,20)); % mean RCS, m^2 
param.numTarget = length(param.dTarget);
%% Target setup
param.target = phased.RadarTarget('MeanRCS',param.RCS,'PropagationSpeed',param.c,...
    'OperatingFrequency',param.fc);
param.targetMotion = phased.Platform('InitialPosition',...
    [param.dTarget.*cosd(param.azTarget);param.dTarget.*sind(param.azTarget);zeros(1,length(param.dTarget))],...
    'Velocity',[-param.vTarget;zeros(1,length(param.dTarget));zeros(1,length(param.dTarget))]);
%% Victim radar RF setup
% Adapted from: AutomatedDrivingRadarSimulationExample. 
% Web: https://www.mathworks.com/help/radar/ug/radar-signal-simulation-and-processing-for-automated-driving.html
% Model the antenna element
param.antElmnt = phased.IsotropicAntennaElement('BackBaffled',true);
% Construct the receive array
param.Nr = 8; % Number of receive elements
param.rxEleSpacing = param.lambda/2; % Receive element spacing
param.rxArray = phased.ULA('Element',param.antElmnt,'NumElements',param.Nr,...
    'ElementSpacing',param.rxEleSpacing);
% Construct the transmit array
param.Nt = 4; % Number of transmit elements
param.txEleSpacing = param.Nr * param.lambda/2; % Transmit element spacing
param.txArray = phased.ULA('Element',param.antElmnt,'NumElements',param.Nt,...
    'ElementSpacing',param.txEleSpacing);
% Model radar transmitter for a single transmit channel, and model a receiver preamplifier for each receive channel
param.ant_aperture = 6.06e-4;                         % antenna aperture in m^2
param.ant_gain = aperture2gain(param.ant_aperture,param.lambda);  % antenna gain in dB
param.tx_ppower = db2pow(5)*1e-3;                     % tx power in watts
param.tx_gain = 9+param.ant_gain;                           % tx gain in dB
param.rx_gain = 15+param.ant_gain;                          % rx gain in dB
param.rx_nf = 70;                                    % noise figure in dB
param.referenceTemperature = 290;                     % reference temperature of receiver in kelvin
%  System parameters            Value
%  ----------------------------------
%  Transmit power (dBm)         5
%  Transmit antenna gain (dBi)  36.0042
%  Receive antenna gain (dBi)   42.0042
%  Noise figure (dB)            4.5
%% Victim radar parameters
% Adapted from: FMCWExample. 
% Web: https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
% Parameter setup: TI AWR1243 + TDA3x processor medium-range radar
param.Bc = 460*1e6; % chirp bandwidth: 460MHz
param.hc = 15e12; % chirp slope: 15MHz/us
param.Tc = 30.7/1e6; % chirp duration (sweep time): 30.7 us
param.Ti = 7/1e6; % idle inter-chirp duration
param.Tg = param.Tc + param.Ti; % pulse repetition interval (PRI): 37.7 us
param.fPRF = ceil(1/param.Tg); % pulse repetition frequency (PRF), Hz
param.Fs = 2*param.Bc; % RF sample rate
param.fH = 15e6; % IF/LPF bandwidth (maximum beat frequency): 15MHz
param.Nsweep = 256; % number of chirps per CPI
param.dmax = param.c/2/param.hc*param.fH;  % max detection distance, m
param.vmax = param.lambda/4/param.Tg; % max velocity, m/s
param.FsADC = 2*16.7e6; % ADC sample frequency: 16.7 MSPS (complex) 
param.rangeRes = param.c/2/param.Bc; % Range resolution, m
param.velRes = param.c/2/param.fc/param.Tg/param.Nsweep; % Velocity resolution, m/s
param.Dopplermax = 0.5*param.fPRF; % max Doppler, Hz
%  System parameters            Value
%  ----------------------------------
%  Operating frequency (GHz)    77
%  Maximum target range (m)     150
%  Range resolution (m)         0.3261
%  Maximum target speed (m/s)   25.8591
%  Sweep time (us)              30.7
%  Sweep bandwidth (MHz)        460
%  Maximum beat frequency (MHz) 15
%% Victim radar system setup
% Waveform transmitter
param.transmitter = phased.Transmitter('PeakPower',param.tx_ppower,'Gain',param.tx_gain);
% Radiator for single transmit array
param.radiator = phased.Radiator('Sensor',param.txArray,'OperatingFrequency',param.fc,'WeightsInputPort',true);
% Collector for receive array
param.collector = phased.Collector('Sensor',param.rxArray,'OperatingFrequency',param.fc);
% Receiver preamplifier
param.receiver = phased.ReceiverPreamp('Gain',param.rx_gain,'NoiseFigure',param.rx_nf,...
    'SampleRate',param.Fs,'ReferenceTemperature',param.referenceTemperature);
%% Victim radar movtion and channel setup
victim_speed = 0; % suppose victim radar speed is 0 for simplicity
param.victimMotion = phased.Platform('InitialPosition',[0;0;0],...
    'Velocity',[victim_speed;0;0]);
param.targetChannel = phased.FreeSpace('PropagationSpeed',param.c,...
    'OperatingFrequency',param.fc,'SampleRate',param.Fs,'TwoWayPropagation',true);
%% Victim radar reference signal
% Construct victim radar's reference data 
% with fast time (time within each pulse) along each column and slow time (time between pulses) along each row
refWav = phased.FMCWWaveform('SampleRate',param.Fs,'SweepTime',param.Tc,'SweepBandwidth',param.Bc);
param.Nc = ceil(param.Fs*param.Tc); % number of samples per chirp
param.Ng = ceil(param.Fs*param.Tg); % number of samples per repetition period
param.victRef = complex(zeros(param.Ng,param.Nsweep)); % victim radar's reference waveform
for l = 1:param.Nsweep
    % Victim radar's reference FMCW waveform
    xVictRef = [refWav();zeros(param.Ng-param.Nc,1)]; % victim reference signal in a chirp cycle
    param.victRef(:,l) = xVictRef;
end
%% TX-target-RX-dechirping process at victim radar
% Adapted from: FMCWExample and MIMORadarVirtualArrayExample
% Web: https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
% Web: https://www.mathworks.com/help/phased/ug/increasing-angular-resolution-with-mimo-radars.html
param.dechirpTar = complex(zeros(param.Ng,param.Nr,param.Nsweep)); % dechipred target's data cube
w = complex(zeros(param.Nt,param.Nsweep));
param.slowTimeCodeType = "Chu";
for l = 1:param.Nsweep
    % Update victim radar, target, and interferer positions
    [victim_pos,victim_vel] = param.victimMotion(param.Tg);
    [tgt_pos,tgt_vel] = param.targetMotion(param.Tg);
    [~,tgt_ang] = rangeangle(tgt_pos,victim_pos);
    % Victim radar's reference FMCW waveform
    victRefSig = param.victRef(:,l);
    % Victim radar's transmitted FMCW waveform
    txVictSig = param.transmitter(victRefSig); % transmitted signal in a chirp cycle
    % Phase code added on each element
    wl = zeros(param.Nt,1);
    for mt = 1: param.Nt
        switch param.slowTimeCodeType
            case "Chu" % Chu sequence
                wl(mt) = exp(1j*pi*(l-1)*l*mt/257); % hard coded to be a prime number larger than param.Nsweep = 256
            case "P" % P sequence
                wl(mt) = pSequence(mt,l,257); % hard coded to be a prime number larger than param.Nsweep = 256
            case "Binary" % Binary random code (+1 or -1)   
                wl(mt) = randi([0,1])*2-1; 
            case "TDM" % TDM, time division multiplexing
                wl(mt) = int8(mt == mod(l-1,param.Nt)+1); 
        end
    end
    w(:,l) = wl;
    % Victim radar's radiator that converts signals into radiated
    % wavefields from arrays and individual sensor elements 
    txVictSig = param.radiator(txVictSig,tgt_ang,wl); 
    % Propagate the signal and reflect from the target
    txVictSig = param.targetChannel(txVictSig,victim_pos,tgt_pos,victim_vel,tgt_vel);
    txVictSig = param.target(txVictSig); % target echo
    % Dechirp the received target return
    tarRxSig = param.collector(txVictSig,tgt_ang);
    tarRxSig = param.receiver(tarRxSig);  
    dechirpTarSig = dechirp(tarRxSig,victRefSig); % dechirping conjugates the target signal
    for m = 1:param.Nr
        param.dechirpTar(:,m,l) = dechirpTarSig(:,m); % dechirped target echo in a chirp cycle
    end
end
param.w = w;
%% Interfering radar channel setup
param.intChannel = phased.FreeSpace('PropagationSpeed',param.c,...
    'OperatingFrequency',param.fc,'SampleRate',param.Fs,'TwoWayPropagation',false);
%% Interfering radar RF setup
% Adapted from: AutomatedDrivingRadarSimulationExample. 
% Web: https://www.mathworks.com/help/radar/ug/radar-signal-simulation-and-processing-for-automated-driving.html
% Construct the receive array
param.NrInt = 2; % Number of receive elements
param.rxEleSpacingInt = param.lambda/2; % Receive element spacing
% Construct the transmit array
param.NtInt = 8; % Number of transmit elements
param.txEleSpacingInt = param.NrInt * param.lambda/2; % Transmit element spacing
param.txArrayInt = phased.ULA('Element',param.antElmnt,'NumElements',param.NtInt,...
    'ElementSpacing',param.txEleSpacingInt);
%% Interfering radar system setup
% Waveform transmitter
param.transmitterInt = phased.Transmitter('PeakPower',param.tx_ppower,'Gain',param.tx_gain);
% Radiator for single transmit array
param.radiatorInt = phased.Radiator('Sensor',param.txArrayInt,'OperatingFrequency',param.fc,'WeightsInputPort',true);
%% Simulation under different interference parameters
param.MinIntDist = 2; % Minimum interference range
param.numInt = 1; % Number of interference
param.targetIdx = 2;
param.numIterSim = 1000;
param = getDetectResult(param);
%% RA Interference power CDF
fftIntPowSet = param.fftIntPowSet;
fastClairvoyantIntPowSet = param.fastClairvoyantIntPowSet;
fastRsIntPowSet = param.fastRsIntPowSet;
fastGsIntPowSet = param.fastGsIntPowSet;
fastLcmvIntPowSet = param.fastLcmvIntPowSet;
agsIntPowSet = param.agsIntPowSet;
save cdfv.mat fftIntPowSet fastClairvoyantIntPowSet fastRsIntPowSet fastGsIntPowSet fastLcmvIntPowSet agsIntPowSet
figure
[fftCDFy,fftCDFx] = ecdf(param.fftIntPowSet);
[clairvoyantCDFy,clairvoyantCDFx] = ecdf(param.fastClairvoyantIntPowSet);
[rsCDFy,rsCDFx] = ecdf(param.fastRsIntPowSet);
[gsCDFy,gsCDFx] = ecdf(param.fastGsIntPowSet);
[lcmvCDFy,lcmvCDFx] = ecdf(param.fastLcmvIntPowSet);
[agsCDFy,agsCDFx] = ecdf(param.agsIntPowSet);
plot(fftCDFx,fftCDFy)
hold on
plot(clairvoyantCDFx,clairvoyantCDFy)
hold on
plot(rsCDFx,rsCDFy)
hold on
plot(gsCDFx,gsCDFy)
hold on
plot(lcmvCDFx,lcmvCDFy)
hold on
plot(agsCDFx,agsCDFy)
xlabel('Output Interference Power (dB)')
ylabel('CDF')
legend('Angle FFT','Clairvoyant','RS','GS','LCMV-SMI','AGS')