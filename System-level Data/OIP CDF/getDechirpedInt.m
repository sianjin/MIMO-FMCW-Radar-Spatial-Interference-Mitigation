function dechirpInt = getDechirpedInt(param)
%% Interfering radar reference signal
% Construct interferer's reference data
% with fast time (time within each pulse) along each column and slow time (time between pulses) along each row
intWav = phased.FMCWWaveform('SampleRate',param.Fs,'SweepTime',param.TcInt,'SweepBandwidth',param.Bc);
NgInt = ceil(param.Fs*param.TgInt); % number of samples per repetition period of interference
NoInt = ceil(param.Fs*param.toffInt); % number of samples corresponding to interference
NsweepInt = ceil(param.Tg*param.Nsweep/param.TgInt); % max number of interference chirps in a victim radar's PRI
intRefv = zeros(NoInt + NsweepInt*NgInt,1); % interference reference waveform with multi chirps in a vector
for iter = 1:NsweepInt
    intRefv(NoInt+(iter-1)*NgInt+1:NoInt+iter*NgInt) = [intWav();zeros(NgInt-length(intWav()),1)]; % interference reference signal in a chirp cycle
end
% Rearrange interferer's reference data using victim radar's reference fast
% time and slow time frame
intRefv = intRefv(1:param.Ng*param.Nsweep); % cut the length of interference signals into the length of victim radar signals
intRef = reshape(intRefv,param.Ng,[]); % reshape into matrix: row length = param.Ng, col length = param.Nsweep
%% Interfer TX-victim RX-victim dechirping process
% Adapted from: FMCWExample.
% Web: https://www.mathworks.com/help/radar/ug/automotive-adaptive-cruise-control-using-fmcw-technology.html
dechirpInt = complex(zeros(param.Ng,param.Nr,param.Nsweep)); % dechipred interferer's data cube
for l = 1:param.Nsweep
    % Update victim radar, target, and interferer positions
    [victim_pos,victim_vel] = param.victimMotion(param.Tg);
    [int_pos,int_vel] = param.intMotion(param.Tg);
    [~,int_ang] = rangeangle(int_pos,victim_pos);
    % Victim radar's reference FMCW waveform
    victRefSig = param.victRef(:,l);
    % Interfering radar's reference FMCW waveform
    intRefSig = intRef(:,l);
    % Interfering radar's transmitted FMCW waveform
    txIntSig = param.transmitterInt(intRefSig); % transmitted signal in a chirp cycle
    % Phase code added on each element
    wl = zeros(param.NtInt,1);
    mtSet = randperm(257,param.NtInt);
    for mt = 1: param.NtInt
        switch param.slowTimeCodeType
            case "Chu" % Chu sequence
                wl(mt) = exp(1j*pi*(l-1)*l*mtSet(mt)/257); % hard coded to be a prime number larger than param.Nsweep = 256
            case "P" % P sequence
                wl(mt) = pSequence(mtSet(mt),l,257); % hard coded to be a prime number larger than param.Nsweep = 256
            case "Binary" % Binary random code (+1 or -1)   
                wl(mt) = randi([0,1])*2-1; 
            case "TDM" % TDM, time division multiplexing
                wl(mt) = int8(mt == mod(l-1,param.NtInt)+1);
        end
    end
    % Victim radar's radiator that converts signals into radiated
    % wavefields from arrays and individual sensor elements 
    txIntSig = param.radiatorInt(txIntSig,int_ang,wl);
    % Propagate the signal and reflect off the target
    txIntSig = param.intChannel(txIntSig,int_pos,victim_pos,int_vel,victim_vel);
    % Dechirp the received radar return
    intRxSig = param.collector(txIntSig,int_ang);
    intRxSig = param.receiver(intRxSig);
    dechirpIntSig = dechirp(intRxSig,victRefSig);
    for m = 1:param.Nr
        dechirpInt(:,m,l) = dechirpIntSig(:,m); % dechirped interference in a chirp cycle
    end
end
end