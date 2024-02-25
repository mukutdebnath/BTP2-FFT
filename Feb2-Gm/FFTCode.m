close all; clear;

project_title = "Gm";
datasetVIN = readtable("GMVINdata.csv");
dataVINTime = datasetVIN.VIN_DIFFX;
dataVIN = datasetVIN.VIN_DIFFY;

datasetIOUT = readtable("GMVOUTdata.csv");
dataIOUTTime = datasetIOUT.IOUT_DIFFX;
dataIOUT = datasetIOUT.IOUT_DIFFY;

fs = 100e3;
scount = 1024; %samples/freq
startTime = 10e-3;

ifsavefig = true;
[fftin, frequenciesin, fftin_measurements] = find_fft(fs, scount, startTime, dataVINTime, ...
    dataVIN, 'FFT V_{in}', 'GMVIN_MATLAB.png', ifsavefig);
[fftout, frequenciesout, fftout_measurements] = find_fft(fs, scount, startTime, dataIOUTTime, ...
    dataIOUT, 'FFT I_{out}', 'GMIOUT_MATLAB.png', ifsavefig);

fileID = fopen("log_FFTCode.txt",'w');
fprintf(fileID, '%s\n', 'Generated from FFTCode.m');
fprintf(fileID, '%s\n', 'Title: '+project_title);
fprintf(fileID, '\n%s\n', 'Input FFT Measurements');
dict_print(fileID, fftin_measurements);
fprintf(fileID, '\n%s\n', 'Output FFT Measurements');
dict_print(fileID, fftout_measurements);
fclose(fileID);

function [] = dict_print(fileID, dict)
    dict_keys = keys(dict);
    fprintf(fileID,'%25s \t %12s\r\n','Specification','Value');
    for i = 1:length(dict_keys)
        fprintf(fileID, '%25s \t %12.8f\r\n', dict_keys(i), dict(dict_keys(i)));
    end
end

function [fftmag_se, frequencies, output_measurements] = find_fft(fs, scount, startTime, time, ...
    Y, plotTitle, saveFileName, ifsave)
    endTime = startTime + scount/fs;
    indices = time >= startTime & time <= endTime;
    YSelected = Y(indices);
    L = length(YSelected);

    fftres = fft(YSelected);
    fftmag = abs(fftres/L);
    fftmag_se = fftmag(1:L/2+1);
    fftmag_se(2:end-1) = 2*fftmag_se(2:end-1);
    fftmag_se_dB = 20*log10(abs(fftmag_se));

    deltaf = fs/scount;
    frequencies = 0:deltaf:fs/2;

    figure('Position', [0, 0, 800, 400]);
    plot(frequencies, fftmag_se_dB);
    title(plotTitle);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    grid on;
    if ifsave
        saveas(gcf, saveFileName)
    end
    
    % fft measurements
    fundamental_index = fftmag_se == max(fftmag_se);
    fundamental_frequency = frequencies(fundamental_index);
    tones = fundamental_frequency:fundamental_frequency:max(frequencies);
    harmonic_frequencies = tones(2:end);
    harmonic_indices = ismember(frequencies, harmonic_frequencies);
    harmonic_ffts = fftmag_se(harmonic_indices);
    max_spur = max(harmonic_ffts) ^ 2;
    harmonic_power = sum(harmonic_ffts .^2);
    signalPower = max(fftmag_se) ^ 2;
    totalPower = sum(fftmag_se .^ 2);
    noisePower = totalPower - signalPower - harmonic_power;
    SFDR = 10*log10(signalPower/max_spur);
    SNR = 10*log10(signalPower/noisePower);
    SINAD = 10*log10((signalPower + noisePower + harmonic_power)/(noisePower + harmonic_power));
    THD = sqrt(harmonic_power/signalPower) * 100;
    THDdB = 10*log10(harmonic_power/signalPower);
    ENOB = (SINAD - 10*log10(3/2)) / (20*log10(2));
    dict_keys = ["Fundamental Frequency" "DC Power" "Signal Power" "SFDR" "SNR" "THD" "THDdB" "SINAD" "ENOB"];
    dict_values = [fundamental_frequency fftmag_se_dB(1) max(fftmag_se_dB) SFDR SNR THD THDdB SINAD ENOB];
    output_measurements = dictionary(dict_keys, dict_values);
end