%% Library

function bits = int_to_bit(int, bit_num) % array -> rowwise bits
    int = int';
    bits = zeros(length(int), bit_num);
    for i_regnum = 1:bit_num
        bits(:, i_regnum) = rem(int, 2);
        int = floor(int/2);
    end
end

function ints = bit_to_int(bits) % rowwise bits -> array (rowwise)
    bits_length = length(bits);
    bits_weight = 2.^(0:(bits_length-1));
    ints = bits_weight * bits';
end

function [gm, constellation, bit_to_symbols_lookup_idxp] = graymap(M)  
    M_sqrt = sqrt(M);
    g1d = bitxor(0:(M_sqrt-1), bitshift(0:(M_sqrt-1), -1)); % 1d graymap
    gm = g1d + g1d'*M_sqrt; % 2d graymap(i, j) = cat[1d graymap(i), 1d graymap(j)]

    min_dis = sqrt(6*log2(M)/(M-1));

    constellation1d = ((1:M_sqrt)-(1+M_sqrt)/2)*min_dis;
    constellation = constellation1d + 1j * constellation1d';

    bit_to_symbols_lookup_idxp = zeros(M,1);
    for i = 1:M_sqrt
        for j = 1:M_sqrt
            bit_to_symbols_lookup_idxp(gm(i, j)+1) = constellation(i, j);
        end
    end
    
    % disp(constellation)
    % disp(gm)
    % disp(bit_to_symbols_lookup_idxp)
end
%% Encoders

function encoded_bit = Convolution_Encoder(input_bit, convolution_rules)
    % code : 1x~

    % 필요 수치 정리
    input_length = length(input_bit); % k
    n_over_k = length(convolution_rules); % n/k, codeword length
    
    % register 정의
    reg_num = floor(log2(max(convolution_rules)));
    reg = zeros(reg_num,1); % columnwise
    disp([num2str(convolution_rules) ' convolution, # of registers:' num2str(reg_num)]);
    
    % convolution 연산을 위한 rule bit화
    convolution_rules_bits = zeros([length(convolution_rules),reg_num+1]); % rowwise
    for i_conv = 1:length(convolution_rules)
        rule = convolution_rules(i_conv);
        for i_regnum = 1:reg_num+1
            convolution_rules_bits(i_conv, i_regnum) = rem(rule, 2); 
            rule = floor(rule/2);
        end
    end

    % convolution encoding
    encoded_bit = zeros(1, input_length*n_over_k);
    for i_bit = 1:input_length
        bit_sequence = [input_bit(i_bit); reg];
        encoded_bit((i_bit-1)*n_over_k+1:i_bit*n_over_k) = rem(convolution_rules_bits*bit_sequence, 2);
    
        reg = bit_sequence(1:reg_num);
    end

    disp([num2str(input_length) ' bits encoded (' num2str(n_over_k*input_length) ' bits)']);
end

%% Modulator

function modulated_symbols = moudulator_M_ary_QAM(encoded_bits, M)

    encoded_bits_length = length(encoded_bits); % n
    bit_per_symbols = floor(log2(M));
    symbols_length = encoded_bits_length/bit_per_symbols;
    if floor(bit_per_symbols/2) ~= log2(M)/2
        error('M must be power of power of 2^2');
    end

    [~, ~, bit_to_symbols_lookup_idxp] = graymap(M);
    modulated_symbols = zeros(1,symbols_length);
    bit_sampled = reshape(encoded_bits, bit_per_symbols,symbols_length)';
    
    for i_sym = 1:symbols_length
        symbol_bits = bit_to_int(bit_sampled(i_sym,:));
        modulated_symbols(i_sym) = bit_to_symbols_lookup_idxp(symbol_bits+1);
    end

    disp([num2str(symbols_length) ' symbols modulated (' num2str(M) '-ary QAM)']);
end

%% AWGN Channel

function noised_symbols = AWGN_channel(symbols, SNR, n_over_k)
    disp('AWGN channel added')

    SNR_linear = 10^(SNR/10);

    sigma = sqrt(0.5/SNR_linear*n_over_k);
    noise_bit = randn([1, length(symbols)]) * sigma + 1j*randn([1, length(symbols)]) * sigma;

    noised_symbols = symbols + noise_bit;
end

%% Demodulator

function demodulated_bits = hard_demoudulator_M_ary_QAM(symbols, M)
    sym_length = length(symbols);
    symbol_bit_num = log2(M);

    disp([num2str(sym_length) ' symbols demodulated (' num2str(sym_length*symbol_bit_num) 'bits)']);

    min_dis = sqrt(6*sqrt(M)/(M-1));
    [gray_map, constellation, bit_to_symbols_lookup_idxp] = graymap(M);
    symbols = symbols/min_dis;
    symbols = symbols+0.5+0.5j;
    symbols = round(symbols);
    symbols = symbols-0.5-0.5j;
    symbols = symbols-constellation(1,1)/min_dis+1+1j;

    
    M_sqrt = sqrt(M);
    demodulated_bits = zeros(1, sym_length*symbol_bit_num);
    for idx_sym = 1:sym_length
        sym = symbols(idx_sym);
        i = max(1, min(sqrt(M), round(imag(sym))));
        j = max(1, min(sqrt(M), round(real(sym))));
        

        demodulated_bits((idx_sym-1)*symbol_bit_num+1:idx_sym*symbol_bit_num) = int_to_bit(gray_map(i, j), symbol_bit_num);
    end
end

function LLRs = soft_demoudulator_M_ary_QAM(symbols, M, SNR, n_over_k)
    sym_length = length(symbols);
    SNR_linear = 10^(SNR/10);
    sigma = sqrt(0.5/SNR_linear*n_over_k);
    M_log2 = log2(M);

    %min_dis = sqrt(6*sqrt(M)/(M-1));
    [gray_map, constellation, ~] = graymap(M);
    
    LLRs = zeros(1, sym_length * M_log2);
    for i_sym = 1:sym_length
        %symbol = symbols(i_sym);
        symbols_distances = abs(constellation - symbols(i_sym));
        symbols_distance_square = symbols_distances .* conj(symbols_distances);
        
        sb0 = zeros(1, M_log2);
        sb1 = zeros(1, M_log2);

        for i_gm = 1:M
            sym = gray_map(i_gm);
            symbol_bits = int_to_bit(sym, log2(M));
            for i_sb = 1:M_log2
                if symbol_bits(i_sb) == 0
                    sb0(i_sb) = sb0(i_sb) + exp(-symbols_distance_square(i_gm)/(2*sigma^2));
                else
                    sb1(i_sb) = sb1(i_sb) + exp(-symbols_distance_square(i_gm)/(2*sigma^2));
                end
            end
        end
        LLRs(((i_sym-1)*M_log2+1):i_sym*M_log2) = log(sb1./sb0);
    end
                
end

%% Decoder

function received_coded_bits = Viterbi_Decoder(LLRs, convolution_rules, mode)
    reg_num = floor(log2(max(convolution_rules)));
    N = length(convolution_rules); %codeword size
    branches = zeros(2^reg_num, 2);

    disp([num2str(length(LLRs)/N) ' bits decoded (' mode ' decoding)']);

    % convolution 연산을 위한 rule bit화
    convolution_rules_bits = zeros([length(convolution_rules),reg_num+1]); % rowwise
    for i_conv = 1:length(convolution_rules)
        rule = convolution_rules(i_conv);
        for i_regnum = 1:reg_num+1
            convolution_rules_bits(i_conv, i_regnum) = rem(rule, 2); 
            rule = floor(rule/2);
        end
    end

    for state = 0:2^reg_num-1
        for bit = 0:1
            bsq = [bit, int_to_bit(state,reg_num)];
            rem(convolution_rules_bits*bsq', 2)';
            codeword = bit_to_int(rem(convolution_rules_bits*bsq', 2)');
            branches(state+1, bit+1) = (bit_to_int(bsq(1:reg_num))+1)*2^N+codeword;
        end
    end

    % branches

    % branches = [
    %     1*4+0,2*4+3;
    %     3*4+1,4*4+2;
    %     1*4+3,2*4+0;
    %     3*4+2,4*4+1
    %     ]; %branches(state+1, input+1) = (next_state+1)*4+codeword

    state_num = size(branches, 1);
    PM = zeros(state_num, length(LLRs)/N+1)+inf;
    PM(1,1) = 0;
    PM_input = zeros(state_num, length(LLRs)/N);
    PM_pre_state = zeros(state_num, length(LLRs)/N);

    for idx = 1:(length(LLRs)/N)
        LLR = LLRs((idx-1)*N+1:idx*N);
        for state = 1:state_num
            for bit = 0:1
                temp = branches(state, bit+1);
                next_state = (temp-mod(temp, (2^N)))/(2^N);
                if idx == (length(LLRs)/N) && next_state ~= 1 % 마지막 state 0으로 끝내기 위함
                    continue
                end
                
                codeword = int_to_bit(mod(temp, (2^N)),N);
                distance = 0;
                if mode == 'hard'
                    distance = sum(abs((LLR>0)*1-codeword));
                elseif mode == 'soft'
                    distance = sum((0.5-codeword).*LLR);
                end

                PM_next = PM(state, idx) + distance;
                if PM(next_state, idx+1) > PM_next
                    PM(next_state, idx+1) = PM_next;
                    PM_input(next_state, idx) = bit;
                    PM_pre_state(next_state, idx) = state;
                end
            end
        end
    end
    
    received_coded_bits = zeros(1, length(LLRs)/N);
    min_state = 1; %% 마지막 state 00
    for idx = (length(LLRs)/N):-1:1
        received_coded_bits(idx) = PM_input(min_state, idx);
        min_state = PM_pre_state(min_state, idx);
    end

    
end

%% Simulation BER
function Simulation_BER(bits_num, convolution_rules, QAM, SNRs, mode)
    reg_num = floor(log2(max(convolution_rules)));

    BERs = zeros(1, length(SNRs));
    for i = 1:length(SNRs)
        clc;
        SNR = SNRs(i);
        disp(['SNR:' num2str(SNRs(i)) 'dB'])
    
        bits = randi(2, [1,bits_num])-1;
        bits((bits_num-reg_num+1):bits_num) = 0; % 마지막 state 0으로 끝내기 위함
        if mode == 'unco'
            encoded_bits = bits;
        else
            encoded_bits = Convolution_Encoder(bits, convolution_rules);
        end
        modulated_symbols = moudulator_M_ary_QAM(encoded_bits, QAM);
        noised_symbols = AWGN_channel(modulated_symbols, SNR, length(convolution_rules));

        if mode == 'unco'
            demodulated_bits = hard_demoudulator_M_ary_QAM(noised_symbols, QAM);
            decoded_bits = demodulated_bits;
        else
            LLRs = soft_demoudulator_M_ary_QAM(noised_symbols, QAM, SNR, length(convolution_rules));
            decoded_bits = Viterbi_Decoder(LLRs, convolution_rules, mode);
        end
        
        
        
        BERs(i) = sum(abs(bits - decoded_bits))/bits_num;
    end
    
    semilogy(SNRs, BERs)
end

%% Simulation FER
function FERs = Simulation_FER(frames_num, bits_num, convolution_rules, QAM, SNRs, mode)
    reg_num = floor(log2(max(convolution_rules)));

    FERs = zeros(1, length(SNRs));
    for i = 1:length(SNRs)
        for i_fra = 1:frames_num
            clc;
            SNR = SNRs(i);
            disp(['SNR:' num2str(SNRs(i)) 'dB | ' num2str(i_fra) '/' num2str(frames_num) ' frames'])
        
            bits = randi(2, [1,bits_num])-1;
            bits((bits_num-reg_num+1):bits_num) = 0; % 마지막 state 0으로 끝내기 위함
            if mode == 'unco'
                encoded_bits = bits;
            else
                encoded_bits = Convolution_Encoder(bits, convolution_rules);
            end
            modulated_symbols = moudulator_M_ary_QAM(encoded_bits, QAM);
            noised_symbols = AWGN_channel(modulated_symbols, SNR, length(convolution_rules));
    
            if mode == 'unco'
                demodulated_bits = hard_demoudulator_M_ary_QAM(noised_symbols, QAM);
                decoded_bits = demodulated_bits;
            else
                LLRs = soft_demoudulator_M_ary_QAM(noised_symbols, QAM, SNR, length(convolution_rules));
                decoded_bits = Viterbi_Decoder(LLRs, convolution_rules, mode);
            end
            
            
            if sum(abs(bits - decoded_bits)) ~= 0
                FERs(i) = 1+FERs(i);
            end
        end
        FERs(i) = FERs(i) / frames_num;
    end

    R = 1/length(convolution_rules);

    throughputs = R .* log2(QAM) .* (1 - FERs);
    
    semilogy(SNRs, throughputs)
end

%% Problem 1, BER

clear; clc
D=2;

figure();
hold on

% uncoded Viterbi Decoding
bits_num = 10^6;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 4;
SNRs = 0:7;
mode = 'unco';

Simulation_BER(bits_num, convolution_rule, QAM, SNRs, mode);

% hard Viterbi Decoding
bits_num = 10^5;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 4;
SNRs = 0:7;
mode = 'hard';

Simulation_BER(bits_num, convolution_rule, QAM, SNRs, mode);

% soft Viterbi Decoding
bits_num = 10^5;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 4;
SNRs = 0:7;
mode = 'soft';

Simulation_BER(bits_num, convolution_rule, QAM, SNRs, mode);

legend('uncoded', 'hard 79 109 convolution', 'soft 79 109 convolution')
xlabel('Eb/N0')
ylabel('BER')
set(gca, 'YScale', 'log');
grid on
hold off

%% Problem 2, FER

clear; clc
D=2;

figure();
hold on

frames_num = 10^4;
bits_num = 1024;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 4;
SNRs = 0:5;
mode = 'soft';


FERs = Simulation_FER(frames_num, bits_num, convolution_rule, QAM, SNRs, mode);

legend('soft 79 109 convolution')
xlabel('Eb/N0')
ylabel('η_e_f_f')
set(gca, 'YScale', 'linear');
grid on
hold off

%% Problem 3, Adaptive Modulation and Coding (AMC) Efficiency

clear; clc
D=2;

figure();
hold on

% QPSK
frames_num = 10^2;
bits_num = 512*3;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 4;
SNRs = 0:0.5:10;
mode = 'soft';

FERs_QPSK = Simulation_FER(frames_num, bits_num, convolution_rule, QAM, SNRs, mode);

% 16-QAM
frames_num = 10^2;
bits_num = 512*3;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 16;
mode = 'soft';

FERs_16QAM = Simulation_FER(frames_num, bits_num, convolution_rule, QAM, SNRs, mode);

% 64-QAM
frames_num = 10^2;
bits_num = 512*3;
convolution_rule = [1+D+D^2+D^3+D^6, 1+D^2+D^3+D^5+D^6];
QAM = 64;
mode = 'soft';

FERs_64QAM = Simulation_FER(frames_num, bits_num, convolution_rule, QAM, SNRs, mode);

title('soft 79 109 convolution')
legend('QPSK', '16-QAM', '64-QAM')
xlabel('Eb/N0')
ylabel('η_e_f_f')
xlim([0,10])
set(gca, 'YScale', 'linear');
grid on
hold off

% FER plot
figure();
hold on
semilogy(SNRs, FERs_QPSK, SNRs, FERs_16QAM, SNRs, FERs_64QAM);
title('soft 79 109 convolution')
legend('QPSK', '16-QAM', '64-QAM')
xlabel('Eb/N0')
ylabel('FER')
xlim([0,10])
set(gca, 'YScale', 'log');
grid on
hold off
