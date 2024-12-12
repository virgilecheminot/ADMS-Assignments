function [FRFopt, freqsOpt, xSol] = FRF_fit(freqs, FRF, nModes, freqRange, errType, plotBool)
    if (~exist('errType', 'var'))
        errType = 'complexDiff';
    end
    if (~exist('plotBool', 'var'))
        plotBool = false;
    end

    df = freqs(2) - freqs(1);
    windowSize = round(freqRange / df);

    FRFopt = zeros(windowSize, size(FRF, 2), nModes);
    freqsOpt = zeros(windowSize, nModes);
    xSol = zeros(3*size(FRF, 2)+2, nModes);

    f = freqs;
    nFRF = size(FRF, 2);

    % Find the peaks of the FRFs
    [pks, locs] = findpeaks(abs(max(FRF, [], 2)), 'SortStr', 'descend', 'NPeaks', nModes);
    % sort peaks according from lower to higher frequency
    [~, idx] = sort(f(locs));
    locs = locs(idx);
    omega_i = 2 * pi * f(locs);
    pks = pks(idx);

    % Plot the experimental FRFs
    if (plotBool)
        figure;
        subplot(2, 1, 1);
        semilogy(f, abs(FRF));
        hold on;
        plot(f(locs), pks, 'ro', 'LineWidth', 1.5);
        title(['Number of FRFs = ', num2str(size(FRF, 2))]);
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        grid on;
        subplot(2, 1, 2);
        plot(f, angle(FRF));
        xlabel('Frequency (Hz)');
        ylabel('Phase');
        grid on;
    end

    disp('Identified mode frequencies:');
    for i = 1:nModes
        disp(['Mode ', num2str(i), ': ', num2str(f(locs(i))), ' Hz']);
    end

    % Identification of the modal parameters
    function FRFs = numericalFRFs(freqs, nFRF, x)
        FRFs = zeros(length(freqs), nFRF);
        omega_vals = 2 * pi * freqs;
        omega_i = x(1);
        df_i = x(2);
        A_i = x(3:nFRF+2);
        RL_i = x(nFRF+3:2*nFRF+2);
        RH_i = x(2*nFRF+3:3*nFRF+2);
        % Main loop
        for k = 1:length(freqs)
            omega = omega_vals(k);
            for l = 1:nFRF
                den = -omega^2 + 2j * df_i * omega_i * omega + omega_i^2;
                FRFs(k, l) = FRFs(k, l) + A_i(l) / den + RL_i(l) / omega^2 + RH_i(l);
            end
        end
    end

    for modeNo = 1:nModes
        % Frequency range
        modeIndex = locs(modeNo);
        fCenter = f(modeIndex);
        iMin = modeIndex - round(windowSize/2);
        if (iMin < 1)
            iMin = 1;
            if f(1) == 0
                iMin = 2;
            end
        end
        iMax = iMin + windowSize-1;
        if (iMax > length(f))
            iMax = length(f);
        end
        indexRange = iMin:iMax;
        freqsInRange = f(indexRange);
        freqsOpt(:, modeNo) = freqsInRange;
        disp(['--- Mode number ', num2str(modeNo), ':']);
        disp(['Frequency range: ', num2str(f(iMin)), ' Hz to ', num2str(f(iMax)), ' Hz']);

        % First guess of the modal parameters
        omegaGuess = 2 * pi * fCenter;

        [~, r] = max(abs(FRF(modeIndex, :)));
        phaseDeriv = (angle(FRF(modeIndex+1, r)) - angle(FRF(modeIndex-1, r)))/((f(2)-f(1))*4*pi);
        dampFactorGuess = -1/(omegaGuess * phaseDeriv);

        ArGuess = zeros(nFRF, 1);
        for i = 1:nFRF
            ArGuess(i) = real(FRF(modeIndex, i) * (2j * dampFactorGuess * omegaGuess^2));
        end

        RLGuess = zeros(nFRF, 1);
        RHGuess = zeros(nFRF, 1);

        x0 = [omegaGuess; dampFactorGuess; ArGuess; RLGuess; RHGuess];

        % Optimization
        if strcmp(errType, 'complexDiff')
            err = @(x) abs(reshape(FRF(indexRange, :) - numericalFRFs(freqsInRange, nFRF, x), [], 1));
        elseif strcmp(errType, 'absDiff')
            err = @(x) (reshape(abs(FRF(indexRange, :) - numericalFRFs(freqsInRange, nFRF, x)), [], 1));
        elseif strcmp(errType, 'conjDiff')
            err = @(x) (reshape((FRF(indexRange, :) - numericalFRFs(freqsInRange, nFRF, x)) .* conj(FRF(indexRange, :) - numericalFRFs(freqsInRange, nFRF, x)), [], 1));
        end
        % err = @(x) norm(G(index_range, :) - numericalFRFs(f_range, nFRF, x), 'fro')^2;
        

        options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
        options.Display = "final-detailed";
        xSol(:, modeNo) = lsqnonlin(@(x) err(x), x0, [], [], options);

        FRFopt(:, :, modeNo) = numericalFRFs(freqsInRange, nFRF, xSol(:, modeNo));

        disp('Optimized parameters:');
        disp(['omega_i = ', num2str(real(xSol(1, modeNo)) / (2 * pi)), ' Hz']);
        disp(['damping factor = ', num2str(real(xSol(2, modeNo)) * 100), '%']);
    end

    if (plotBool)
        figure;
        subplot(2, 1, 1);
        semilogy(f, abs(FRF), "LineWidth", 1.5);
        hold on;
        for i = 1:nModes
            semilogy(freqsOpt(:, i), abs(FRFopt(:, :, i)), '--', 'LineWidth', 1.5, 'Color', 'r');
        end
        title('Optimized FRFs');
        xlabel('Frequency (Hz)');
        ylabel('Magnitude');
        
        subplot(2, 1, 2);
        plot(f, angle(FRF), "LineWidth", 1.5);
        hold on;
        for i = 1:nModes
            plot(freqsOpt(:, i), angle(FRFopt(:, :, i)), '--', 'LineWidth', 1.5, 'Color', 'r');
        end
        xlabel('Frequency (Hz)');
        ylabel('Phase');
    end
end