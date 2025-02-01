function [strfitStats, dChbvInterpVector] = checkFitChbvPoly(ui32PolyDeg, ...
    dInterpDomain, ...
    dChbvCoeffs, ...
    dDataMatrix, ...
    dDomainLB, ...
    dDomainUB, ...
    bIS_ATT_QUAT, ...
    dSwitchIntervals)
arguments
    ui32PolyDeg         (1, 1) uint32
    dInterpDomain       (:, 1) double
    dChbvCoeffs         (:, 1) double
    dDataMatrix         (:, :) double
    dDomainLB           (1, 1) double
    dDomainUB           (1, 1) double
    bIS_ATT_QUAT        (1, 1) logical
    dSwitchIntervals    (:, :) double = []
end
%% PROTOTYPE
% [o_strfitStats] = checkFitChbvPoly(i_ui8PolyDeg, i_dInterpDomain, ...
%     i_dDataMatrix, i_dDomainLB, i_dDomainUB)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_ui8PolyDeg    (1, 1) uint8
% i_dInterpDomain (:, 1) double
% i_dChbvCoeffs   (:, 1) double
% i_dDataMatrix   (:, :) double
% i_dDomainLB     (1, 1) double
% i_dDomainUB     (1, 1) double
% i_bIS_ATT_QUAT  (1, 1) logical
% i_dswitchIntervals (:, :) double = []
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_strfitStats
% o_dChbvInterpVector
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 08-05-2024        Pietro Califano         Function adapted from testing script.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code
if bIS_ATT_QUAT == true
    i_ui8OutputSize = 4; % HARDCODED for specialization
    assert( size(dDataMatrix, 1) == i_ui8OutputSize );
else
    i_ui8OutputSize = size(dDataMatrix, 1);
    dSwitchIntervals = [];
end

assert( size(dDataMatrix, 2) == length(dInterpDomain) );

% Evaluation at test points
Npoints = 5000;
Npoints = Npoints - 2;

testpointsIDs = sort( randi( length(dInterpDomain), Npoints, 1 ), 'ascend' );

TestPoints_Time = [dInterpDomain(1); dInterpDomain(testpointsIDs); dInterpDomain(end)];
TestPoints_Labels = [dDataMatrix(:, 1),...
    dDataMatrix(:, testpointsIDs), ...
    dDataMatrix(:, end)];

dChbvInterpVector = zeros(i_ui8OutputSize, length(TestPoints_Time));

evalRunTime = zeros(length(TestPoints_Time), 1);

for idP = 1:length(TestPoints_Time)
    i_dEvalPoint = TestPoints_Time(idP);

    if bIS_ATT_QUAT == true 
        tic
        dChbvInterpVector(:, idP) = evalAttQuatChbvPolyWithCoeffs(ui32PolyDeg, i_ui8OutputSize, ...
            i_dEvalPoint, dChbvCoeffs, dSwitchIntervals, dDomainLB, dDomainUB);

    else
        tic
        dChbvInterpVector(:, idP) = evalChbvPolyWithCoeffs(ui32PolyDeg, i_ui8OutputSize, ...
            i_dEvalPoint, dChbvCoeffs, dDomainLB, dDomainUB);
    end

    evalRunTime(idP) = toc;
end
fprintf("\nAverage interpolant evaluation time: %4.4g [s]\n", mean(evalRunTime))

% Error evaluation
strfitStats = struct();


if bIS_ATT_QUAT
    strfitStats.dAbsErrVec = abs(abs(dot(dChbvInterpVector, TestPoints_Labels, 1)) - 1);

    strfitStats.dMaxAbsErr = max(strfitStats.dAbsErrVec, [], 'all');
    strfitStats.dAvgAbsErr = mean(strfitStats.dAbsErrVec, 2);

    fprintf('Max absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strfitStats.dMaxAbsErr);
    fprintf('Average absolute difference of (q1-dot-q2 - 1): %4.4g [-]\n', strfitStats.dAvgAbsErr);

else
    strfitStats.absErrVec = abs(dChbvInterpVector - TestPoints_Labels);
    strfitStats.relErrVec = strfitStats.absErrVec./vecnorm(TestPoints_Labels, 2, 1);

    strfitStats.maxAbsErr = max(strfitStats.absErrVec, [], 'all');
    strfitStats.avgAbsErr = mean(strfitStats.absErrVec, 2);

    strfitStats.maxRelErr = 100*max(strfitStats.relErrVec, [], 'all');
    strfitStats.avgRelErr = 100*mean(strfitStats.relErrVec, 2);

    % Printing
    fprintf('Max absolute error: %4.4g [-]\n', strfitStats.maxAbsErr);
    fprintf('Average absolute error: %4.4g, %4.4g, %4.4g [-]\n', strfitStats.avgAbsErr(1), ...
        strfitStats.avgAbsErr(2), strfitStats.avgAbsErr(3));

    fprintf('\nMax relative error: %4.4g [%%]\n', strfitStats.maxRelErr);
    fprintf('Average relative error: %4.4g, %4.4g, %4.4g [%%]\n', strfitStats.avgRelErr(1), ...
        strfitStats.avgRelErr(2), strfitStats.avgRelErr(3));
end


end
