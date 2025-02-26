function [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
                                                            bEliminateInPlace) %#codegen
arguments
    dTargetMatrix      
    bEliminateInPlace (1,1) logical {islogical, isscalar} = true;
end
%% SIGNATURE
% [dTargetMatrix, dOrthogonalQ] = GivensEliminateQR(dTargetMatrix, ...
%                                                   bEliminateInPlace) %#codegen
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function performing QR decomposition using Givens Rotations.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% dTargetMatrix         (:,:) double
% bEliminateInPlace     (1,1) logical
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% dTargetMatrix % Target matrix transformed in place
% dOrthogonalQ TODO
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 26-02-2025    Pietro Califano     First version implemented and validated.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------

if coder.target('MATLAB') || coder.target("MEX")
    % Asserts

end

% TODO add capability of computing Q = G1 * ... * Gt
dOrthogonalQ = 0; 

% NOTE: these branches are intended as separate functions, thus self-contained. The coder should generate
% two versions according to the flag.
if bEliminateInPlace
    ui16NumOfCols = uint16(size(dTargetMatrix, 2));
    
    for ui16IdElimCol = uint16(1:ui16NumOfCols)

        for ui16IdElimRow = uint16(size(dTargetMatrix, 1):-1:ui16IdElimCol+1)

            % Run GivensEliminateRow
            ui16TargetSubscript = [ui16IdElimRow, ui16IdElimCol];
            ui16AuxRowId        = ui16TargetSubscript(1) - 1;
            
            % Use "inlined" function instead of call (reduces copy of arrays)
            % dTargetMatrix = GivensEliminateRow(dTargetMatrix, ...
            %                                    ui16TargetIds, ...
            %                                    ui16AuxRowId);

            % Get (a,b) entries
            dVal1 = dTargetMatrix(ui16AuxRowId,           ui16TargetSubscript(2));
            dVal2 = dTargetMatrix(ui16TargetSubscript(1), ui16TargetSubscript(2));

            % Compute Givens rotation values
            [dCos, dSin] = GivensRotVals([dVal1; dVal2]);

            % Transform rows in place
            for ui16IdCol = 1:ui16NumOfCols

                dTmp1 = dTargetMatrix(ui16AuxRowId, ui16IdCol);
                dTmp2 = dTargetMatrix(ui16TargetSubscript(1), ui16IdCol);

                % Rotate dTargetMatrix, equivalent to Arot = G^T * A
                dTargetMatrix(ui16AuxRowId, ui16IdCol) = dCos * dTmp1 - dSin * dTmp2;
                dTargetMatrix(ui16TargetSubscript(1), ui16IdCol) = dSin * dTmp1 + dCos * dTmp2;

            end

        end % Slide across rows
    end % Slide across columns
    
    % Remove numerical zeros
    dTargetMatrix(abs(dTargetMatrix) < 1.5*eps ) = 0.0;

    % TODO add capability of computing Q incrementally
    return
else
    error('Not implemented yet. Use in-place rotations')
    % TODO TBD, not sure this makes any sense in MATLAB tbh
    return
end




end

