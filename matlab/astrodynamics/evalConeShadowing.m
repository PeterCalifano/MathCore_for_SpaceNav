function o_ui8ShadowType = evalConeShadowing(i_dSCpos_IN, ...
    i_dSunPos_IN, ...
    i_dRbody, ...
    i_dAlphaPenumbra, ...
    i_dAlphaUmbra) %#codegen
%% PROTOTYPE
% o_ui8ShadowType = evalConeShadowing(i_dSCpos_IN, i_dSunPos_IN, i_dRbody, i_dAlphaPenumbra, i_dAlphaUmbra)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function checking shadowing condition for a spherical planet/body orbiting the Sun. The reference frame is
% assumed centred in the body CoM; the input values (i_dAlphaPenumbra, i_dAlphaUmbra) depends on the
% distance from the Sun, but are relatively constant for circular orbits. See [1].
% REFERENCE
% [1] Fundamentals of Astrodynamics and Applications - D. Vallado, 4th Edition, page 301, Algorithm 34.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% i_dSCpos_IN:       [3, 1]   Spacecraft position in IN frame centred in the attractor
% i_dSunPos_IN:      [3, 1]   Position vector to the Sun in IN frame
% i_dRbody:          [1]      Radius of the spherical attractor (CoM := IN origin)
% i_dAlphaPenumbra:  [1]      Complementary half-angle of Penumbra shadow cone. See [1] for details
% i_dAlphaUmbra:     [1]      Complementary half-angle of Umbra shadow cone. See [1] for details
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% o_ui8ShadowType:   [1]      Integer indicating shadow type: 0: no shadow; 1: Penumbra; 2: Umbra
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 14-02-2024        Pietro Califano         Function coded. Not verified.
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades
% [-]
% -------------------------------------------------------------------------------------------------------------
%% Function code

% Default value assignment
o_ui8ShadowType = uint8(0); % No shadowing case

% Compute auxiliary variables
dSCposDotSunPos = dot(i_dSCpos_IN, i_dSunPos_IN);
dSCposNorm = norm(i_dSCpos_IN);

if dSCposDotSunPos < 0 % If true --> there could be shadowing
    % Check for penumbra case

    % Compute angle between directions
    dCosTheta = dSCposDotSunPos / ( dSCposNorm * norm(i_dSunPos_IN);
    dthetaAngle = acos( dCosTheta );

    % Decompose SC position vector
    dhorizSatDist = dSCposNorm * dCosTheta;
    dVertSatDist = dSCposNorm * sin(dthetaAngle);

    % Compute distance from Planet centre to intersection point
    dxPoint = i_dRbody/( sin(i_dAlphaPenumbra) );
    dPenumbraVert = tan(i_dAlphaPenumbra) * (dxPoint + dhorizSatDist);

    if dVertSatDist <= dPenumbraVert
        o_ui8ShadowType = uint8(1); % Penumbra shadow

        % Check for umbra case
        dyPoint = i_dRbody/( sin(i_dAlphaUmbra) );
        dUmbraVert = tan(i_dAlphaUmbra) * (dyPoint - dhorizSatDist);

        if dVertSatDist <= dUmbraVert
            o_ui8ShadowType = uint8(2); % Umbra shadow
        end
    end
end

end



