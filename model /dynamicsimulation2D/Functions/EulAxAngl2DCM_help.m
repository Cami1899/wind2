%function [DCM_BA] = EulAxAngl2DCM(EulerAxis, EulerAngle)

% OROLAB, 2023
% This function implements the Euler formula for a DCM
% (c) Dr. Marcello Romano: written on 2019-04-11
%
% INPUTS
% EulerAxis = [3,1 matrix of real numbers] = COMPONENTS of the
% Euler axis unit vector of relative orientation of CCS B w.r.t. CCS A ALONG 
% either CCS A or CCS B (they are the same!)
% EulerAngle = [1 real number] = Euler angle in radians of the relative 
% orientation of CCS B w.r.t. CCS A  
% 
%
% OUTPUTS
% DCM_BA = [3,3 matrix in SO(3)]= Direction Cosine Matrix from CCS A to
% CCS B = Rotation Matrix from CCS B to CCS A
