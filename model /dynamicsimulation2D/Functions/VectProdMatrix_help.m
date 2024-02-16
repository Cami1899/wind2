
%function [a_cross] = VectProdMatrix(a)

% OROLAB 2023
% (c) Dr. Marcello Romano on 2019-04-11
% Ref: 3815 Unit 1 Slides, Marcello Romano, 2019
%
% INPUT
% a = [3,1 matrix of real numbers] = components of a Euclidean vector a_vect in a
% CCS K 
%
% OUTPUT
% a_cross = [3,3 antisyimmetric matrix of real numbers] matrix which right
% multiplied by a 3,1 matrix b (of components in K of Euclidean vector b_vect) 
% gives a 3,1 matrix c ofcomponents in K of Euclidean vector
% c_vect =a_vect X b_vector
%
% In other words a_cross b = c is the matrix-algebra equivalent of the vector product. 