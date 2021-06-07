%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines a pointer class, which enables to pass the Components structure 
% by reference to other functions.
%
% USAGE:
%{
    compPtr = ComponentsPtr(Components);
%}
%
% Authors:
% Ido Marcus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ComponentsPtr < handle
   properties
      comp; % short for 'components'. This field is a struct returned by initializeComponents.m
   end
 
   methods
      function compPtr=ComponentsPtr(Components) % (constructor)
         compPtr.comp=Components;
      end
   end
end