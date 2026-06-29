function A = n2A(nin, nout)
% A = n2A(nin, nout)
%   Inputs:
%       nin  - Index of refraction inside medium
%       nout - Index of refraction outside medium
%   Output:
%       A    - Index of refraction mismatch parameter

    if nargin==0
        dan12=1.4;
    else
        dan12=nin/nout;
    end

    if dan12>1
        A=504.332889-2641.00214*dan12+...
            5923.699064*dan12.^2-7376.355814*dan12^3+...
            5507.53041*dan12^4-2463.357945*dan12^5+...
            610.956547*dan12^6-64.8047*dan12^7;
    elseif dan12<1
        A=3.084635-6.531194*dan12+...
            8.357854*dan12^2-5.082751*dan12^3+1.171382*dan12^4;
    else
        A=1;
    end

end