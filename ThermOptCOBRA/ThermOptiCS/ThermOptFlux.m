function [flux] = ThermOptFlux(model,flux,TICs,Dir)
% Gets the TIC free flux values for the given model and a feasible flux
%
% USAGE: 
%   [flux] = ThermOptFlux(model,flux,TICs,Dir)
%
% INPUTS:
%     model:     COBRA model structure for from which the flux values are
%                obtained
%     flux:      Flux distribution obtained from any flux analysis methods 
%                (FBA or flux sampling) 
%     TICs:      List of all the Thermodynamically infeasible cycles in
%                the given input model
%     Dir:       The flux directions for reactions in the corresponding
%                TICs
% OUTPUTS:
%     flux:      TIC free flux
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


for i =1:numel(TICs)
    t = TICs{i}; c = Dir{i};
    ids = find(ismember(model.rxns,t));
    f=flux(ids);
    if sum(f>0&c>0)+sum(f<0&c<0)==numel(t)
        [~,mID] = min(abs(f));
        c = (c/abs(c(mID)))*abs(f(mID));
        flux(ids) = flux(ids)-c;
    end
end

end