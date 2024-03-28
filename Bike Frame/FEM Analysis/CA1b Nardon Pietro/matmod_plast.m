function [stress,dstress_dstrain,state] = matmod_plast(strain,state_old,mpar)
% INPUT: strain
%        state_old(2), old state variables (plastic strain, hardening)
%        mpar, material parameters mpar.Emod, mpar:Sy, mpar.Hmod
% OUTPUT: stress
%        dstress_dstrain
%        state(2), updated state variables

plastic_strain_old=state_old(1); kappa_old=state_old(2);
stress_trial=mpar.Emod*(strain-plastic_strain_old);
phi_trial=abs(stress_trial)-(mpar.Sy+kappa_old);
if phi_trial<0 %elastic response 
    stress=stress_trial; dstress_dstrain=mpar.Emod;
    state=state_old;
else %plastic repsonse
    dlambda=phi_trial/(mpar.Emod+mpar.Hmod);
    stress=stress_trial-mpar.Emod*dlambda*sign(stress_trial);
    dstress_dstrain=mpar.Emod*mpar.Hmod/(mpar.Emod+mpar.Hmod);
    state(1)=state_old(1)+dlambda*sign(stress_trial);
    state(2)=state_old(2)+mpar.Hmod*dlambda;
end

