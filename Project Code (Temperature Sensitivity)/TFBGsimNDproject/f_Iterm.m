%{

J. Albert

This script is a function file which calculates the integral "I" needed to
calculate the coupling coefficients
Information comes from Ergogan and Sipe and the mask period used here is
the "effective" mask period along the fiber axis.

@param n: Constant, index of refraction of the core
@param r: array of radial coordinate values (over the core only)
@param MaskPer: Period of the mask projected on the z axis
@param theta_int: tilt angle inside the fiber
@param E1: array of Electric field values, first mode to be coupled
@param E2: array of Electric field values, first mode to be coupled
@param k: order of Bessel function to be used inside integral for Imn,pq,k

@return Iterm: value of integral Eqn(116) Michel
%}

function Iterm = f_Iterm(n, r, MaskPer, theta_int, E1, E2, k)
	
    arg=r*4*pi*tand(theta_int)/MaskPer;
	Integrand = (pi*n*besselj(k,arg)).*E1.*E2.*r;
	Iterm = trapz(r, Integrand);
    

end