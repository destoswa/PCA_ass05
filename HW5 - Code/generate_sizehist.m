function [Dp, dDp, n] = generate_sizehist
  %% no input arguments

  %% generate points evenly in log10 space
  Dp_min = .003;	% um
  Dp_max = 1;		% um
  npts = 200;
  %%
  logDp = linspace(log10(Dp_min), log10(Dp_max), npts)';
  dlogDp = diff(logDp([1,end]))/(npts-1);
  Dp = 10.^logDp;
  dDp = Dp.*(10^dlogDp-1);

  %% logarithmic size distribution (Seinfeld and Pandis, 2006, Ch. 8 and 13 for parameters)
  Dp_bar = 0.2;		% um
  sigma_g = 1.5;	% um
  %%
  dn = 1./(sqrt(2*pi).*Dp.*log(sigma_g)) .* ...
       exp(-log(Dp./Dp_bar).^2./(2.*log(sigma_g).^2)); % cm^-3 cm^-1
  n = dn .* dDp;	% cm^-3 (integral)
  n = n./sum(n);	% normalize

end