function idt = IDTKriging(T,P,phi,comp,model)
% compute the IDT of a 6-component mixture with order of components
% [isocetane, n-decane, n-dodecane, n-heptane, trimethylbenzene, xylene]
% using Kriging Models
% Inputs can have multiple rows - to compute many IDTs
% T in [K], P in [Pa], phi [-], comp [-]
% 'model' is a struct containing the whole Kriging Model
  npoints = size(T,1);
  idt = zeros(npoints,1);
  for i = 1:npoints
    if P(i) < min(model.P) || P(i) > max(model.P)
        error('The pressure you want is out of bounds: %f',P(i));
    end
    idxMLow = 1;
    idxMUp  = 2;
    while P(i) >= model.P(idxMUp) && idxMUp < length(model.P)
            idxMLow = idxMLow + 1;
            idxMUp  = idxMUp + 1;
    end
    xin = [1000/T(i), phi(i), comp(i,:)];
    % Estimate Low-bound
    idtLow = krigingsingleP(xin,...
        reshape( model.X(idxMLow,1:model.Npt(idxMLow),:), [model.Npt(idxMLow),size(model.X,3)] ), ...
        reshape( model.Y(idxMLow,1:model.Npt(idxMLow)), [model.Npt(idxMLow),1]),...
        model.thetas(idxMLow,:),  model.pkriging(idxMLow),...
        reshape( model.U(idxMLow,1:model.Npt(idxMLow),1:model.Npt(idxMLow)),[model.Npt(idxMLow),model.Npt(idxMLow)]));
    % Estimate Up-bound
    idtUp = krigingsingleP(xin,...
        reshape( model.X(idxMUp,1:model.Npt(idxMUp),:), [model.Npt(idxMUp),size(model.X,3)] ), ...
        reshape( model.Y(idxMUp,1:model.Npt(idxMUp)), [model.Npt(idxMUp),1]),...
        model.thetas(idxMUp,:),  model.pkriging(idxMUp),...
        reshape( model.U(idxMUp,1:model.Npt(idxMUp),1:model.Npt(idxMUp)),[model.Npt(idxMUp),model.Npt(idxMUp)]));
    % Linear interpolation of pressure
    idt(i) = idtLow + (P(i)-model.P(idxMLow)) * (idtUp-idtLow)/(model.P(idxMUp)-model.P(idxMLow));
  end
  idt = exp(idt); % convert ln(tau) to tau
end

function f = krigingsingleP(xtoeval,X,y,theta,p,U)
n=size(X,1);
one=ones(n,1);
mu=(one'*(U\(U'\y)))/(one'*(U\(U'\one)));
psi=ones(n,1);
for i=1:n
	psi(i)=exp(-sum(theta.*abs(X(i,:)-xtoeval).^p));
end
f=mu+psi'*(U\(U'\(y-one*mu)));
end