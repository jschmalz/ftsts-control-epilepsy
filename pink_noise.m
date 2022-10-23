function I = pink_noise(sigma,tau_x,tau_t,dt,Iz,n)

[x2,x1] = meshgrid(1:n(2),1:n(1));
I = 0;


% Random part - Ornstein–Uhlenbeck process
% Space
z = randn(n);
if ~isempty(tau_x)
    z = SpaceFilter(z,tau_x);
end
    function z = SpaceFilter(z,tau_x)
        dx = 1;
        for i = 1:numel(tau_x)
            if tau_x(i)>0
                z = permute(z,[i sort(setdiff(1:numel(size(z)),i))]);
                nz = size(z);
                z = reshape(z,size(z,1),[]);
                z = filter(sqrt(2*dx/tau_x(i)),[1, (dx-tau_x(i))/tau_x(i)],z, ...
                    (1-sqrt(2*dx/tau_x(i)))*z(1,:));
                % Initial condition needs to be set to avoid 'warm up' phenomenon
                z = reshape(z, nz);
                z = ipermute(z,[i sort(setdiff(1:numel(size(z)),i))]);
            end
        end
    end

% Time - Ornstein–Uhlenbeck process
if ~isempty(tau_t) && tau_t > 0
    
    if all(Iz(:)) == 0
        z0 = randn(n);
        if ~isempty(tau_x)
            z0 = SpaceFilter(z0,tau_x);
        end
        Iz = z0;
    end
    Iz = (Iz./sigma) .* exp(-dt ./ tau_t) + ... % The autoregression part
        sqrt(2*dt/tau_t) .* z; % Innovation part
else
    Iz = z;
end

% Calculate result
I = I + sigma*Iz;
end
