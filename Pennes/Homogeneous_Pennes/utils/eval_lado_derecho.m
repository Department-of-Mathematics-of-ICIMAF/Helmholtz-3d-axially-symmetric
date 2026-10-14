function d_alpha_dt=eval_lado_derecho(t, alpha0, tstop, K_0, K_g, alphag, q)
%     fprintf('t= %12.8f; ', t);
    epsilon = 0.1;
    tini=tstop + epsilon;
    tfin=tstop + 2*epsilon;

    %% contruir suavizante
    if t<=tini
%         fprintf('I(t) = 1\n')
        qon=q;
    elseif t>tini && t< tfin
        % tes- tescalon, valor de t donde hay que evaluar la funcion escalon
        delta=tfin-tini;
        tes=-0.0025+(0.005/delta)*(t-tini);
        sp=stepf;
        tstep=fnval(sp,tes);
%         fprintf('I(t) = %e \n', tstep)
        qon=tstep*q;
    else
%         fprintf('I(t) = 0 \n')
        n0=size(K_0,1);
        qon=zeros(n0,1);  
    end


    z=-K_g*alphag + qon;
    d_alpha_dt = -K_0*alpha0 + z;