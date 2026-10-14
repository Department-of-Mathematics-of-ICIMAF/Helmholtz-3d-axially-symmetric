function d_alpha_dt=eval_lado_derecho(t, alpha0, tstop, K_0, K_g, alphag, q)

    tini=tstop-0.6;
    tfin=tstop+0.3;

    %% contruir suabizante
    if t<=tini
        qon=q;
    elseif t>tini && t< tfin 
        % tes- tescalon, valor de t donde hay que evaluar la funcion escalon
        delta=tfin-tini;
        tes=-0.0025+(0.005/delta)*(t-tini);
        sp=stepf;
        tstep=fnval(sp,tes);
        qon=tstep*q;
    else
        n0=size(K_0,1);
        qon=zeros(n0,1);  
    end


    z=-K_g*alphag + qon;
    d_alpha_dt = -K_0*alpha0 + z;