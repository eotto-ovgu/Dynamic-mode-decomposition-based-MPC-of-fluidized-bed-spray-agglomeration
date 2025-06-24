function n_dot_enuc = generateFeed(m_dot_enuc,L_enuc,sig_enuc,xgrid,dx,dv,rho_p)
    
    q0_enuc_x     =       normpdf(xgrid,L_enuc,sig_enuc);
    q0_enuc       =       q0_enuc_x.*dx./dv/(dv*(q0_enuc_x'.*dx'./dv'));
    N_dot_enuc    =       m_dot_enuc/rho_p/pi*6/((xgrid.^3.*dx)*q0_enuc_x');  
    n_dot_enuc    =       N_dot_enuc*q0_enuc;
end