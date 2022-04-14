function state_dot = orbdyn(t, state, mu)
    
    x = state(1);
    y = state(2);
    z = state(3);
	r = [state(1),state(2), state(3)];
	
    %% Sun's gravity

	accx = -(mu/((norm(r))^3))*x;
	accy = -(mu/((norm(r))^3))*y;
	accz = -(mu/((norm(r))^3))*z;

    state_dot(1) = state(4);
    state_dot(2) = state(5);
    state_dot(3) = state(6);
    state_dot(4) = accx; 
	state_dot(5) = accy; 
	state_dot(6) = accz; 
	state_dot= state_dot';
end