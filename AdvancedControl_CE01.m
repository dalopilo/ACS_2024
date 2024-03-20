%% 1.1 Norms of SISO systems

%% 1.1.1
%question 1)
G = tf([1, -1], [1, 2, 10]);
delta_w = 0.01;
w = 0.01:delta_w:1000; %frequency steps    last value is 1/100 so it is small enough, i am ignoring the rest :)
[magnitude_G, phase] = bode(G,w);
bode(G,w) %bode plots

%question 2) approximation Gtwo_norm ~= 1/pi * sum(mag_G)^2 3 delta_w)
approx_2norm = sqrt(1/pi * sum(magnitude_G.^2)*delta_w);

% question 3)
t = 0:0.01:10;
g = impulse(G,t);
impulse_2norm = sqrt(trapz(t,g.^2));

%question 4)
[A,B,C,D] = ssdata(G);
L = are(A',zeros(2,2),B*B');

ARE_norm2_G = sqrt(trace(C*L*C'));

%check question 5)
norm(G,2);

%% 1.1.2
%question 1) freq response inf_norm = sup_w(|G(jw)|)

freq_resp_inf_norm = max(magnitude_G);

%question 2) Bounded Real Lemma
lower_b = 0;
upper_b = 100;
tol = 10^-4;

while 1
    if (((upper_b - lower_b)/lower_b)<tol)
        break;
    end

    gamma = (upper_b + lower_b)/2;

    H = [A, B*B'/(gamma^2) ; -C'*C, -A'];
    eig_H = eig(H);

    if any(real(eig_H)== 0) 
        lower_b = gamma;
    else 
        upper_b = gamma;
    end 
end

inf_norm_BRL = (upper_b + lower_b)/2;

%question 3) check with norm fct
check = norm(G,inf);

%% 1.2 Norms of MIMO systems

%% 1.2.1
A = [20, -27, 7;
    53, -63, 13;
    -5, 12, -8];

B = [1, -1;
    -2, -1;
    -3, 0];

C = [0, 0, -2;
    1, -1, -1];

D = [0, 0; 
    0, 0];

%question 1) 
[num_coeff,denom_coeff] = ss2tf(A,B,C,D,1);
G = tf(num_coeff(1,:), denom_coeff);
delta_w = 0.01;
w = 0.01:delta_w:1000;
[magnitude_G, ~] = bode(G,w);
approx_2norm = sqrt(1/pi * sum(magnitude_G.^2)*delta_w);
%slide 15

%question 2)
L = are(A',zeros(3,3),B*B');
ARE_norm2 = sqrt(trace(C*L*C'));

%question 3)
norm2_G = norm(G,2);

%% 1.2.2 
%question 1)
g = ss(A, B, C, D);
g11 = g(1,1);
g12 = g(1,2);
g21 = g(2,1);
g22 = g(2,2);

freq_resp_11 = freqresp(g11, w);
freq_resp_12 = freqresp(g12, w);
freq_resp_21 = freqresp(g21, w);
freq_resp_22 = freqresp(g22, w);

g_freqresp_inf_norm = 0;
for i=1:length(w)
    g_mat = [freq_resp_11(:, :, i), freq_resp_12(:, :, i);
             freq_resp_21(:, :, i), freq_resp_22(:, :, i)];
    g_svd = svd(g_mat);
    if g_svd(1) > g_freqresp_inf_norm
        g_freqresp_inf_norm = g_svd(1);
    end
end

% For recall : svd(A) = eig(AA'), so we could use this relationship also

%%
%question 2) Bounded Real Lemma
lower_b = 0;
upper_b = 100;
tol = 10^-4;

while 1
    if (((upper_b - lower_b)/lower_b)<tol)
        break;
    end

    gamma = (upper_b + lower_b)/2;

    H = [A, B*B'/(gamma^2);
        -C'*C, -A'];
    eig_H = eig(H);

    if any(real(eig_H) >= 0) 
        lower_b = gamma;
    else 
        upper_b = gamma;
    end 
end

inf_norm_BRL = (upper_b + lower_b)/2;

%question 3)
check = norm(G,inf);






