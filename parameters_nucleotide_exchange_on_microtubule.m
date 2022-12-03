N_MT = 1;

%length distribution of GMPCPP seeds
%We sample the seed lengths from a exponential distribution
mu = 1.0; %micron, mean for exponential distribution

num_sampling = 1E+5;
R_exp = exprnd(mu,1,num_sampling);
seed_length_ll = 8; %micron 
seed_length_ul = 10; %micron
Rbounded = R_exp((R_exp>=seed_length_ll) & (R_exp<=seed_length_ul));
num_sampling = N_MT;
half_seed_length = 0.5*Rbounded(randperm(numel(Rbounded),num_sampling));
half_seed_length = half_seed_length';

lattice_unit = 0.008; %micron
Npf = 14; %Number of protofilaments

delta_T = 0.01; %time step in sec
T_tot = 33000; % Total simulation time (T_tot*delta_T/60.0 min)
frame_interval = 1500; %frame_interval*delta_T sec

k_base_species1_plus = 0.017; %/sec
k_base_species1_minus = 0.014; %/sec
k_base_species1_lattice = 0.0*k_base_species1_minus; %/sec

k_exchange_species1_plus = 0.00678; %/sec/nM
k_exchange_species1_minus = 0.00167; %/sec/nM
k_exchange_species1_lattice = 0*k_exchange_species1_minus; %/sec/nM
    
k_base_species2_plus = 0.0;
k_base_species2_minus = 0.0;

k_exchange_species2_plus = 0.0;
k_exchange_species2_minus = 0.0;

clasp_conc = 61.5; %nM

N_thres_plus = 6;
N_thres_minus = 6;

pixel_resolution = 0.16; %micron

flag_lattice_exchange = 0; %if it is 1, exchange at the lattice is 'on'



