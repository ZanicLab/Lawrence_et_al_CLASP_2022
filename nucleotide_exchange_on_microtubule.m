%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%          Modeling CLASP mediated Nucleotide Exchange Dependent          %
%                    Microtubule Depolymerization                         %         
%                                                                         %   		                                          	      								
%	  THIS CODE UTILIZES STOCHASTIC MONTE CARLO ALGORITHM TO GENERATE,    %
%     REPLICATE AND PREDICT EVENTS OF CLASP-MEDIATED NUCLEOTIDE-DEPENDENT %
%     MICROTUBULE DEPOLYMERIZATION AS OBSERVED IN THE IN VITRO TIRF       %
%       RECONSTITUTION EXPERIMENTS. THE MODELING PARAMETERS ARE           % 			      
%         AIMED TO BE GROUNDED IN EXPERIMENTAL QUANTIFICATIONS            %      			                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			                                                              %
%                      DESIGN & DEVELOPMENTS 				              %
%									                                      %
%			     Saptarshi Chatterjee and Marija Zanic                    %       
%	Dept. of Cell and Developmental Biology, Vanderbilt University, USA   %    
%      			            (c) 2022                                      %				                                
%                                                                         %
%                 last updated on November, 2022                          %
%                                                                         %
%This code generates data set for number of nucleotide exchange events on %
%microtubules over time followed by microtubule depolymeriztion.          %
%The rates of exchange at both microtubule ends and lattice are such that %
%the average microtubule depolymerization rates obtained from simulations %
%match with the experimental depolymerization rates at both minus and plus%
%ends                                                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

delete *.txt %This removes previosuly stored .txt files in the directory

clear all  %clear memory
close all  %close all MATLAB popup windows (figure windows mostly)

N_ensemble = 100; %Number of independent simulation runs

%In this code each ensemble has N_MT = 1; So, N_ensemble denotes the number
%of independent microtubules

for ens_counter = 1:N_ensemble %loop for the number of ensembles

%clear memory 
clearvars -except ens_counter N_ensemble %deletes all variables except 
                  %ens_counter & N_ensemble in workspace
                    
close all %close all MATLAB popup windows (figure windows mostly)

rng (ens_counter) % rng => random number generator, ens_counter => rng seed

%path to the directory where parameters_exchange.m file exists
%addpath ('/Users/saptarshi/Desktop/Nucleotide_Exchange')

%parameters values are stored in this file
parameters_nucleotide_exchange_on_microtubule 

%fileIDs for data files storing simulation data output

%This file stores the coordinates of microtubule ends and length at certain
%time intervals

str = sprintf('data_microtubule_length_vs_time_seed_%d_ens_%d.txt',...
        N_MT,ens_counter);
fileID = fopen(str,'a');

%This file stores the depolymerization rates for microtubule plus and 
%minus ends evaluated for each microtubule (in this case each ensemble)

str = sprintf('data_rates_%d_ens_%d.txt',...
        N_MT,ens_counter);
fileID1 = fopen(str,'a');

%The function discretize_segment discretizes the entire microtubule into 
%segments of length lattice_unit. num_microtubule_segment stores the 
%rounded up value of (microtubule length/lattice unit). The array named 
%array_seed converts the discretized microtubules into an array of dimension
%N_MT X num_microtubule_segment X Npf and stores index of each discretized
%segment for Npf number of protofilaments.

[num_microtubule_segment,array_seed] = ...
    discretize_segment(lattice_unit,2.0*half_seed_length,N_MT,Npf);

%This function gives the visual rendition of a microtubules at a resolution
%of lattice_unit. If the lattice_unit is set at 8 nm (0.008 micron), each
%circle in the rendered image denotes a tubulin dimer or the absence of it
%depending on the context.

%render_Visual(array_seed,num_microtubule_segment,N_MT,Npf,t);

%These two arrays store the number of how many microtubule rings have fallen
%off from the microtubule ends due to nucleotide exchange over time
%since the onset of simulation. At the onset of simulation, the arrays are 
%initialized with zeros as no ring has fallen off at t = 0.

empty_ring_plus = zeros(N_MT,1);
empty_ring_minus = zeros(N_MT,1);

%The duration of the simulation run is T_tot*delta_T/60.0 min.

for t=0:T_tot % This loop facilitates the duration of the simulation run

%In the following loop (or just a single iteration as N_MT = 1), the 
%decision that whether a microtubule ring would fall off or not is taken. 
%If the number of exchanged sites in a ring equals or exceeds a 
%predetermined thereshold, we decide to depolymerize that ring entirely.

% 1 2 3 4 5 6 7 8 9 .....num_microtubule_segment(i) lets say (100)
%   2 3 4 5 6 7 8 9 .....99
%     3 4 5 6 7 8 9 ...98
%       4 5 6 7 8 9 ..97 so on and so forth

%At the beginning the minus end index starts from 1 and the plus end from 
%num_microtubule_segment(i). Now as time progresses, rings fall off. As a
%consequence of that, minus ends take index values greater than 1 and plus
%ends less than num_microtubule_segment(i). The numbers jplus and jminus
%give the indices of the plus and minus ends at a particular time instant.

%In the beginning the array_seed elements are ones where microtubule exists. 
%At plus end sites, when the nucleotide exchange occurs, the ones are 
%changed to -1's. Similarly, nucleotide exchange occuring at minus ends
%changes the array_seed element value to -2 and nucleotide exchange 
%occuring at the microtubule lattice changes the array_seed element value 
%to -3. In the simulations, we impose the following rule that if the
%total number of exchanged sites at the plus end ring equals or exceeds
%the value of N_thres_plus, the concerned ring falls off.
%Similarly, at the minus ends, the ring falls off if the total
%number of exchanged sites equals or exceeds the value of N_thres_minus.

 flag_lattice_plus_fall_off = 0;
 flag_lattice_minus_fall_off = 0;

for i=1:N_MT
    
    while 1<2
        
    jplus = num_microtubule_segment(i) - empty_ring_plus(i);
    jminus = 1 + empty_ring_minus(i);
        
    flag_plus_end = sum(array_seed(i,jplus,:) == -1) + ...
                    sum(array_seed(i,jplus,:) == -3);
    flag_minus_end = sum(array_seed(i,jminus,:) == -2) + ...
                    sum(array_seed(i,jminus,:) == -3);
    
    flag_lattice_adjacent_plus = sum(array_seed(i,jplus-1,:) == -3);
    flag_lattice_adjacent_minus = sum(array_seed(i,jminus+1,:) == -3);
    
    flag_lattice_plus_fall_off = 0;
    flag_lattice_minus_fall_off = 0;
    
    if flag_plus_end >= N_thres_plus
       empty_ring_plus(i) = empty_ring_plus(i) + 1;
       array_seed(i,jplus,:) = 0;
       if flag_lattice_adjacent_plus >= N_thres_plus
          empty_ring_plus(i) = empty_ring_plus(i) + 1;
          array_seed(i,jplus-1,:) = 0;
          flag_lattice_plus_fall_off = 1;
       end
    end
       
    if flag_minus_end >= N_thres_minus
       empty_ring_minus(i) = empty_ring_minus(i) + 1;
       array_seed(i,jminus,:) = 0;
       if flag_lattice_adjacent_minus >= N_thres_minus
          empty_ring_minus(i) = empty_ring_minus(i) + 1;
          array_seed(i,jminus+1,:) = 0;
          flag_lattice_minus_fall_off = 1;
       end
    end
    
    if ((flag_lattice_plus_fall_off + flag_lattice_minus_fall_off) == 0)
       break; 
    else 
        fprintf('Lattice ring %d %d had fallen off at %.2f min...\n',...
            jplus,jminus,t*delta_T/60)
    end
    
    end
                   
end
    
%This function facilitates the exchange at both microtubule plus
%and minus ends. The code is generalized for 2 species of nucleotides.
%Each species has 2 distinct exchange rates for plus and minus ends. These
%rates are considered to follow exponential distribution. The
%probability of exchange depends on the rate times clasp concentration.
%Since in experiments, without the clasp, a slow depolymerization is
%observed, we also have considered base rates (clasp independent) for both
%the species at both ends. The probability of exchange depends on the
%following factor: 1 - exp(-(k_base + k_exchange_species*clasp_conc)*delta_T)
% At each iteration, we randomly visit all Npf number of sites at the 
%microtubule ends and check if those sites are un-exchanged (i.e the 
%element holds a value of 1). At an un-exchanged site, the exchange happens 
%with the previously mentioned probability.

[exchanged_array] =...
    microtubule_lattice_exchange(array_seed,num_microtubule_segment,...
    N_MT,Npf,k_base_species1_plus,k_base_species1_minus,k_base_species1_lattice,...
    k_exchange_species1_plus,k_exchange_species1_minus,k_exchange_species1_lattice,...
    k_base_species2_plus,k_base_species2_minus,k_exchange_species2_plus,...
    k_exchange_species2_minus,clasp_conc,delta_T,...
    empty_ring_plus,empty_ring_minus,flag_lattice_exchange);

%records the new state of the array after exchange

array_seed = exchanged_array; 

%This function calculates and provides the microtubule_length and the 
%coordinates of microtubules plus and minus ends at each iteration as the output.
   
[microtubule_length,microtubule_end_plus,microtubule_end_minus]= ...
Measure_microtubule_length (N_MT,lattice_unit,num_microtubule_segment,...
empty_ring_plus,empty_ring_minus);


%If the microtubule length falls below the pixel resolution, the simulation 
%for the ensemble under consideration stops. It moves onto the next ensemble. 
   
if microtubule_length < pixel_resolution
   flag_pixel_resolution = 1;
end

if (mod(t,frame_interval) == 0 || flag_pixel_resolution == 1)
    
   flag_pixel_resolution = 0; 
   
   %t*delta_T/60;
   
   %render_Visual(array_seed,num_microtubule_segment,N_MT,Npf,t,frame_interval);
   
   %fprintf(fileID,'%f %f\n',(t*delta_T/60.0),microtubule_length);
   
   fprintf(fileID,'%f %f\n',microtubule_end_plus,t*delta_T/60.0);
   fprintf(fileID,'%f %f\n',microtubule_end_minus,t*delta_T/60.0);
   fprintf(fileID,'\n\n');
   
   if microtubule_length < pixel_resolution
      fprintf('The microtubule %d entirely depolymerized at %.2f min...\n',ens_counter,t*delta_T/60)
      break; 
   end
   
      
end

end 

%The function Measure_Depolymerization_Rate calculates the depolymerization 
%rates at both microtubule ends. It calculates the rate by evaluating the 
%difference between end positions at t_initial and t_final and dividing it
%by the time duration (t_final - t_initial).

t_initial = 0; %sec
t_final = t*delta_T; %sec

[rate_depoly_plus,rate_depoly_minus] = Measure_Depolymerization_Rate (...
    N_MT,lattice_unit,num_microtubule_segment,empty_ring_plus,...
    empty_ring_minus,t_initial,t_final);

%The multiplication by a factor of 1000 changes the unit to nm/s from um/s

fprintf(fileID1,'%f %f\n',1000*rate_depoly_plus,1000*rate_depoly_minus); 

%close the output data file when the simulation run for each individula 
%ensemble ends. 

fclose(fileID); 
fclose(fileID1); 

fprintf('%dth ensemble ends...\n',ens_counter);

end

function [num_microtubule_segment,array_seed] = discretize_segment(lattice_unit,segment_length,N_MT,Npf)

array_seed = zeros(N_MT,100,Npf); %Initialize array_seed with zeros
num_microtubule_segment = zeros(N_MT,1); %Initialize array storing discretized segment info with zeros

for i=1:N_MT
    num_microtubule_segment(i) = round(segment_length(i)/lattice_unit);
    for j = 1:num_microtubule_segment(i)
        array_seed(i,j,:) = 1;
    end
end 
                  
end

%This function named microtubule_lattice_exchange facilitates the lattice 
%exchange at both microtubule plus and minus ends. The code is generalized
%for 2 species of nucleotides. Each species has 2 distinct exchange rates 
%for plus and minus ends. These rates are considered to follow exponential 
%distribution. The probability of exchange depends on the rate times clasp 
%concentration. Since in experiments, without the clasp, a slow 
%depolymerization is observed, we also have considered base rates 
%(clasp independent) for both the species at both ends. The probability of 
%exchange depends on the following factor: 
% (1 - exp(-(k_base + k_exchange_species*clasp_conc)*delta_T))
% At each iteration, we randomly visit Npf number of sites at the microtubule
%ends and check if those sites are un-exchanged (i.e the element holds a
%value of 1). At an un-exchanged site, the exchange happens with the
%previously mentioned probability.

%At plus end sites, when the nucleotide exchange occurs, the ones are 
%changed to -1's. Similarly, nucleotide exchange occuring at minus ends
%change the array_seed element value to -2. 

function [exchanged_array] =...
    microtubule_lattice_exchange(array_seed,num_microtubule_segment,...
    N_MT,Npf,k_base_species1_plus,k_base_species1_minus,k_base_species1_lattice,...
    k_exchange_species1_plus,k_exchange_species1_minus,k_exchange_species1_lattice,...
    k_base_species2_plus,k_base_species2_minus,k_exchange_species2_plus,...
    k_exchange_species2_minus,clasp_conc,delta_T,...
    empty_ring_plus,empty_ring_minus,flag_lattice_exchange)

exchanged_array = array_seed;

p_exchange_species1_plus = 1.0 - ...
 exp(-(k_base_species1_plus+k_exchange_species1_plus*clasp_conc)*delta_T);
p_exchange_species1_minus = 1.0 - ...
 exp(-(k_base_species1_minus+k_exchange_species1_minus*clasp_conc)*delta_T);

p_exchange_species1_lattice = 1.0 - ...
 exp(-(k_base_species1_lattice+k_exchange_species1_lattice*clasp_conc)*delta_T);


p_exchange_species2_plus = 1.0 - ...
 exp(-(k_base_species2_plus+k_exchange_species2_plus*clasp_conc)*delta_T);
p_exchange_species2_minus = 1.0 - ...
 exp(-(k_base_species2_minus+k_exchange_species2_minus*clasp_conc)*delta_T);


for i=1:N_MT
    
    jplus = num_microtubule_segment(i) - empty_ring_plus(i);
    jminus = 1 + empty_ring_minus(i);
    
    kpf = 1;       
    pf_list = randperm(Npf);
    while kpf<=Npf
          %k = randi(Npf);
          k = pf_list(kpf);
          
          if (exchanged_array(i,jplus,k) == 1)
             p = rand;
             if p <= p_exchange_species1_plus
                exchanged_array(i,jplus,k) = -1;
             end
          end
              
          if (exchanged_array(i,jminus,k) == 1)
             p = rand;
             if p <= p_exchange_species1_minus
                exchanged_array(i,jminus,k) = -2;
             end
          end
          
          if (exchanged_array(i,jplus,k) == 1)
             p = rand;
             if p <= p_exchange_species2_plus
                exchanged_array(i,jplus,k) = -1;
             end
          end
              
          if (exchanged_array(i,jminus,k) == 1)
             p = rand;
             if p <= p_exchange_species2_minus
                exchanged_array(i,jminus,k) = -2;
             end
          end
              
          kpf = kpf + 1;
    end
    
  if flag_lattice_exchange == 1
    
    for jlattice = jminus+1:jplus-1
        
        kpf_lattice = 1;
        pf_list = randperm(Npf);
        
        while kpf_lattice<=Npf
          %k = randi(Npf);
          k = pf_list(kpf_lattice);
          
          if (exchanged_array(i,jlattice,k) == 1)
             p = rand;
             if p <= p_exchange_species1_lattice
                exchanged_array(i,jlattice,k) = -3;
             end
          end
              
          kpf_lattice = kpf_lattice + 1;
          
        end
        
    end
    
  end
    
end 

end

%The function Measure_Depolymerization_Rate calculates the depolymerization 
%rates at both microtubule ends. It calculates the rate by evaluating the 
%difference between end positions at t_initial and t_final and dividing it
%by the time duration (t_final - t_initial).

function [rate_depoly_plus,rate_depoly_minus] = ...
    Measure_Depolymerization_Rate (N_MT,lattice_unit,num_microtubule_segment,...
    empty_ring_plus,empty_ring_minus,t_initial,t_final)
    
    rate_depoly_plus = zeros(N_MT,1);
    rate_depoly_minus = zeros(N_MT,1);

    for i = 1:N_MT
        j_plus_initial = num_microtubule_segment(i);
        j_plus_final = num_microtubule_segment(i) - empty_ring_plus(i);
        
        j_minus_initial = 1;
        j_minus_final = 1 + empty_ring_minus(i);
        
        rate_depoly_plus (i) = ...
            (j_plus_initial - j_plus_final)*lattice_unit/(t_final - t_initial);
        
        rate_depoly_minus (i) = ...
            (j_minus_final - j_minus_initial)*lattice_unit/(t_final - t_initial);
        
    end
    
end

%This function named Measure_microtubule_length calculates and provides the 
%microtubule_length and the coordinates of microtubules plus and minus ends 
%at each iteration as the output.

function [microtubule_length,microtubule_end_plus,microtubule_end_minus] = ...
    Measure_microtubule_length (N_MT,lattice_unit,num_microtubule_segment,...
    empty_ring_plus,empty_ring_minus)
    
    microtubule_length = zeros(N_MT,1);

    for i = 1:N_MT

        j_plus = num_microtubule_segment(i) - empty_ring_plus(i);
        j_minus = 1 + empty_ring_minus(i);
    
        microtubule_length(i) = (j_plus - j_minus)*lattice_unit;
        
        microtubule_end_plus = j_plus*lattice_unit;
        microtubule_end_minus = j_minus*lattice_unit;
       
    
    end
    
end















