function [number_of_objectives, number_of_decision_variables, min_range_of_decesion_variable, max_range_of_decesion_variable] = objective_description_function()
number_of_objectives = 2;
% number_of_objectives = 3;
number_of_decision_variables = 10;
%% FDA1,FDA2,FDA3,DF5,DF6,DF8,DF9,JY1,JY2,JY3,JY4,JY5,JY6,JY7,JY8,DIMP1
% min_range_of_decesion_variable = -1*ones(1,number_of_decision_variables);
% min_range_of_decesion_variable(1)=0;
% max_range_of_decesion_variable =ones(1,number_of_decision_variables); 
%% DF3,
min_range_of_decesion_variable = -1*ones(1,number_of_decision_variables);
min_range_of_decesion_variable(1)=0;
max_range_of_decesion_variable =2*ones(1,number_of_decision_variables); 
max_range_of_decesion_variable(1)=1;
%% FDA4(ÈýÎ¬£©,FDA5(ÈýÎ¬£©,DF1,DF2,DF11,DMOP1,DMOP2,DMOP3,
% min_range_of_decesion_variable = zeros(1,number_of_decision_variables);
% max_range_of_decesion_variable =ones(1,number_of_decision_variables); 
%% F5,F6,F7,
% min_range_of_decesion_variable = zeros(1,number_of_decision_variables);
% max_range_of_decesion_variable = 5*ones(1,number_of_decision_variables); 
%% DF4,
% min_range_of_decesion_variable = -2*ones(1,number_of_decision_variables);
% max_range_of_decesion_variable =2*ones(1,number_of_decision_variables); 
%% DF7
% min_range_of_decesion_variable = zeros(1,number_of_decision_variables);
% min_range_of_decesion_variable(1)=1; 
% max_range_of_decesion_variable=ones(1,number_of_decision_variables);
% max_range_of_decesion_variable(1)=4;
%% %% DIMP2
% min_range_of_decesion_variable = -2*ones(1,number_of_decision_variables);
% min_range_of_decesion_variable(1)=0;
% max_range_of_decesion_variable =2*ones(1,number_of_decision_variables); 
% max_range_of_decesion_variable(1)=1;
%% DF10.DF12,DF13,DF14
% min_range_of_decesion_variable = -1*ones(1,number_of_decision_variables);
% min_range_of_decesion_variable(1)=0;
% min_range_of_decesion_variable(2)=0;
% max_range_of_decesion_variable =1*ones(1,number_of_decision_variables); 
g = sprintf('\n Now edit the function named "evaluate_objective" appropriately to match your needs.\n Make sure that the number of objective functions and decision variables match your numerical input. \n Make each objective function as a corresponding array element. \n After editing do not forget to save. \n Press "c" and enter to continue... ');

end

